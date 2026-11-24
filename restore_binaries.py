#!/usr/bin/env python3
"""Reassemble split text parts and base64-encoded binaries. Standard library only."""
from __future__ import annotations

import base64
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[tuple[str, str]]:
    items: list[tuple[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        digest, name = line.split(None, 1)
        items.append((digest, name))
    return items


def assemble_numbered_parts(prefix: Path) -> bytes:
    """Read prefix.part001, prefix.part002, ... in order."""
    chunks: list[bytes] = []
    index = 1
    while True:
        part = prefix.parent / f"{prefix.name}.part{index:03d}"
        if not part.exists():
            break
        chunks.append(part.read_bytes())
        index += 1
    if not chunks:
        raise FileNotFoundError(f"no parts for {prefix}")
    # Do not skip a hole: if part001 and part003 exist but part002 does not, stop at 1.
    leftover = sorted(prefix.parent.glob(prefix.name + ".part*"))
    expected = {f"{prefix.name}.part{i:03d}" for i in range(1, index)}
    unexpected = [p.name for p in leftover if p.name not in expected]
    if unexpected:
        raise SystemExit(f"unexpected part names for {prefix}: {unexpected}")
    return b"".join(chunks)


def restore_one(rel: str) -> bytes:
    dest_hint = ROOT / rel
    text_part1 = ROOT / f"{rel}.part001"
    b64_single = ROOT / f"{rel}.b64"
    b64_part1 = ROOT / f"{rel}.b64.part001"
    if text_part1.exists():
        return assemble_numbered_parts(ROOT / rel)
    if b64_part1.exists():
        raw = assemble_numbered_parts(ROOT / f"{rel}.b64")
        return base64.b64decode(raw, validate=True)
    if b64_single.exists():
        return base64.b64decode(b64_single.read_bytes(), validate=True)
    if dest_hint.exists():
        return dest_hint.read_bytes()
    raise FileNotFoundError(f"no encoded sources for {rel}")


def main() -> int:
    manifest = ROOT / "BINARIES.sha256"
    if not manifest.exists():
        print("missing BINARIES.sha256", file=sys.stderr)
        return 1
    failed = False
    for expected, rel in read_manifest(manifest):
        dest = ROOT / rel
        try:
            data = restore_one(rel)
        except Exception as exc:
            print(f"MISSING {rel}: {exc}")
            failed = True
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        got = sha256_file(dest)
        if got != expected:
            print(f"MISMATCH {rel}")
            print(f"  expected {expected}")
            print(f"  got      {got}")
            failed = True
        else:
            print(f"OK {rel}")
    if failed:
        return 1
    print("All files restored and hashes match.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
