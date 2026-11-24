# riemann-hourglass

Author: Heywood Geblomi. Public snapshot built 24 November 2026.

## What this is

A one-page writeup of an hourglass picture for the zeros of the Riemann zeta function, plus a later phase-2 note that tries to attach standard objects (the completed xi function, the explicit formula, the smooth term of N(T)) to that picture.

This is a vision, not a proof. It does not prove the Riemann Hypothesis, and it does not introduce a new equation that locates the zeros. The original one-pager says so on the first line.

Phase 2 (`phase2/REPORT.md`) uses standard mathematics and a computation on the first 1000 non-trivial zeros (mpmath, 30 decimal digits). Read that report for what was formalized and what was not claimed.

The figure files are listed below. The SVG source is included.

## How to read it

Start with `ONEPAGER_riemann_hourglass.md` and `hourglass_figure.svg`. Then read `phase2/REPORT.md`. The Python under `src/` and `phase2/src/` renders the one-pagers.

A previous wording that included percents is under `prev_with_percent/`.

## Added data

Text under `phase2/out/` is in this repo (`CLAIM.txt`, `results.json`, `zeros_meta.json`). In `results.json` only, the author line was rewritten from a two-name form to `Heywood Geblomi`.

These binaries could not be uploaded: the GitHub connector used here accepts text only, and `$file:` rejection for non-text files blocked them. I did not use `gh` or a raw API token as a workaround.

- `ONEPAGER_riemann_hourglass.pdf`
- `ONEPAGER_riemann_hourglass.png`
- `hourglass_figure.png`
- `phase2/ONEPAGER_phase2.pdf`
- `phase2/ONEPAGER_phase2.png`
- `prev_with_percent/ONEPAGER_riemann_hourglass.pdf`
- `prev_with_percent/ONEPAGER_riemann_hourglass.png`
- `phase2/out/plot_data.npz`
- `phase2/out/zeros_t.npy`

## Restoring large and binary files

The one-pager PDFs and PNGs, plus `phase2/out/plot_data.npz` and `phase2/out/zeros_t.npy`, could not be stored as binary blobs through the connector used to publish this repo. They are here as base64 text in matching `*.b64` files.

`BINARIES.sha256` lists the sha256 of each original file. To restore those originals, run:

```
python3 restore_binaries.py
```

That script uses only the Python standard library. It decodes the base64 files, writes the original paths, and checks each sha256.
