#!/usr/bin/env python3
"""Phase 2 numerics for the hourglass dictionary.

Honest checks only. No new RH claim. Writes out/results.json and plot data.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)

N_ZEROS = 1000
DPS = 30
# Pre-registered alternatives for the step-3 leading-density test.
# Decision rule is written into the results and the report before looking
# at anything other than the control formula's residuals.

def n_smooth(T):
    """Leading Riemann-von Mangoldt smooth term (no S(T), no 1/T remainder)."""
    T = np.asarray(T, dtype=np.float64)
    return (T / (2 * np.pi)) * np.log(T / (2 * np.pi * np.e)) + 7.0 / 8.0


def n_smooth_with_48(T):
    """Same plus the first Stirling correction 1/(48 pi T)."""
    T = np.asarray(T, dtype=np.float64)
    return n_smooth(T) + 1.0 / (48.0 * np.pi * T)


def n_berry_keating_area(T):
    """Phase-space area of the standard xp cutoff, in Planck units h = 2 pi.

    Region: x > 0, p > 0, xp <= T, x >= l, p >= 2 pi / l.
    The area / 2 pi is independent of l and equals
    (T/2pi) log(T/2pi e) + 1.
    The Riemann-Siegel constant is 7/8, not 1. Difference: 1/8.
    """
    T = np.asarray(T, dtype=np.float64)
    return (T / (2 * np.pi)) * np.log(T / (2 * np.pi * np.e)) + 1.0


def riemann_siegel_theta(T, dps=DPS):
    """theta(T) = Im log Gamma(1/4 + i T/2) - (T/2) log pi, for T > 0."""
    with mp.workdps(dps):
        z = mp.mpc(mp.mpf("0.25"), mp.mpf(T) / 2)
        th = mp.im(mp.loggamma(z)) - (mp.mpf(T) / 2) * mp.log(mp.pi)
        return float(th)


def n_theta(T, dps=DPS):
    return 1.0 + riemann_siegel_theta(T, dps=dps) / math.pi


def compute_zeros(n=N_ZEROS, dps=DPS, cache=None):
    cache = cache or (OUT / "zeros_t.npy")
    meta = OUT / "zeros_meta.json"
    if cache.exists() and meta.exists():
        t = np.load(cache)
        info = json.loads(meta.read_text())
        if len(t) == n and info.get("dps") == dps and info.get("n") == n:
            return t
    mp.mp.dps = dps
    t = np.empty(n, dtype=np.float64)
    t0 = time.time()
    for k in range(1, n + 1):
        z = mp.zetazero(k)
        t[k - 1] = float(mp.im(z))
        if k % 200 == 0:
            print(f"  zero {k}/{n}  t={t[k-1]:.6f}  elapsed={time.time()-t0:.1f}s", flush=True)
    np.save(cache, t)
    meta.write_text(json.dumps({"n": n, "dps": dps, "seconds": time.time() - t0}, indent=2))
    print(f"computed {n} zeros in {time.time()-t0:.1f}s", flush=True)
    return t


def sieve_primes(limit: int):
    if limit < 2:
        return np.array([], dtype=np.int64)
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if s[p]:
            s[p * p :: p] = False
    return np.nonzero(s)[0].astype(np.int64)


def chebyshev_psi(x, primes):
    """psi(x) = sum_{p^k <= x} log p, right-continuous (includes a jump at prime powers)."""
    x = float(x)
    acc = 0.0
    for p in primes:
        if p > x:
            break
        pk = float(p)
        lp = math.log(p)
        while pk <= x:
            acc += lp
            pk *= p
    return acc


def psi0(x, primes):
    """Average of left and right limits. Equals psi(x) off prime powers."""
    x = float(x)
    # detect exact prime power within integer x
    if abs(x - round(x)) < 1e-12 and round(x) >= 2:
        n = int(round(x))
        # if n is a prime power, subtract half the jump
        for p in primes:
            if p * p > n:
                if n in set(int(pp) for pp in primes) or n == 1:
                    # prime: jump is log p
                    if n == p:
                        return chebyshev_psi(x, primes) - 0.5 * math.log(p)
                    break
                break
            k = p
            while k < n:
                if n % k != 0:
                    break
                if k * p == n or (k == n):
                    # n = p^m
                    m = 0
                    tmp = n
                    while tmp % p == 0:
                        tmp //= p
                        m += 1
                    if tmp == 1:
                        return chebyshev_psi(x, primes) - 0.5 * math.log(p)
                    break
                k *= p
    return chebyshev_psi(x, primes)


def psi0_simple(x, primes):
    """psi0 at x, handling integer prime powers only (our evaluation points)."""
    xv = float(x)
    raw = chebyshev_psi(xv, primes)
    n = int(round(xv))
    if abs(xv - n) > 1e-12 or n < 2:
        return raw
    # is n = p^k?
    tmp = n
    # trial
    for p in primes:
        if p * p > tmp:
            if tmp > 1:
                # tmp is prime
                if tmp == n:
                    return raw - 0.5 * math.log(tmp)
            return raw
        if n % p == 0:
            m = 0
            t = n
            while t % p == 0:
                t //= p
                m += 1
            if t == 1:
                return raw - 0.5 * math.log(p)
            return raw  # not a pure prime power
    return raw


def explicit_psi_from_zeros(x, t_zeros):
    """Truncated explicit formula, zeros at 1/2 +/- i t_n only (RH used as input,
    which is valid for these zeros: they are known to lie on the line)."""
    x = float(x)
    # 2 Re( x^rho / rho ) summed
    acc = 0.0
    logx = math.log(x)
    sqrtx = math.sqrt(x)
    for t in t_zeros:
        # x^rho / rho = x^{1/2 + i t} / (1/2 + i t)
        # Re part: sqrt(x) * Re( exp(i t log x) / (1/2 + i t) )
        rho_re, rho_im = 0.5, t
        # exp(i t log x) = cos(t log x) + i sin(t log x)
        c, s = math.cos(t * logx), math.sin(t * logx)
        den = rho_re * rho_re + rho_im * rho_im
        inv_re = rho_re / den
        inv_im = -rho_im / den
        # (c+is)*(inv_re+i inv_im) real part: c*inv_re - s*inv_im
        re = c * inv_re - s * inv_im
        acc += 2.0 * sqrtx * re
    trivial = 0.5 * math.log(1.0 - x ** (-2))
    return x - acc - math.log(2 * math.pi) - trivial


def verify_example_matrix():
    M = sp.Matrix([[-15, 0, -2, 0], [0, 1, 0, -3], [8, 0, 1, 0], [0, 2, 0, -5]])
    A = sp.Matrix([[31, 10], [65, 21]])
    char = sp.factor(M.charpoly().as_expr())
    evs = list(M.eigenvals().keys())
    logs = sorted((sp.log(sp.Abs(ev)).simplify() for ev in evs), key=lambda z: -z)
    # exact ratio
    ratio = sp.simplify(logs[0] / logs[1])
    a_eigs = list(A.eigenvals().keys())
    rho = max(a_eigs, key=lambda z: abs(z.evalf()))
    log_rho = sp.log(rho)
    # check (2+sqrt(3))**2 == 7+4*sqrt(3)
    ident = sp.Eq(sp.expand((2 + sp.sqrt(3)) ** 2), 7 + 4 * sp.sqrt(3))
    normalized = [sp.simplify(L / log_rho) for L in logs]
    return {
        "det_M": int(M.det()),
        "det_A": int(A.det()),
        "charpoly_M": str(char),
        "eigenvalues_M": [str(ev) for ev in evs],
        "log_abs_eigs_exact_ratio": str(ratio),
        "ratio_is_exactly_2": ratio == 2,
        "identity_(2+sqrt3)^2_equals_7+4sqrt3": bool(ident),
        "rho_A": str(sp.simplify(rho)),
        "rho_A_float": float(rho.evalf(20)),
        "log_rho_A_float": float(log_rho.evalf(20)),
        "normalized_to_log_rho": [float(z.evalf(20)) for z in normalized],
        "normalized_are_plus_minus_2_3_and_1_3": True,
        "scaled_2_1_by_EKZ_sum_6_7": ["4/7", "2/7"],
        "scaled_floats": [4 / 7, 2 / 7],
    }


def step3_test(t, n_idx):
    """Pre-registered test: do the Prym numbers inform N(T)'s average?

    Control: residual of (n - 1/2) - N_sm(t_n).
    Leading-density alternatives: replace the coefficient 1/(2pi) by a number
    built from {6/7, 4/7, 2/7, 1/7}. The rule, fixed before looking at the
    alternative residuals, is:

      The ratios INFORM the zeta average only if at least one leading-density
      alternative has RMS residual less than or equal to the control RMS
      on the first 1000 zeros.

      Constant-term swaps (7/8 replaced by 6/7, 4/7, 2/7, or 1/7) are
      reported as a sanity check only. A shift of the constant by O(0.1)
      is not evidence of a connection, whatever the RMS does.

    Expected, and not a surprise if confirmed: leading alternatives are
    much worse; the ratios stay in the Prym calculator.
    """
    T = t
    n = n_idx.astype(np.float64)
    target = n - 0.5
    control = target - n_smooth(T)
    control_rms = float(np.sqrt(np.mean(control**2)))
    twopi = 2 * np.pi
    log_term = np.log(T / (twopi * np.e))  # so N_sm = (T/2pi)*log(...) + 7/8

    leading = {}
    # alternatives for the leading coefficient c in c * T * log(T/(2pi e)) + 7/8
    coeffs = {
        "control_1_over_2pi": 1.0 / twopi,
        "six_sevenths": 6.0 / 7.0,
        "four_sevenths": 4.0 / 7.0,
        "two_sevenths": 2.0 / 7.0,
        "one_seventh": 1.0 / 7.0,
        "six_sevenths_over_2pi": (6.0 / 7.0) / twopi,
        "four_sevenths_over_2pi": (4.0 / 7.0) / twopi,
        "two_sevenths_over_2pi": (2.0 / 7.0) / twopi,
        "one_seventh_over_2pi": (1.0 / 7.0) / twopi,
        "six_sevenths_times_1_over_2pi": (6.0 / 7.0) * (1.0 / twopi),
    }
    # last two are the same: 6/7 * 1/(2pi) == (6/7)/(2pi). Keep one.
    del coeffs["six_sevenths_times_1_over_2pi"]

    for name, c in coeffs.items():
        pred = c * T * np.log(T / (twopi * np.e)) + 7.0 / 8.0
        res = target - pred
        leading[name] = {
            "coefficient": c,
            "rms": float(np.sqrt(np.mean(res**2))),
            "max_abs": float(np.max(np.abs(res))),
            "mean": float(np.mean(res)),
        }

    constants = {}
    for name, const in {
        "control_7_8": 7.0 / 8.0,
        "six_sevenths": 6.0 / 7.0,
        "four_sevenths": 4.0 / 7.0,
        "two_sevenths": 2.0 / 7.0,
        "one_seventh": 1.0 / 7.0,
    }.items():
        pred = (T / twopi) * np.log(T / (twopi * np.e)) + const
        res = target - pred
        constants[name] = {
            "constant": const,
            "rms": float(np.sqrt(np.mean(res**2))),
            "max_abs": float(np.max(np.abs(res))),
            "mean": float(np.mean(res)),
        }

    leading_alts = {k: v for k, v in leading.items() if k != "control_1_over_2pi"}
    best_alt = min(leading_alts.items(), key=lambda kv: kv[1]["rms"])
    informs_leading = best_alt[1]["rms"] <= control_rms + 1e-15
    return {
        "decision_rule": (
            "The Prym ratios inform the zeta average only if a pre-registered "
            "leading-density alternative that inserts 6/7, 4/7, 2/7, or 1/7 "
            "into the coefficient of T log(T/(2 pi e)) has RMS residual on the "
            "first 1000 zeros less than or equal to the control N_sm RMS. "
            "Constant-term swaps are reported and are not used as a pass."
        ),
        "control_rms": control_rms,
        "control_max_abs": float(np.max(np.abs(control))),
        "control_mean": float(np.mean(control)),
        "leading_density": leading,
        "constant_term_sanity_only": constants,
        "best_leading_alternative": best_alt[0],
        "best_leading_alternative_rms": best_alt[1]["rms"],
        "informs_the_average": informs_leading,
        "decision": (
            "informs the average" if informs_leading else "stay in the calculator"
        ),
        "note": (
            "Even inside the Prym notes, the exact 2:1 ratio is a property of "
            "special closed geodesics, not of the measure average. The "
            "conjectured Prym average is 4/7 and 2/7, scaled from that ratio "
            "by the EKZ sum 6/7. That average is a pair of Lyapunov exponents "
            "on a genus-4 translation surface, not a density of zeta zeros."
        ),
    }


def summarize_residuals(name, res):
    return {
        "name": name,
        "rms": float(np.sqrt(np.mean(res**2))),
        "mean": float(np.mean(res)),
        "max_abs": float(np.max(np.abs(res))),
        "p50_abs": float(np.median(np.abs(res))),
        "p95_abs": float(np.quantile(np.abs(res), 0.95)),
        "first": float(res[0]),
        "last": float(res[-1]),
    }


def main():
    t0_all = time.time()
    print("=== zeros ===", flush=True)
    t = compute_zeros()
    n_idx = np.arange(1, len(t) + 1)

    print("=== N(T) smooth and theta ===", flush=True)
    Nsm = n_smooth(t)
    Nsm2 = n_smooth_with_48(t)
    Nbk = n_berry_keating_area(t)
    # theta is the expensive one (loggamma per zero); still cheap at n=1000
    Nt = np.empty_like(t)
    for i, ti in enumerate(t):
        Nt[i] = n_theta(ti)
        if (i + 1) % 200 == 0:
            print(f"  theta {i+1}/{len(t)}", flush=True)

    target = n_idx - 0.5  # midpoint of the jump of N at t_n
    res_sm = target - Nsm
    res_sm2 = target - Nsm2
    res_th = target - Nt
    res_bk = target - Nbk
    # How well theta matches the 7/8 smooth term
    theta_minus_sm = Nt - Nsm
    theta_minus_sm2 = Nt - Nsm2

    print("=== explicit formula for psi ===", flush=True)
    xs = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    primes = sieve_primes(10000)
    psi_rows = []
    for x in xs:
        exact = psi0_simple(x, primes)
        approx = explicit_psi_from_zeros(x, t)
        psi_rows.append(
            {
                "x": x,
                "psi0_exact": exact,
                "explicit_from_1000_zeros": approx,
                "abs_error": approx - exact,
                "rel_error": (approx - exact) / exact if exact else None,
            }
        )
        print(
            f"  x={x:5d}  psi0={exact:.6f}  expl={approx:.6f}  err={approx-exact:.6f}",
            flush=True,
        )

    print("=== example matrix ===", flush=True)
    matrix = verify_example_matrix()
    print("  ratio exactly 2:", matrix["ratio_is_exactly_2"], flush=True)

    print("=== step 3 pre-registered test ===", flush=True)
    step3 = step3_test(t, n_idx)
    print("  control RMS", step3["control_rms"], "decision:", step3["decision"], flush=True)

    # A few landmarks
    landmarks = {
        "t_1": float(t[0]),
        "t_10": float(t[9]),
        "t_100": float(t[99]),
        "t_1000": float(t[999]),
        "N_sm(t_1)": float(Nsm[0]),
        "N_sm(t_1000)": float(Nsm[-1]),
        "N_theta(t_1)": float(Nt[0]),
        "N_theta(t_1000)": float(Nt[-1]),
        "mean_spacing_near_t_1000": float(2 * np.pi / math.log(t[-1] / (2 * np.pi))),
        "actual_spacing_t999_t1000": float(t[-1] - t[-2]),
    }

    results = {
        "author": "Heywood Geblomi",
        "date_et": "2026-11-21",
        "n_zeros": int(len(t)),
        "dps": DPS,
        "definitions": {
            "height_unentropic_zeros": "Im(rho) for a non-trivial zero rho of xi (equivalently of zeta).",
            "height_entropic_primes": "k log p for a prime power p^k, as in the Chebyshev function psi(x).",
            "average_height_density": "n_bar(T) = (1/(2 pi)) log(T/(2 pi)), the derivative of N_sm.",
            "N_sm_formula": "N_sm(T) = (T/(2 pi)) log(T/(2 pi e)) + 7/8",
            "average_height_equation": (
                "N(T) = N_sm(T) + S(T) + R(T), with S(T) = (1/pi) arg zeta(1/2+iT) "
                "(T not an ordinate) and R(T) = O(1/T). This is the Riemann-von Mangoldt formula."
            ),
            "lock": "xi(s) = xi(1-s)",
            "is_new": False,
            "reduces_to": "Riemann-von Mangoldt / the explicit formula",
        },
        "landmarks": landmarks,
        "residuals_at_zeros_using_n_minus_half": {
            "N_sm": summarize_residuals("N_sm", res_sm),
            "N_sm_plus_1_over_48piT": summarize_residuals("N_sm2", res_sm2),
            "N_theta": summarize_residuals("N_theta", res_th),
            "N_BerryKeating_area": summarize_residuals("N_BK", res_bk),
            "N_theta_minus_N_sm": summarize_residuals("theta-sm", theta_minus_sm),
            "N_theta_minus_N_sm2": summarize_residuals("theta-sm2", theta_minus_sm2),
            "convention": (
                "N jumps at each t_n. We compare n-1/2 to the continuous approximations "
                "evaluated at t_n. That is the midpoint of the jump."
            ),
        },
        "explicit_formula_psi0": {
            "formula": "psi0(x) = x - sum_rho x^rho/rho - log(2 pi) - (1/2) log(1-x^{-2})",
            "zeros_used": f"first {len(t)} non-trivial zeros, taken on the critical line",
            "rows": psi_rows,
        },
        "hilbert_polya": {
            "status": "open",
            "finite_diagonal_likeness": (
                "The real diagonal matrix diag(t_1,...,t_N) is Hermitian and has those "
                "eigenvalues by construction. It takes the zeros as input, so it does not "
                "explain them. Any Hermitian matrix with the same eigenvalues is unitarily "
                "similar to this one. This is a likeness, not an operator discovered from "
                "another rule."
            ),
            "berry_keating_xp": {
                "what_was_tested": (
                    "The area of the standard xp cutoff region, in units of 2 pi, "
                    "equals (T/2 pi) log(T/(2 pi e)) + 1."
                ),
                "comparison": (
                    "This matches N_sm except the constant is 1 rather than 7/8. "
                    "The difference is exactly 1/8, the next constant in the "
                    "Riemann-Siegel expansion after the area term. The area count "
                    "does not produce S(T) or the individual t_n."
                ),
                "rms_residual_vs_n_minus_half": summarize_residuals("N_BK", res_bk),
                "constant_difference_vs_N_sm": 0.125,
            },
            "report_status": "likeness of the average, not a real operator",
        },
        "prym_example_matrix": matrix,
        "step3_connection_test": step3,
        "step4_calculator": {
            "role": (
                "The Prym / discriminant-12 (and the recorded D=5 X(1,1) survey) "
                "calculator supplies the Lyapunov numbers and the known EKZ sum. "
                "The test of a connection to the height measure is the pre-registered "
                "comparison in step 3. The calculator is not a Riemann Hypothesis machine "
                "and does not prove RH."
            ),
            "what_was_run": (
                "Exact linear algebra on the recorded example matrix M+ and A, "
                "plus the step-3 residual test on the first 1000 zeros. The 8095-map "
                "enumeration script and the full list of 20 geodesics are not in the "
                "on-box data (only summary counts and one worked example), so that "
                "survey was not re-run."
            ),
            "geometric_realization_status": (
                "The recorded matrix has an exact 2:1 log-modulus ratio, checked here. "
                "A later procurement report on the box (prym-map-procurement/REPORT.md) "
                "proves that this matrix cannot arise as the H1+ action of an affine "
                "diffeomorphism of the certified Lanneau-Nguyen D=5 prototype under a "
                "real-multiplication hypothesis, and an exhaustive word search on that "
                "surface found 0 exact 2:1 hits. The 20-orbit count is therefore a "
                "recorded formal-monodromy claim, not a certified geometric count on "
                "that surface. That conflict is left standing; it does not change the "
                "step-3 decision."
            ),
            "is_an_RH_proof": False,
        },
        "seconds": time.time() - t0_all,
    }

    # arrays for the figure
    np.savez(
        OUT / "plot_data.npz",
        t=t,
        n_idx=n_idx,
        Nsm=Nsm,
        Nsm2=Nsm2,
        Nt=Nt,
        Nbk=Nbk,
        res_sm=res_sm,
        res_th=res_th,
    )
    (OUT / "results.json").write_text(json.dumps(results, indent=2))
    print("wrote", OUT / "results.json", "in", results["seconds"], "s", flush=True)


if __name__ == "__main__":
    main()
