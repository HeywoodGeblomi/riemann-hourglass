#!/usr/bin/env python3
"""One-page summary of phase 2. Does not touch the phase-1 one-pager files."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = np.load(ROOT / "out/plot_data.npz")

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Liberation Serif", "DejaVu Serif", "Nimbus Roman"],
        "mathtext.fontset": "dejavuserif",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.6,
    }
)

SERIF = dict(fontfamily="serif")


def wrap(text, width):
    words = text.split()
    lines, cur = [], ""
    for w in words:
        trial = w if not cur else cur + " " + w
        if len(trial) <= width:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


def draw_figure(fig):
    t = DATA["t"]
    n_idx = DATA["n_idx"]
    Nsm = DATA["Nsm"]
    res_sm = DATA["res_sm"]

    # Left: staircase of N(T) vs N_sm for the first 80 zeros.
    m = 80
    T_end = float(t[m - 1]) + 8
    ax = fig.add_axes([0.07, 0.620, 0.46, 0.230])
    # staircase
    Ts = [0.0]
    Ns = [0.0]
    for i in range(m):
        Ts.extend([t[i], t[i]])
        Ns.extend([i, i + 1])
    Ts.append(T_end)
    Ns.append(m)
    grid = np.linspace(max(8.0, 2 * np.pi + 0.2), T_end, 400)
    ax.plot(Ts, Ns, color="#111", lw=0.9, solid_capstyle="butt", label=r"$N(T)$, the zeros")
    ax.plot(
        grid,
        (grid / (2 * np.pi)) * np.log(grid / (2 * np.pi * np.e)) + 7.0 / 8.0,
        color="#1f4e79",
        lw=1.3,
        label=r"$N_{\mathrm{sm}}(T)$, the average",
    )
    ax.plot(t[:m], n_idx[:m] - 0.5, "o", ms=2.0, color="#111", zorder=3)
    ax.set_xlim(0, T_end)
    ax.set_ylim(-2, m + 6)
    ax.set_xlabel("height $T$ (imaginary part of a zero)", fontsize=7.2)
    ax.set_ylabel("$N(T)$", fontsize=7.2)
    ax.tick_params(labelsize=6.5, length=2.5)
    ax.legend(frameon=False, fontsize=6.8, loc="upper left")
    ax.set_title("The waist is the known smooth term of $N(T)$", fontsize=8.2, pad=3, loc="left")
    for s in ax.spines.values():
        s.set_color("#333")
        s.set_linewidth(0.6)

    # Right: residuals for all 1000 zeros, the leftover S(T).
    axr = fig.add_axes([0.62, 0.620, 0.325, 0.230])
    axr.axhline(0, color="#888", lw=0.5)
    axr.plot(n_idx, res_sm, color="#1f4e79", lw=0.55)
    axr.set_xlim(1, 1000)
    axr.set_ylim(-0.95, 0.95)
    axr.set_xlabel("zero index $n$", fontsize=7.2)
    axr.set_ylabel(r"$n-\frac{1}{2}-N_{\mathrm{sm}}(t_n)$", fontsize=7.2)
    axr.tick_params(labelsize=6.5, length=2.5)
    axr.set_title("Leftover on 1000 zeros: usual $S(T)$, RMS $0.26$", fontsize=8.2, pad=3, loc="left")
    for s in axr.spines.values():
        s.set_color("#333")
        s.set_linewidth(0.6)


def render():
    fig = plt.figure(figsize=(8.5, 11), facecolor="white")

    fig.text(0.50, 0.975, "The waist is a known formula", ha="center", va="top", fontsize=17, fontweight="bold", color="#111", **SERIF)
    fig.text(
        0.50,
        0.954,
        "Phase 2 of the hourglass picture: a definition, a check, and a claim that can fail",
        ha="center",
        va="top",
        fontsize=9.0,
        **SERIF,
        color="#222",
    )
    fig.text(
        0.50,
        0.937,
        "Heywood Geblomi",
        ha="center",
        va="top",
        fontsize=9.5,
        fontweight="bold", color="#111",
        **SERIF,
    )
    fig.text(
        0.50,
        0.922,
        "21 November 2026 (America/New_York).  First 1000 zeros computed with mpmath.  Not a proof of the Riemann Hypothesis.",
        ha="center",
        va="top",
        fontsize=7.5,
        **SERIF,
        color="#333",
    )

    lede = (
        "Once height is defined in the usual way, the average-height equation is the Riemann-von Mangoldt formula. It is not new. "
        "The first 1000 zeros follow the smooth term with a leftover of size about 0.26, which is the usual fluctuation S(T)."
    )
    y = 0.902
    for line in wrap(lede, 118):
        fig.text(0.07, y, line, ha="left", va="top", fontsize=8.3, color="#111", **SERIF)
        y -= 0.0138

    draw_figure(fig)

    left = 0.07
    y = 0.568

    def heading(title):
        nonlocal y
        fig.text(left, y, title, ha="left", va="top", fontsize=9.2, fontweight="bold", color="#111", **SERIF)
        y -= 0.0155

    def para(text, size=7.85, width=118, leading=0.0122):
        nonlocal y
        for line in wrap(text, width):
            fig.text(left, y, line, ha="left", va="top", fontsize=size, color="#111", **SERIF)
            y -= leading
        y -= 0.0035

    heading("1.  Definitions, and the equation they give")
    para(
        "Zeros: the height of a non-trivial zero is t = Im(rho). Primes: the height of a prime power p^k is k log p, as in psi(x). The lock is xi(s) = xi(1-s). "
        "The waist is the smooth density (1/2 pi) log(T/2 pi), with running count N_sm(T) = (T/2 pi) log(T/(2 pi e)) + 7/8. The primes gatekeep the waist through the explicit formula. "
        "Those definitions make the average-height equation N(T) = N_sm(T) + S(T) + O(1/T), the Riemann-von Mangoldt formula, not a new equation. On the first 1000 zeros, n minus 1/2 minus N_sm(t_n) has RMS 0.259 and mean nearly zero. The explicit formula for psi_0(x), truncated to these zeros, has relative error below 0.001 at every tested x from 20 to 10,000."
    )

    heading("2.  Hilbert-Polya: still open")
    para(
        "There is no proven self-adjoint operator, built from a rule that does not already use the zeros, whose eigenvalues are the ordinates t_n. "
        "A diagonal matrix of those ordinates is a likeness: it matches them because it takes them as input. A Berry-Keating xp area count recovers the same average as N_sm, except the constant is 1 rather than 7/8. It does not produce S(T) or the individual zeros. Status: likeness of the average, not a real operator."
    )

    heading("3.  The 2:1 ratios stay in the calculator")
    para(
        "In the Prym notes the 2:1 ratios are Lyapunov exponents along special closed Teichmuller geodesics, not a density of zeta zeros. The recorded example matrix has an exact 2:1 log-modulus ratio, because (2+sqrt(3))^2 = 7+4 sqrt(3). The notes record 20 such geodesics, a known sum 6/7, a related plus-or-minus 1/7 piece, and candidates 4/7 and 2/7. Even there, 2:1 is a property of special orbits, not of the average. "
        "A pre-registered test inserted 6/7, 4/7, 2/7, or 1/7 into the leading coefficient of N_sm. The control RMS was 0.259; the least bad leading alternative had RMS 59. The ratios stay in the calculator. No link was forced."
    )

    heading("4.  What the calculator is for, and what it is not")
    para(
        "The Prym work supplies those Lyapunov numbers and runs the test in section 3. It is not a Riemann Hypothesis machine and does not prove the hypothesis. The 8095-map survey was not re-run (only summary counts and one example are on the box). A later procurement report on a certified Lanneau-Nguyen D=5 prototype finds 0 exact 2:1 hits and shows that the recorded matrix cannot arise as an affine map of that surface under its real-multiplication hypothesis. That conflict is left standing. It does not change the step-3 decision."
    )

    heading("5.  A claim that can fail")
    para(
        "The average-height equation is found if, and only if, there is a minimal relation, not equivalent to the Riemann-von Mangoldt formula or the explicit formula, that is built from the heights t_n and k log p defined in this report, that fails if the primes are removed, and that locates the first 1000 non-trivial zero ordinates to a pre-registered error; it is missed if every such candidate either reduces to one of those two known identities or fails that location test."
    )
    para(
        "SHA-256 of out/CLAIM.txt (UTF-8, 472 bytes, one trailing newline), recorded before any further testing: d0a59770f542b002fefb880188f98de44b14f48e9d7eaa8714539899f51aa800. Under this claim, the verification of N(T) in step 1 is a miss for finding a new equation, because that formula is a known identity."
    )

    fig.text(
        0.50,
        0.018,
        "Files in /workspace/riemann_onepager/phase2/.  Existing one-pager files and the Riemann project screenshots were not changed.",
        ha="center",
        va="bottom",
        fontsize=7.0,
        color="#444",
        **SERIF,
    )

    pdf = ROOT / "ONEPAGER_phase2.pdf"
    png = ROOT / "ONEPAGER_phase2.png"
    fig.savefig(pdf, dpi=200, facecolor="white")
    fig.savefig(png, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", pdf, png)


if __name__ == "__main__":
    render()
