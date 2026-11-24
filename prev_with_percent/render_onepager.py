#!/usr/bin/env python3
"""Render the Riemann hourglass one-pager and the standalone figure."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.lines import Line2D
import numpy as np

OUT = Path(__file__).resolve().parent.parent

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Liberation Serif", "DejaVu Serif", "Nimbus Roman"],
        "mathtext.fontset": "dejavuserif",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
    }
)

# Hourglass geometry in data coordinates.
# Symmetric polygonal hourglass, two trapezoids locked at a short waist.
Y_TOP, Y_UP_BULGE, Y_WAIST = 3.15, 1.72, 0.0
Y_LO_BULGE, Y_BOT = -1.72, -3.15
X_TOP = 1.55
X_BULGE = 2.15
X_WAIST = 0.64

UPPER = np.array(
    [
        [-X_TOP, Y_TOP],
        [X_TOP, Y_TOP],
        [X_BULGE, Y_UP_BULGE],
        [X_WAIST, Y_WAIST],
        [-X_WAIST, Y_WAIST],
        [-X_BULGE, Y_UP_BULGE],
    ]
)
LOWER = np.array(
    [
        [-X_WAIST, Y_WAIST],
        [X_WAIST, Y_WAIST],
        [X_BULGE, Y_LO_BULGE],
        [X_TOP, Y_BOT],
        [-X_TOP, Y_BOT],
        [-X_BULGE, Y_LO_BULGE],
    ]
)

# First few imaginary parts of zeta zeros, used only as spacing for the
# eigenvalue dots (a likeness, not a computation).
ZERO_T = np.array([14.135, 21.022, 25.011, 30.425, 32.935, 37.586, 40.919, 43.327])

PRIMES = [2, 3, 5, 7, 11]


def _serif(size, weight="normal", style="normal", color="#111"):
    return dict(
        fontsize=size,
        fontweight=weight,
        fontstyle=style,
        color=color,
        fontfamily="serif",
    )


def draw_hourglass(ax, standalone=False):
    """Draw the hourglass into ax. Leaves a translation-surface feel, not clipart."""
    ax.set_aspect("equal")
    ax.axis("off")

    fill_up = "#f3efe6"
    fill_lo = "#e7eef3"
    edge = "#111111"
    lw = 1.15

    ax.add_patch(
        Polygon(
            UPPER,
            closed=True,
            facecolor=fill_up,
            edgecolor=edge,
            linewidth=lw,
            joinstyle="miter",
            zorder=2,
        )
    )
    ax.add_patch(
        Polygon(
            LOWER,
            closed=True,
            facecolor=fill_lo,
            edgecolor=edge,
            linewidth=lw,
            joinstyle="miter",
            zorder=2,
        )
    )

    # Interior dashed lines, like cylinder cuts on a glued polygon.
    def hcut(y, xspan, lo=0.12):
        ax.plot(
            [-xspan, xspan],
            [y, y],
            ls=(0, (2.2, 2.0)),
            lw=0.7,
            color="#444",
            zorder=3,
            solid_capstyle="butt",
        )

    hcut(2.48, 1.68)
    hcut(-2.48, 1.68)

    # Vertices as small filled discs.
    verts = np.vstack([UPPER, LOWER[2:]])
    ax.scatter(
        verts[:, 0],
        verts[:, 1],
        s=22,
        c="#111",
        zorder=5,
        linewidths=0,
        clip_on=False,
    )

    # Side-height dimension on the left: the two heights whose average is the waist.
    dim_x = -3.05
    tick = 0.10
    ax.plot([dim_x, dim_x], [Y_BOT, Y_TOP], color="#222", lw=0.6, zorder=4, clip_on=False)
    for y in (Y_TOP, Y_WAIST, Y_BOT):
        ax.plot(
            [dim_x - tick, dim_x + tick],
            [y, y],
            color="#222",
            lw=0.6,
            zorder=4,
            clip_on=False,
        )
    ax.text(
        dim_x - 0.18,
        (Y_TOP + Y_WAIST) / 2,
        "height",
        ha="right",
        va="center",
        rotation=90,
        **_serif(8.0, color="#222"),
    )
    ax.text(
        dim_x - 0.18,
        (Y_BOT + Y_WAIST) / 2,
        "height",
        ha="right",
        va="center",
        rotation=90,
        **_serif(8.0, color="#222"),
    )

    # Average-height line: the locked waist, extended to the right.
    x_line_l, x_line_r = -2.55, 6.15
    ax.plot(
        [x_line_l, x_line_r],
        [Y_WAIST, Y_WAIST],
        ls=(0, (3.4, 2.2)),
        lw=0.95,
        color="#111",
        zorder=4,
        clip_on=False,
        solid_capstyle="butt",
    )

    # Prime gatekeepers sitting on the waist segment itself.
    pxs = np.linspace(-X_WAIST + 0.08, X_WAIST - 0.16, len(PRIMES))
    for x, p in zip(pxs, PRIMES):
        ax.plot(
            [x, x],
            [-0.18, 0.18],
            color="#111",
            lw=1.2,
            solid_capstyle="butt",
            zorder=6,
        )
        ax.plot(x, 0, "o", ms=3.4, color="#111", zorder=7, clip_on=False)
        ax.text(
            x,
            0.28,
            str(p),
            ha="center",
            va="bottom",
            **_serif(7.0, color="#111"),
        )
    ax.text(
        pxs[-1] + 0.18,
        0.28,
        "...",
        ha="left",
        va="bottom",
        **_serif(7.0, color="#111"),
    )
    ax.text(
        0.0,
        -0.62,
        "primes as gatekeepers",
        ha="center",
        va="top",
        **_serif(8.0, color="#111"),
    )

    # Eigenvalue dots on the extended waist line (to the right of the glass).
    xs = np.linspace(2.55, 5.85, len(ZERO_T))
    # Irregular spacing, scaled from the first zeta ordinates, then recentered.
    t = ZERO_T - ZERO_T[0]
    t = t / t[-1]
    eq = np.linspace(0, 1, len(t))
    t = 0.7 * t + 0.3 * eq  # keep irregular spacing, avoid a pile-up at the end
    xs = 2.55 + t * (5.85 - 2.55)
    ax.scatter(xs, np.zeros_like(xs), s=22, c="#111", zorder=7, linewidths=0, clip_on=False)
    ax.annotate(
        "eigenvalues of a self-adjoint\noperator  (Hilbert-Polya)",
        xy=(xs[4], 0.0),
        xytext=(4.25, 1.15),
        ha="center",
        va="bottom",
        arrowprops=dict(
            arrowstyle="-|>",
            color="#333",
            lw=0.6,
            shrinkA=0,
            shrinkB=3,
        ),
        **_serif(8.0, color="#111"),
    )

    # Labels for the two sides.
    ax.text(
        0.0,
        Y_TOP + 0.38,
        "entropic side",
        ha="center",
        va="bottom",
        **_serif(11.0, weight="bold", color="#111"),
    )
    ax.text(
        0.0,
        Y_BOT - 0.38,
        "unentropic side",
        ha="center",
        va="top",
        **_serif(11.0, weight="bold", color="#111"),
    )

    # Average-height callouts, to the right but above/below the eigenvalue note.
    ax.text(
        6.28,
        0.22,
        "average height",
        ha="left",
        va="bottom",
        **_serif(9.5, weight="bold", color="#111"),
    )
    ax.text(
        6.28,
        -0.12,
        "hoped-for minimal\nnon-reducible equation\ncomes from here",
        ha="left",
        va="top",
        **_serif(8.0, color="#111"),
        linespacing=1.25,
    )
    ax.text(
        6.28,
        -1.55,
        "the same equation\nwould also locate the\nnon-trivial zeros",
        ha="left",
        va="top",
        **_serif(8.0, color="#111"),
        linespacing=1.25,
    )

    # Locked-at-waist, quiet and small, left of the neck just outside the glass.
    ax.text(
        -1.55,
        -0.28,
        "locked at the waist",
        ha="right",
        va="top",
        **_serif(7.4, style="italic", color="#333"),
    )

    # Vision badge under the lower half, centered.
    ax.text(
        0.0,
        Y_BOT - 1.05,
        "A vision, not a proof.",
        ha="center",
        va="top",
        **_serif(10.0, style="italic", color="#111"),
    )

    # Limits leave room for labels.
    ax.set_xlim(-3.85, 9.15)
    ax.set_ylim(-4.35, 3.85)
    ax.set_anchor("C")


def render_standalone():
    fig = plt.figure(figsize=(11.2, 6.4), facecolor="white")
    ax = fig.add_axes([0.02, 0.10, 0.96, 0.80])
    draw_hourglass(ax, standalone=True)
    fig.text(
        0.50,
        0.965,
        "The hourglass picture",
        ha="center",
        va="top",
        **_serif(16, weight="bold"),
    )
    fig.text(
        0.50,
        0.925,
        "A vision, not a proof.  Primes sit at the waist as gatekeepers.  "
        "The average height is where a hoped-for equation, and the zeros, would live.",
        ha="center",
        va="top",
        **_serif(9.0, color="#222"),
    )
    fig.text(
        0.50,
        0.022,
        "Heywood Geblomi  |  from the conversation "
        "Riemann Hypothesis Algebraic Solution, opened August 19, 2026",
        ha="center",
        va="bottom",
        **_serif(8.0, color="#333"),
    )
    png = OUT / "hourglass_figure.png"
    svg = OUT / "hourglass_figure.svg"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    return png, svg


def wrap_words(text, width):
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
    return "\n".join(lines)


def render_page():
    fig = plt.figure(figsize=(8.5, 11), facecolor="white")

    # Title block.
    fig.text(
        0.50,
        0.975,
        "The hourglass picture",
        ha="center",
        va="top",
        **_serif(17.5, weight="bold"),
    )
    fig.text(
        0.50,
        0.954,
        "A vision, not a proof, for a minimal non-reducible equation that would locate the zeros of zeta",
        ha="center",
        va="top",
        **_serif(9.0, color="#222"),
    )
    fig.text(
        0.50,
        0.937,
        "Heywood Geblomi",
        ha="center",
        va="top",
        **_serif(9.5, weight="bold"),
    )
    fig.text(
        0.50,
        0.922,
        "From the conversation Riemann Hypothesis Algebraic Solution, opened August 19, 2026, in the Grok project Riemann.",
        ha="center",
        va="top",
        **_serif(7.8, color="#333"),
    )

    # The figure carries the page, placed near the top.
    ax = fig.add_axes([0.02, 0.480, 0.96, 0.430])
    draw_hourglass(ax, standalone=False)

    # Body text under the figure.
    left = 0.075
    width_chars = 118
    y = 0.465

    def heading(title):
        nonlocal y
        fig.text(left, y, title, ha="left", va="top", **_serif(9.6, weight="bold"))
        y -= 0.018

    def para(text, size=8.4, leading=0.0138, color="#111"):
        nonlocal y
        lines = wrap_words(text, width_chars).split("\n")
        for line in lines:
            fig.text(left, y, line, ha="left", va="top", **_serif(size, color=color))
            y -= leading
        y -= 0.004

    def bullets(items, size=8.2, leading=0.0134):
        nonlocal y
        for item in items:
            lines = wrap_words(item, width_chars - 3).split("\n")
            fig.text(left, y, "•", ha="left", va="top", **_serif(size))
            for i, line in enumerate(lines):
                fig.text(
                    left + 0.018,
                    y,
                    line,
                    ha="left",
                    va="top",
                    **_serif(size),
                )
                y -= leading
        y -= 0.004

    heading("The idea")
    para(
        "Two sides lock together in an hourglass: an entropic side above and an unentropic side below, "
        "joined at a narrow waist. The primes sit there as gatekeepers. The hope is a minimal non-reducible "
        "equation, derived from the average height of the two sides, that would also locate the non-trivial "
        "zeros of the zeta function. In the conversation, Grok treated this picture as a vision, not a proof."
    )

    heading("What the conversation stated")
    # Formula as its own line so the typography stays clean.
    bullets(
        [
            "Grok's opening finding: there is no algebraic solution, and the Riemann Hypothesis was still unsolved as of August 19, 2026.",
            "The only literature echo named was Hilbert-Polya: the zeros as eigenvalues of a self-adjoint operator.",
        ]
    )
    # Insert the completed-zeta line with mathtext, then the share-of-zeros bullet.
    fig.text(left, y, "•", ha="left", va="top", **_serif(8.2))
    fig.text(
        left + 0.018,
        y,
        r"The completed zeta function was written  "
        r"$\xi(s)=\frac{1}{2}\,s(s-1)\,\pi^{-s/2}\,\Gamma(s/2)\,\zeta(s)$.",
        ha="left",
        va="top",
        **_serif(8.2),
    )
    y -= 0.0134
    bullets(
        [
            "The rigorously known share of non-trivial zeros on the critical line was stated as having risen from about 41.6% to about 67.2%.",
        ]
    )
    y -= 0.002

    heading("A later turn, kept separate")
    para(
        "The thread later turned to Prym and discriminant-12 work: McMullen's S-shaped translation surface, "
        "residual-zero real multiplication, and Lyapunov growth along a Teichmuller stretch. That work was "
        "treated as a specialist calculator, not a Riemann Hypothesis machine."
    )

    heading("What to do next")
    bullets(
        [
            "Formalize the average-height equation, with a precise meaning of height on each side and of the primes as gatekeepers.",
            "Check whether the Hilbert-Polya operator can be written down explicitly, so the dots on the waist are eigenvalues of a real operator, not only a likeness.",
            "Connect the 2:1 geodesic ratios (20 exact, known sum 6/7 with a related plus-or-minus 1/7 piece, and candidates 4/7 and 2/7) to the height measure, and say whether they inform the average or stay in the calculator.",
            "Keep the discriminant-12 Prym calculator in its own lane: use it to test that connection, and do not treat it as a proof of the hypothesis.",
            "Write a one-line claim that can fail: what would count as finding the equation, and what would count as missing it.",
        ]
    )

    fig.text(
        0.50,
        0.018,
        "One page.  The figure is a picture of the vision, not a derivation.",
        ha="center",
        va="bottom",
        **_serif(7.2, color="#444"),
    )

    pdf = OUT / "ONEPAGER_riemann_hourglass.pdf"
    png = OUT / "ONEPAGER_riemann_hourglass.png"
    fig.savefig(pdf, dpi=200, facecolor="white")
    fig.savefig(png, dpi=200, facecolor="white")
    plt.close(fig)
    return pdf, png


def main():
    png_s, svg_s = render_standalone()
    pdf, png = render_page()
    print("standalone", png_s, svg_s)
    print("page", pdf, png)


if __name__ == "__main__":
    main()
