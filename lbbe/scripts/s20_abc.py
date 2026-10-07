"""
Slide 20 -- when the likelihood is unavailable: approximate Bayesian
computation (Nuismer & Week 2019, PLoS Comput. Biol. 15:e1006988).

Three things on the slide:
  left    the three assumptions of the 2019 likelihood that ABC relaxes:
          selection of any strength, gene flow among populations, abiotic
          optima that differ among populations;
  middle  the idea of ABC in three steps: draw parameters from their priors,
          simulate data with them, keep the draw if the simulated data are
          close to the observed (the paper's five summary statistics: the
          mean of each species' population means, the spread of each, and
          their correlation, each within a tolerance);
  right   the published result for Camellia japonica and its weevil
          Curculio camelliae: the posteriors of the two strengths of biotic
          selection, redrawn from the paper's Fig. 5 (bar heights read off the
          figure; modes 2.37 and 0.21 and 95 % intervals 0.60-2.94 and 0-2.40
          from its text).  Selection on the weevil is clearly above zero;
          selection on the camellia piles up toward zero, so the analysis
          supports coevolution but cannot exclude one-sided evolution.

scripts/abc_sim.py (a small port of the paper's C++ simulator) is kept in the
bundle but no longer used by this figure.
"""
import numpy as np
from matplotlib.patches import Circle, Ellipse

from glyphs import Pen, arrow, species
from semstyle import *  # noqa: F403

def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)

    # ---- what is relaxed ------------------------------------------------------
    cv.text(30, 508, "No longer assumed", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    rows = [(404, "weak selection", "any strength now"),
            (286, "isolated populations", "gene flow among them"),
            (168, "one optimum everywhere", "optima differ by place")]
    for y, was, now in rows:
        cv.text(150, y + 14, was, ha="left", va="center", fontsize=fs(FS_SMALL),
                color=MUTED)
        cv.plot([150, 150 + 9.3 * len(was)], [y + 14, y + 14], color=MUTED, lw=lw(1.4))
        cv.text(150, y - 14, now, ha="left", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold")
    ax = px_axes(fig, 36, 372, 92, 66)                     # selection of any strength
    d = np.linspace(-3, 3, 100)
    ax.plot(d, 0.5 + 0.12 * d, color=FAINT, lw=lw(2.0), ls=(0, (3, 2)))
    ax.plot(d, 1 / (1 + np.exp(-2.6 * d)), color=INK, lw=lw(2.8))
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-3, 3)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    pts = [(48, 300), (84, 318), (120, 296), (62, 262), (104, 262)]   # gene flow
    for i in range(5):
        for j in range(i + 1, 5):
            if (i + j) % 2 == 1 or j == i + 1:
                arrow(pen, pts[i], pts[j], px=1.2, color=MUTED, head=5, shrink=(8, 8),
                      both=True)
    for p in pts:
        species(pen, *p, 7, "N", edge_px=1.2)
    ax = px_axes(fig, 36, 138, 92, 62)                     # optima vary
    th = np.array([0.50, 0.72, 0.38, 0.62, 0.30, 0.80])
    ax.axhline(0.5, color=FAINT, lw=lw(2.0), ls=(0, (3, 2)))
    ax.plot(range(6), th, ls="none", marker="_", ms=13, mew=lw(3.0), color=INK)
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.6, 5.6)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    cv.text(30, 96, "dashed: what the likelihood assumed", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)

    # ---- the idea of ABC, in three steps ------------------------------------------
    X0 = 440
    cv.plot([X0 - 34, X0 - 34], [60, 520], color=LINE, lw=lw(1.2))
    cv.text(X0 - 10, 508, "Simulate instead", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(X0 - 10, 482, "approximate Bayesian computation (ABC)", ha="left",
            va="baseline", fontsize=fs(FS_TINY), color=NAVY, fontweight="bold")
    steps = [(404, "draw parameters", "from their priors"),
             (286, "simulate the data", "with those parameters"),
             (168, "keep the draw", "if its data are close to ours")]
    for k, (y, head, line) in enumerate(steps):
        cv.text(X0 + 118, y + 14, head, ha="left", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold")
        cv.text(X0 + 118, y - 14, line, ha="left", va="center", fontsize=fs(FS_SMALL),
                color=MUTED)
        if k < 2:
            arrow(pen, (X0 + 46, y - 44), (X0 + 46, y - 74), px=2.0, color=NAVY, head=9)
    rng = np.random.default_rng(8)
    ax = px_axes(fig, X0, 374, 92, 62)                      # priors
    xx = np.linspace(0, 3, 100)
    ax.plot([0, 0, 3, 3], [0, 0.62, 0.62, 0], color=INDIGO, lw=lw(2.4))
    ax.plot([0, 0, 3, 3], [0, 0.5, 0.5, 0], color=OCHRE, lw=lw(2.4), ls=(0, (3.2, 1.6)))
    for xv, col in ((0.9, INDIGO), (2.1, OCHRE)):
        ax.plot([xv], [0.18], marker="v", ms=8, mfc=col, mec=PAPER, mew=lw(0.8))
    ax.set_xlim(-0.2, 3.2)
    ax.set_ylim(0, 0.9)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    ax = px_axes(fig, X0, 254, 92, 66)                      # one simulated data set
    sim = rng.multivariate_normal([0, 0], [[1, 0.72], [0.72, 1]], 13)
    ax.plot(sim[:, 0], sim[:, 1], ls="none", marker="o", ms=4.6, mfc=MUTED, mec=PAPER,
            mew=lw(0.6))
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax = px_axes(fig, X0 + 8, 130, 76, 76)                  # close enough?
    pts = np.clip(rng.normal(0, 1.25, (26, 2)), -2.7, 2.7)
    rr = np.hypot(pts[:, 0], pts[:, 1])
    ax.add_patch(Circle((0, 0), 1.15, facecolor=SOFT, edgecolor=INK, linewidth=lw(1.3),
                        linestyle=(0, (3, 2))))
    ax.plot(pts[rr > 1.15, 0], pts[rr > 1.15, 1], ls="none", marker="o", ms=4.6,
            mfc=PAPER, mec=FAINT, mew=lw(1.1))
    ax.plot(pts[rr <= 1.15, 0], pts[rr <= 1.15, 1], ls="none", marker="o", ms=4.6,
            mfc=INK, mec=PAPER, mew=lw(0.6))
    ax.plot([0], [0], marker="*", ms=13, mfc=PAPER, mec=INK, mew=lw(1.4))
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_aspect("equal")
    blank(ax)
    cv.text(X0 - 10, 100, "close: five summaries within tolerance", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(X0 - 10, 78, "(two means, two spreads, one correlation)", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)
    arrow(pen, (X0 + 356, 200), (X0 + 386, 200), px=2.2, color=NAVY, head=10)

    # ---- the published result: camellia and its weevil ------------------------------------
    # Bar heights read from Fig. 5 of Nuismer & Week (2019); modes and 95 % intervals
    # from its text.  A redrawing of the published posteriors, not a new analysis.
    weevil = [0.00, 0.005, 0.03, 0.08, 0.12, 0.16, 0.26, 0.28, 0.30, 0.33, 0.39, 0.36, 0.35,
              0.41, 0.37, 0.41, 0.43, 0.41, 0.41, 0.43, 0.40, 0.43, 0.42, 0.46, 0.42, 0.45,
              0.47, 0.45, 0.47, 0.45]
    camellia = [0.61, 0.69, 0.74, 0.65, 0.64, 0.66, 0.54, 0.55, 0.48, 0.46, 0.43, 0.38, 0.36,
                0.28, 0.29, 0.28, 0.21, 0.22, 0.18, 0.18, 0.17, 0.15, 0.13, 0.14, 0.11, 0.10,
                0.09, 0.09, 0.08, 0.07]
    X1 = 838
    cv.plot([X1 - 34, X1 - 34], [60, 520], color=LINE, lw=lw(1.2))
    cv.text(X1 - 10, 508, "A camellia and its weevil", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    # the interaction, drawn after Toju (2008, Fig. 1) and Iseki et al. (2011, Fig. 1):
    # a camellia fruit in section (thick pericarp around the seeds) and a weevil
    # boring through it with its long rostrum
    fx, fy, ro, ri = X1 + 50, 426, 56, 27
    cv.add_patch(Circle((fx, fy), ro, facecolor=OCHRE_T, edgecolor=OCHRE, linewidth=lw(2.4),
                        zorder=2))
    cv.add_patch(Circle((fx, fy), ri, facecolor=PAPER, edgecolor=OCHRE, linewidth=lw(1.4),
                        zorder=3))
    for dx, dy, ang in ((-9, 5, -30), (9, 4, 25), (0, -10, 0)):
        cv.add_patch(Ellipse((fx + dx, fy + dy), 16, 12, angle=ang, facecolor="#4A3A2A",
                             edgecolor="none", zorder=4))
    wx, wy = fx + 104, fy + 46                                 # the weevil, standing on the fruit
    for lx in (-14, 2, 18):                                    # legs, down to the fruit
        cv.plot([wx + lx, wx + lx - 8, wx + lx - 16], [wy - 6, wy - 22, wy - 30], color=INDIGO,
                lw=lw(1.8), solid_capstyle="round", zorder=5)
    cv.add_patch(Ellipse((wx + 4, wy), 58, 32, angle=-14, facecolor=INDIGO, edgecolor="none",
                         zorder=6))
    cv.add_patch(Ellipse((wx - 20, wy + 6), 22, 20, angle=-14, facecolor="#3D2FA8",
                         edgecolor="none", zorder=7))
    hd = np.array([wx - 36, wy + 9])
    cv.add_patch(Circle(hd, 9, facecolor="#3D2FA8", edgecolor="none", zorder=7))
    tt = np.linspace(0, 1, 50)                                 # the rostrum: a long arc into the pericarp
    p0, pc, p2_ = hd + (-6, 2), np.array([fx + 46, fy + 62]), np.array([fx + 19, fy + 18])
    arc_ = ((1 - tt) ** 2)[:, None] * p0 + (2 * tt * (1 - tt))[:, None] * pc + (tt ** 2)[:, None] * p2_
    cv.plot(arc_[:, 0], arc_[:, 1], color="#3D2FA8", lw=lw(3.0), solid_capstyle="round", zorder=8)
    k_ = 16
    cv.plot([arc_[k_, 0], arc_[k_, 0] - 4, arc_[k_, 0] - 16], [arc_[k_, 1], arc_[k_, 1] + 14, arc_[k_, 1] + 18],
            color="#3D2FA8", lw=lw(1.3), zorder=7)              # antenna
    cv.annotate("", xy=(fx - ri, fy), xytext=(fx - ro, fy),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=lw(1.5), shrinkA=0, shrinkB=0),
                zorder=9)
    cv.text(fx - ro + 2, fy - ro - 12, "pericarp thickness", ha="left", va="center",
            fontsize=fs(FS_TINY), color=OCHRE, fontweight="bold")
    cv.text(wx + 40, wy + 6, "rostrum length", ha="left", va="center", fontsize=fs(FS_TINY),
            color=INDIGO, fontweight="bold")
    edges = np.linspace(0, 3, 31)
    for k, (h, col, hatch, who, mode, ci, say, sx, sha) in enumerate([
            (weevil, INDIGO, "////", "selection on the weevil", 2.37, (0.60, 2.94),
             "above zero", 0.06, "left"),
            (camellia, OCHRE, "....", "selection on the camellia", 0.21, (0.0, 2.40),
             "piles up at zero", 2.94, "right")]):
        y0 = 240 - k * 126
        ax = px_axes(fig, X1 + 24, y0, 290, 60)
        ax.bar(edges[:-1], h, width=0.1, align="edge", facecolor=PAPER, edgecolor=col,
               hatch=hatch, linewidth=lw(1.2), zorder=2)
        ax.plot(ci, [-0.10, -0.10], color=INK, lw=lw(3.0), solid_capstyle="butt",
                clip_on=False, zorder=5)
        ax.plot([mode], [-0.10], marker="o", ms=8, mfc=PAPER, mec=INK, mew=lw(1.8),
                clip_on=False, zorder=6)
        ax.text(sx, 0.74, say, ha=sha, va="center", fontsize=fs(FS_TINY), color=INK)
        ax.set_xlim(0, 3)
        ax.set_ylim(0, 0.86)
        ax.set_xticks([0, 1, 2, 3])
        ax.tick_params(axis="x", pad=9)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        cv.text(X1 + 24, y0 + 74, who, ha="left", va="center", fontsize=fs(FS_TINY),
                color=col, fontweight="bold")
    cv.text(X1 - 10, 62, "Support for coevolution; but one-sided", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold")
    cv.text(X1 - 10, 40, "evolution cannot be ruled out", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold")
    status(cv, "published",
           "Nuismer & Week 2019, PLoS Comput. Biol.; posteriors (mode, 95 % interval) redrawn from their Fig. 5; drawing after Toju 2008")
    save(fig, "s20_abc")


if __name__ == "__main__":
    main()
