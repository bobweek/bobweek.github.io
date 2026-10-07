"""
Slide 07 -- different interfaces imply different selection.

Week & Nuismer (2021) Am. Nat. 198:195-205, Fig. 1, from Eqs 1-2:
  trait differences   W_X ~ exp( B_X (x - y)),        W_Y ~ exp( B_Y (y - x))
                      B_X = 1.5, B_Y = 1
  offset matching     W_X ~ exp(-B_X/2 (y + d - x)^2), W_Y ~ exp(-B_Y/2 (x + d - y)^2)
                      B_X = 15, B_Y = 5, d = 0.5
Curves are scaled to a maximum of one in the window, as in the paper.
"""
import numpy as np

from semstyle import *  # noqa: F403

DASH_Y = (0, (3.2, 1.6))


def panel(fig, cv, x0, title, sub, kind):
    w, y0, h = 440, 110, 300
    cv.text(x0, 492, title, ha="left", va="baseline", fontsize=fs(FS_TITLE),
            fontweight="bold")
    cv.text(x0, 464, sub, ha="left", va="baseline", fontsize=fs(FS_SMALL),
            color=MUTED)
    ax = px_axes(fig, x0, y0, w, h)
    d = np.linspace(-1, 1, 400)
    if kind == "diff":
        wx, wy = np.exp(1.5 * d), np.exp(-1.0 * d)
    else:
        wx, wy = np.exp(-7.5 * (d - 0.5) ** 2), np.exp(-2.5 * (d + 0.5) ** 2)
    wx, wy = wx / wx.max(), wy / wy.max()
    ax.axvline(0, color=LINE, lw=lw(1.4), zorder=0)
    ax.plot(d, wx, color=INDIGO, lw=lw(LW_DATA + 0.4))
    ax.plot(d, wy, color=OCHRE, lw=lw(LW_DATA + 0.4), ls=DASH_Y)
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1.12)
    ax.set_xticks([-1, 0, 1])
    ax.set_xticklabels(["Y ahead", "equal", "X ahead"])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "max"])
    ax.set_xlabel(r"trait difference $x - y$")
    return ax


def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    a1 = panel(fig, cv, 112, "Trait differences", "a bigger lead is always better", "diff")
    a1.set_ylabel("fitness")
    a1.text(0.66, 0.92, "X", ha="right", va="center", color=INDIGO, fontweight="bold",
            fontsize=fs(FS_LABEL))
    a1.text(-0.66, 0.80, "Y", ha="left", va="center", color=OCHRE, fontweight="bold",
            fontsize=fs(FS_LABEL))
    a1.text(0.0, 0.66, "no best lead:\nselection never relaxes", ha="center",
            va="center", fontsize=fs(FS_SMALL), color=MUTED, linespacing=1.2,
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))

    a2 = panel(fig, cv, 690, "Offset matching", "fitness peaks at a particular lead", "off")
    for xo, col in ((0.5, INDIGO), (-0.5, OCHRE)):
        a2.plot([xo, xo], [0, 1.0], color=col, lw=lw(1.4), ls=(0, (1, 2.2)), zorder=1)
    a2.text(0.5, 1.035, r"X is best ahead by $\delta$", ha="center", va="bottom",
            color=INDIGO, fontsize=fs(FS_SMALL))
    a2.text(-0.5, 1.035, r"Y is best ahead by $\delta$", ha="center", va="bottom",
            color=OCHRE, fontsize=fs(FS_SMALL))
    a2.text(0.74, 0.86, "X", ha="left", va="center", color=INDIGO, fontweight="bold",
            fontsize=fs(FS_LABEL))
    a2.text(-0.93, 0.86, "Y", ha="left", va="center", color=OCHRE, fontweight="bold",
            fontsize=fs(FS_LABEL))
    status(cv, "published",
           "Week & Nuismer 2021, Am. Nat. 198:195 \u2014 Fig. 1 (Eqs 1\u20132, their parameters)")
    save(fig, "s07_interfaces")


if __name__ == "__main__":
    main()
