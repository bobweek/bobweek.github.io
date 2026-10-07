"""
Slide 08 -- Same arms race, different ecological fate.

Computed from Week & Nuismer (2021) Am. Nat. 198:195-205.

Trait differences: the parameters of their Fig. 2 (left)
    d xbar/dt = G_X B_X,  d ybar/dt = G_Y B_Y                       (Eq. 5a)
    d I_Y/dt = B_Y (G_Y B_Y - G_X B_X)                               (Eq. 7b)
    B_X = 5e-3, B_Y = 1.5e-2, G_X = 10, G_Y = 1, I_Y(0) = 3

Offset matching:
    d xbar/dt = G_X B_X (ybar + delta - xbar),
    d ybar/dt = G_Y B_Y (xbar + delta - ybar)                        (Eq. 5b)
    I_Y = e_Y - B_Y/2 [(delta + D)^2 + s2_X + s2_Y],  D = xbar - ybar (Eq. 6b)
    D0 = 0, B_X = 5e-5, B_Y = 1e-5, G_X = G_Y = s2_X = s2_Y = 10

Offset matching uses the statement's selection strengths and variances
(B_X = 5e-5, B_Y = 1e-5, G = sigma^2 = 10) with a smaller offset, delta = 474
instead of 550, and both baselines set so that each partner's benefit starts
where it starts under trait differences (I_X(0) = 1, I_Y(0) = 3).  With
delta = 550 Y's benefit fell to 0.30, close enough to zero to invite the
question whether it is still a mutualism; with delta = 474 it levels off at
1.0 (Eq. 11b), a third of where it started, and X's rises to 6.0.
"""
import numpy as np

from glyphs import *  # noqa: F403
from semstyle import *  # noqa: F403

TIME_LABEL = "time (thousands of years)"
T = np.linspace(0, 10_000, 801)


def trait_differences():
    BX, BY, GX, GY, IX0, IY0 = 5e-3, 1.5e-2, 10.0, 1.0, 1.0, 3.0
    x = GX * BX * T
    y = GY * BY * T
    IX = IX0 + BX * (GX * BX - GY * BY) * T
    IY = IY0 + BY * (GY * BY - GX * BX) * T
    t_cross = IY0 / (BY * (GX * BX - GY * BY))
    assert IY[-1] < 0 < IX.min(), "Y should end up exploited, X should not"
    return dict(x=x, y=y, IX=IX, IY=IY, t_cross=t_cross)


DELTA = 474.0     # the offset; chosen so that Y's benefit levels off at 1 (see docstring)
E_X = 1.0 + 5e-5 / 2 * (DELTA ** 2 + 20.0)   # so that I_X(0) = 1, as on the left
E_Y = 3.0 + 1e-5 / 2 * (DELTA ** 2 + 20.0)   # so that I_Y(0) = 3, as on the left


def offset_matching():
    D0, delta, BX, BY, eY = 0.0, DELTA, 5e-5, 1e-5, E_Y
    GX = GY = 10.0
    s2 = 20.0                                   # s2_X + s2_Y
    a, b = GX * BX, GY * BY
    Dhat = (a - b) / (a + b) * delta            # Eq. 10
    D = Dhat + (D0 - Dhat) * np.exp(-(a + b) * T)          # solves Eq. 9
    intD = Dhat * T + (D0 - Dhat) * (1 - np.exp(-(a + b) * T)) / (a + b)
    x = a * (delta * T - intD)                  # Eq. 5b integrated
    y = b * (delta * T + intD)
    IY = eY - BY / 2 * ((delta + D) ** 2 + s2)
    IYhat = eY - BY / 2 * ((2 * a * delta / (a + b)) ** 2 + s2)   # Eq. 11b
    assert IY.min() > 0, "mutualism should persist"
    assert abs((E_X - BX / 2 * (delta ** 2 + s2)) - 1.0) < 1e-9
    assert abs(IY[-1] - IYhat) < 0.02 and abs(IYhat - 1.0) < 0.02
    assert np.allclose(x - y, D, atol=1e-6)
    IX = None if E_X is None else E_X - BX / 2 * ((delta - D) ** 2 + s2)
    return dict(x=x, y=y, IX=IX, IY=IY, Dhat=Dhat, IYhat=IYhat)


def fitness_icon(fig, x, y, kind):
    """Individual fitness vs. x - y (Eqs 1-2, parameters of their Fig. 1)."""
    ax = px_axes(fig, x, y, 132, 58)
    d = np.linspace(-1, 1, 200)
    if kind == "diff":
        wx, wy = np.exp(1.5 * d), np.exp(-1.0 * d)
    else:
        wx, wy = np.exp(-7.5 * (d - 0.5) ** 2), np.exp(-2.5 * (d + 0.5) ** 2)
    wx, wy = wx / wx.max(), wy / wy.max()
    ax.plot(d, wx, color=INDIGO, lw=lw(2.6))
    ax.plot(d, wy, color=OCHRE, lw=lw(2.6), ls=(0, (3.2, 1.6)))
    ax.axvline(0, color=LINE, lw=lw(1.2), zorder=0)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-0.04, 1.12)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    return ax


def column(fig, cv, x0, title, sub, kind, res):
    w = 440
    cv.text(x0, 526, title, ha="left", va="baseline", fontsize=fs(FS_TITLE),
            fontweight="bold", color=INK)
    cv.text(x0, 500, sub, ha="left", va="baseline", fontsize=fs(FS_SMALL),
            color=MUTED)
    fitness_icon(fig, x0 + w - 132, 494, kind)
    t = T / 1000.0

    # ---- mean traits ---------------------------------------------------
    at = px_axes(fig, x0, 348, w, 122)
    at.fill_between(t, res["x"], res["y"], color=SOFT, lw=0, zorder=0)
    at.plot(t, res["x"], color=INDIGO, lw=lw(LW_DATA))
    at.plot(t, res["y"], color=OCHRE, lw=lw(LW_DATA), ls=(0, (3.2, 1.6)))
    at.set_xlim(0, 10)
    at.set_xticks([0, 5, 10])
    at.set_xticklabels([])
    top = max(res["x"].max(), res["y"].max())
    at.set_ylim(0, top * 1.08)
    step = 500 if top > 700 else (250 if top > 300 else 100)
    at.set_yticks(np.arange(0, top * 1.05, step))

    # ---- interaction effects on both partners ---------------------------
    ab = px_axes(fig, x0, 100, w, 208)
    ab.axhspan(-3.2, 0, color=SOFT, zorder=0, lw=0)
    ab.axhline(0, color=FAINT, lw=lw(1.4), zorder=1)
    ab.plot(t, res["IX"], color=INDIGO, lw=lw(LW_DATA), zorder=3)
    ab.plot(t, res["IY"], color=OCHRE, lw=lw(LW_DATA), ls=(0, (3.2, 1.6)), zorder=3)
    ab.set_xlim(0, 10)
    ab.set_ylim(-2.8, 7.2)
    ab.set_xticks([0, 5, 10])
    ab.set_yticks([-2, 0, 2, 4, 6])
    ab.set_xlabel(TIME_LABEL)
    return at, ab


def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    D, O = trait_differences(), offset_matching()
    assert D["IX"].min() > 0 and O["IX"].min() > 0, "X benefits throughout in both"

    at, ab = column(fig, cv, 112, "Trait differences",
                    "a bigger lead is always better", "diff", D)
    at.set_ylabel("mean trait\n(change)")
    ab.set_ylabel("benefit from\nthe interaction")
    j = 200
    at.text(T[j] / 1000 - 0.15, D["x"][j] + 30, "X", ha="right", va="bottom",
            color=INDIGO, fontweight="bold", fontsize=fs(FS_LABEL))
    at.text(3.55, D["y"][284] + 14, "Y", ha="center", va="bottom",
            color=OCHRE, fontweight="bold", fontsize=fs(FS_LABEL))
    at.text(9.85, 238, "the lead keeps growing", ha="right", va="center",
            fontsize=fs(FS_SMALL), color=MUTED)
    tc = D["t_cross"] / 1000
    ab.plot([tc], [0], marker="o", ms=9, mfc=PAPER, mec=INK, mew=lw(2.2), zorder=5)
    ab.annotate("Y is now exploited:\nX has become a parasite", xy=(tc, 0),
                xytext=(5.2, -1.75), ha="right", va="center",
                fontsize=fs(FS_SMALL), color=INK,
                arrowprops=dict(arrowstyle="-", color=INK, lw=lw(1.4),
                                shrinkA=3, shrinkB=7))
    ab.text(9.8, D["IX"][-1] + 0.45, "X benefits throughout", ha="right", va="bottom",
            fontsize=fs(FS_SMALL), color=INDIGO, fontweight="bold")
    ab.text(0.25, D["IY"][0] + 0.35, "Y", ha="left", va="bottom", color=OCHRE,
            fontweight="bold", fontsize=fs(FS_LABEL))
    ab.text(9.8, 0.2, "mutualism", ha="right", va="bottom", fontsize=fs(FS_TINY),
            color=MUTED, style="italic")
    ab.text(9.8, -0.2, "antagonism", ha="right", va="top", fontsize=fs(FS_TINY),
            color=MUTED, style="italic")

    at2, ab2 = column(fig, cv, 690, "Offset matching",
                      "fitness peaks at a particular lead", "off", O)
    i = 480
    at2.text(8.7, O["x"][696] + 60, "X", ha="right", va="bottom",
             color=INDIGO, fontweight="bold", fontsize=fs(FS_LABEL))
    at2.text(T[i] / 1000 + 0.35, O["y"][i] - 50, "Y", ha="left", va="top",
             color=OCHRE, fontweight="bold", fontsize=fs(FS_LABEL))
    at2.text(0.35, 0.80 * at2.get_ylim()[1], "both escalate; the lead settles",
             ha="left", va="center", fontsize=fs(FS_SMALL), color=MUTED)
    ab2.plot([10], [O["IYhat"]], marker="o", ms=9, mfc=PAPER, mec=INK,
             mew=lw(2.2), zorder=5, clip_on=False)
    ab2.annotate("Y's benefit shrinks, then holds\nwell above zero: still a mutualism",
                 xy=(10, O["IYhat"]), xytext=(9.7, 3.2), ha="right", va="center",
                 fontsize=fs(FS_SMALL), color=INK,
                 arrowprops=dict(arrowstyle="-", color=INK, lw=lw(1.4),
                                 shrinkA=3, shrinkB=7))
    ab2.text(9.8, O["IX"][-1] - 0.5, "X benefits throughout", ha="right", va="top",
             fontsize=fs(FS_SMALL), color=INDIGO, fontweight="bold")
    ab2.text(2.0, 1.7, "Y", ha="left", va="bottom", color=OCHRE, fontweight="bold",
             fontsize=fs(FS_LABEL))
    ab2.text(0.22, -0.2, "antagonism", ha="left", va="top", fontsize=fs(FS_TINY),
             color=MUTED, style="italic")
    ab2.text(0.22, 0.2, "mutualism", ha="left", va="bottom", fontsize=fs(FS_TINY),
             color=MUTED, style="italic")

    status(cv, "published",
           "Week & Nuismer 2021, Am. Nat. 198:195 \u2014 Eqs 5\u20139; left: their Fig. 2 parameters; right: \u03b4 = 474, so Y's benefit levels off at 1")
    save(fig, "s08_arms_race")


if __name__ == "__main__":
    main()
