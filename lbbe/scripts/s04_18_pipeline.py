"""
Slides 04 and 18 -- how I approach these questions, and the way back.

One picture, drawn twice with different emphasis:

  events among individuals (a close-up)
    -> a stochastic process for the whole populations
         -> analytical approximation -> foundational insights
         -> measurable patterns
  and two arrows back to the stochastic process: from the patterns (linking
  data to mechanistic models) and from the insights (approximations that keep
  that inference tractable).  They point at the process, not at individual
  events: what is inferred is the evolutionary process.

s04_approach   forward arrows prominent
s18_backward   the two arrows back prominent, the process panel highlighted

What is computed and what is drawn
----------------------------------
* The population panel is `toy_ibm.simulate`: a toy individual-based
  simulation of a predator and its prey on a plane (twenty-odd thousand
  individuals) whose habitat quality keeps shifting, not a model from any of
  the papers.  Trait matching with weak stabilising selection gives Red Queen
  cycles in the two trait means, and the abundances cycle with them.  Four snapshots show the two populations
  later and later; a dot is an individual, its hue its species, its darkness
  its trait.  The strips below are the two abundances and the two trait means
  +- s.d. from the same run.  Dense on purpose: the approach takes the limit
  of many individuals.
* The left panel is a drawn close-up of a few such individuals.
* "Analytical approximation" shows one generic example of the kind of
  relation such limits yield.  "Foundational insights" shows G of one
  simulated population at three times under drift alone, against the
  deterministic expectation (dashed), from s05_drift_G (Week 2026).
* The "measurable patterns" scatter and histogram are placeholders for
  empirical data; they are not results.
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle

import s05_drift_G as g5
import toy_ibm
from glyphs import Pen, arrow, recip
from semstyle import *  # noqa: F403

# two hues that separate at a glance (the two ends of viridis); darkness = trait
PREY = LinearSegmentedColormap.from_list("prey", ["#F4F1A6", "#C6DC1C", "#7E9C00"])
PRED = LinearSegmentedColormap.from_list("pred", ["#C9BEF4", "#6A45CF", "#32127A"])
LIME, LIME_T, VIOLET, VIOLET_T = "#93B200", "#EDF4B8", "#5B36C4", "#E3DCF8"
P1 = (24, 110, 210, 360)         # x, y, w, h
P2 = (262, 110, 378, 360)
R1A = (672, 306, 224, 164)       # analytical approximation
R1B = (928, 306, 228, 164)       # foundational insights
R2 = (672, 110, 484, 164)        # measurable patterns


def _sim():
    o = toy_ibm.simulate()
    out = {k: o[k] for k in ("t", "NH", "NP", "mH", "mP", "sH", "sP")}
    fr = o["frames"]
    ks = [int((len(fr) - 1) * f) for f in (0.03, 0.34, 0.63, 0.97)]
    for n, k in enumerate(ks):
        t_, hxy, hz, pxy, pz = fr[k]
        out[f"hxy{n}"], out[f"hz{n}"], out[f"pxy{n}"], out[f"pz{n}"] = hxy, hz, pxy, pz
    out["snap_t"] = np.array([fr[k][0] for k in ks])
    return out


def box(cv, P, title, sub="", hot=False):
    x, y, w, h = P
    cv.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=14",
                                facecolor=PAPER, edgecolor=NAVY if hot else LINE,
                                linewidth=lw(2.6 if hot else 1.6), zorder=0))
    lines = title.split("\n")
    for n, line in enumerate(lines):
        cv.text(x + 16, y + h - 32 - 22 * n, line, ha="left", va="baseline",
                fontsize=fs(FS_SMALL), fontweight="bold")
    for n, line in enumerate(sub.split("\n") if sub else []):
        cv.text(x + 16, y + h - 34 - 22 * len(lines) - 19 * n, line, ha="left",
                va="baseline", fontsize=fs(FS_TINY), color=NAVY if hot else MUTED,
                fontweight="bold" if hot else "normal")


def scene(fig, cv, pen, backward=False):
    S = cached("toy_ibm_run3", _sim)
    allh = np.concatenate([S[f"hz{k}"] for k in range(4)])
    allp = np.concatenate([S[f"pz{k}"] for k in range(4)])
    lo_h, hi_h = np.percentile(allh, [2, 98])
    lo_p, hi_p = np.percentile(allp, [2, 98])
    shade_h = lambda z: PREY(np.clip((z - lo_h) / (hi_h - lo_h), 0, 1))
    shade_p = lambda z: PRED(np.clip((z - lo_p) / (hi_p - lo_p), 0, 1))

    # ---- 1  a close-up: events among individuals --------------------------------
    box(cv, P1, "Events among\nindividuals", "births, deaths,\nencounters")
    x0, y0, w, h = P1
    rng = np.random.default_rng(12)
    pts = []
    while len(pts) < 10:
        p = np.array([rng.uniform(x0 + 30, x0 + w - 34), rng.uniform(y0 + 34, y0 + h - 128)])
        if all(np.hypot(*(p - q)) > 46 for q in pts):
            pts.append(p)
    pts = np.array(pts)
    is_pred = np.zeros(len(pts), bool)
    is_pred[[1, 4, 8]] = True
    tz = rng.uniform(0.15, 0.9, len(pts))
    for i, (p, v) in enumerate(zip(pts, tz)):
        dead = i == 3
        fc = PAPER if dead else (PRED(v) if is_pred[i] else PREY(v))
        cv.add_patch(Circle(p, 11, facecolor=fc, edgecolor=FAINT if dead else INK,
                            linewidth=lw(1.2), linestyle=(0, (2, 1.6)) if dead else "-",
                            zorder=5))
        if dead:                                              # a death
            for sgn in (1, -1):
                cv.plot([p[0] - 6, p[0] + 6], [p[1] - 6 * sgn, p[1] + 6 * sgn],
                        color=MUTED, lw=lw(1.6), zorder=6)
    par = pts[6]                                              # a birth
    kid = par + np.array([24, 22])
    cv.add_patch(Circle(kid, 6.5, facecolor=PREY(min(tz[6] + 0.12, 1)), edgecolor=INK,
                        linewidth=lw(1.0), zorder=5))
    arrow(pen, par, kid, px=1.5, color=INK, head=6, shrink=(12, 8), z=6)
    dd = np.hypot(pts[:, None, 0] - pts[None, :, 0], pts[:, None, 1] - pts[None, :, 1])
    dd[~(is_pred[:, None] ^ is_pred[None, :])] = 1e9          # an encounter between species
    dd[[3, 6], :] = 1e9
    dd[:, [3, 6]] = 1e9
    i, j = np.unravel_index(np.argmin(dd), dd.shape)
    recip(pen, pts[i], pts[j], gap=13, rad=0.5, px=1.5, colors=(INK, INK), head=6)

    # ---- 2  the populations as a stochastic process --------------------------------
    box(cv, P2, "A stochastic process for the populations",
        "which process, and how strong?" if backward else "abundances and trait variation",
        hot=backward)
    x0, y0, w, h = P2
    SN, gap = 82, 10
    sx0 = x0 + (w - 4 * SN - 3 * gap) / 2
    sy0 = y0 + 168
    for k in range(4):
        ax = px_axes(fig, sx0 + k * (SN + gap), sy0, SN, SN)
        ax.scatter(S[f"hxy{k}"][:, 0], S[f"hxy{k}"][:, 1], s=0.34, c=shade_h(S[f"hz{k}"]),
                   linewidths=0, zorder=3, rasterized=True)
        ax.scatter(S[f"pxy{k}"][:, 0], S[f"pxy{k}"][:, 1], s=0.34, c=shade_p(S[f"pz{k}"]),
                   linewidths=0, zorder=4, rasterized=True)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        blank(ax)
        for sp_ in ax.spines.values():
            sp_.set_visible(True)
            sp_.set_color(LINE)
            sp_.set_linewidth(lw(1.2))
    arrow(pen, (sx0, sy0 - 11), (sx0 + 4 * SN + 3 * gap, sy0 - 11), px=1.4, color=MUTED,
          head=8)
    cv.text(sx0 + 4 * SN + 3 * gap, sy0 - 25, "time", ha="right", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(sx0, sy0 - 25, "hue: species; darkness: trait", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    t = S["t"]
    wS = 4 * SN + 3 * gap
    ax = px_axes(fig, sx0, y0 + 76, wS, 46)                       # abundances
    ax.plot(t, S["NH"], color=LIME, lw=lw(2.2), ls=(0, (3.2, 1.6)))
    ax.plot(t, S["NP"], color=VIOLET, lw=lw(2.2))
    ax.set_xlim(0, t[-1])
    ax.set_ylim(0, 1.3 * S["NH"].max())
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    ax.text(0.99, 0.97, "abundances", transform=ax.transAxes, ha="right", va="top",
            fontsize=fs(FS_TINY), color=MUTED)
    ax = px_axes(fig, sx0, y0 + 16, wS, 46)                       # traits
    for mean, sd, col, tint, ls in ((S["mH"], S["sH"], LIME, LIME_T, (0, (3.2, 1.6))),
                                    (S["mP"], S["sP"], VIOLET, VIOLET_T, "-")):
        ax.fill_between(t, mean - sd, mean + sd, color=tint, lw=0, alpha=0.8)
        ax.plot(t, mean, color=col, lw=lw(2.0), ls=ls)
    ax.axhline(0, color=FAINT, lw=lw(1.0), zorder=0)
    ax.set_xlim(0, t[-1])
    lo = min((S["mH"] - S["sH"]).min(), (S["mP"] - S["sP"]).min())
    hi = max((S["mH"] + S["sH"]).max(), (S["mP"] + S["sP"]).max())
    ax.set_ylim(lo - 0.05, hi + 0.38 * (hi - lo))
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    ax.text(0.99, 0.97, "trait means and spread", transform=ax.transAxes, ha="right",
            va="top", fontsize=fs(FS_TINY), color=MUTED)
    for k in range(4):                                  # when the four snapshots were taken
        xs = sx0 + wS * S["snap_t"][k] / t[-1]
        cv.plot([xs, xs], [y0 + 66, y0 + 72], color=MUTED, lw=lw(1.4), zorder=3)
    # zoom lines: panel 1 is a close-up inside the first snapshot
    zx, zy = sx0 + 0.56 * SN, sy0 + 0.47 * SN
    cv.add_patch(Rectangle((zx - 6, zy - 6), 12, 12, facecolor="none", edgecolor=INK,
                           linewidth=lw(1.2), zorder=6))
    for (px_, py_), (qx, qy) in [((P1[0] + P1[2], P1[1] + P1[3] - 14), (zx - 6, zy + 6)),
                                 ((P1[0] + P1[2], P1[1] + 14), (zx - 6, zy - 6))]:
        cv.plot([px_, qx], [py_, qy], color=FAINT, lw=lw(1.1), ls=(0, (2, 2)), zorder=1)

    # ---- 3a  analytical approximation -> foundational insights ---------------------------
    box(cv, R1A, "Analytical\napproximation")
    x0, y0, w, h = R1A
    cv.text(x0 + w / 2, y0 + 56, r"$D \approx e\,u\,v\,(1-u)(1-v)$", ha="center",
            va="center", fontsize=fs(FS_LABEL))
    cv.text(x0 + w / 2, y0 + 24, "for example", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    box(cv, R1B, "Foundational\ninsights")
    x0, y0, w, h = R1B
    tau, Gs, rho = g5.simulate_G()
    dtau = tau[1] - tau[0]
    k15, k3 = int(round(1.5 / dtau)), int(round(3.0 / dtau))
    hi_ = next(j for j in range(rho.shape[1])
               if rho[k15, j] > -0.1 and -0.75 < rho[k3, j] < -0.3 and rho[-1, j] < -0.9)
    for n_, ts_ in enumerate((0.0, 3.0, 7.0)):
        ea = px_axes(fig, x0 + 22 + n_ * 64, y0 + 46, 56, 56)
        ea.set_xlim(-1, 1)
        ea.set_ylim(-1, 1)
        ea.set_aspect("equal")
        blank(ea)
        size = np.exp(-g5.SHRINK * ts_)
        if ts_ > 0:
            g5.draw_cov(ea, Gs[0, hi_], scale=size, facecolor="none", edgecolor=INK,
                        lw=lw(1.4), ls=(0, (3.0, 1.6)))
        g5.draw_cov(ea, Gs[int(round(ts_ / dtau)), hi_], scale=size, facecolor=CRIMSON_T,
                    edgecolor=CRIMSON, lw=lw(2.0), alpha=0.95)
    cv.text(x0 + 16, y0 + 32, "new consequences of drift", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(x0 + 16, y0 + 14, "for genetic correlations", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)

    # ---- 3b  measurable patterns ---------------------------------------------------------
    box(cv, R2, "Measurable patterns", "where theory meets observation")
    x0, y0, w, h = R2
    rg = np.random.default_rng(4)
    ax = px_axes(fig, x0 + 286, y0 + 18, 176, 92)
    xy = rg.multivariate_normal([0, 0], [[1, 0.72], [0.72, 1]], 24)
    ax.plot(xy[:, 0] * 1.25, xy[:, 1] * 0.62, ls="none", marker="o", ms=5, mfc=INK,
            mec=PAPER, mew=lw(0.8))
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-2.0, 2.0)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax = px_axes(fig, x0 + 40, y0 + 18, 190, 72)
    ax.hist(rg.gamma(3.0, 1.0, 160), bins=13, facecolor=SOFT, edgecolor=MUTED,
            linewidth=lw(1.2))
    blank(ax)
    ax.spines["bottom"].set_visible(True)

    # ---- arrows ---------------------------------------------------------------------------
    col_f = FAINT if backward else NAVY
    arrow(pen, (P1[0] + P1[2] + 4, 180), (P2[0] - 4, 180), px=2.6, color=col_f, head=11)
    xm = P2[0] + P2[2] + 4
    ya, yb = R1A[1] + R1A[3] / 2, R2[1] + R2[3] / 2
    arrow(pen, (xm, ya), (R1A[0] - 4, ya), px=2.6, color=col_f, head=11)
    arrow(pen, (R1A[0] + R1A[2] + 4, ya), (R1B[0] - 4, ya), px=2.6, color=col_f, head=11)
    arrow(pen, (xm, yb), (R2[0] - 4, yb), px=2.6, color=col_f, head=11)
    cv.text((R2[0] + R2[0] + R2[2]) / 2, (R2[1] + R2[3] + R1A[1]) / 2, "worked on in parallel",
            ha="center", va="center", fontsize=fs(FS_TINY), color=MUTED, style="italic")
    # the two ways back, to the process
    col_b, px_b = (INK, 2.8) if backward else (MUTED, 2.0)
    wt = "bold" if backward else "normal"
    x2 = P2[0] + P2[2] / 2
    arrow(pen, (R1A[0] + R1A[2] / 2, R1A[1] + R1A[3] + 6), (x2 + 60, P2[1] + P2[3] + 6),
          px=px_b, color=col_b, rad=0.26, head=11)
    cv.text(648, 538, "approximations that keep inference tractable", ha="center",
            va="center", fontsize=fs(FS_SMALL), fontweight=wt, color=col_b)
    arrow(pen, (R2[0] + R2[2] / 2, R2[1] - 6), (x2 + 60, P2[1] - 6), px=px_b, color=col_b,
          rad=-0.16, head=11)
    cv.text(760, 46, "link data to mechanistic models", ha="center", va="center",
            fontsize=fs(FS_SMALL), fontweight=wt, color=col_b)


def forward():
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    scene(fig, cv, pen)
    status(cv, "schematic", "the approach; the population panel is a toy predator\u2013prey simulation, the ellipses are from Week 2026")
    save(fig, "s04_approach")


def backward():
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    scene(fig, cv, pen, backward=True)
    status(cv, "schematic", "slide 4 again, with the two ways back to the process brought forward")
    save(fig, "s18_backward")


def main():
    forward()
    backward()


if __name__ == "__main__":
    main()
