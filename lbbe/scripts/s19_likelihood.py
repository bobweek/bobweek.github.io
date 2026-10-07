"""
Slide 19 -- from a correlation to a mechanistic likelihood.

Week & Nuismer (2019) Ecology Letters 22:717-725, recomputed.  Data and
background parameters are those of estimate.R in the measuring.coevolution
repository; the formulas are scripts/coev_likelihood.py (functions.R,
transliterated).

s19a_pattern   the pattern: the map of Pauw et al. (2009, Fig. 2) and the
               eight site means taken from it, with their correlation.  Two
               lines put the classical reading beside the point of Nuismer,
               Gomulkiewicz & Ridenhour (2010, Am. Nat. 175:525): coevolution
               is neither necessary nor sufficient for spatially correlated
               traits (one-sided selection and correlated abiotic optima
               produce them too, and reciprocal selection often does not).
s19b_model     the model and how it gives a likelihood, without the algebra:
               pairwise coevolution (offset matching, shown as an inset) at
               many locations, each also under abiotic selection and drift;
               the model therefore predicts the spatial pattern of mean-trait
               pairs for given parameters, P(data | theta), which read as a
               function of the parameters is the likelihood L(theta | data);
               its maximum is the estimate and its curvature the uncertainty;
               coevolution is supported when both one-sided models (B_1 = 0,
               B_2 = 0) are rejected by likelihood-ratio tests (their Fig. 1).
               The scatter, the two candidate distributions and the
               likelihood curve in the middle block are schematic.
s19c_result    the fly and flower: left, the two strengths of biotic selection
               with 95 % bootstrap intervals and p-values; right, the effect
               size, the fitted distribution against the one the model predicts
               without biotic selection (their Fig. 4), with the fly and flower
               drawn to the two pairs of mean lengths.  The build asserts the
               published Table 2 values.
s19x_both_systems  the earlier three-panel figure with the camellia weevil as
               well; kept as a backup.
"""
import matplotlib.image as mpimg
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle
from scipy import stats

from glyphs import Pen, arrow, flower, fly, landscape, on_tile, species

import coev_likelihood as cl
from semstyle import *  # noqa: F403


def _fly_flower():
    snout = np.array([79.3, 54.8, 57.0, 72.8, 42.9, 58.4, 85.8, 75.6])
    tubes = np.array([77.0, (50.0 + 27.5) / 2, 49.1, 52.6, 41.1, 58.6, 60.8, 68.0])
    G1 = np.mean(np.array([4.6, 4.4, 5.9, 5.8, 1.6, 6.8, 6.5, 7.4]) ** 2)
    G2 = np.mean(np.array([5.9, 4.9, 4.3, 4.7, 2.5, 7.5, 4.6, 5.5]) ** 2)
    th1 = np.mean([11.5, 32, (38 + 44) / 2])      # sister species of the fly
    th2 = np.mean([31.2, 51.7])                   # L. anceps without the fly
    return dict(d=np.c_[snout, tubes], n1=100.0, n2=1000.0, th1=th1, th2=th2,
                G1=G1, G2=G2)


def _weevil_camellia():
    fs_ = np.array([9.63, 7.98, 9.13, 10.42, 10.79, 9.89, 10.31, 10.05, 11.91, 10.66,
                    9.12, 9.21, 9.61, 13.63, 10.06, 12.98, 11.68, 11.48, 14.54, 16.92,
                    12.99, 19.48, 21.11, 18.28, 20.59])
    pt = np.array([5.42, 4.32, 3.88, 6.07, 4.90, 6.30, 6.13, 6.66, 7.95, 7.73, 6.76,
                   6.42, 7.52, 11.65, 7.77, 12.80, 11.89, 11.13, 12.49, 17.83, 12.97,
                   20.41, 19.69, 21.21, 19.35])
    sv = np.mean(np.array([0.72, 0.61, 1.02, 0.83, 0.87, 1.17, 0.74, 0.69, 0.59, 0.91,
                           0.66, 0.93, 0.94, 0.70, 1.87, 1.114, 1.87, 1.85, 1.55, 1.61,
                           0.46]) ** 2)
    pv = np.mean(np.array([0.85, 1.19, 0.59, 0.87, 0.71, 1.65, 1.20, 1.38, 2.09, 2.47,
                           1.14, 1.08, 7.52, 2.82, 1.69, 2.45, 2.61, 2.68, 3.39, 1.8,
                           3.6, 3.99, 3.39, 2.38, 2.68]) ** 2)
    th1 = np.mean([6.19, 5.57, 5.5, 5.45, 5.52, 5.71, 6.01, 5.95, 6.33, 6.42, 6.38])
    th2 = np.mean([5.22, 5.27, 5.48, 7.54, 4.68, 5.19, 4.91, 4.99, 5.55, 7.65, 6.82,
                   6.42, 7.36, 7.69, 7.28, 7.56, 6.61])
    n1 = 1 / np.mean(1 / np.array([26389, 41944, 19167, 66389, 44167, 23333, 30278]))
    n2 = 1 / np.mean(1 / np.array([1178, 1334, 2155, 1841, 3112, 2389, 1770]))
    return dict(d=np.c_[fs_, pt], n1=n1, n2=n2, th1=th1, th2=th2,
                G1=sv * np.mean([0.496, 0.353, 0.206, 0.595]),
                G2=pv * np.mean([0.8, 0.82, 0.7, 0.63]))


def fit(sysd, n_boot=1000, seed=1):
    d = sysd["d"]
    N = len(d)
    bg = (sysd["n1"], sysd["n2"], sysd["th1"], sysd["th2"], sysd["G1"], sysd["G2"])
    m, S = d.mean(0), np.cov(d.T, bias=True)
    ML = cl.ml_sol(m[0], m[1], S[0, 0], S[1, 1], S[0, 1], *bg)
    # the model reproduces the sample moments exactly
    assert np.allclose(cl.mu(ML["A1"], ML["A2"], ML["B1"], ML["B2"], ML["k"], *bg[:4]), m)
    assert np.allclose(cl.SIGMA(ML["A1"], ML["A2"], ML["B1"], ML["B2"], ML["k"], *bg), S)
    # the same populations without coevolution (estimate.R, "effect size")
    A1n, A2n = 1 / (2 * S[0, 0] * bg[0]), 1 / (2 * S[1, 1] * bg[1])
    m0 = cl.mu(A1n, A2n, 0.0, 0.0, ML["k"], *bg[:4])
    S0 = cl.SIGMA(A1n, A2n, 0.0, 0.0, ML["k"], *bg)
    rng = np.random.default_rng(seed)
    B = np.full((n_boot, 2), np.nan)
    with np.errstate(all="ignore"):
        for i in range(n_boot):
            f = rng.multivariate_normal(m, S, N)
            fm, fS = f.mean(0), np.cov(f.T, bias=True)
            b = cl.ml_sol(fm[0], fm[1], fS[0, 0], fS[1, 1], fS[0, 1], *bg)
            B[i] = b["B1"], b["B2"]
    ci = np.nanquantile(B, [0.025, 0.975], axis=0)
    # likelihood-ratio tests against the two one-sided models (estimate.R):
    # with B_i = 0 the mean of species i sits at its abiotic optimum
    ll = lambda mean: stats.multivariate_normal(mean, S).logpdf(d).sum()
    lnL = ll(m)
    p = [stats.chi2.sf(2 * (lnL - ll(np.array([bg[2], m[1]]))), 1),
         stats.chi2.sf(2 * (lnL - ll(np.array([m[0], bg[3]]))), 1)]
    b1, b2 = abs(ML["B1"]), abs(ML["B2"])
    w1, w2 = b1 / (b1 + b2), b2 / (b1 + b2)
    return dict(m=m, S=S, ML=ML, m0=m0, S0=S0, ci=ci, p=p, C=np.sqrt(b1 * b2),
                balance=(w1 * np.log(w1) + w2 * np.log(w2)) / np.log(0.5),
                theta=np.array([bg[2], bg[3]]))


def cov_ellipse(ax, m, S, nsd, **kw):
    vals, vecs = np.linalg.eigh(S)
    ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    ax.add_patch(Ellipse(m, 2 * nsd * np.sqrt(vals[1]), 2 * nsd * np.sqrt(vals[0]),
                         angle=ang, **kw))


def scatter_panel(fig, cv, x0, title, sub, sysd, F, xlab, ylab, pad):
    cv.text(x0, 496, title, ha="left", va="baseline", fontsize=fs(FS_LABEL),
            fontweight="bold")
    cv.text(x0, 471, sub, ha="left", va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    ax = px_axes(fig, x0, 112, 300, 330)
    for nsd in (2, 1):
        cov_ellipse(ax, F["m0"], F["S0"], nsd, facecolor="none", edgecolor=MUTED,
                    linewidth=lw(1.8), linestyle=(0, (3.2, 2.0)), zorder=2)
        cov_ellipse(ax, F["m"], F["S"], nsd, facecolor=INDIGO_T if nsd == 2 else "none",
                    edgecolor=INDIGO, linewidth=lw(2.2), zorder=3, alpha=0.9)
    d = sysd["d"]
    ax.plot(d[:, 0], d[:, 1], ls="none", marker="o", ms=7, mfc=INK, mec=PAPER,
            mew=lw(1.0), zorder=5)
    ax.plot(*F["m0"], marker="+", ms=11, mec=MUTED, mew=lw(2.0), zorder=4)
    sd0, sd1 = np.sqrt(np.diag(F["S0"])), np.sqrt(np.diag(F["S"]))
    lo = np.minimum(F["m0"] - 2.3 * sd0, F["m"] - 2.3 * sd1)
    hi = np.maximum(F["m0"] + 2.3 * sd0, F["m"] + 2.3 * sd1)
    ax.set_xlim(lo[0] - pad, hi[0] + pad)
    ax.set_ylim(lo[1] - pad, hi[1] + pad)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax


MAP = HERE.parent / "assets" / "pauw2009_fig2_map.jpg"


def pattern():
    """The pattern, and why a correlation alone does not settle it."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    ff = _fly_flower()
    if MAP.exists():
        img = mpimg.imread(MAP)
        h_px = 322
        w_px = h_px * img.shape[1] / img.shape[0]
        kw = dict(cmap="gray", vmin=0, vmax=255) if img.ndim == 2 else {}
        cv.imshow(img, extent=[40, 40 + w_px, 186, 186 + h_px], zorder=3,
                  aspect="auto", interpolation="lanczos", **kw)
        cv.add_patch(Rectangle((40, 186), w_px, h_px, facecolor="none", edgecolor=INK,
                               linewidth=lw(1.2), zorder=4))
        cv.set_xlim(0, fig._px[0])
        cv.set_ylim(0, fig._px[1])
    cv.text(40, 164, "bars: mean tube depth of the flower (open) and mean proboscis",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(40, 144, "length of the fly (filled) at each site", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    arrow(pen, (486, 348), (590, 348), px=2.6, color=NAVY, head=12)
    cv.text(538, 378, "one point\nper site", ha="center", va="bottom",
            fontsize=fs(FS_SMALL), color=NAVY, linespacing=1.15)
    d = ff["d"]
    r, pr = stats.pearsonr(d[:, 0], d[:, 1])
    cv.text(704, 508, "Longer flies where tubes are deeper", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    ax = px_axes(fig, 704, 190, 320, 290)
    lo, hi = 35, 92
    ax.plot([lo, hi], [lo, hi], color=FAINT, lw=lw(1.4), ls=(0, (1, 2.2)), zorder=1)
    ax.plot(d[:, 0], d[:, 1], ls="none", marker="o", ms=10, mfc=INK, mec=PAPER,
            mew=lw(1.2), zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(22, hi)
    ax.set_xticks([40, 60, 80])
    ax.set_yticks([40, 60, 80])
    ax.set_xlabel("mean proboscis length (mm)")
    ax.set_ylabel("mean tube depth (mm)")
    ax.text(88, 84.5, "equal", ha="right", va="top", fontsize=fs(FS_TINY), color=MUTED,
            rotation=38, rotation_mode="anchor")
    ax.text(0.96, 0.07, f"r = {r:.2f} across {len(d)} sites", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=fs(FS_SMALL), color=INK)
    cv.plot([40, 1140], [112, 112], color=LINE, lw=lw(1.2))
    cv.text(40, 86, "Classically: a correlation like this is evidence of coevolution.",
            ha="left", va="center", fontsize=fs(FS_SMALL), color=MUTED)
    cv.text(40, 56, "But one-sided selection or covarying optima give the same pattern, and coevolution need not give it.",
            ha="left", va="center", fontsize=fs(FS_SMALL), fontweight="bold", color=INK)
    status(cv, "published",
           "map: Pauw et al. 2009, Evolution 63:268, Fig. 2; site means as in Week & Nuismer 2019; Nuismer, Gomulkiewicz & Ridenhour 2010, Am. Nat. 175:525")
    save(fig, "s19a_pattern")


def _mini_pair(pen, cx, cy, u, on=(True, True)):
    """Two species; an arrow from i to j means selection on j caused by i."""
    p, q = (cx - 34 * u, cy), (cx + 34 * u, cy)
    if on[1]:                                        # species 2 responds
        arrow(pen, p, q, px=2.0 * u, color=INDIGO, rad=0.5, head=8 * u,
              shrink=(13 * u, 13 * u), z=4)
    if on[0]:                                        # species 1 responds
        arrow(pen, q, p, px=2.0 * u, color=OCHRE, rad=0.5, head=8 * u,
              shrink=(13 * u, 13 * u), z=4)
    species(pen, *p, 9.5 * u, "A")
    species(pen, *q, 9.5 * u, "B")


def model():
    """The model in one population, the distribution it predicts across
    populations, and the test that distribution allows."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(2)

    # ---- A: pairwise coevolution at many locations ---------------------------------
    cv.text(30, 508, "Coevolution, place by place", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    ax = px_axes(fig, 36, 372, 150, 84)                    # offset matching, as on slide 7
    dd = np.linspace(-1.15, 1.15, 200)
    ax.plot(dd, np.exp(-7.5 * (dd - 0.5) ** 2), color=INDIGO, lw=lw(2.6))
    ax.plot(dd, np.exp(-2.5 * (dd + 0.5) ** 2), color=OCHRE, lw=lw(2.6), ls=(0, (3.2, 1.6)))
    ax.axvline(0, color=LINE, lw=lw(1.2), zorder=0)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(0, 1.12)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    cv.text(36, 356, "offset matching", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    for y, head, line in [(446, "coevolution", r"strengths $B_1$, $B_2$"),
                          (400, "abiotic selection", r"strengths $A_1$, $A_2$"),
                          (354, "drift", "finite populations")]:
        cv.text(214, y, head, ha="left", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold")
        cv.text(214, y - 20, line, ha="left", va="center", fontsize=fs(FS_TINY),
                color=MUTED)
    quad = landscape(pen, 36, 124, 330, 176, skew=0.22, n_lines=3, seed=5, px=1.6)
    sites = [np.array(on_tile(quad, u, v)) for u, v in
             [(0.15, 0.70), (0.40, 0.26), (0.53, 0.76), (0.76, 0.34), (0.90, 0.78)]]
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 4)):          # no gene flow between populations
        p, q = sites[i], sites[j]
        u_ = (q - p) / np.hypot(*(q - p))
        cv.plot(*zip(p + 44 * u_, q - 44 * u_), color=FAINT, lw=lw(1.4), ls=(0, (2.4, 2.4)),
                zorder=3)
        mid = (p + q) / 2
        for sgn in (1, -1):
            cv.plot([mid[0] - 5, mid[0] + 5], [mid[1] - 5 * sgn, mid[1] + 5 * sgn], color=CRIMSON,
                    lw=lw(2.0), zorder=5)
    for k, p in enumerate(sites):
        cv.add_patch(Ellipse(p, 78, 40, facecolor=PAPER, edgecolor=MUTED, linewidth=lw(1.3),
                             zorder=4))
        _mini_pair(pen, p[0], p[1], 0.50)
    cv.text(36, 100, "separate populations, no gene flow:", ha="left", va="center",
            fontsize=fs(FS_TINY), color=CRIMSON, fontweight="bold")
    cv.text(36, 80, "each one an independent run of the same process", ha="left",
            va="center", fontsize=fs(FS_TINY), color=INK)

    # ---- B: the likelihood ---------------------------------------------------------
    X0 = 446
    cv.plot([X0 - 34, X0 - 34], [60, 520], color=LINE, lw=lw(1.2))
    arrow(pen, (X0 - 50, 326), (X0 - 18, 326), px=2.6, color=NAVY, head=11)
    cv.text(X0 - 6, 508, "The likelihood", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(X0 - 6, 482, "predicted pattern against observed pattern", ha="left",
            va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(X0 + 99, 286, "one point per population", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    ax = px_axes(fig, X0 + 4, 318, 190, 148)
    S = np.array([[1.0, 0.72], [0.72, 1.0]])
    xy = rng.multivariate_normal([0, 0], S, 11)
    for nsd in (2, 1):
        cov_ellipse(ax, (0, 0), S, nsd, facecolor=INDIGO_T if nsd == 2 else "none",
                    edgecolor=INDIGO, linewidth=lw(2.0), zorder=2)
    cov_ellipse(ax, (-1.5, 0.5), np.array([[0.8, -0.1], [-0.1, 0.6]]), 2, facecolor="none",
                edgecolor=MUTED, linewidth=lw(1.8), linestyle=(0, (3.2, 2.0)), zorder=1)
    ax.plot(xy[:, 0], xy[:, 1], ls="none", marker="o", ms=6, mfc=INK, mec=PAPER,
            mew=lw(0.9), zorder=3)
    ax.set_xlim(-3.6, 3.0)
    ax.set_ylim(-3.0, 3.0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("mean trait, species 1", fontsize=fs(FS_TINY))
    ax.set_ylabel("species 2", fontsize=fs(FS_TINY))
    ax.text(0.97, 0.05, r"$\hat{\theta}$", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=fs(FS_SMALL), color=INDIGO)
    ax.text(0.05, 0.95, r"$\theta'$", transform=ax.transAxes, ha="left", va="top",
            fontsize=fs(FS_SMALL), color=MUTED)
    xe = X0 + 206
    cv.text(xe, 452, r"$P(\mathrm{data} \mid \theta)$", ha="left", va="center",
            fontsize=fs(FS_LABEL))
    cv.text(xe, 420, r"$=\ \prod_i\ P(\bar{x}_i,\ \bar{y}_i \mid \theta)$", ha="left",
            va="center", fontsize=fs(FS_SMALL))
    cv.text(xe, 392, "a product over populations,", ha="left", va="center",
            fontsize=fs(FS_TINY), color=CRIMSON, fontweight="bold")
    cv.text(xe, 373, "because they are independent", ha="left", va="center",
            fontsize=fs(FS_TINY), color=CRIMSON, fontweight="bold")
    arrow(pen, (xe + 16, 358), (xe + 16, 338), px=1.8, color=MUTED, head=8)
    cv.text(xe + 28, 348, "read as", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    cv.text(xe, 318, r"$L(\theta \mid \mathrm{data})$", ha="left", va="center",
            fontsize=fs(FS_LABEL))
    cv.text(xe + 124, 318, "the likelihood", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(xe, 286, r"$\theta$: all the parameters,", ha="left", va="center",
            fontsize=fs(FS_TINY))
    cv.text(xe, 266, r"$B_1,\ B_2,\ A_1,\ A_2,\ \delta$", ha="left", va="center",
            fontsize=fs(FS_TINY))
    ax = px_axes(fig, X0 + 4, 110, 330, 128)                 # a likelihood curve
    th = np.linspace(-3, 3, 300)
    L = np.exp(-0.5 * (th / 0.8) ** 2)
    ax.fill_between(th, 0, L, color=INDIGO_T, lw=0)
    ax.plot(th, L, color=INDIGO, lw=lw(2.6))
    ax.plot([0, 0], [0, 1], color=INDIGO, lw=lw(1.4), ls=(0, (1, 2.2)))
    ax.plot([0], [1], marker="o", ms=9, mfc=INDIGO, mec=PAPER, mew=lw(1.4), zorder=5)
    hw = 0.8 * 1.2
    ax.plot([-hw, hw], [0.42, 0.42], color=INK, lw=lw(2.0), zorder=4)
    for xx in (-hw, hw):
        ax.plot([xx, xx], [0.36, 0.48], color=INK, lw=lw(2.0), zorder=4)
    ax.plot([-1.9], [np.exp(-0.5 * (1.9 / 0.8) ** 2)], marker="o", ms=8, mfc=PAPER,
            mec=MUTED, mew=lw(1.8), zorder=5)
    ax.text(-2.05, 0.16, r"$\theta'$", ha="right", va="center", fontsize=fs(FS_SMALL),
            color=MUTED)
    ax.text(0.22, 1.08, r"estimate $\hat{\theta}$", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INDIGO)
    ax.text(1.22, 0.80, "uncertainty,\nfrom the curvature", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INK, linespacing=1.15)
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, 1.2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(r"one of the parameters in $\theta$",
                  fontsize=fs(FS_TINY))

    # ---- C: a likelihood and a test -----------------------------------------------
    X1 = 834
    cv.plot([X1 - 34, X1 - 34], [60, 520], color=LINE, lw=lw(1.2))
    arrow(pen, (X1 - 50, 228), (X1 - 18, 228), px=2.6, color=NAVY, head=11)
    cv.text(X1 - 6, 508, "A likelihood, and a test", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    top, left, right = (X1 + 150, 396), (X1 + 56, 236), (X1 + 244, 236)
    for nd, lab, col in [(left, r"$p_1$", MUTED), (right, r"$p_2$", MUTED)]:
        cv.plot([top[0], nd[0]], [top[1] - 34, nd[1] + 34], color=MUTED, lw=lw(1.6),
                zorder=1)
        cv.text((top[0] + nd[0]) / 2 + (-16 if nd is left else 16), (top[1] + nd[1]) / 2,
                lab, ha="center", va="center", fontsize=fs(FS_SMALL))
    for nd, on, hot in [(top, (True, True), True), (left, (False, True), False),
                        (right, (True, False), False)]:
        cv.add_patch(FancyBboxPatch((nd[0] - 62, nd[1] - 30), 124, 60,
                                    boxstyle="round,pad=0,rounding_size=12",
                                    facecolor=INDIGO_T if hot else PAPER,
                                    edgecolor=INDIGO if hot else LINE,
                                    linewidth=lw(2.0 if hot else 1.4), zorder=2))
        _mini_pair(pen, nd[0], nd[1], 0.9, on=on)
    cv.text(top[0], top[1] + 54, r"coevolution:  $B_1,\ B_2 \neq 0$", ha="center",
            va="center", fontsize=fs(FS_SMALL), fontweight="bold")
    cv.text(left[0], left[1] - 50, r"only 2 responds" + "\n" + r"$B_1 = 0$", ha="center",
            va="center", fontsize=fs(FS_TINY), color=MUTED, linespacing=1.2)
    cv.text(right[0], right[1] - 50, r"only 1 responds" + "\n" + r"$B_2 = 0$", ha="center",
            va="center", fontsize=fs(FS_TINY), color=MUTED, linespacing=1.2)
    cv.text(X1 - 6, 122, "likelihood-ratio tests: coevolution is supported",
            ha="left", va="center", fontsize=fs(FS_TINY), color=INK)
    cv.text(X1 - 6, 100, "only if both one-sided models are rejected",
            ha="left", va="center", fontsize=fs(FS_TINY), color=INK)
    status(cv, "published",
           "Week & Nuismer 2019, Ecol. Lett. 22:717 \u2014 discrete populations without gene flow; the tests are their Fig. 1; the middle block is schematic")
    save(fig, "s19b_model")


def result():
    """The fly and the flower.  Left: the two strengths of biotic selection
    with 95 % bootstrap intervals and likelihood-ratio p-values (Table 2 of
    the paper, recomputed).  Right: the effect size, the fitted distribution
    against the one the model predicts without biotic selection (their
    Fig. 4), with a cartoon of the fly and flower: short parts in grey without
    coevolution, long ones in purple with it (exaggerated, not to scale)."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    ff = _fly_flower()
    F = fit(ff, seed=1)
    B = [F["ML"]["B1"], F["ML"]["B2"]]
    gain = (F["m"] - F["m0"]) / F["m0"] * 100
    assert abs(B[0] - 6.40e-5) / 6.40e-5 < 0.02 and abs(B[1] - 1.84e-6) / 1.84e-6 < 0.02
    assert F["p"][0] < 1e-15 and abs(np.log10(F["p"][1]) - np.log10(1.19e-7)) < 0.05
    assert abs(gain[0] - 134) < 1.5 and abs(gain[1] - 34.5) < 1.0
    print(f"  [check] FF: p1 = {F['p'][0]:.1e}, p2 = {F['p'][1]:.2e}, gains = {gain[0]:.0f} %, {gain[1]:.1f} %")

    # ---- what was inferred ----------------------------------------------------------
    cv.text(40, 496, "Biotic selection on each partner", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(40, 471, "strength, 95 % interval, and test against none", ha="left",
            va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    bx = px_axes(fig, 150, 232, 300, 170)
    for i, (v, col, mk) in enumerate(zip(B, [INDIGO, OCHRE], ["o", "s"])):
        y = 1 - i
        lo_, hi_ = max(F["ci"][0, i], 1e-8), F["ci"][1, i]
        bx.plot([lo_, hi_], [y, y], color=col, lw=lw(3.0), zorder=2)
        for xx in (lo_, hi_):
            bx.plot([xx, xx], [y - 0.13, y + 0.13], color=col, lw=lw(3.0), zorder=2)
        bx.plot([v], [y], marker=mk, ms=15, mfc=col, mec=PAPER, mew=lw(1.6), zorder=3)
    bx.set_xscale("log")
    bx.set_xlim(1e-7, 1e-3)
    bx.set_ylim(-0.6, 1.6)
    bx.set_xticks([1e-7, 1e-5, 1e-3])
    bx.set_yticks([])
    bx.spines["left"].set_visible(False)
    bx.set_xlabel(r"strength of biotic selection, $B$ (mm$^{-2}$)", fontsize=fs(FS_SMALL))
    yrow = lambda v: 232 + 170 * (v + 0.6) / 2.2
    fly(pen, 46, yrow(1) + 12, s=16, proboscis=34, color=INDIGO)
    cv.text(100, yrow(1), "fly", ha="left", va="center", fontsize=fs(FS_SMALL),
            color=INDIGO, fontweight="bold")
    flower(pen, 60, yrow(0) + 16, s=18, tube=30)
    cv.text(84, yrow(0), "flower", ha="left", va="center", fontsize=fs(FS_SMALL),
            color=OCHRE, fontweight="bold")
    e2 = int(np.floor(np.log10(F["p"][1])))
    cv.text(468, yrow(1), r"$p < 10^{-15}$", ha="left", va="center", fontsize=fs(FS_SMALL))
    cv.text(468, yrow(0), rf"$p = {F['p'][1] / 10 ** e2:.1f} \times 10^{{{e2}}}$", ha="left",
            va="center", fontsize=fs(FS_SMALL))
    cv.text(40, 112, "Both differ from zero:", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold")
    cv.text(40, 84, "reciprocal selection, not one-sided", ha="left", va="center",
            fontsize=fs(FS_LABEL), fontweight="bold")

    # ---- what it has done: the effect size ---------------------------------------------
    X = 670
    cv.plot([X - 46, X - 46], [50, 520], color=LINE, lw=lw(1.2))
    ax = scatter_panel(fig, cv, X + 16, "Effect size: what coevolution has done",
                       "8 sites; with coevolution (fitted) and without (predicted)",
                       ff, F, "proboscis length (mm)", "tube depth (mm)", 2)
    ax.annotate("", xy=F["m"], xytext=F["m0"],
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=lw(2.0), shrinkA=8,
                                shrinkB=4), zorder=6)
    # the fly and the flower as a cartoon of the elongation: long parts with
    # coevolution, short ones without (exaggerated, not to scale)
    xc = X + 338
    for (pro, tub, col, tint, yb, lab) in [(84, 70, INDIGO, INDIGO_T, 300, "with coevolution"),
                                           (16, 18, MUTED, SOFT, 84, "without")]:
        fly(pen, xc + 6, yb + 118, s=12, proboscis=pro, color=col)
        flower(pen, xc + 84, yb + 96, s=14, tube=tub, color=col, tint=tint)
        cv.text(xc - 6, yb + 152, lab, ha="left", va="center", fontsize=fs(FS_TINY),
                color=col, fontweight="bold")
    cv.text(xc - 6, 304, f"proboscis +{gain[0]:.0f} %", ha="left", va="center",
            fontsize=fs(FS_SMALL), color=INDIGO, fontweight="bold")
    cv.text(xc - 6, 280, f"tube +{gain[1]:.1f} %", ha="left", va="center",
            fontsize=fs(FS_SMALL), color=INDIGO, fontweight="bold")
    status(cv, "published",
           "Week & Nuismer 2019, Ecol. Lett. 22:717 \u2014 their Table 2 and Fig. 4 for the fly and flower, recomputed; the drawings are a cartoon, not to scale")
    save(fig, "s19c_result")


def likelihood():
    fig = figure_px(FULL)
    cv = canvas(fig)
    ff, cw = _fly_flower(), _weevil_camellia()
    Fff, Fcw = fit(ff, seed=1), fit(cw, seed=2)
    est = [Fff["ML"]["B1"], Fff["ML"]["B2"], Fcw["ML"]["B1"], Fcw["ML"]["B2"]]
    for got, want in zip(est, [6.40e-5, 1.84e-6, 7.17e-4, 5.00e-6]):
        assert abs(got - want) / want < 0.02, (got, want)
    print("  [check] B estimates:", ", ".join(f"{v:.3g}" for v in est))

    a1 = scatter_panel(fig, cv, 92, "Fly and flower", "8 sites, Pauw et al. 2009",
                       ff, Fff, "proboscis length (mm)", "tube depth (mm)", 2)
    a1.text(0.97, 0.04, "without\ncoevolution\n(dashed)", transform=a1.transAxes,
            ha="right", va="bottom", fontsize=fs(FS_TINY), color=MUTED,
            linespacing=1.15)
    a1.text(0.04, 0.96, "with coevolution\n(fitted)", transform=a1.transAxes,
            ha="left", va="top", fontsize=fs(FS_TINY), color=INDIGO,
            fontweight="bold", linespacing=1.15)
    scatter_panel(fig, cv, 478, "Weevil and camellia", "25 sites, Toju & Sota",
                  cw, Fcw, "rostrum length (mm)", "pericarp thickness (mm)", 0.5)

    X0 = 884
    cv.text(X0 - 36, 496, "Strength of biotic selection", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(X0 - 36, 471, "on each partner; 95 % intervals",
            ha="left", va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    ax = px_axes(fig, X0, 112, 262, 330)
    ci = np.c_[Fff["ci"], Fcw["ci"]]
    for i, (v, col, mk) in enumerate(zip(est, [INDIGO, OCHRE, INDIGO, OCHRE],
                                         ["o", "s", "o", "s"])):
        lo_, hi_ = max(ci[0, i], 1e-8), ci[1, i]
        ax.plot([i, i], [lo_, hi_], color=col, lw=lw(2.2), zorder=2)
        for yy in (lo_, hi_):
            ax.plot([i - 0.16, i + 0.16], [yy, yy], color=col, lw=lw(2.2), zorder=2)
        ax.plot([i], [v], marker=mk, ms=11, mfc=col, mec=PAPER, mew=lw(1.4), zorder=3)
    ax.axvline(1.5, color=LINE, lw=lw(1.2))
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 1e-2)
    ax.set_xlim(-0.6, 3.6)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["fly", "flower", "weevil", "camellia"], fontsize=fs(FS_TINY))
    ax.set_yticks([1e-7, 1e-5, 1e-3])
    ax.set_ylabel(r"$B$ (mm$^{-2}$)")
    status(cv, "published",
           "Week & Nuismer 2019, Ecol. Lett. 22:717 \u2014 Figs 3\u20134 recomputed from the repository's data and formulas")
    save(fig, "s19x_both_systems")


def main():
    pattern()
    model()
    result()
    likelihood()


if __name__ == "__main__":
    main()
