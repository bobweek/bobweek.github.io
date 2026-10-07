"""
Slides 14-17 -- space.

s14_mosaic        the geographic mosaic in one picture: the same community at
                  several places, with reciprocal selection strong in some
                  (hot spots) and weak in others (cold spots), linked by gene
                  flow.  Schematic.  Ideas: Thompson (2005); Nuismer, Thompson
                  & Gomulkiewicz (2000) Evolution 54:1102 (selection mosaics
                  and gene flow give clines); Gomulkiewicz, Nuismer & Thompson
                  (2003) Am. Nat. 162:S80 (an interaction that is mutualistic
                  in some places and antagonistic in others); Nuismer,
                  Thompson & Gomulkiewicz (2003) J. Evol. Biol. 16 (hosts and
                  parasites, hot spots, local adaptation).
s15_pair_in_space one host and one parasite in continuous space, and what
                  "local adaptation at distance d" means: the parasite's
                  success on local hosts minus its success on hosts from d
                  away (Week & Bradburd 2024, Eq. 7b).  Schematic.
s16_model         the four forces of the model and what it returns.  The two
                  surfaces are illustrations (random fields with the model's
                  covariance form and a five-fold difference in spatial scale,
                  as in the paper's Fig. 1), not a run of the model.  The
                  covariance curves use the form in hpla/fig3.py.
s17_scale         Week & Bradburd (2024) Am. Nat. 203(1).
                  Left, their Fig. 3: parasite local adaptation against
                  distance for nine combinations of the two species' spatial
                  scales of trait autocorrelation (lambda = 1, 10, 20), all
                  else equal (V_H = V_P = 1, B/A = 1/250).  Recomputed here
                  from Eq. 7b, l_P(d) ~ C_HP(0) - C_HP(d), with C_HP the
                  order-0 Hankel transform of the cross-spectrum (S11c) in
                  closed form (`cross_cov`).  Signs and shapes match the
                  published panels; the vertical scale differs from the
                  published one by a constant factor of about 1.6 (a scaling
                  convention not recoverable from the text), so the axis
                  marks zero only.  In this grid one species is ahead at every
                  distance, or neither is.
                  Right, their Fig. 4: two parameter sets with further
                  asymmetries, where the locally adapted species changes with
                  distance.  Read unchanged from data/hpla/ (the repository's
                  computed values; those sets have strong biotic selection,
                  so they are not recomputed from the weak-selection formula).
                  Vertical lines: lambda = sqrt(sigma^2 / (G A)), as in fig4.R.
"""
import numpy as np
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle
from scipy.special import kn, kv

from glyphs import (Pen, arc, arrow, landscape, link, on_tile, species)
from semstyle import *  # noqa: F403

DATA = HERE.parent / "data" / "hpla"
TILE = dict(x=60, y=64, w=900, h=300, skew=0.26)


# =================================================================== s14 ===
WIDE = dict(x=44, y=44, w=980, h=300, skew=0.26)       # the landscape, nearly the whole slide


def mini_net(pen, cx, cy, seed, hot, missing=()):
    rng = np.random.default_rng(seed)
    top = [(cx + dx, cy + 26) for dx in (-34, 0, 34)]
    bot = [(cx + dx, cy - 26) for dx in (-48, -16, 16, 48)]
    A = (rng.uniform(size=(3, 4)) < 0.55).astype(int)
    for i in range(3):
        A[i, rng.integers(0, 4)] = 1
    for j in range(4):
        A[rng.integers(0, 3), j] = 1
    for i, p in enumerate(top):
        for j, q in enumerate(bot):
            if A[i, j] and ("t", i) not in missing and ("b", j) not in missing:
                if hot:
                    link(pen, p, q, px=1.0 + 2.2 * rng.uniform(0.4, 1), color=INK, z=6)
                else:
                    link(pen, p, q, px=1.1, color=FAINT, ls=(0, (2.4, 2)), z=6)
    for i, p in enumerate(top):
        species(pen, *p, 8, "A", "ghost" if ("t", i) in missing else "on", z=7,
                edge_px=1.4)
    for j, q in enumerate(bot):
        species(pen, *q, 7.5, "B", "ghost" if ("b", j) in missing else "on", z=7,
                edge_px=1.4)
    return top, bot


def _callout(cv, x, head, line, col, target, y=520):
    """A two-line label in the band above the landscape, with a leader down to
    the thing it names."""
    cv.text(x, y, head, ha="center", va="center", fontsize=fs(FS_SMALL),
            fontweight="bold", color=col)
    cv.text(x, y - 22, line, ha="center", va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.plot([x, target[0]], [y - 36, target[1]], color=col, lw=lw(1.2), zorder=9)
    cv.plot([target[0]], [target[1]], marker="o", ms=5, mfc=col, mec=PAPER, mew=lw(0.8),
            zorder=10)


def mosaic():
    """The geographic mosaic: the landscape is the slide.  Four labels sit in a
    band above it, each with a leader to one thing it names.  A faint bracket
    marks a distance d between two places, for the next slide."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    quad = landscape(pen, **WIDE, n_lines=4, seed=5, px=1.8)
    sites = [((0.12, 0.64), True, ()), ((0.33, 0.24), False, ()),
             ((0.53, 0.66), True, ()), ((0.72, 0.26), True, (("t", 2), ("b", 0))),
             ((0.90, 0.68), False, (("t", 0),))]
    pts = [on_tile(quad, *uv) for uv, _, _ in sites]
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (0, 2), (2, 4)]:
        arc(pen, pts[i], pts[j], bulge=0.10, px=1.8, color=NAVY, ls=(0, (3, 2.4)), z=3)
    nets = []
    for k, ((uv, hot, missing), p) in enumerate(zip(sites, pts)):
        cv.add_patch(Ellipse(p, 150, 62, facecolor=INDIGO_T if hot else PAPER,
                             edgecolor=INDIGO if hot else FAINT, linewidth=lw(1.4),
                             linestyle="-" if hot else (0, (3, 2)), zorder=2))
        cv.plot([p[0], p[0]], [p[1], p[1] + 44], color=MUTED, lw=lw(1.2), zorder=4)
        cv.add_patch(Circle(p, 4, facecolor=INK, edgecolor="none", zorder=5))
        nets.append(mini_net(pen, p[0], p[1] + 82, seed=20 + k, hot=hot, missing=missing))
    top = lambda k: (pts[k][0], pts[k][1] + 122)
    _callout(cv, 196, "hot spot", "strong reciprocal selection", INDIGO, top(0))
    _callout(cv, 410, "cold spot", "weak reciprocal selection", MUTED, top(1))
    ghost = nets[3][0][2]
    _callout(cv, 640, "partners change", "who is there differs by place", INK,
             (ghost[0], ghost[1] + 12))
    mid = (pts[3] + pts[4]) / 2
    _callout(cv, 880, "gene flow", "connects local histories", NAVY, (mid[0] - 6, mid[1] + 14))
    # a distance between two places, for the next slide
    p, q = pts[1], pts[3]
    yb = min(p[1], q[1]) - 46
    cv.plot([p[0], q[0]], [yb, yb], color=MUTED, lw=lw(1.2), zorder=4)
    for x_, y_ in ((p[0], p[1]), (q[0], q[1])):
        cv.plot([x_, x_], [yb - 5, yb + 5], color=MUTED, lw=lw(1.2), zorder=4)
        cv.plot([x_, x_], [yb + 5, y_ - 33], color=FAINT, lw=lw(1.0), ls=(0, (1.5, 2)),
                zorder=3)
    cv.text((p[0] + q[0]) / 2, yb, r"a distance $d$", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED, zorder=8,
            bbox=dict(boxstyle="square,pad=0.25", fc=SOFT, ec="none"))
    status(cv, "schematic",
           "geographic mosaic theory: Thompson 2005; Nuismer, Thompson & Gomulkiewicz 2000, 2003; Gomulkiewicz, Nuismer & Thompson 2003")
    save(fig, "s14_mosaic")


# =================================================================== s15 ===
def _trait_fields(seed=11, n=180, lH=0.25, lP=0.05):
    """The pair of surfaces used on the three spatial slides (one realisation)."""
    rng = np.random.default_rng(seed)
    zH = field(n, lH, rng)
    zP = 0.55 * zH + 0.84 * field(n, lP, rng)
    return zH, zP


def pair_in_space():
    """One complication at a time.  Many species was one; space is another, so
    take it with a single host and a single parasite.  The landscape carries
    pairs of populations whose mean traits differ from place to place (shade),
    more alike nearby than far apart.  Three double arrows mark reciprocal
    transplants at three distances: would they give the same answer?"""
    from matplotlib.colors import LinearSegmentedColormap
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    quad = landscape(pen, **WIDE, n_lines=4, seed=5, px=1.8)
    rng = np.random.default_rng(2)
    zH, zP = _trait_fields(lH=0.45, lP=0.22)
    n = zH.shape[0]
    cmH = LinearSegmentedColormap.from_list("h", ["#F6E3C2", OCHRE, "#5A3306"])
    cmP = LinearSegmentedColormap.from_list("p", ["#DDD9F8", INDIGO, "#1E1670"])
    shade = lambda cm, v: cm(float(np.clip((v + 2.0) / 4.0, 0, 1)))
    # the three comparisons, at three distances, in three parts of the landscape
    comps = [((0.10, 0.30), (0.19, 0.30), "near"),
             ((0.34, 0.72), (0.60, 0.72), "further"),
             ((0.30, 0.20), (0.92, 0.42), "far")]
    special = [np.array(on_tile(quad, *uv)) for c_ in comps for uv in c_[:2]]
    uvs = [uv for c_ in comps for uv in c_[:2]]
    placed, puv = list(special), list(uvs)
    while len(placed) < 58:
        u, v = rng.uniform(0.04, 0.96), rng.uniform(0.08, 0.92)
        p = np.array(on_tile(quad, u, v))
        if all(np.hypot(*(p - q)) > 46 for q in placed):
            placed.append(p)
            puv.append((u, v))
    for k, (p, (u, v)) in enumerate(zip(placed, puv)):
        i, j = int(u * (n - 1)), int(v * (n - 1))
        big = k < len(special)
        r = 8.5 if big else 6.0
        for dx, z, cm, sq in ((-r - 1.5, zH[i, j], cmH, True), (r + 1.5, zP[i, j], cmP, False)):
            col = shade(cm, z)
            if sq:
                cv.add_patch(Rectangle((p[0] + dx - r, p[1] - r), 2 * r, 2 * r, facecolor=col,
                                       edgecolor=INK if big else "none", linewidth=lw(1.4),
                                       zorder=6 if big else 4))
            else:
                cv.add_patch(Circle((p[0] + dx, p[1]), r, facecolor=col,
                                    edgecolor=INK if big else "none", linewidth=lw(1.4),
                                    zorder=6 if big else 4))
    box = dict(boxstyle="square,pad=0.2", fc=SOFT, ec="none")
    for k, (ua, ub, lab) in enumerate(comps):
        p, q = special[2 * k], special[2 * k + 1]
        rad = (-0.55, -0.30, -0.16)[k]
        arrow(pen, tuple(p), tuple(q), px=2.6, color=INK, rad=rad, head=11, shrink=(22, 22),
              both=True, z=7)
        mid = (p + q) / 2
        lift = abs(rad) * np.hypot(*(q - p)) / 2 + 16
        cv.text(mid[0], mid[1] + lift, lab, ha="center", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", zorder=9, bbox=box)
    # ---- the band above the landscape ---------------------------------------------
    cv.text(44, 528, "Many species was one complication. Space is another.", ha="left",
            va="center", fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(44, 500, "So take space with one host and one parasite, in many places at once.",
            ha="left", va="center", fontsize=fs(FS_SMALL))
    cv.add_patch(Rectangle((48, 461), 15, 15, facecolor=OCHRE, edgecolor="none"))
    cv.text(70, 469, "host population", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    cv.add_patch(Circle((204, 469), 8, facecolor=INDIGO, edgecolor="none"))
    cv.text(218, 469, "parasite population", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    cv.text(372, 469, "darker: larger mean trait", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(1136, 528, "Swap parasites and hosts between two places:", ha="right",
            va="center", fontsize=fs(FS_SMALL))
    cv.text(1136, 500, "does the answer depend on how far apart they are?", ha="right",
            va="center", fontsize=fs(FS_SMALL), fontweight="bold", color=INDIGO)
    status(cv, "schematic", "after Week & Bradburd 2024, Am. Nat. 203; trait values are an illustration")
    save(fig, "s15_pair_in_space")


def _row(cv, fig, x0, y, pieces, size, widths=None):
    """Lay maths pieces out left to right; return the x-range of each.  With
    `widths`, each piece is centred in a column of that width."""
    r = fig.canvas.get_renderer()
    k = fig._px[0] / fig.bbox.width
    out, x = [], x0
    for i, txt in enumerate(pieces):
        t = cv.text(x, y, txt, ha="left", va="center", fontsize=fs(size))
        w = t.get_window_extent(r).width * k
        if widths is not None:
            t.set_x(x + (widths[i] - w) / 2)
            w = widths[i]
        out.append((x, x + w))
        x += w
    return out


def _marks(ax, n, kind, p, filled=True, r=0.034):
    """A host population (square) or a parasite population (circle) drawn on a
    trait map; ghosted when it stands for the other species' place."""
    x, y = p[0] * n, p[1] * n
    kw = dict(facecolor="white" if filled else "none", edgecolor=INK if filled else "white",
              linewidth=lw(2.0), linestyle="-" if filled else (0, (2, 1.4)), zorder=8)
    if kind == "H":
        ax.add_patch(Rectangle((x - r * n, y - r * n), 2 * r * n, 2 * r * n, **kw))
    else:
        ax.add_patch(Circle((x, y), 1.12 * r * n, **kw))


def spde():
    """The model in one slide (Week & Bradburd 2024, Eqs 1a-b, main-text
    parameterisation).  Top: the two equations with each term named.  Bottom:
    what they produce, a pair of random surfaces (a bivariate Gaussian random
    field), each with its own characteristic length, lambda = sigma / sqrt(A)
    (their Eq. 4), the distance over which that species' mean trait turns
    over.  Marked on the surfaces: a pair of places within one species, a
    distance lambda apart (what the autocovariance is about), and a host
    place paired with a parasite place a distance d away (what the
    cross-covariance is about; the other species' place is ghosted on each
    map)."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    fig.canvas.draw()
    rows = [(446, [r"$\frac{\partial \bar{z}_H}{\partial t}\ =$", r"$-\,A_H\,\bar{z}_H$",
                   r"$-\,B_H\,(\bar{z}_P - \bar{z}_H)$",
                   r"$+\ \frac{\sigma_H^2}{2}\,\nabla^2 \bar{z}_H$", r"$+\ D_H\,\xi_H$"], OCHRE,
             "host"),
            (378, [r"$\frac{\partial \bar{z}_P}{\partial t}\ =$", r"$-\,A_P\,\bar{z}_P$",
                   r"$+\,B_P\,(\bar{z}_H - \bar{z}_P)$",
                   r"$+\ \frac{\sigma_P^2}{2}\,\nabla^2 \bar{z}_P$", r"$+\ D_P\,\xi_P$"], INDIGO,
             "parasite")]
    widths = [150, 150, 250, 210, 150]
    x0 = 150
    heads = [None, ("abiotic selection", "toward a local optimum"),
             ("biotic selection", "host evades, parasite follows"),
             ("dispersal", "offspring settle nearby"),
             ("random genetic drift", "space-time white noise")]
    xs = np.cumsum([x0] + widths)
    for i, h in enumerate(heads):
        if h is None:
            continue
        cv.add_patch(Rectangle((xs[i] + 6, 340), widths[i] - 12, 208, facecolor=SOFT,
                               edgecolor="none", zorder=0))
        cv.text(xs[i] + widths[i] / 2, 528, h[0], ha="center", va="center",
                fontsize=fs(FS_SMALL), fontweight="bold")
        cv.text(xs[i] + widths[i] / 2, 504, h[1], ha="center", va="center",
                fontsize=fs(FS_TINY), color=MUTED)
    for y, pieces, col, who in rows:
        _row(cv, fig, x0, y, pieces, 25, widths)
        cv.text(40, y, who, ha="left", va="center", fontsize=fs(FS_LABEL), fontweight="bold",
                color=col)

    # ---- what comes out: two surfaces, each with its own scale -----------------------
    lH, lP = 0.25, 0.10
    zH, zP = _trait_fields(lH=lH, lP=lP)
    n = zH.shape[0]
    S = 232
    home, away = (0.30, 0.34), (0.72, 0.60)                 # a host place, a parasite place
    for x_, z, lab, col, kind, ell, sym, cov in (
            (318, zH, "host mean trait", OCHRE, "H", lH, r"$\lambda_H$", r"$C_H(d)$"),
            (590, zP, "parasite mean trait", INDIGO, "P", lP, r"$\lambda_P$", r"$C_P(d)$")):
        ax = px_axes(fig, x_, 46, S, S)
        ax.imshow(z, cmap="cividis", vmin=-2.6, vmax=2.6, origin="lower",
                  interpolation="bilinear")
        blank(ax)
        cv.text(x_, 294, lab, ha="left", va="center", fontsize=fs(FS_SMALL), color=col,
                fontweight="bold")
        # within the species: two places a distance lambda apart
        p1 = (0.16, 0.84)
        p2 = (p1[0] + ell, p1[1])
        ax.plot([p1[0] * n, p2[0] * n], [p1[1] * n] * 2, color="white", lw=lw(4.2), zorder=6)
        ax.plot([p1[0] * n, p2[0] * n], [p1[1] * n] * 2, color=INK, lw=lw(2.0), zorder=7)
        _marks(ax, n, kind, p1, r=0.028)
        _marks(ax, n, kind, p2, r=0.028)
        ax.text((p2[0] + 0.06) * n, p1[1] * n, cov, ha="left", va="center",
                fontsize=fs(FS_SMALL), color=INK, zorder=9,
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.92))
        ax.text((p1[0] + p2[0]) / 2 * n, (p1[1] - 0.085) * n, sym, ha="center", va="center",
                fontsize=fs(FS_TINY), color=INK, zorder=9,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.92))
        # between the species: a host place and a parasite place, a distance d apart
        ax.plot([home[0] * n, away[0] * n], [home[1] * n, away[1] * n], color="white",
                lw=lw(1.8), ls=(0, (2.4, 1.8)), zorder=6)
        _marks(ax, n, "H", home, filled=kind == "H")
        _marks(ax, n, "P", away, filled=kind == "P")
        ax.text((home[0] + away[0]) / 2 * n + 0.03 * n, (home[1] + away[1]) / 2 * n - 0.075 * n,
                r"$C_{HP}(d)$", ha="left", va="center", fontsize=fs(FS_SMALL), color=INK, zorder=9,
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.92))
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(-0.5, n - 0.5)
    arrow(pen, (150, 322), (150, 296), px=2.6, color=NAVY, head=11)
    cv.text(40, 262, "Its solutions:", ha="left", va="center", fontsize=fs(FS_SMALL))
    cv.text(40, 236, "two random surfaces", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold")
    cv.text(40, 210, "a bivariate Gaussian random field", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(40, 150, "each with its own spatial scale", ha="left", va="center",
            fontsize=fs(FS_SMALL))
    cv.text(40, 112, r"$\lambda_H = \sigma_H \,/\, \sqrt{A_H}$", ha="left", va="center",
            fontsize=fs(FS_LABEL), color=OCHRE)
    cv.text(40, 72, r"$\lambda_P = \sigma_P \,/\, \sqrt{A_P}$", ha="left", va="center",
            fontsize=fs(FS_LABEL), color=INDIGO)
    # what the marked pairs are
    X = 850
    cv.text(X, 236, r"$C_H(d),\ C_P(d)$", ha="left", va="center", fontsize=fs(FS_LABEL))
    cv.text(X, 204, "intraspecific spatial\nautocovariance functions", ha="left", va="center",
            fontsize=fs(FS_SMALL), linespacing=1.2)
    cv.text(X, 142, r"$C_{HP}(d)$", ha="left", va="center", fontsize=fs(FS_LABEL))
    cv.text(X, 110, "interspecific spatial\ncross-covariance function", ha="left", va="center",
            fontsize=fs(FS_SMALL), linespacing=1.2)
    cv.text(X, 58, "dashed outline: the other species", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    status(cv, "published", "Week & Bradburd 2024, Am. Nat. 203(1), Eqs 1a\u2013b and 4; the surfaces are illustrations with the model's covariance form")
    save(fig, "s15b_spde")


def la_index():
    """Local adaptation as a function of distance (Week & Bradburd 2024).
    Definition (their Eq. 2): expected success of a population on the partners
    it lives with, minus on partners a distance d away.  In the model (their
    Eqs 7a-b) it is the host-parasite cross-covariance at the same place minus
    at distance d, times the strength of biotic selection, and the host's is
    the mirror image.  Left: the two surfaces, with a parasite population at
    home and hosts at three distances.  Right: the cross-covariance and the
    local adaptation it gives, from `cross_cov` (one case of their Fig. 3;
    vertical scale as on the next slide)."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    cv.text(40, 522, "Local adaptation at distance $d$", ha="left", va="center",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(40, 494, "how well a population does at home, minus a distance $d$ away",
            ha="left", va="center", fontsize=fs(FS_SMALL))
    cv.text(640, 522, r"$\ell_S(d)$ := $\mathbb{E}\,[\,\bar{m}_S(0) - \bar{m}_S(d)\,], \quad S = H,\,P$",
            ha="left", va="center", fontsize=fs(FS_LABEL))
    cv.text(640, 494, r"$\bar{m}_S(d)$: mean growth rate of species $S$ when transplanted a distance $d$ away",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    span = 64.0                                   # the maps are this many distance units wide
    ds = (7.0, 17.0, 36.0)
    zH, zP = _trait_fields(lH=10.0 / span * 1.6, lP=0.05)
    n = zH.shape[0]
    home = np.array([0.22, 0.46])
    S = 236
    for x_, z, lab, col, kind in ((40, zP, "the parasite, at home", INDIGO, "P"),
                                  (300, zH, "hosts: at home, and further and further away",
                                   OCHRE, "H")):
        ax = px_axes(fig, x_, 190, S, S)
        ax.imshow(z, cmap="cividis", vmin=-2.6, vmax=2.6, origin="lower",
                  interpolation="bilinear")
        blank(ax)
        cv.text(x_, 446, lab, ha="left", va="center", fontsize=fs(FS_TINY), color=col,
                fontweight="bold")
        if kind == "P":
            _marks(ax, n, "P", home, r=0.04)
        else:
            _marks(ax, n, "H", home, r=0.04)
            for k, d_ in enumerate(ds):
                th = np.radians((-38.0, 24.0, 8.0)[k])
                q = home + d_ / span * np.array([np.cos(th), np.sin(th)])
                ax.add_patch(Circle(home * n, d_ / span * n, facecolor="none", edgecolor="white",
                                    linewidth=lw(1.1), linestyle=(0, (2, 2)), zorder=5))
                ax.plot([home[0] * n, q[0] * n], [home[1] * n, q[1] * n], color="white",
                        lw=lw(1.4), zorder=6)
                _marks(ax, n, "H", q, r=0.034)
                ax.text(q[0] * n, (q[1] + 0.075) * n, str(k + 1), ha="center", va="center",
                        fontsize=fs(FS_TINY), color=INK, zorder=9,
                        bbox=dict(boxstyle="circle,pad=0.12", fc="white", ec="none", alpha=0.95))
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(-0.5, n - 0.5)

    # ---- the cross-covariance, and the local adaptation it gives ---------------------------
    X = 640
    d = np.linspace(0, 48, 300)
    C = cross_cov(np.maximum(d, 1e-9), 10.0, 1.0) * 250.0
    C0 = float(cross_cov(1e-9, 10.0, 1.0) * 250.0)
    ell = C0 - C
    assert C0 > 0 and np.all(np.diff(ell) > -1e-9) and abs(ell[-1] - C0) < 0.05 * C0
    ax = px_axes(fig, X + 30, 190, 440, 236)
    ax.axhline(0, color=FAINT, lw=lw(1.2), zorder=1)
    ax.axhline(C0, color=FAINT, lw=lw(1.2), ls=(0, (2, 2)), zorder=1)
    ax.plot(d, C, color=MUTED, lw=lw(LW_DATA), ls=(0, (3.2, 1.6)), zorder=3)
    ax.plot(d, ell, color=INDIGO, lw=lw(LW_DATA + 0.4), zorder=4)
    for k, d_ in enumerate(ds):
        y_ = float(np.interp(d_, d, ell))
        ax.plot([d_, d_], [0, y_], color=FAINT, lw=lw(1.2), zorder=2)
        ax.plot([d_], [y_], marker="o", ms=10, mfc=PAPER, mec=INDIGO, mew=lw(2.2), zorder=5)
        ax.text(d_, -0.11 * C0, str(k + 1), ha="center", va="center", fontsize=fs(FS_TINY))
    ax.text(47, 0.80 * C0, r"local adaptation $\ell_P(d)$", ha="right", va="center",
            fontsize=fs(FS_SMALL), color=INDIGO, fontweight="bold")
    ax.text(47, 0.13 * C0, r"cross-covariance $C_{HP}(d)$", ha="right", va="center",
            fontsize=fs(FS_SMALL), color=MUTED, fontweight="bold")
    ax.set_xlim(0, 48)
    ax.set_ylim(-0.2 * C0, 1.14 * C0)
    ax.set_xticks([])
    ax.set_yticks([0, C0])
    ax.set_yticklabels(["0", r"$C_{HP}(0)$"])
    ax.set_xlabel(r"distance $d$")
    cv.text(X + 30, 446, "in the model", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)

    cv.plot([40, 1140], [150, 150], color=LINE, lw=lw(1.2))
    cv.text(40, 112, r"$\ell_P(d) \;=\; B_P\,[\,C_{HP}(0) - C_{HP}(d)\,]$", ha="left",
            va="center", fontsize=fs(FS_LABEL), color=INDIGO)
    cv.text(40, 76, "in the model: how well host and parasite match in the same place,", ha="left",
            va="center", fontsize=fs(FS_TINY), color=INK)
    cv.text(40, 56, "minus how well they match a distance $d$ apart", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INK)
    cv.text(X + 30, 112, r"$\ell_H(d) \;=\; B_H\,[\,C_{HP}(d) - C_{HP}(0)\,]$", ha="left",
            va="center", fontsize=fs(FS_LABEL), color=OCHRE)
    cv.text(X + 30, 76, "for the host, the mirror image", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INK)
    cv.text(X + 30, 56, r"$B$: strength of biotic selection", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    status(cv, "published", "Week & Bradburd 2024, Am. Nat. 203(1), Eqs 2 and 7; curves from the model's cross-covariance (one case of their Fig. 3); surfaces are illustrations")
    save(fig, "s16_la_index")


# =================================================================== s16 ===
def field(n, ell, rng):
    """A random surface with the model's covariance form: white noise passed
    through (kappa^2 - Laplacian)^-1, kappa = sqrt(2)/ell (box of side 1)."""
    k = 2 * np.pi * np.fft.fftfreq(n, d=1.0 / n)
    kx, ky = np.meshgrid(k, k)
    z = np.real(np.fft.ifft2(np.fft.fft2(rng.normal(size=(n, n)))
                             / (2.0 / ell ** 2 + kx ** 2 + ky ** 2)))
    return (z - z.mean()) / z.std()


def model():
    """Each species' mean trait is a random surface with its own spatial scale
    (lambda: the distance over which mean traits stay alike), set by dispersal
    against selection; reciprocal selection ties the two surfaces together."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    lH, lP = 0.25, 0.05
    zH, zP = _trait_fields(lH=lH, lP=lP)
    for x0, z, title, col, sub, ell, sym in [
            (60, zH, "Host mean trait", OCHRE, "disperses far: varies over a long distance", lH,
             r"$\lambda_H$"),
            (390, zP, "Parasite mean trait", INDIGO, "disperses little: varies over a short distance",
             lP, r"$\lambda_P$")]:
        ax = px_axes(fig, x0, 170, 300, 300)
        ax.imshow(z, cmap="cividis", vmin=-2.6, vmax=2.6, origin="lower",
                  interpolation="bilinear")
        blank(ax)
        cv.text(x0, 508, title, ha="left", va="baseline", fontsize=fs(FS_LABEL),
                fontweight="bold", color=col)
        cv.text(x0, 484, sub, ha="left", va="baseline", fontsize=fs(FS_TINY),
                color=MUTED)
        # the spatial scale, to scale, on the map itself
        n = z.shape[0]
        ax.add_patch(Rectangle((0.03 * n, 0.03 * n), (ell + 0.22) * n, 0.13 * n,
                               facecolor="white", edgecolor="none", alpha=0.92, zorder=5))
        ax.plot([0.06 * n, (0.06 + ell) * n], [0.095 * n] * 2, color=INK, lw=lw(3.4),
                solid_capstyle="butt", zorder=6)
        for xx in (0.06 * n, (0.06 + ell) * n):
            ax.plot([xx, xx], [0.07 * n, 0.12 * n], color=INK, lw=lw(1.8), zorder=6)
        ax.text((0.085 + ell) * n, 0.095 * n, sym, ha="left", va="center",
                fontsize=fs(FS_SMALL), zorder=7)
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(-0.5, n - 0.5)
    # how the model's parameters set the scale and the amount of variation
    for y, head, formula, words in [
            (130, "spatial scale", r"$\lambda = \sigma \,/\, \sqrt{G A}$",
             "dispersal distance, over the root of genetic\nvariance times abiotic selection"),
            (72, "variation among places", r"$V = 1 \,/\, (N A\, \sigma^2)$",
             "more in small populations (drift), less with\nstrong selection and wide dispersal")]:
        cv.text(60, y, head, ha="left", va="center", fontsize=fs(FS_SMALL), fontweight="bold")
        cv.text(262, y, formula, ha="left", va="center", fontsize=fs(FS_LABEL))
        cv.text(424, y, words, ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED,
                linespacing=1.2)

    X = 800
    cv.text(X - 40, 508, "Each species has its own spatial scale", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(X - 40, 484, r"$\lambda$: the distance over which mean traits stay alike",
            ha="left", va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    ax = px_axes(fig, X, 170, 340, 290)
    d = np.linspace(1e-6, 0.6, 300)
    for ell, col, ls in [(lH, OCHRE, (0, (3.2, 1.6))), (lP, INDIGO, "-")]:
        r = np.sqrt(2) * d / ell
        ax.plot(d, r * kn(1, r), color=col, lw=lw(LW_DATA), ls=ls)
    for ell, col, y, sym in [(lH, OCHRE, 0.09, r"$\lambda_H$"), (lP, INDIGO, 0.21, r"$\lambda_P$")]:
        ax.plot([0, ell], [y, y], color=col, lw=lw(3.0), solid_capstyle="butt", zorder=5)
        ax.plot([ell, ell], [y - 0.03, y + 0.03], color=col, lw=lw(1.8), zorder=5)
        ax.text(ell + 0.012, y, sym, ha="left", va="center", fontsize=fs(FS_SMALL),
                color=col, zorder=6,
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.9))
    ax.text(0.35, 0.47, "host", ha="left", va="center", fontsize=fs(FS_SMALL),
            color=OCHRE, fontweight="bold")
    ax.text(0.10, 0.40, "parasite", ha="left", va="center", fontsize=fs(FS_SMALL),
            color=INDIGO, fontweight="bold")
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, 0.3, 0.6])
    ax.set_xticklabels(["0", "", "far"])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "same"])
    ax.set_xlabel("distance between two places")
    ax.set_ylabel("similarity of mean trait", labelpad=-14)
    cv.text(X - 40, 112, "one scale for each species; and between them,\na cross-covariance at every distance,\nmade by biotic selection",
            ha="left", va="top", fontsize=fs(FS_TINY), color=MUTED, linespacing=1.25)
    status(cv, "schematic",
           "model: Week & Bradburd 2024, Am. Nat. 203; the surfaces are illustrations with the model's covariance form, not a model run")
    save(fig, "s16_model")


# ================================================================== s17a ===
def _layer(cv, z, x0, y0, w, h, skew, cmap, vlim=2.6, z_=3):
    """A surface drawn as a tilted sheet: the unit square mapped to a
    parallelogram with lower-left corner (x0, y0)."""
    from matplotlib.transforms import Affine2D
    tr = Affine2D.from_values(w, 0, skew * h, h, x0, y0) + cv.transData
    cv.imshow(z, extent=[0, 1, 0, 1], origin="lower", cmap=cmap, vmin=-vlim, vmax=vlim,
              transform=tr, interpolation="bilinear", aspect="auto", zorder=z_)
    quad = np.array([(x0, y0), (x0 + w, y0), (x0 + w + skew * h, y0 + h), (x0 + skew * h, y0 + h)])
    cv.add_patch(Polygon(quad, closed=True, facecolor="none", edgecolor=INK,
                         linewidth=lw(1.4), zorder=z_ + 1))
    return quad


def optima():
    """A direction for this model: let the abiotic optima vary across the
    landscape too, as Gaussian random fields of their own.  Four layers
    instead of two: the two mean-trait surfaces, each pulled toward its own
    optimum surface, and tied to each other by biotic selection.  Planned
    work: the sheets are illustrations, nothing here is a result."""
    from matplotlib.colors import LinearSegmentedColormap
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(21)
    env = LinearSegmentedColormap.from_list("env", ["#FFFFFF", "#C9D1D9", "#4A5560"])
    tH = field(120, 0.30, rng)
    tP = 0.75 * tH + 0.66 * field(120, 0.30, rng)            # optima that partly covary
    zH = 0.7 * tH + 0.7 * field(120, 0.22, rng)
    zP = 0.5 * tP + 0.4 * zH + 0.75 * field(120, 0.07, rng)
    W, H, SK = 270, 96, 0.55
    cols = [(70, zH, tH, "host", OCHRE), (470, zP, tP, "parasite", INDIGO)]
    tops = []
    for x0, z, th, who, col in cols:
        q_env = _layer(cv, th, x0, 96, W, H, SK, env, z_=3)
        q_tr = _layer(cv, z, x0, 300, W, H, SK, "cividis", z_=6)
        tops.append(q_tr)
        for u in (0.2, 0.5, 0.8):                             # each place pulled to its optimum
            xa = x0 + W * u + SK * H * 0.5
            arrow(pen, (xa, 96 + H * 0.5 + 8), (xa, 300 + H * 0.5 - 8), px=1.6, color=MUTED,
                  head=8, z=5)
        cv.text(x0 + SK * H + 4, 300 + H + 18, f"{who} mean trait", ha="left", va="center",
                fontsize=fs(FS_SMALL), color=col, fontweight="bold")
        cv.text(x0 + 4, 78, f"{who}'s abiotic optimum", ha="left", va="center",
                fontsize=fs(FS_SMALL), color=INK, fontweight="bold")
    cv.set_xlim(0, fig._px[0])
    cv.set_ylim(0, fig._px[1])
    # biotic selection ties the two trait surfaces; the two optima may covary
    xa, xb = 70 + W + SK * H * 0.5 + 10, 470 + SK * H * 0.5 - 10
    arrow(pen, (xa, 358), (xb, 358), px=2.4, color=INDIGO, head=10)
    arrow(pen, (xb, 338), (xa, 338), px=2.4, color=OCHRE, head=10)
    cv.text((xa + xb) / 2, 382, "coevolution", ha="center", va="center", fontsize=fs(FS_TINY),
            color=INK)
    cv.plot([xa, xb], [144, 144], color=MUTED, lw=lw(1.6), ls=(0, (3, 2.4)))
    cv.text((xa + xb) / 2, 162, "may covary", ha="center", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    cv.text(70 + W * 0.5 + SK * H * 0.5 + 12, 246, "abiotic\nselection", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED, linespacing=1.15)
    cv.text(40, 520, "Next: let the optima vary across the landscape too", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(40, 494, "so far each species had one optimum everywhere; now each optimum is a random surface of its own",
            ha="left", va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    # what it is for
    X = 870
    cv.plot([X - 30, X - 30], [60, 470], color=LINE, lw=lw(1.2))
    cv.text(X, 430, "Four random fields", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold")
    cv.text(X, 404, "instead of two", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold")
    cv.text(X, 330, "How does adaptation to the\nbackground interact with\ncoevolution to shape\nspatial patterns?",
            ha="left", va="center", fontsize=fs(FS_SMALL), linespacing=1.3)
    cv.text(X, 196, "It matters for any method\nthat looks for the signature\nof coevolution in data.",
            ha="left", va="center", fontsize=fs(FS_SMALL), color=NAVY, fontweight="bold",
            linespacing=1.3)
    status(cv, "planned", "a planned extension of Week & Bradburd 2024 (see their Discussion); the surfaces are illustrations")
    save(fig, "s17a_optima")


# =================================================================== s17 ===
def _hankel(d, lam_outer, lam_inner):
    """Order-0 Hankel transform of 1 / [(1 + a k^2)(1 + b k^2)^2] with
    a = lam_outer^2 / 2 and b = lam_inner^2 / 2, by partial fractions."""
    a, b = lam_outer ** 2 / 2, lam_inner ** 2 / 2
    if abs(a - b) < 1e-9 * a:
        a *= 1 + 1e-6
    A, B, C = a * a / (a - b) ** 2, -a * b / (b - a) ** 2, b / (b - a)
    d = np.maximum(d, 1e-9)
    return (A / a * kv(0, d / np.sqrt(a)) + B / b * kv(0, d / np.sqrt(b))
            + C * d / (2 * b ** 1.5) * kv(1, d / np.sqrt(b)))


def cross_cov(d, lH, lP, VH=1.0, VP=1.0, bH=1 / 250, bP=1 / 250):
    """Host-parasite cross-covariance at distance d under weak biotic selection
    (supplement, Eq. S11c); bS = B_S / A_S."""
    return bP * lH ** 2 * VH * _hankel(d, lP, lH) - bH * lP ** 2 * VP * _hankel(d, lH, lP)


def local_adaptation(d, lH, lP, **kw):
    """Parasite local adaptation (Eq. 7b), scaled by (B/A) V."""
    return (cross_cov(1e-9, lH, lP, **kw) - cross_cov(d, lH, lP, **kw)) * 250.0


def _shade(ax, top):
    ax.axhspan(0, top, color=INDIGO_T, lw=0, zorder=0)
    ax.axhspan(-top, 0, facecolor=OCHRE_T, edgecolor=PAPER, hatch="///", lw=0, zorder=0)
    ax.axhline(0, color=FAINT, lw=lw(1.4), zorder=1)


def scale():
    """Three statements, left to right.
    (a) Host and parasite local adaptation mirror each other (Eqs 7a-b), drawn
        for one case of the paper's Fig. 3 with equal biotic selection.
    (b) Which one is ahead is set by the two spatial scales: the paper's Fig. 3
        grid, recomputed (parasite local adaptation in each panel).
    (c) With unequal abiotic selection and population sizes as well, which one
        is ahead can change with distance: the paper's Fig. 4, from the
        repository's values."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    d = np.linspace(0, 80, 300)

    # ---- (a) mirror images ---------------------------------------------------------
    cv.text(44, 524, "Host and parasite", ha="left", va="baseline", fontsize=fs(FS_SMALL),
            fontweight="bold")
    cv.text(44, 502, "mirror each other", ha="left", va="baseline", fontsize=fs(FS_SMALL),
            fontweight="bold")
    ax = px_axes(fig, 70, 150, 250, 310)
    y = local_adaptation(d, 10.0, 1.0)
    ax.axhline(0, color=FAINT, lw=lw(1.4), zorder=1)
    ax.plot(d, y, color=INDIGO, lw=lw(LW_DATA + 0.4), zorder=4)
    ax.plot(d, -y, color=OCHRE, lw=lw(LW_DATA + 0.4), ls=(0, (3.2, 1.6)), zorder=4)
    ax.text(78, y[-1] + 0.14, "parasite", ha="right", va="bottom", fontsize=fs(FS_SMALL),
            color=INDIGO, fontweight="bold")
    ax.text(78, -y[-1] - 0.14, "host", ha="right", va="top", fontsize=fs(FS_SMALL),
            color=OCHRE, fontweight="bold")
    ax.text(2, 1.22, "locally adapted", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    ax.text(2, -1.22, "locally maladapted", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    ax.set_xlim(0, 80)
    ax.set_ylim(-1.35, 1.35)
    ax.set_xticks([0, 40, 80])
    ax.set_yticks([0])
    ax.set_xlabel(r"distance $d$")
    ax.set_ylabel("local adaptation", labelpad=4)
    cv.text(44, 76, "here the parasite's traits vary", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(44, 56, "over the shorter distance", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)

    # ---- (b) who is ahead: the two spatial scales -------------------------------------
    XB = 392
    cv.plot([XB - 36, XB - 36], [40, 540], color=LINE, lw=lw(1.2))
    cv.text(XB - 12, 524, "Who is ahead is set by", ha="left", va="baseline",
            fontsize=fs(FS_SMALL), fontweight="bold")
    cv.text(XB - 12, 502, "the two spatial scales", ha="left", va="baseline",
            fontsize=fs(FS_SMALL), fontweight="bold")
    lams = (1.0, 10.0, 20.0)
    words = ("short", "medium", "long")
    gx, gy, pw, ph, sx, sy = XB + 78, 150, 86, 82, 96, 96
    cv.text(gx + sx + pw / 2, 466, r"parasite's scale $\lambda_P$", ha="center", va="center",
            fontsize=fs(FS_TINY), color=INDIGO)
    cv.text(XB - 10, gy + sy + ph / 2, r"host's scale $\lambda_H$", ha="center", va="center",
            rotation=90, fontsize=fs(FS_TINY), color=OCHRE)
    ends = {}
    for i, lH in enumerate(lams):
        y0 = gy + (2 - i) * sy
        cv.text(gx - 8, y0 + ph / 2, words[i], ha="right", va="center", fontsize=fs(FS_TINY),
                color=OCHRE)
        for j, lP in enumerate(lams):
            if i == 0:
                cv.text(gx + j * sx + pw / 2, 444, words[j], ha="center", va="center",
                        fontsize=fs(FS_TINY), color=INDIGO)
            ax = px_axes(fig, gx + j * sx, y0, pw, ph)
            yy = local_adaptation(d, lH, lP)
            _shade(ax, 1.35)
            ax.plot(d, yy, color=INK, lw=lw(LW_DATA), zorder=4)
            ax.set_xlim(0, 80)
            ax.set_ylim(-1.35, 1.35)
            ax.set_xticks([])
            ax.set_yticks([])
            ends[(i, j)] = yy[-1]
    assert all(abs(ends[(k, k)]) < 1e-3 for k in range(3))
    assert ends[(0, 1)] < -0.5 and ends[(1, 0)] > 0.5 and ends[(1, 2)] < -0.2 and ends[(2, 1)] > 0.2
    cv.add_patch(Rectangle((XB - 12, 96), 16, 12, facecolor=INDIGO_T, edgecolor="none"))
    cv.text(XB + 10, 102, "parasite locally adapted", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INDIGO)
    cv.add_patch(Rectangle((XB - 12, 74), 16, 12, facecolor=OCHRE_T, edgecolor=PAPER, hatch="///",
                           linewidth=0))
    cv.text(XB + 10, 80, "host locally adapted", ha="left", va="center", fontsize=fs(FS_TINY),
            color=OCHRE)
    cv.text(XB - 12, 54, "shorter scale: ahead; equal scales: neither", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)

    # ---- (c) and it can change with distance ---------------------------------------------
    XC = 812
    cv.plot([XC - 36, XC - 36], [40, 540], color=LINE, lw=lw(1.2))
    cv.text(XC - 12, 524, "With unequal selection and population", ha="left",
            va="baseline", fontsize=fs(FS_SMALL), fontweight="bold")
    cv.text(XC - 12, 502, "sizes too, it can change with distance", ha="left",
            va="baseline", fontsize=fs(FS_SMALL), fontweight="bold")
    xr = np.loadtxt(DATA / "xcov-xrng.csv")[:128]
    for k, case in enumerate("ab"):
        yv = np.loadtxt(DATA / f"fig4{case}-vals.csv") * 1e3
        par = np.loadtxt(DATA / f"fig4{case}-pars.csv")
        Gh, Gp, Ah, Ap, sh, sp_ = par[0], par[1], par[4], par[5], par[8], par[9]
        lh, lp = np.sqrt(sh / (Gh * Ah)), np.sqrt(sp_ / (Gp * Ap))
        m = xr <= 100
        x, y = xr[m], yv[:len(xr)][m]
        ax = px_axes(fig, XC + 6, 316 - k * 166, 322, 142)
        top = 1.3 * np.abs(y).max()
        _shade(ax, top)
        first = "parasite" if y[3] > 0 else "host"
        second = "host" if first == "parasite" else "parasite"
        side = 1 if second == "parasite" else -1
        for ell, col, ls, sym in ((lh, OCHRE, (0, (3.2, 1.6)), r"$\lambda_H$"),
                                  (lp, INDIGO, (0, (1, 1.8)), r"$\lambda_P$")):
            ax.axvline(ell, color=col, lw=lw(1.8), ls=ls, zorder=2)
            ax.text(ell + 1.5, -0.82 * top * side, sym, ha="left", va="center",
                    fontsize=fs(FS_TINY), color=col, zorder=6)
        ax.plot(x, y, color=INK, lw=lw(LW_DATA + 0.4), zorder=4)
        flips = np.where(np.diff(np.sign(y[1:])) != 0)[0]
        assert len(flips) >= 1, "the published curves change sign"
        i = flips[0] + 1
        xc = x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
        ax.plot([xc], [0], marker="o", ms=10, mfc=PAPER, mec=INK, mew=lw(2.2), zorder=5)
        ax.set_xlim(0, 100)
        ax.set_ylim(-top, top)
        ax.set_xticks([0, 50, 100])
        ax.set_yticks([])
        if k == 0:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(r"distance $d$")
        ax.text(98, 0.80 * top * side, f"near: {first} ahead; far: {second}", ha="right",
                va="center", fontsize=fs(FS_TINY), color=INK, zorder=6,
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85))
    status(cv, "published",
           "Week & Bradburd 2024, Am. Nat. 203(1) \u2014 their Fig. 3 recomputed (left, middle); their Fig. 4 from the repository's values (right)")
    save(fig, "s17_scale")


def grid_backup():
    fig = figure_px(FULL)
    cv = canvas(fig)

    # ---- the grid: nine combinations of the two spatial scales ---------------------
    lams = (1.0, 10.0, 20.0)
    d = np.linspace(0, 80, 240)
    gx, gy, pw, ph, sx, sy = 168, 138, 132, 92, 148, 108
    cv.text(60, 528, "Usually one species is ahead at every distance", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(gx + sx + pw / 2, 482, "distance over which the parasite's traits stay alike",
            ha="center", va="center", fontsize=fs(FS_TINY), color=INDIGO)
    cv.text(34, gy + sy + ph / 2, "the same, for the host", ha="center", va="center",
            rotation=90, fontsize=fs(FS_TINY), color=OCHRE)
    words = ("short", "medium", "long")
    ends = {}
    for i, lH in enumerate(lams):                 # rows, top to bottom
        y0 = gy + (2 - i) * sy
        cv.text(gx - 44, y0 + ph / 2, words[i], ha="right", va="center",
                fontsize=fs(FS_TINY), color=OCHRE)
        for j, lP in enumerate(lams):
            x0 = gx + j * sx
            if i == 0:
                cv.text(x0 + pw / 2, 458, words[j], ha="center", va="center",
                        fontsize=fs(FS_TINY), color=INDIGO)
            ax = px_axes(fig, x0, y0, pw, ph)
            y = local_adaptation(d, lH, lP)
            _shade(ax, 1.35)
            ax.plot(d, y, color=INK, lw=lw(LW_DATA), zorder=4)
            ax.set_xlim(0, 80)
            ax.set_ylim(-1.35, 1.35)
            ax.set_xticks([0, 40, 80])
            ax.set_yticks([0])
            if i < 2:
                ax.set_xticklabels([])
            if j > 0:
                ax.set_yticklabels([])
            ends[(i, j)] = y[-1]
    # the paper's pattern: the species whose traits vary over the shorter distances
    # is the locally adapted one; equal scales give none
    assert all(abs(ends[(k, k)]) < 1e-3 for k in range(3))
    assert ends[(0, 1)] < -0.5 and ends[(0, 2)] < -0.5 and ends[(1, 2)] < -0.2
    assert ends[(1, 0)] > 0.5 and ends[(2, 0)] > 0.5 and ends[(2, 1)] > 0.2
    print("  [check] Fig. 3 plateaus (row 1):", ", ".join(f"{ends[(0, j)]:+.2f}" for j in range(3)))
    cv.text(gx + sx + pw / 2, 100, "distance to the hosts the parasite is tested on",
            ha="center", va="center", fontsize=fs(FS_SMALL))
    cv.text(60, 62, "The species whose traits vary over the shorter distances is the one locally adapted.",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(60, 42, "Equal scales: neither.", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)

    # ---- two cases where the answer flips ---------------------------------------------
    X = 730
    cv.plot([X - 62, X - 62], [40, 540], color=LINE, lw=lw(1.2))
    cv.text(X - 36, 528, "With other asymmetries, it can flip", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    xr = np.loadtxt(DATA / "xcov-xrng.csv")[:128]
    for k, case in enumerate("ab"):
        yv = np.loadtxt(DATA / f"fig4{case}-vals.csv") * 1e3
        par = np.loadtxt(DATA / f"fig4{case}-pars.csv")
        Gh, Gp, Ah, Ap, sh, sp_ = par[0], par[1], par[4], par[5], par[8], par[9]
        lh, lp = np.sqrt(sh / (Gh * Ah)), np.sqrt(sp_ / (Gp * Ap))
        m = xr <= 100
        x, y = xr[m], yv[:len(xr)][m]
        ax = px_axes(fig, X, 318 - k * 190, 410, 150)
        top = 1.3 * np.abs(y).max()
        _shade(ax, top)
        ax.axvline(lh, color=OCHRE, lw=lw(1.8), ls=(0, (3.2, 1.6)), zorder=2)
        ax.axvline(lp, color=INDIGO, lw=lw(1.8), ls=(0, (1, 1.8)), zorder=2)
        ax.plot(x, y, color=INK, lw=lw(LW_DATA + 0.4), zorder=4)
        s_ = np.sign(y[1:])
        flips = np.where(np.diff(s_) != 0)[0]
        assert len(flips) >= 1, "the published curves change sign"
        i = flips[0] + 1
        xc = x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
        ax.plot([xc], [0], marker="o", ms=10, mfc=PAPER, mec=INK, mew=lw(2.2), zorder=5)
        ax.set_xlim(0, 100)
        ax.set_ylim(-top, top)
        ax.set_xticks([0, 50, 100])
        ax.set_yticks([0])
        if k == 0:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("distance to the hosts the parasite is tested on",
                          fontsize=fs(FS_SMALL))
        ax.text(98, 0.80 * top, "parasite locally adapted", ha="right", va="center",
                fontsize=fs(FS_TINY), color=INDIGO, fontweight="bold")
        ax.text(98, -0.80 * top, "host locally adapted", ha="right", va="center",
                fontsize=fs(FS_TINY), color=OCHRE, fontweight="bold")
    cv.text(X - 36, 496, "two parameter sets; dashed, dotted: the two spatial scales",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    status(cv, "published",
           "Week & Bradburd 2024, Am. Nat. 203(1) \u2014 left: their Fig. 3, recomputed; right: their Fig. 4, from the repository's values")
    save(fig, "s17x_grid")


def main():
    mosaic()
    pair_in_space()
    spde()
    la_index()
    model()
    scale()
    optima()
    grid_backup()


if __name__ == "__main__":
    main()
