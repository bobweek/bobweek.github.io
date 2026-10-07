"""
Slide 22 -- host-microbiome evolution combines all three, in one nested picture.

A direct descendant of the geographic-mosaic slide (s14).  The nesting is

    landscape -> locality -> social network of hosts -> microbes within a host

At each locality a patchy environmental pool of microbes lies on the
landscape, and a social network of hosts floats above it, each host carrying
its own microbial community.  Three kinds of movement are drawn differently:
social contact between hosts (edges of the network), uptake from the local
environmental pool (dotted, upward), and hosts moving between localities
(dashed arcs on the landscape).  For one microbial taxon, hosts plus pools
form a metapopulation; over many taxa, a metacommunity.

s22a_hinge       the nested picture, with the three lessons pointed at it;
                 one host is enlarged to show genes and microbes together
s22b_hinge_new   the same picture plus two insets: inheritance and timescales

Schematic.  After Week et al. (2025) Evolution 79:2487 and Week et al. (2025)
Nat. Ecol. Evol. 9:1769.
"""
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch

from glyphs import (Pen, arc, arrow, host, landscape, link, microbe,
                    microbe_positions, on_tile, sys_inherit, tick_timeline)
from semstyle import *  # noqa: F403

TILE = dict(x=56, y=56, w=840, h=276, skew=0.26)
SITES = [((0.16, 0.52), 0), ((0.50, 0.34), 1), ((0.83, 0.60), 2)]     # (u, v), dominant taxon
HOSTS = [(-60, 6), (-24, 44), (24, 40), (58, 4), (0, -6)]
EDGES = [[(0, 1), (1, 2), (2, 3), (1, 4), (3, 4)],
         [(0, 1), (1, 2), (2, 3), (0, 4), (2, 4), (3, 4)],
         [(0, 1), (1, 2), (2, 3), (1, 4), (0, 4)]]
R_HOST = 17


def _host(pen, x, y, r, taxa, seed, edge=1.8, ms=3.0):
    host(pen, x, y, r, edge_px=edge)
    for (mx, my), t in zip(microbe_positions(x, y, r, len(taxa), seed=seed), taxa):
        microbe(pen, mx, my, ms, int(t), z=8)


def scene(cv, pen):
    rng = np.random.default_rng(11)
    quad = landscape(pen, **TILE, n_lines=4, seed=5, px=1.8)
    pts = [on_tile(quad, *uv) for uv, _ in SITES]
    # hosts move between localities, on the ground
    for i, j in ((0, 1), (1, 2)):
        arc(pen, pts[i], pts[j], bulge=-0.16, px=1.8, color=NAVY, ls=(0, (4, 2.6)), z=3)
    nets = []
    for k, ((uv, dom), p) in enumerate(zip(SITES, pts)):
        # the local environmental pool: patchy, mostly one taxon
        cv.add_patch(Ellipse(p, 176, 60, facecolor=PAPER, edgecolor=FAINT,
                             linewidth=lw(1.3), zorder=2))
        for _ in range(13):
            a, rr = rng.uniform(0, 2 * np.pi), np.sqrt(rng.uniform(0, 1))
            t = dom if rng.uniform() < 0.75 else int(rng.integers(0, 3))
            microbe(pen, p[0] + 74 * rr * np.cos(a), p[1] + 22 * rr * np.sin(a), 3.0,
                    taxon=t, state="muted", z=3)
        # the social network of hosts above it
        c = (p[0], p[1] + 112)
        hs = [(c[0] + dx, c[1] + dy) for dx, dy in HOSTS]
        for i, j in EDGES[k]:
            link(pen, hs[i], hs[j], px=2.0, color=MUTED, z=5)
        for n, h in enumerate(hs):
            taxa = [dom if rng.uniform() < 0.7 else int(rng.integers(0, 3)) for _ in range(4)]
            _host(pen, h[0], h[1], R_HOST, taxa, seed=40 + 7 * k + n)
        # uptake from the pool
        for n in (0, 3):
            arrow(pen, (hs[n][0], p[1] + 20), (hs[n][0], hs[n][1] - R_HOST - 3), px=1.5,
                  color=TEAL, head=7, ls=(0, (1.2, 2.2)), z=4)
        nets.append(hs)
    return quad, pts, nets


def lessons(cv, pen, nets):
    cv.text(56, 540, "three lessons from coevolution, in one system", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)
    # one host, enlarged: genes and microbes together
    zx, zy, zr = 112, 458, 52
    src = nets[0][1]
    for ang in (35, 215):
        a = np.radians(ang)
        cv.plot([src[0] + R_HOST * np.cos(a), zx + zr * np.cos(a)],
                [src[1] + R_HOST * np.sin(a), zy + zr * np.sin(a)], color=FAINT,
                lw=lw(1.1), ls=(0, (2, 2)), zorder=1)
    host(pen, zx, zy, zr, edge_px=2.8)
    mp = microbe_positions(zx + 12, zy - 2, zr - 14, 8, seed=3)
    for i, (mx, my) in enumerate(mp):
        microbe(pen, mx, my, 5.0, i % 3, z=8)
    g = (zx - 26, zy + 14)
    for i in (0, 3, 5):
        cv.plot([g[0] + 10, mp[i][0]], [g[1], mp[i][1]], color=TEAL, lw=lw(1.2),
                ls=(0, (1.2, 1.8)), zorder=7)
    cv.add_patch(FancyBboxPatch((g[0] - 13, g[1] - 10), 26, 20,
                                boxstyle="round,pad=0,rounding_size=5", facecolor=PAPER,
                                edgecolor=TEAL, linewidth=lw(1.8), zorder=9))
    t = cv.text(g[0], g[1], r"$g$", ha="center", va="center", fontsize=fs(FS_TINY),
                color=TEAL, zorder=10)
    t._allow_overlap = True
    for y, head, line in [(496, "interface", "genes and microbes build the host trait"),
                          (440, "community", "the partner is many microbial taxa")]:
        cv.text(184, y, head, ha="left", va="center", fontsize=fs(FS_LABEL),
                fontweight="bold", color=INDIGO)
        cv.text(184, y - 23, line, ha="left", va="center", fontsize=fs(FS_TINY),
                color=MUTED)
    # space, three ways
    X = 986
    cv.text(X, 404, "space", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold", color=INDIGO)
    rows = [(356, "social contacts", "microbes pass\nbetween hosts", MUTED, "-"),
            (276, "the environment", "local, patchy pools", TEAL, (0, (1.2, 2.2))),
            (212, "localities", "hosts move between\nthem, with their microbes", NAVY,
             (0, (4, 2.6)))]
    for y, head, line, col, ls in rows:
        cv.plot([X, X + 30], [y, y], color=col, lw=lw(2.2), ls=ls)
        cv.text(X + 40, y, head, ha="left", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold")
        cv.text(X, y - 16, line, ha="left", va="top", fontsize=fs(FS_TINY), color=MUTED,
                linespacing=1.2)


def new_problems(cv, pen):
    X = 548
    cv.text(X - 28, 540, "and two things the pair never had", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    sys_inherit(pen, X, 488, 0.62)
    cv.text(X + 40, 498, "inheritance", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold", color=TEAL)
    cv.text(X + 40, 475, "which microbial ancestry stays with the host lineage?",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    tick_timeline(pen, X - 26, 436, 52, 9, color=TEAL, tick=7, px=1.6, arrow_head=6)
    tick_timeline(pen, X - 26, 420, 52, 2, color=INK, tick=7, px=1.6, arrow_head=6)
    cv.text(X + 40, 440, "timescales", ha="left", va="center", fontsize=fs(FS_LABEL),
            fontweight="bold", color=TEAL)
    cv.text(X + 40, 417, "microbes turn over fast, hosts evolve slowly", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)


def _frame(name, new):
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    quad, pts, nets = scene(cv, pen)
    lessons(cv, pen, nets)
    if new:
        new_problems(cv, pen)
    status(cv, "schematic", "after Week et al. 2025, Evolution 79:2487; Week et al. 2025, Nat. Ecol. Evol. 9:1769")
    save(fig, name)


def main():
    _frame("s22a_hinge", False)
    _frame("s22b_hinge_new", True)


if __name__ == "__main__":
    main()
