"""
glyphs.py -- the seminar's shared visual vocabulary.

The notation follows the notebook sketch the talk's map grew out of:

    small filled marker      a species (population) at the "ordinary" scale
        circle  = guild A (pollinators, parasites, competitors, ...)
        square  = guild B (plants, hosts in a host-parasite pair, ...)
    large open circle        a host: a unit that *contains* its partners
    tiny markers inside it   microbes (shape = taxon)
    two curved arrows        reciprocal selection between a pair
    straight edge            an interaction in a network
    arc                      an interaction within one guild (unipartite)
    tilted tile              a landscape
    soft band / dashed line  host lineage / microbial ancestry

Every glyph takes a `Pen`: the axes plus `k`, the number of slide pixels per
scene unit.  Sizes are in scene units (so they zoom with the camera); stroke
widths are in slide pixels at k = 1 and are scaled by k, so a zoomed frame
looks like the same drawing seen closer.
"""
from __future__ import annotations

import numpy as np
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, PathPatch,
                                Polygon, Rectangle, RegularPolygon)
from matplotlib.path import Path as MPath

from semstyle import *  # noqa: F401,F403


class Pen:
    def __init__(self, ax, k: float = 1.0):
        self.ax = ax
        self.k = k

    def lw(self, px: float) -> float:
        return px * PT * self.k

    def pts(self, units: float) -> float:
        """Scene units -> points (for arrow shrink / mutation_scale)."""
        return units * PT * self.k


# colour + shape per guild: redundancy for greyscale
GUILD = {
    "A": dict(color=INDIGO, tint=INDIGO_T, shape="circle"),
    "B": dict(color=OCHRE, tint=OCHRE_T, shape="square"),
    "N": dict(color=INK, tint=SOFT, shape="circle"),       # neutral species
}


def _state(color, tint, state):
    """face, edge, alpha for on / open / muted / ghost."""
    if state == "on":
        return color, PAPER, 1.0
    if state == "open":
        return tint, color, 1.0
    if state == "muted":
        return FAINT, PAPER, 1.0
    if state == "ghost":
        return LINE, PAPER, 1.0
    raise KeyError(state)


def species(pen, x, y, r=11, guild="A", state="on", z=5, edge_px=2.0):
    g = GUILD[guild]
    face, edge, alpha = _state(g["color"], g["tint"], state)
    if g["shape"] == "circle":
        p = Circle((x, y), r, facecolor=face, edgecolor=edge,
                   linewidth=pen.lw(edge_px), zorder=z, alpha=alpha)
    else:
        s = r * 1.78
        p = Rectangle((x - s / 2, y - s / 2), s, s, facecolor=face,
                      edgecolor=edge, linewidth=pen.lw(edge_px), zorder=z,
                      alpha=alpha, joinstyle="round")
    pen.ax.add_patch(p)
    return p


def link(pen, p, q, px=1.6, color=FAINT, ls="-", z=2, alpha=1.0):
    return pen.ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=pen.lw(px),
                       ls=ls, zorder=z, alpha=alpha, solid_capstyle="round")[0]


def arc(pen, p, q, bulge=0.35, px=1.6, color=FAINT, ls="-", z=2, alpha=1.0):
    """Quadratic arc from p to q; bulge is the sagitta as a fraction of |pq|,
    positive to the left of the direction p -> q."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    d = q - p
    n = np.array([-d[1], d[0]])
    c = (p + q) / 2 + 2 * bulge * n          # control point (2x sagitta)
    path = MPath([p, c, q], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
    patch = PathPatch(path, facecolor="none", edgecolor=color,
                      linewidth=pen.lw(px), linestyle=ls, zorder=z, alpha=alpha,
                      capstyle="round")
    pen.ax.add_patch(patch)
    return patch


def arrow(pen, p, q, px=2.2, color=INK, rad=0.0, head=11, shrink=(0, 0), z=4,
          ls="-", both=False):
    """Arrow p -> q.  `head` = arrowhead length in slide px at k = 1;
    `shrink` in scene units."""
    hl = pen.pts(head)
    kind = "<|-|>" if both else "-|>"
    a = FancyArrowPatch(p, q, connectionstyle=f"arc3,rad={rad}",
                        arrowstyle=f"{kind},head_length={hl},head_width={0.40 * hl}",
                        mutation_scale=1.0, linewidth=pen.lw(px), color=color,
                        zorder=z, shrinkA=pen.pts(shrink[0]),
                        shrinkB=pen.pts(shrink[1]), linestyle=ls,
                        capstyle="round", joinstyle="miter")
    pen.ax.add_patch(a)
    return a


def recip(pen, p, q, gap=16, rad=0.42, px=2.4, colors=(INK, INK), head=11, z=4):
    """Reciprocal selection: two curved arrows, p -> q and q -> p."""
    arrow(pen, p, q, px=px, color=colors[0], rad=rad, head=head,
          shrink=(gap, gap), z=z)
    arrow(pen, q, p, px=px, color=colors[1], rad=rad, head=head,
          shrink=(gap, gap), z=z)


# ------------------------------------------------------------ host, microbes --
_TAXA = [("o", None), ("s", None), ("^", None)]


def microbe_positions(x, y, r, n, seed=0, margin=0.80):
    """n points inside a disc, spread out (simple dart throwing)."""
    rng = np.random.default_rng(seed)
    pts = []
    dmin = 1.55 * r * margin / np.sqrt(max(n, 1))
    tries = 0
    while len(pts) < n and tries < 20000:
        tries += 1
        a, rr = rng.uniform(0, 2 * np.pi), r * margin * np.sqrt(rng.uniform())
        c = np.array([x + rr * np.cos(a), y + rr * np.sin(a)])
        if all(np.hypot(*(c - p)) > dmin for p in pts):
            pts.append(c)
        if tries % 4000 == 0:
            dmin *= 0.9
    return np.array(pts)


def microbe(pen, x, y, s=4.2, taxon=0, color=None, z=7, state="on"):
    """One microbe.  Shape encodes taxon; colour is redundant."""
    col = color or [TEAL, INK, OCHRE][taxon % 3]
    if state == "muted":
        col = FAINT
    if taxon % 3 == 0:
        p = Circle((x, y), s, facecolor=col, edgecolor="none", zorder=z)
    elif taxon % 3 == 1:
        p = Rectangle((x - s * 0.88, y - s * 0.88), 1.76 * s, 1.76 * s,
                      facecolor=col, edgecolor="none", zorder=z)
    else:
        p = RegularPolygon((x, y - 0.15 * s), 3, radius=1.35 * s, facecolor=col,
                           edgecolor="none", zorder=z)
    pen.ax.add_patch(p)
    return p


def host(pen, x, y, r=48, state="on", n=0, seed=0, taxa=3, links=0, edge_px=3.0,
         micro_s=None, z=4, color=None, tint=None):
    """A host (large open circle), optionally with a microbiome inside.

    n      number of microbes drawn
    taxa   number of taxa (marker shapes) among them
    links  number of microbe-microbe interaction edges to draw (community)
    """
    col, tnt = color or TEAL, tint or TEAL_T
    face, edge = (tnt, col)
    if state == "muted":
        face, edge = SOFT, FAINT
    elif state == "ghost":
        face, edge = PAPER, LINE
    pen.ax.add_patch(Circle((x, y), r, facecolor=face, edgecolor=edge,
                            linewidth=pen.lw(edge_px), zorder=z))
    pts = np.zeros((0, 2))
    if n:
        pts = microbe_positions(x, y, r, n, seed)
        rng = np.random.default_rng(seed + 101)
        tx = rng.integers(0, taxa, size=len(pts))
        ms = micro_s or max(2.6, r * 0.085)
        if links:
            # nearest-neighbour style edges: a small community network
            d = np.hypot(pts[:, None, 0] - pts[None, :, 0],
                         pts[:, None, 1] - pts[None, :, 1])
            iu = np.triu_indices(len(pts), 1)
            order = np.argsort(d[iu])[:links]
            for o in order:
                i, j = iu[0][o], iu[1][o]
                link(pen, pts[i], pts[j], px=1.3,
                     color=MUTED if state == "on" else LINE, z=z + 1)
        for (px_, py_), t in zip(pts, tx):
            microbe(pen, px_, py_, ms, int(t), z=z + 2,
                    state="on" if state == "on" else "muted")
    return pts


# ------------------------------------------------------------------ networks --
def bipartite(pen, top, bot, A, r=11, edge_px=1.5, edge_color=FAINT,
              top_state="on", bot_state="on", z=2):
    """Two rows of species and the links between them.  A[i, j] = 1 links
    top[i] to bot[j]."""
    for i, p in enumerate(top):
        for j, q in enumerate(bot):
            if A[i][j]:
                link(pen, p, q, px=edge_px, color=edge_color, z=z)
    for p in top:
        species(pen, p[0], p[1], r, "A", top_state)
    for q in bot:
        species(pen, q[0], q[1], r, "B", bot_state)


def unipartite(pen, pts, W, r=11, guild="A", px=1.6, color=FAINT, bulge=0.0,
               state="on"):
    n = len(pts)
    for i in range(n):
        for j in range(i + 1, n):
            if W[i][j]:
                if bulge:
                    arc(pen, pts[i], pts[j], bulge=bulge, px=px, color=color)
                else:
                    link(pen, pts[i], pts[j], px=px, color=color)
    for p in pts:
        species(pen, p[0], p[1], r, guild, state)


# ------------------------------------------------------------------ landscape --
def landscape(pen, x, y, w, h, skew=0.32, color=MUTED, fill=SOFT, n_lines=3,
              seed=0, z=1, px=1.6):
    """A tilted tile standing for geographic space, with a few contour lines.
    (x, y) is the lower-left corner of the *base*; the top edge is shifted
    right by skew * h."""
    s = skew * h
    quad = np.array([[x, y], [x + w, y], [x + w + s, y + h], [x + s, y + h]])
    pen.ax.add_patch(Polygon(quad, closed=True, facecolor=fill, edgecolor=color,
                             linewidth=pen.lw(px), zorder=z, joinstyle="round"))
    rng = np.random.default_rng(seed)
    u = np.linspace(0.04, 0.96, 60)
    for i in range(n_lines):
        v0 = (i + 1) / (n_lines + 1)
        ph, amp = rng.uniform(0, 2 * np.pi), rng.uniform(0.05, 0.11)
        v = v0 + amp * np.sin(2 * np.pi * u * rng.uniform(0.8, 1.4) + ph)
        v = np.clip(v, 0.05, 0.95)
        px_ = x + u * w + s * v
        py_ = y + v * h
        pen.ax.plot(px_, py_, color=LINE if fill != PAPER else LINE,
                    lw=pen.lw(1.2), zorder=z + 0.1)
    return quad


def on_tile(quad, u, v):
    """Point at fractional position (u, v) on a landscape tile."""
    a, b, c, d = quad
    return (1 - v) * ((1 - u) * a + u * b) + v * ((1 - u) * d + u * c)


# ------------------------------------------------------- lineages, timescales --
def lineage_band(pen, pts, px=11, color=LINE, z=1):
    pts = np.asarray(pts, float)
    return pen.ax.plot(pts[:, 0], pts[:, 1], color=color, lw=pen.lw(px),
                       solid_capstyle="round", solid_joinstyle="round",
                       zorder=z)[0]


def ancestry(pen, pts, px=2.0, color=TEAL, z=6, marker=True, ms=3.8):
    pts = np.asarray(pts, float)
    ln = pen.ax.plot(pts[:, 0], pts[:, 1], color=color, lw=pen.lw(px),
                     ls=(0, (3.2, 2.2)), zorder=z)[0]
    if marker:
        for px_, py_ in pts:
            pen.ax.add_patch(Rectangle((px_ - ms, py_ - ms), 2 * ms, 2 * ms,
                                       facecolor=color, edgecolor="none",
                                       zorder=z + 1))
    return ln


def tick_timeline(pen, x, y, w, n, color=INK, tick=9, px=2.0, arrow_head=8):
    """A time axis with n event ticks: dense = fast process, sparse = slow."""
    arrow(pen, (x, y), (x + w, y), px=1.4, color=FAINT, head=arrow_head, z=3)
    for t in np.linspace(0.04, 0.88, n):
        pen.ax.plot([x + t * w] * 2, [y - tick, y + tick], color=color,
                    lw=pen.lw(px), zorder=4, solid_capstyle="butt")


# ------------------------------------------------------------------ niche --
def niche(pen, mu, base, sd, height, color=INDIGO, fill=None, px=2.4, z=3,
          hatch=None, xs=None):
    xs = np.linspace(mu - 3.2 * sd, mu + 3.2 * sd, 120) if xs is None else xs
    ys = base + height * np.exp(-0.5 * ((xs - mu) / sd) ** 2)
    if fill or hatch:
        pen.ax.fill_between(xs, base, ys, facecolor=fill or "none",
                            edgecolor=color if hatch else "none", hatch=hatch,
                            linewidth=0, zorder=z - 0.5, alpha=0.9)
    pen.ax.plot(xs, ys, color=color, lw=pen.lw(px), zorder=z)
    return xs, ys


# ------------------------------------------------------------------ organisms --
def fly(pen, hx, hy, s=22, proboscis=38, color=INK, z=6):
    """Moegistorhynchus longirostris in side view, hovering, proboscis hanging
    straight down (as nemestrinids carry it).  (hx, hy) is where the proboscis
    leaves the head; the body lies to the right.  Head + body is 1.6 s long,
    so s = 10 mm gives the 16 mm fly of Pauw et al. (2009).
    Returns the (top, tip) points of the proboscis."""
    ax = pen.ax
    ax.add_patch(Ellipse((hx + 0.98 * s, hy + 0.66 * s), 1.10 * s, 0.30 * s,
                         angle=32, facecolor=SOFT, edgecolor=MUTED,
                         linewidth=pen.lw(1.2), zorder=z))
    for i, (x0, dx) in enumerate([(0.38, 0.10), (0.56, 0.20), (0.72, 0.30)]):
        ax.plot([hx + x0 * s, hx + (x0 + dx) * s, hx + (x0 + dx + 0.10) * s],
                [hy - 0.02 * s, hy - 0.24 * s, hy - 0.46 * s], color=color,
                lw=pen.lw(1.1), zorder=z, solid_capstyle="round")
    ax.add_patch(Ellipse((hx + 1.06 * s, hy + 0.13 * s), 0.80 * s, 0.44 * s,
                         angle=-8, facecolor=color, edgecolor="none", zorder=z + 1))
    ax.add_patch(Ellipse((hx + 0.52 * s, hy + 0.18 * s), 0.62 * s, 0.52 * s,
                         facecolor=color, edgecolor="none", zorder=z + 1))
    ax.add_patch(Circle((hx + 0.13 * s, hy + 0.11 * s), 0.18 * s, facecolor=color,
                        edgecolor="none", zorder=z + 1))
    ax.add_patch(Ellipse((hx + 0.11 * s, hy + 0.13 * s), 0.15 * s, 0.20 * s,
                         angle=-15, facecolor=FAINT, edgecolor="none", zorder=z + 2))
    top, tip = (hx, hy), (hx, hy - proboscis)
    ax.plot([top[0], tip[0]], [top[1], tip[1]], color=color, lw=pen.lw(1.7),
            zorder=z + 3, solid_capstyle="round")
    return top, tip


def flower(pen, cx, cy, s=20, tube=34, color=OCHRE, tint=OCHRE_T, z=6):
    """A long-tubed iris (Lapeirousia anceps) in side view: a slender tube with
    narrow tepals at its mouth.  (cx, cy) is the mouth; the tube runs straight
    down for `tube` units.  Returns (mouth, base)."""
    ax = pen.ax
    for ang in (158, 122, 90, 58, 22, -24, 204):
        a = np.radians(ang)
        L = (0.62 if 0 < ang < 180 else 0.50) * s
        ax.add_patch(Ellipse((cx + 0.5 * L * np.cos(a), cy + 0.5 * L * np.sin(a)),
                             L, 0.20 * s, angle=ang, facecolor=tint,
                             edgecolor=color, linewidth=pen.lw(1.3), zorder=z))
    tw = 0.075 * s
    ax.add_patch(Polygon([[cx - tw, cy], [cx + tw, cy],
                          [cx + 0.6 * tw, cy - tube], [cx - 0.6 * tw, cy - tube]],
                         closed=True, facecolor=tint, edgecolor=color,
                         linewidth=pen.lw(1.3), zorder=z + 1, joinstyle="round"))
    ax.add_patch(Ellipse((cx, cy - tube - 0.10 * s), 0.22 * s, 0.26 * s,
                         facecolor=color, edgecolor="none", zorder=z + 1))
    return (cx, cy), (cx, cy - tube)


def dim_bracket(pen, x, y0, y1, color=MUTED, px=1.3, cap=3.5):
    """Vertical dimension line with end caps."""
    ax = pen.ax
    ax.plot([x, x], [y0, y1], color=color, lw=pen.lw(px), zorder=6)
    for yy in (y0, y1):
        ax.plot([x - cap, x + cap], [yy, yy], color=color, lw=pen.lw(px),
                zorder=6)


# ------------------------------------------------------------------ pipeline --
def chevron(pen, x, y, w, h, color=NAVY, fill=None, z=3, px=1.6, tip=0.22,
            flip=False):
    """A pipeline stage: a block with an arrow point.  flip = points left."""
    t = tip * h
    if not flip:
        pts = [[x, y], [x + w - t, y], [x + w, y + h / 2], [x + w - t, y + h],
               [x, y + h], [x + t, y + h / 2]]
    else:
        pts = [[x + t, y], [x + w, y], [x + w - t, y + h / 2], [x + w, y + h],
               [x + t, y + h], [x, y + h / 2]]
    p = Polygon(pts, closed=True, facecolor=fill or PAPER, edgecolor=color,
                linewidth=pen.lw(px), zorder=z, joinstyle="round")
    pen.ax.add_patch(p)
    return p


# ===================================================================== systems
# The systems of the notebook sketch, each drawn about a centre (cx, cy) at a
# unit size u (u = 1 gives a glyph roughly 150 px wide).  These are what the
# map (slides 3, 31), the hinge (21, 22) and the title card are built from.

def sys_pair(pen, cx, cy, u=1.0, guilds=("A", "B")):
    r = 11 * u
    p, q = (cx - 46 * u, cy), (cx + 46 * u, cy)
    recip(pen, p, q, gap=r + 6 * u, rad=0.48, px=2.2 * u,
          colors=(GUILD[guilds[0]]["color"], GUILD[guilds[1]]["color"]),
          head=9 * u)
    species(pen, *p, r, guilds[0])
    species(pen, *q, r, guilds[1])


def sys_star(pen, cx, cy, u=1.0, n=5, focal=2):
    r = 10 * u
    top = (cx, cy + 40 * u)
    bots = [(cx + dx * u, cy - 40 * u) for dx in np.linspace(-84, 84, n)]
    for i, q in enumerate(bots):
        link(pen, top, q, px=1.5 * u)
    species(pen, *top, r * 1.1, "A")
    for q in bots:
        species(pen, *q, r, "B")


_BIP = [[1, 1, 0, 0, 0], [1, 1, 1, 0, 0], [0, 1, 1, 1, 0], [0, 0, 1, 1, 1]]


def sys_bip(pen, cx, cy, u=1.0):
    r = 10 * u
    top = [(cx + dx * u, cy + 40 * u) for dx in (-66, -22, 22, 66)]
    bot = [(cx + dx * u, cy - 40 * u) for dx in (-84, -42, 0, 42, 84)]
    bipartite(pen, top, bot, _BIP, r=r, edge_px=1.5 * u)
    return top, bot


def sys_uni(pen, cx, cy, u=1.0, guild="A"):
    """One guild: the projection of sys_bip onto its top row (species linked
    when they share a partner), drawn as arcs along a row."""
    r = 10 * u
    pts = [(cx + dx * u, cy - 14 * u) for dx in (-66, -22, 22, 66)]
    M = np.array(_BIP)
    for i in range(4):
        for j in range(i + 1, 4):
            w = int((M[i] & M[j]).sum())
            if w:
                arc(pen, pts[i], pts[j], bulge=0.34 + 0.10 * (j - i - 1),
                    px=(1.0 + 0.9 * w) * u, color=GUILD[guild]["color"], z=3)
    for p in pts:
        species(pen, *p, r, guild)
    return pts


def _taxon_row(pen, xs, y, s, taxa=None):
    for i, x in enumerate(xs):
        microbe(pen, x, y, s, taxon=(taxa[i] if taxa else i))


def sys_host_microbe(pen, cx, cy, u=1.0):
    """A host and one microbe: the pair, with one partner scaled up."""
    R = 27 * u
    p, q = (cx - 40 * u, cy), (cx + 56 * u, cy)
    arrow(pen, p, q, px=2.2 * u, color=TEAL, rad=0.48, head=9 * u,
          shrink=(R + 5 * u, 13 * u))
    arrow(pen, q, p, px=2.2 * u, color=INK, rad=0.48, head=9 * u,
          shrink=(13 * u, R + 5 * u))
    host(pen, *p, R, edge_px=2.6 * u)
    microbe(pen, q[0], q[1], 6.5 * u, taxon=1)


def sys_host_microbiome(pen, cx, cy, u=1.0, n=5):
    """One host and a community of microbes: the star, scaled up on one side."""
    R = 25 * u
    top = (cx, cy + 30 * u)
    xs = [cx + dx * u for dx in np.linspace(-84, 84, n)]
    for x in xs:
        link(pen, top, (x, cy - 44 * u), px=1.5 * u)
    host(pen, *top, R, edge_px=2.6 * u)
    _taxon_row(pen, xs, cy - 44 * u, 6.0 * u)


_HM = [[1, 1, 1, 0, 0, 0, 0], [0, 1, 1, 1, 1, 1, 0], [0, 0, 0, 1, 0, 1, 1]]


def sys_hosts_microbiomes(pen, cx, cy, u=1.0):
    """Many hosts and many microbes: the bipartite web, scaled up on one side."""
    R = 22 * u
    tops = [(cx + dx * u, cy + 32 * u) for dx in (-74, 0, 74)]
    xs = [cx + dx * u for dx in np.linspace(-102, 102, 7)]
    for i, t in enumerate(tops):
        for j, x in enumerate(xs):
            if _HM[i][j]:
                link(pen, t, (x, cy - 44 * u), px=1.5 * u)
    for t in tops:
        host(pen, *t, R, edge_px=2.4 * u)
    _taxon_row(pen, xs, cy - 44 * u, 5.6 * u)


def sys_space(pen, cx, cy, u=1.0, seed=4):
    """A pair spread over a landscape."""
    quad = landscape(pen, cx - 78 * u, cy - 34 * u, 130 * u, 68 * u, seed=seed,
                     px=1.5 * u)
    for (a, b), g in [((0.16, 0.30), "A"), ((0.30, 0.62), "B"),
                      ((0.52, 0.28), "B"), ((0.66, 0.70), "A"),
                      ((0.86, 0.36), "A"), ((0.90, 0.74), "B")]:
        p = on_tile(quad, a, b)
        species(pen, p[0], p[1], 6.5 * u, g, edge_px=1.4 * u)
    return quad


def sys_inherit(pen, cx, cy, u=1.0):
    """A host lineage and a microbe that does, or does not, follow it."""
    R = 17 * u
    par, off = (cx - 10 * u, cy + 34 * u), (cx - 10 * u, cy - 34 * u)
    lineage_band(pen, [par, off], px=10 * u, z=2.2)
    host(pen, *par, R, edge_px=2.2 * u)
    host(pen, *off, R, edge_px=2.2 * u)
    microbe(pen, par[0], par[1], 5.2 * u, taxon=0)
    microbe(pen, off[0], off[1], 5.2 * u, taxon=0)
    arrow(pen, (par[0] + R + 3 * u, par[1] - 6 * u),
          (off[0] + R + 3 * u, off[1] + 6 * u), px=1.8 * u, color=TEAL,
          rad=-0.5, head=8 * u, ls=(0, (3, 2)))
