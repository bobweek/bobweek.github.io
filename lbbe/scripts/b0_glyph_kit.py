"""Reference sheet for the glyph kit (not a slide)."""
import numpy as np

from glyphs import *  # noqa: F403
from semstyle import *  # noqa: F403

W, H = 1180, 700
COLS, ROWS = 4, 3
CW, CH = W / COLS, (H - 40) / ROWS


def cell(c, r):
    """Centre of the drawing area of cell (c, r) and its label baseline."""
    x0, y0 = c * CW, H - (r + 1) * CH
    return x0 + CW / 2, y0 + CH / 2 + 18, y0 + 22


def label(cv, c, r, text, sub=""):
    cx, _, yb = cell(c, r)
    cv.text(cx, yb + (20 if sub else 0), text, ha="center", va="baseline",
            fontsize=fs(FS_SMALL), fontweight="bold")
    if sub:
        cv.text(cx, yb - 1, sub, ha="center", va="baseline",
                fontsize=fs(FS_TINY), color=MUTED)


def main():
    fig = figure_px((W, H))
    cv = canvas(fig)
    pen = Pen(cv)
    for c in range(1, COLS):
        cv.plot([c * CW, c * CW], [44, H - 10], color=LINE, lw=lw(1))
    for r in range(1, ROWS):
        cv.plot([14, W - 14], [H - r * CH, H - r * CH], color=LINE, lw=lw(1))

    # 0,0 species
    cx, cy, _ = cell(0, 0)
    species(pen, cx - 60, cy + 8, 13, "A")
    species(pen, cx + 60, cy + 8, 13, "B")
    cv.text(cx - 60, cy - 28, "guild A", ha="center", fontsize=fs(FS_TINY), color=INDIGO)
    cv.text(cx + 60, cy - 28, "guild B", ha="center", fontsize=fs(FS_TINY), color=OCHRE)
    label(cv, 0, 0, "species", "shape + colour = guild")

    # 1,0 pair
    cx, cy, _ = cell(1, 0)
    species(pen, cx - 52, cy, 13, "A")
    species(pen, cx + 52, cy, 13, "B")
    recip(pen, (cx - 52, cy), (cx + 52, cy), gap=19, colors=(INDIGO, OCHRE))
    label(cv, 1, 0, "pairwise coevolution", "reciprocal selection")

    # 2,0 star
    cx, cy, _ = cell(2, 0)
    top = (cx, cy + 42)
    bots = [(cx + dx, cy - 42) for dx in (-90, -45, 0, 45, 90)]
    for q in bots:
        link(pen, top, q)
    link(pen, top, bots[2], px=3.0, color=INK)
    species(pen, *top, 12, "A")
    for q in bots:
        species(pen, *q, 11, "B")
    label(cv, 2, 0, "one species, its guild", "focal link in ink")

    # 3,0 bipartite
    cx, cy, _ = cell(3, 0)
    top = [(cx + dx, cy + 42) for dx in (-78, -26, 26, 78)]
    bot = [(cx + dx, cy - 42) for dx in (-100, -50, 0, 50, 100)]
    A = [[1, 1, 0, 0, 0], [1, 1, 1, 0, 0], [0, 1, 1, 1, 0], [0, 0, 1, 1, 1]]
    bipartite(pen, top, bot, A)
    label(cv, 3, 0, "bipartite community", "links only across guilds")

    # 0,1 unipartite
    cx, cy, _ = cell(0, 1)
    ang = np.linspace(0, 2 * np.pi, 7)[:-1] + 0.3
    pts = [(cx + 62 * np.cos(a), cy + 50 * np.sin(a)) for a in ang]
    Wm = np.zeros((6, 6))
    for i, j in [(0, 1), (0, 2), (1, 3), (2, 4), (3, 5), (4, 5), (0, 4), (1, 5), (2, 3)]:
        Wm[i, j] = 1
    unipartite(pen, pts, Wm)
    label(cv, 0, 1, "one guild (diffuse)", "everyone can meet everyone")

    # 1,1 niche overlap
    cx, cy, _ = cell(1, 1)
    base = cy - 46
    cv.plot([cx - 120, cx + 120], [base, base], color=FAINT, lw=lw(1.4))
    xs = np.linspace(cx - 120, cx + 120, 300)
    g1 = base + 78 * np.exp(-0.5 * ((xs - (cx - 30)) / 30) ** 2)
    g2 = base + 78 * np.exp(-0.5 * ((xs - (cx + 34)) / 30) ** 2)
    cv.fill_between(xs, base, np.minimum(g1, g2), facecolor=SOFT, edgecolor=MUTED,
                    hatch="////", linewidth=0, zorder=2)
    cv.plot(xs, g1, color=INDIGO, lw=lw(2.6))
    cv.plot(xs, g2, color=INDIGO, lw=lw(2.6), ls=(0, (4, 2)))
    label(cv, 1, 1, "niche overlap", "hatched = competition")

    # 2,1 landscape
    cx, cy, _ = cell(2, 1)
    quad = landscape(pen, cx - 120, cy - 44, 200, 92, seed=3)
    for (u, v), g in [((0.2, 0.3), "A"), ((0.55, 0.65), "A"), ((0.8, 0.25), "A")]:
        p = on_tile(quad, u, v)
        species(pen, p[0], p[1], 8, g)
    label(cv, 2, 1, "landscape", "the same system, many places")

    # 3,1 host
    cx, cy, _ = cell(3, 1)
    host(pen, cx - 62, cy, 40)
    host(pen, cx + 50, cy, 40, n=1, seed=2, taxa=1, micro_s=5)
    label(cv, 3, 1, "host  ·  host + one microbe", "a unit that contains its partner")

    # 0,2 host + microbiome
    cx, cy, _ = cell(0, 2)
    host(pen, cx, cy, 56, n=16, seed=5, links=9)
    label(cv, 0, 2, "host + microbiome", "marker shape = taxon")

    # 1,2 hosts-microbiome
    cx, cy, _ = cell(1, 2)
    hs = [(cx - 84, cy + 6), (cx, cy - 12), (cx + 84, cy + 6)]
    arc(pen, hs[0], hs[1], bulge=-0.18, color=FAINT, px=1.8)
    arc(pen, hs[1], hs[2], bulge=-0.18, color=FAINT, px=1.8)
    arc(pen, hs[0], hs[2], bulge=0.26, color=FAINT, px=1.8)
    for i, (hx, hy) in enumerate(hs):
        host(pen, hx, hy, 33, n=8, seed=10 + i, micro_s=3.2)
    label(cv, 1, 2, "hosts–microbiome", "hosts as patches; contacts as edges")

    # 2,2 lineage
    cx, cy, _ = cell(2, 2)
    ys = np.linspace(cy - 52, cy + 56, 5)
    hx = [cx - 40, cx - 40, cx - 12, cx - 12, cx - 40]
    lineage_band(pen, list(zip(hx, ys)))
    for x_, y_ in zip(hx, ys):
        cv.add_patch(Circle((x_, y_), 5, facecolor=MUTED, edgecolor="none", zorder=5))
    mx = [cx - 40, cx + 14, cx - 12, cx + 44, cx + 44]
    ancestry(pen, list(zip(mx, ys)))
    label(cv, 2, 2, "lineages", "band = host · dashed = microbe")

    # 3,2 timescales + status
    cx, cy, _ = cell(3, 2)
    tick_timeline(pen, cx - 100, cy + 38, 200, 16, color=TEAL)
    tick_timeline(pen, cx - 100, cy + 2, 200, 3, color=INK)
    for i, (kind, word) in enumerate([("published", "published"), ("inprep", "in prep."),
                                      ("planned", "planned"), ("schematic", "schematic")]):
        xx = cx - 96 + (i % 2) * 112
        yy = cy - 30 - (i // 2) * 24
        status_mark(cv, kind, xx, yy)
        cv.text(xx + 13, yy, word, ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    label(cv, 3, 2, "fast / slow  ·  status marks")

    status(cv, "schematic", "glyph kit, batch 0 — reference sheet, not a slide")
    save(fig, "b0_glyph_kit")


if __name__ == "__main__":
    main()
