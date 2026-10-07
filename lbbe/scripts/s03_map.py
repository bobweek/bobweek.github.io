"""
Slide 03 -- a map of the interaction systems in the talk.  Two candidates.

A  `s03_mapA`  a typology read straight off the notebook sketch.
   The six cartoons differ in two ways only: whether the partners sit at the
   same scale (dots with dots) or one is a host to the other (circle with
   dots), and how many partners there are on each side.  The unipartite case
   is the many-with-many cell seen from inside one guild -- the projection
   that slide 11 performs.  `s03_mapA_route` adds the route of the talk.

B  `s03_mapB`  a chart with bearings instead of axes, after latest/ (the
   tangent chart).  The pair is the origin; each bearing is one way a system
   departs from it.  Systems that depart in two respects sit between bearings.

Neither is data: positions record kind, not measurement, and both say so.
"""
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

from glyphs import *  # noqa: F403
from semstyle import *  # noqa: F403

# =========================================================== version A =====
COLX = [(170, 440), (440, 730), (730, 1170)]       # column spans
ROWY = [(292, 498), (52, 270)]                     # row spans (top, bottom)
CX = [305, 585, None]
CY = [412, 182]


def _name(cv, x, y, title, sub, muted=False):
    cv.text(x, y, title, ha="center", va="baseline", fontsize=fs(FS_SMALL),
            fontweight="bold", color=MUTED if muted else INK)
    cv.text(x, y - 21, sub, ha="center", va="baseline", fontsize=fs(FS_TINY),
            color=MUTED)


def map_A(route=False):
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)

    # frame: column heads, row heads, hairlines
    for (x0, x1), head in zip(COLX, ["one with one", "one with many",
                                     "many with many"]):
        cv.text((x0 + x1) / 2, 528, head, ha="center", va="baseline",
                fontsize=fs(FS_LABEL), color=MUTED, fontweight="bold")
    cv.plot([24, 1170], [512, 512], color=LINE, lw=lw(1.2))
    cv.plot([24, 1170], [281, 281], color=LINE, lw=lw(1.2))
    for x in (COLX[0][1], COLX[1][1]):
        cv.plot([x, x], [52, 505], color=LINE, lw=lw(1.2))
    cv.text(24, 432, "same scale", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold", color=INDIGO)
    cv.text(24, 408, "species with\nspecies", ha="left", va="top",
            fontsize=fs(FS_TINY), color=MUTED, linespacing=1.25)
    cv.text(24, 202, "nested", ha="left", va="baseline", fontsize=fs(FS_LABEL),
            fontweight="bold", color=TEAL)
    cv.text(24, 178, "a host with\nits microbes", ha="left", va="top",
            fontsize=fs(FS_TINY), color=MUTED, linespacing=1.25)

    # ---- row 1: same scale ------------------------------------------------
    sys_pair(pen, CX[0], CY[0] + 6, 1.15)
    _name(cv, CX[0], 336, "pairwise coevolution", "fly and flower, host and parasite")
    sys_star(pen, CX[1], CY[0] + 4, 1.05)
    _name(cv, CX[1], 336, "a species and its guild", "one fly, twenty flowers")
    xb, xu = 838, 1068
    sys_bip(pen, xb, CY[0] + 4, 1.0)
    _name(cv, xb, 336, "between two guilds", "mutualistic networks")
    sys_uni(pen, xu, CY[0] + 4, 1.0)
    _name(cv, xu, 336, "within one guild", "diffuse competition")
    arrow(pen, (xb + 104, CY[0] + 4), (xu - 92, CY[0] + 4), px=1.4, color=MUTED,
          head=8)
    cv.text((xb + 104 + xu - 92) / 2, CY[0] + 13, "seen from\none side",
            ha="center", va="bottom", fontsize=fs(FS_TINY), color=MUTED,
            linespacing=1.15)

    # ---- row 2: nested ----------------------------------------------------
    sys_host_microbe(pen, CX[0] - 8, CY[1] + 8, 1.15)
    _name(cv, CX[0], 98, "host and microbe", "one symbiont")
    sys_host_microbiome(pen, CX[1], CY[1] + 14, 1.05)
    _name(cv, CX[1], 98, "host and microbiome", "one host, a community")
    sys_hosts_microbiomes(pen, 950, CY[1] + 14, 1.12)
    _name(cv, 950, 98, "hosts and microbiomes", "many hosts, shared microbes")

    if route:
        stops = [(196, 486, 1), (466, 486, 2), (756, 486, 3), (986, 486, 4),
                 (466, 256, 5), (196, 256, 6), (756, 256, 7)]
        for x, y, n in stops:
            cv.add_patch(Circle((x, y), 13, facecolor=NAVY, edgecolor=PAPER,
                                linewidth=lw(1.6), zorder=9))
            t = cv.text(x, y - 0.5, str(n), ha="center", va="center",
                        fontsize=fs(FS_TINY), color=PAPER, fontweight="bold",
                        zorder=10)
            t._allow_overlap = True
        note = cv.text(968, 281, "between 4 and 5:  the pair again, in space  \u00B7  then inference",
                       ha="center", va="center", fontsize=fs(FS_TINY), color=NAVY,
                       bbox=dict(boxstyle="round,pad=0.35", facecolor=PAPER,
                                 edgecolor=NAVY, linewidth=lw(1.0)), zorder=9)
    status(cv, "schematic", "the systems of the notebook sketch; position records kind, not measurement")
    save(fig, "s03_mapA_route" if route else "s03_mapA")


# =========================================================== version B =====
SERIF = dict(family="serif")


def _polar(c, r, bearing):
    """Bearing clockwise from north, as on a chart."""
    a = np.radians(90 - bearing)
    return c[0] + r * np.cos(a), c[1] + r * np.sin(a)


def _rose(cv, c, R):
    for i in range(8):
        b = i * 45
        tip = _polar(c, R if i % 2 == 0 else 0.55 * R, b)
        l = _polar(c, 0.16 * R, b - 45)
        r_ = _polar(c, 0.16 * R, b + 45)
        cv.add_patch(Polygon([c, l, tip], closed=True, facecolor=INK,
                             edgecolor=INK, linewidth=lw(0.6), zorder=1.6))
        cv.add_patch(Polygon([c, r_, tip], closed=True, facecolor=PAPER,
                             edgecolor=INK, linewidth=lw(0.6), zorder=1.6))


def _token(cv, x, y, text):
    t = cv.text(x, y, text, ha="center", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", color=INK, zorder=12, **SERIF,
                bbox=dict(boxstyle="square,pad=0.32", facecolor=PAPER,
                          edgecolor=INK, linewidth=lw(1.4)))
    return t


def _site(cv, c, R, frac, bearing, label, dx=0, dy=0, ha="left", mark=True):
    x, y = _polar(c, frac * R, bearing)
    if mark:
        cv.add_patch(Circle((x, y), 5.2, facecolor=PAPER, edgecolor=INK,
                            linewidth=lw(1.3), zorder=10))
        cv.add_patch(Circle((x, y), 1.6, facecolor=INK, edgecolor="none",
                            zorder=11))
    cv.text(x + dx, y + dy, label, ha=ha, va="center", fontsize=fs(FS_TINY),
            style="italic", color=INK, zorder=12, linespacing=1.1, **SERIF,
            path_effects=[pe.withStroke(linewidth=lw(6), foreground=PAPER)])
    return x, y


HALO = [pe.withStroke(linewidth=lw(6), foreground=PAPER)]


def _label(cv, x, y, name, example, ha="center"):
    cv.text(x, y, name, ha=ha, va="baseline", fontsize=fs(FS_SMALL),
            fontweight="bold", color=INK, zorder=12, path_effects=HALO, **SERIF)
    cv.text(x, y - 21, example, ha=ha, va="baseline", fontsize=fs(FS_TINY),
            style="italic", color=MUTED, zorder=12, path_effects=HALO, **SERIF)


ENTRIES = [
    # glyph, (x, y, size), (label x, label y, align), name, example, lesson, stop
    (sys_pair, (215, 287, 0.95), (215, 240, "center"), "pairwise coevolution",
     "fly and flower, predator and prey", "the interface sets the outcome", 1),
    (sys_star, (385, 368, 0.62), (385, 322, "center"), "a species and its guild",
     "one fly, its flowers", "one link sits among many", 2),
    (sys_bip, (545, 455, 0.70), (545, 401, "center"), "between two guilds",
     "plants and their pollinators", "the web around the link", 3),
    (sys_uni, (236, 436, 0.85), (236, 390, "center"), "within one guild",
     "competitors for a resource", "strong competitors select strongly", 4),
    (sys_host_microbiome, (975, 392, 0.85), (975, 330, "center"), "host and microbiome",
     "one host, a community", "genes and microbes build the trait", 5),
    (sys_inherit, (852, 150, 0.9), (884, 158, "left"), "host and inherited symbiont",
     "passed from parent to offspring", "what selection reaches is what responds", 6),
    (sys_host_microbe, (748, 287, 0.85), (756, 236, "center"), "host and microbe",
     "host and parasite, one symbiont", "a microbe can buffer, or speed adaptation", 7),
    (sys_hosts_microbiomes, (752, 455, 0.70), (752, 401, "center"),
     "hosts and microbiomes", "many hosts, shared microbes",
     "microbes have a geography of their own", 8),
]


# a soft white plate behind each system and its label (x0, y0, x1, y1)
VEIL = {1: (104, 204, 326, 324), 2: (292, 290, 480, 408), 3: (458, 366, 634, 496),
        4: (134, 358, 340, 488), 5: (874, 298, 1078, 454), 6: (818, 102, 1150, 198),
        7: (644, 202, 868, 324), 8: (646, 366, 858, 496)}


def _veil(cv, box, strength=0.36):
    """Three nested translucent plates: whitest in the middle, fading out, so
    the chart stays visible but recedes behind the glyph and its words."""
    x0, y0, x1, y1 = box
    for grow in (18, 9, 0):
        cv.add_patch(FancyBboxPatch((x0 - grow, y0 - grow), x1 - x0 + 2 * grow,
                                    y1 - y0 + 2 * grow,
                                    boxstyle=f"round,pad=0,rounding_size={16 + grow}",
                                    facecolor=PAPER, edgecolor="none", alpha=strength,
                                    zorder=1.8))


BADGE = {1: (136, 318), 2: (385, 414), 3: (468, 458), 4: (150, 452),
         5: (878, 418), 6: (822, 172), 7: (686, 320), 8: (656, 458)}


def map_B(mode="map"):
    """The chart.  Across: same scale -> nested.  Up: more partners (one guild,
    then two).  Down: a more intimate partnership, kept across generations.

    The compass ring, spokes, rose and fray are the chart's background; the
    systems are laid over it and are free to cross the ring or sit outside it.

    mode = "map"    slide 3: each system with a name and an example
           "close"  slide 31: the stops of the talk numbered, each with what
                    it taught in place of the example
           "title"  slide 1: no words, veiled, to sit behind the title
    """
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    c, R = (590, 287), 205
    rng = np.random.default_rng(7)
    words = mode != "title"

    # ---- the chart (all of it below the glyphs) -----------------------------
    n = 3400
    rr = R * (1.0 + rng.exponential(0.26, n))
    th = rng.uniform(0, 2 * np.pi, n)
    sx, sy = c[0] + rr * np.cos(th), c[1] + rr * np.sin(th)
    ok = (sy > 30) & (sy < 556) & (sx > 6) & (sx < 1174)
    cv.scatter(sx[ok], sy[ok], s=1.6, color=MUTED, linewidths=0, zorder=1.0,
               alpha=0.75)
    for b_ in range(0, 360, 30):
        p0, p1 = _polar(c, 0.17 * R, b_), _polar(c, R, b_)
        cv.plot([p0[0], p1[0]], [p0[1], p1[1]], color=LINE, lw=lw(1.0), zorder=1.1)
        q0, q1 = _polar(c, R, b_), _polar(c, 1.16 * R, b_)
        if b_ % 90:
            cv.plot([q0[0], q1[0]], [q0[1], q1[1]], color=FAINT, lw=lw(0.9),
                    ls=(0, (1.5, 2.5)), zorder=1.1)
    cv.add_patch(Circle(c, 0.5 * R, facecolor="none", edgecolor=FAINT,
                        linewidth=lw(1.0), linestyle=(0, (4, 3)), zorder=1.1))
    cv.add_patch(Circle(c, R, facecolor="none", edgecolor=INK, linewidth=lw(2.6),
                        zorder=1.5))
    cv.add_patch(Circle(c, R - 9, facecolor="none", edgecolor=INK,
                        linewidth=lw(0.9), zorder=1.5))
    for b_ in range(0, 360, 10):
        p0, p1 = _polar(c, R - (15 if b_ % 30 == 0 else 9), b_), _polar(c, R, b_)
        cv.plot([p0[0], p1[0]], [p0[1], p1[1]], color=INK, lw=lw(0.9), zorder=1.5)
    for p0, p1 in [((c[0] - R, c[1]), (92, c[1])), ((c[0] + R, c[1]), (1088, c[1])),
                   ((c[0], c[1] + R), (c[0], 520)), ((c[0], c[1] - R), (c[0], 58))]:
        cv.plot([p0[0], p1[0]], [p0[1], p1[1]], color=INK, lw=lw(1.1),
                ls=(0, (5, 3)), zorder=1.4)
    _rose(cv, c, 28)

    if words:
        sub = dict(fontsize=fs(FS_TINY), style="italic", color=MUTED, zorder=12,
                   path_effects=HALO, **SERIF)
        _token(cv, c[0], 536, "DIFFUSE")
        cv.text(c[0] + 62, 536, "many partners", ha="left", va="center", **sub)
        _token(cv, c[0], 42, "INTIMATE")
        cv.text(c[0] + 66, 42, "one partner, often inherited", ha="left",
                va="center", **sub)
        _token(cv, 50, c[1], "SAME SCALE").set_rotation(90)
        cv.text(77, c[1], "species with species", ha="center", va="center",
                rotation=90, **sub)
        _token(cv, 1130, c[1], "NESTED").set_rotation(270)
        cv.text(1103, c[1], "a host with its microbes", ha="center", va="center",
                rotation=270, **sub)

    for glyph, (gx, gy, u), (lx, ly, ha), name, example, lesson, stop in ENTRIES:
        if words:
            x0, y0, x1, y1 = VEIL[stop]
            _veil(cv, (x0, y0 - (24 if mode == "close" and stop == 1 else 0), x1, y1))
        glyph(pen, gx, gy, u)
        if not words:
            continue
        _label(cv, lx, ly, name, lesson if mode == "close" else example, ha)
        if mode == "close":
            bx, by = BADGE[stop]
            cv.add_patch(Circle((bx, by), 11, facecolor=NAVY, edgecolor=PAPER,
                                linewidth=lw(1.4), zorder=13))
            t = cv.text(bx, by - 0.5, str(stop), ha="center", va="center",
                        fontsize=fs(FS_TINY), color=PAPER, fontweight="bold",
                        zorder=14)
            t._allow_overlap = True
    if mode == "close":
        cv.text(215, 198, "across space: who is ahead depends on distance",
                ha="center", va="baseline", fontsize=fs(FS_TINY), style="italic",
                color=MUTED, zorder=12, path_effects=HALO, **SERIF)
    if mode == "title":
        # a light overall veil; the deck's CSS whitens further just under the words
        cv.add_patch(Rectangle((0, 0), 1180, 560, facecolor=PAPER, edgecolor="none",
                               alpha=0.42, zorder=30))
        save(fig, "s01_title_bg")
        return
    status(cv, "schematic", "position records kind, not measurement")
    save(fig, "s31_closing_map" if mode == "close" else "s03_mapB")


def main():
    map_A()
    map_A(route=True)
    map_B()
    map_B("close")
    map_B("title")


if __name__ == "__main__":
    main()
