"""
Slide 5 -- the route of the talk, drawn as a route.

Six stops on one line: the three things learned from coevolution (interface,
community context, space), a pit stop where they are put together for
inference, then host-microbiome evolution, where the same three return with
inheritance and timescales, and last the prospects for inference from spatial
molecular data.  Each stop reuses the glyph that stands for it elsewhere in
the deck.  Part I is indigo, Part II teal, the pit stop navy.
"""
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

from glyphs import (Pen, arrow, host, microbe, microbe_positions, species, sys_bip,
                    sys_hosts_microbiomes, sys_pair, sys_space)
from semstyle import *  # noqa: F403

XS = [104, 298, 492, 686, 882, 1076]
Y0 = 250                         # the route


def _wave(x):
    return Y0 + 10 * np.sin((x - 60) / 1080 * 2 * np.pi * 1.5)


def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(3)
    stops = [("1", "the interface", "what happens when two\nindividuals meet", INDIGO),
             ("2", "community context", "the others each interaction\nis embedded in", INDIGO),
             ("3", "space", "the patterns left\nacross a landscape", INDIGO),
             ("", "pit stop: inference", "putting the pieces together:\nprocess from pattern", NAVY),
             ("4", "hosts and microbiomes", "the same lessons, plus\ninheritance and timescales", TEAL),
             ("5", "prospects", "inference from spatial\nmolecular data", TEAL)]
    # the route: one dashed line, coloured by part
    xx = np.linspace(40, 1140, 500)
    for lo, hi, col in ((40, 589, INDIGO), (589, 784, NAVY), (784, 1140, TEAL)):
        m = (xx >= lo) & (xx <= hi)
        cv.plot(xx[m], _wave(xx[m]), color=col, lw=lw(2.6), ls=(0, (4, 2.6)), zorder=2)
    arrow(pen, (1130, _wave(1130)), (1150, _wave(1150)), px=2.6, color=TEAL, head=12)
    # the two parts
    for x0, x1, lab, col in ((40, 589, "Part I  \u00b7  coevolution", INDIGO),
                             (784, 1140, "Part II  \u00b7  hosts and their microbiomes", TEAL)):
        cv.plot([x0, x1], [516, 516], color=col, lw=lw(3.0), solid_capstyle="butt")
        cv.text(x0, 534, lab, ha="left", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", color=col)
    for k, (x, (num, head, line, col)) in enumerate(zip(XS, stops)):
        y = _wave(x)
        # the glyph on a soft plate above the route
        cv.add_patch(FancyBboxPatch((x - 82, 322), 164, 150,
                                    boxstyle="round,pad=0,rounding_size=14", facecolor=PAPER,
                                    edgecolor=LINE, linewidth=lw(1.4), zorder=1))
        cv.plot([x, x], [y + 16, 322], color=col, lw=lw(1.4), zorder=1)
        gy = 398
        if k == 0:
            sys_pair(pen, x, gy, 0.95)
        elif k == 1:
            sys_bip(pen, x, gy, 0.66)
        elif k == 2:
            sys_space(pen, x, gy - 4, 0.86)
        elif k == 3:                                   # data read backwards to process
            pts = rng.multivariate_normal([0, 0], [[420, 300], [300, 420]], 14)
            cv.plot(x + 22 + pts[:, 0] * 0.9, gy + 6 + pts[:, 1] * 0.9, ls="none", marker="o",
                    ms=4.4, mfc=INK, mec=PAPER, mew=lw(0.6), zorder=5)
            arrow(pen, (x + 6, gy - 44), (x - 56, gy - 44), px=2.4, color=NAVY, head=10)
            cv.text(x - 60, gy + 30, "process", ha="left", va="center", fontsize=fs(FS_TINY),
                    color=NAVY)
            cv.text(x - 60, gy + 8, "?", ha="left", va="center", fontsize=fs(FS_TITLE),
                    color=NAVY, fontweight="bold")
        elif k == 4:
            sys_hosts_microbiomes(pen, x, gy + 4, 0.56)
        else:                                          # genomes and microbiomes across places
            for i, (dx, dy) in enumerate([(-44, 26), (10, 40), (40, -8), (-18, -30)]):
                px_, py_ = x + dx, gy + dy
                for j in range(6):
                    on = rng.uniform() < 0.5
                    cv.add_patch(Rectangle((px_ - 20 + j * 6.6, py_ + 14), 5.4, 8,
                                           facecolor=TEAL if on else PAPER, edgecolor=TEAL,
                                           linewidth=lw(0.9), zorder=6))
                host(pen, px_, py_, 11, edge_px=1.3)
                for (mx, my), t in zip(microbe_positions(px_, py_, 11, 3, seed=i),
                                       rng.integers(0, 3, 3)):
                    microbe(pen, mx, my, 2.5, int(t), z=8)
        # the stop itself
        if num:
            cv.add_patch(Circle((x, y), 17, facecolor=col, edgecolor=PAPER, linewidth=lw(2.4),
                                zorder=6))
            t = cv.text(x, y, num, ha="center", va="center", fontsize=fs(FS_SMALL),
                        color=PAPER, fontweight="bold", zorder=7)
        else:
            cv.add_patch(Polygon([(x, y + 21), (x + 21, y), (x, y - 21), (x - 21, y)],
                                 closed=True, facecolor=PAPER, edgecolor=col,
                                 linewidth=lw(2.8), zorder=6))
        cv.text(x, 204, head, ha="center", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", color=col)
        cv.text(x, 184, line, ha="center", va="top", fontsize=fs(FS_TINY), color=MUTED,
                linespacing=1.2)
    cv.text(590, 78, "Three lessons from coevolution, carried into host\u2013microbiome evolution, with inference as the thread.",
            ha="center", va="center", fontsize=fs(FS_SMALL), color=INK)
    status(cv, "schematic", "the route of the talk")
    save(fig, "s05_roadmap")


STOPS = [("1", "interface", INDIGO), ("2", "community", INDIGO), ("3", "space", INDIGO),
         ("", "inference", NAVY), ("4", "hosts and microbiomes", TEAL), ("5", "prospects", TEAL)]


def _strip(cv, pen, y=84, active=(), done=()):
    """The route, small, along the bottom of a divider slide: where we are."""
    xs = np.linspace(150, 1030, 6)
    cv.plot([xs[0] - 60, xs[-1] + 60], [y, y], color=FAINT, lw=lw(2.0), ls=(0, (4, 2.6)),
            zorder=1)
    for k, (x, (num, name, col)) in enumerate(zip(xs, STOPS)):
        on, was = k in active, k in done
        c = col if (on or was) else FAINT
        if num:
            cv.add_patch(Circle((x, y), 15 if on else 11, facecolor=c if on else PAPER,
                                edgecolor=c, linewidth=lw(2.4), zorder=3))
            if on:
                t = cv.text(x, y, num, ha="center", va="center", fontsize=fs(FS_TINY),
                            color=PAPER, fontweight="bold", zorder=4)
                t._allow_overlap = True
            elif was:
                cv.plot([x - 5, x - 1.5, x + 6], [y, y - 4.5, y + 5], color=c, lw=lw(2.4),
                        solid_capstyle="round", solid_joinstyle="round", zorder=4)
        else:
            s_ = 19 if on else 13
            cv.add_patch(Polygon([(x, y + s_), (x + s_, y), (x, y - s_), (x - s_, y)],
                                 closed=True, facecolor=c if on else PAPER, edgecolor=c,
                                 linewidth=lw(2.6), zorder=3))
        cv.text(x, y - 34, name, ha="center", va="center", fontsize=fs(FS_TINY),
                color=INK if on else (MUTED if was else FAINT),
                fontweight="bold" if on else "normal")


def part_one():
    """A pause: Part I, coevolution.  Mostly empty on purpose."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    cv.text(120, 372, "Part I", ha="left", va="center", fontsize=fs(30), color=MUTED)
    cv.text(116, 296, "Coevolution", ha="left", va="center", fontsize=fs(76),
            fontweight="bold", color=INDIGO)
    sys_pair(pen, 930, 330, 2.3)
    _strip(cv, pen, active=(0, 1, 2))
    status(cv, "schematic", "the route of the talk: the next three stops")
    save(fig, "s05c_part1")


def pit_stop():
    """A pause before inference: three things were shown to shape what
    coevolution does and the patterns it leaves; can the spatial patterns be
    read back to the process?"""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(3)
    cv.text(80, 500, "So far", ha="left", va="center", fontsize=fs(FS_LABEL), color=MUTED)
    for x, name, draw in [(170, "interface", lambda: sys_pair(pen, 170, 392, 0.95)),
                          (360, "community", lambda: sys_bip(pen, 360, 392, 0.66)),
                          (550, "space", lambda: sys_space(pen, 550, 388, 0.86))]:
        cv.add_patch(FancyBboxPatch((x - 82, 322), 164, 140,
                                    boxstyle="round,pad=0,rounding_size=14", facecolor=PAPER,
                                    edgecolor=LINE, linewidth=lw(1.4), zorder=1))
        draw()
        cv.text(x, 300, name, ha="center", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", color=INDIGO)
    arrow(pen, (660, 392), (800, 392), px=3.0, color=INDIGO, head=13)
    cv.text(730, 414, "shape", ha="center", va="center", fontsize=fs(FS_TINY), color=INDIGO)
    ax = px_axes(fig, 840, 330, 190, 130)
    xy = rng.multivariate_normal([0, 0], [[1, 0.75], [0.75, 1]], 26)
    ax.plot(xy[:, 0], xy[:, 1], ls="none", marker="o", ms=6.5, mfc=INK, mec=PAPER,
            mew=lw(1.0))
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    cv.text(935, 300, "outcomes, and the patterns they leave", ha="center", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold")
    arrow(pen, (935, 276), (360, 276), px=3.2, color=NAVY, rad=-0.14, head=14,
          ls=(0, (3.4, 2.4)))
    cv.text(648, 196, "Now the other way: from spatial patterns back to process?",
            ha="center", va="center", fontsize=fs(FS_TITLE), fontweight="bold", color=NAVY)
    _strip(cv, pen, active=(3,), done=(0, 1, 2))
    status(cv, "schematic", "the route of the talk: a pit stop")
    save(fig, "s17b_pitstop")


def main_all():
    main()
    part_one()
    pit_stop()


if __name__ == "__main__":
    main_all()
