"""
Slides 02, 21, 23, 28, 30 -- first drafts, kept deliberately simple.

s02_question   the two questions of the talk in one picture
s20b_generative  a light bridge from ABC to the three lessons: once the model
               is simulated rather than solved, almost anything could go in;
               which ingredients does theory say matter for interacting
               species?  The slide leaves three slots empty for the next one.
s21_lessons    the three lessons of the coevolution half, each with a
               thumbnail of the result it came from (read from figures/, so
               this script must run after those are built)
s23_qgmmt      what builds a microbiome-mediated trait, and what selection on
               hosts can act on (framework of Week et al. 2025, Evolution
               79:2487; drawn without values)
s28_ecology    hosts as patches in a microbial metacommunity (ongoing work:
               schematic only)
s30_sbi        the planned inference from joint host-genomic and microbiome
               data (planned work: schematic only, no outcomes drawn)
"""
import matplotlib.image as mpimg
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Polygon, Rectangle

from glyphs import (Pen, arc, arrow, host, landscape, link, microbe,
                    microbe_positions, on_tile, species, sys_bip, sys_host_microbiome,
                    sys_pair, sys_space)
from variance_blocks import variance_block
from semstyle import *  # noqa: F403


def panel(cv, x0, y0, w, h, title, sub="", dashed=False):
    cv.add_patch(FancyBboxPatch((x0, y0), w, h,
                                boxstyle="round,pad=0,rounding_size=14",
                                facecolor=PAPER, edgecolor=FAINT if dashed else LINE,
                                linewidth=lw(1.6),
                                linestyle=(0, (4, 3)) if dashed else "-", zorder=0))
    cv.text(x0 + 16, y0 + h - 34, title, ha="left", va="baseline",
            fontsize=fs(FS_SMALL), fontweight="bold")
    if sub:
        cv.text(x0 + 16, y0 + h - 58, sub, ha="left", va="baseline",
                fontsize=fs(FS_TINY), color=MUTED)


# =================================================================== s02 ===
def _badge(cv, x, y, num, col, r=19):
    cv.add_patch(Circle((x, y), r, facecolor=col, edgecolor="none", zorder=5))
    t = cv.text(x, y, num, ha="center", va="center", fontsize=fs(FS_LABEL), color=PAPER,
                fontweight="bold", zorder=6)
    t._allow_overlap = True


def question():
    """The two questions, each sitting on its own arrow: the first runs forward
    over three small scenes (interactions, evolution, the variation we
    observe), the second runs back underneath them."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(3)
    cards = [(70, "interactions between species"), (450, "evolution"),
             (830, "the variation we observe")]
    W, Y0, H = 280, 184, 190
    for x0, lab in cards:
        cv.add_patch(FancyBboxPatch((x0, Y0), W, H, boxstyle="round,pad=0,rounding_size=16",
                                    facecolor=PAPER, edgecolor=LINE, linewidth=lw(1.6),
                                    zorder=1))
        cv.text(x0 + W / 2, Y0 + 24, lab, ha="center", va="center", fontsize=fs(FS_SMALL),
                color=MUTED)
    # scene 1: interactions, of more than one kind
    sys_pair(pen, 210, 316, 1.25)
    sys_bip(pen, 130, 250, 0.34)
    sys_host_microbiome(pen, 292, 250, 0.36)
    # scene 2: traits shift
    ax = px_axes(fig, 478, 236, 224, 116)
    x = np.linspace(-4, 4, 300)
    for mu0, mu1, col, ls in [(-1.6, -0.5, INDIGO, "-"), (1.7, 0.9, OCHRE, (0, (3.2, 1.6)))]:
        ax.plot(x, np.exp(-0.5 * ((x - mu0) / 0.7) ** 2), color=FAINT, lw=lw(1.6), ls=ls)
        ax.plot(x, np.exp(-0.5 * ((x - mu1) / 0.7) ** 2), color=col, lw=lw(LW_DATA), ls=ls)
        ax.annotate("", xy=(mu1, 1.14), xytext=(mu0, 1.14),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=lw(1.8)))
    ax.set_xlim(-4, 4)
    ax.set_ylim(0, 1.3)
    blank(ax)
    ax.spines["bottom"].set_visible(True)
    # scene 3: variation, among populations and across places
    bx = px_axes(fig, 852, 232, 118, 124)
    xy = rng.multivariate_normal([0, 0], [[1, 0.75], [0.75, 1]], 26)
    bx.plot(xy[:, 0], xy[:, 1], ls="none", marker="o", ms=5.6, mfc=INK, mec=PAPER,
            mew=lw(0.9))
    bx.set_xlim(-3, 3)
    bx.set_ylim(-3, 3)
    blank(bx)
    bx.spines["bottom"].set_visible(True)
    bx.spines["left"].set_visible(True)
    sys_space(pen, 1036, 296, 0.80)
    # between the scenes
    for xa in (350, 730):
        arrow(pen, (xa + 18, Y0 + H / 2 + 14), (xa + 82, Y0 + H / 2 + 14), px=3.0,
              color=INDIGO, head=13)
        arrow(pen, (xa + 82, Y0 + H / 2 - 14), (xa + 18, Y0 + H / 2 - 14), px=2.4,
              color=NAVY, head=11, ls=(0, (3.2, 2.2)))
    # question 1 runs forward, above
    _badge(cv, 90, 486, "1", INDIGO)
    cv.text(122, 486, "What are the evolutionary consequences of interspecific interactions?",
            ha="left", va="center", fontsize=fs(FS_TITLE), fontweight="bold")
    arrow(pen, (70, 428), (1110, 428), px=3.4, color=INDIGO, head=15)
    # question 2 runs back, below
    arrow(pen, (1110, 132), (70, 132), px=3.0, color=NAVY, head=15, ls=(0, (3.4, 2.4)))
    _badge(cv, 90, 74, "2", NAVY)
    cv.text(122, 74, "How can we recognize these consequences from observed variation?",
            ha="left", va="center", fontsize=fs(FS_TITLE), fontweight="bold")
    status(cv, "schematic", "the two questions of the talk")
    save(fig, "s02_question")


# ================================================================== s20b ===
def generative():
    """A light bridge.  Simulation frees inference from the likelihood, and then
    almost anything could go into the simulator.  For interacting species,
    theory can say which ingredients matter.  The slide does not name them:
    the three slots are empty, and the next slide fills them."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(5)
    cv.text(40, 520, "If we simulate, a lot could go into the simulator.", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")

    # ---- the loop of the previous slide --------------------------------------------
    nodes = {"draw": (170, 392), "simulate": (256, 248), "keep if close": (84, 248)}
    order = ["draw", "simulate", "keep if close"]
    for k, name in enumerate(order):
        p, q = nodes[name], nodes[order[(k + 1) % 3]]
        arrow(pen, p, q, px=2.2, color=NAVY, rad=-0.28, head=10, shrink=(44, 44))
    for name, (x, y) in nodes.items():
        cv.add_patch(Circle((x, y), 36, facecolor=PAPER, edgecolor=NAVY if name == "simulate" else LINE,
                            linewidth=lw(2.4 if name == "simulate" else 1.6), zorder=3))
    x, y = nodes["draw"]
    xs = np.linspace(-22, 22, 40)
    cv.plot(x + xs, y - 12 + 26 * np.exp(-0.5 * ((xs + 4) / 9) ** 2), color=INDIGO,
            lw=lw(2.2), zorder=4)
    cv.plot(x + xs, y - 12 + 20 * np.exp(-0.5 * ((xs - 6) / 11) ** 2), color=OCHRE,
            lw=lw(2.2), ls=(0, (3.2, 1.6)), zorder=4)
    x, y = nodes["simulate"]
    pts = rng.multivariate_normal([0, 0], [[90, 60], [60, 90]], 12)
    cv.plot(x + pts[:, 0], y + pts[:, 1], ls="none", marker="o", ms=4.2, mfc=MUTED,
            mec=PAPER, mew=lw(0.5), zorder=4)
    x, y = nodes["keep if close"]
    cv.add_patch(Circle((x, y), 15, facecolor=SOFT, edgecolor=INK, linewidth=lw(1.2),
                        linestyle=(0, (3, 2)), zorder=4))
    cv.plot([x], [y], marker="*", ms=12, mfc=PAPER, mec=INK, mew=lw(1.2), zorder=5)
    for dx, dy, inside in [(-7, 6, True), (8, -4, True), (22, 14, False), (-20, -16, False),
                           (4, 24, False)]:
        cv.plot([x + dx], [y + dy], marker="o", ms=4, mfc=INK if inside else PAPER,
                mec=PAPER if inside else FAINT, mew=lw(0.9), zorder=5)
    cv.text(170, 442, "draw", ha="center", va="center", fontsize=fs(FS_SMALL))
    cv.text(256, 198, "simulate", ha="center", va="center", fontsize=fs(FS_SMALL),
            fontweight="bold", color=NAVY)
    cv.text(84, 198, "keep if close", ha="center", va="center", fontsize=fs(FS_SMALL))

    # ---- everything one could put in ------------------------------------------------------
    words = ["demography", "mating system", "life history", "genetic architecture",
             "plasticity", "mutation", "phenology", "environmental change", "behaviour",
             "development", "physiology", "sampling design"]
    spots = [(392, 430), (540, 452), (690, 424), (402, 372), (606, 386), (760, 364),
             (430, 310), (566, 322), (744, 296), (396, 252), (540, 262), (680, 236)]
    for w_, (x, y) in zip(words, spots):
        cv.text(x, y, w_, ha="left", va="center", fontsize=fs(FS_SMALL), color=FAINT if (x + y) % 3 else MUTED)
    cv.text(840, 452, "\u2026", ha="left", va="center", fontsize=fs(FS_LABEL), color=FAINT)
    arrow(pen, (296, 262), (372, 300), px=2.2, color=NAVY, head=10, ls=(0, (3, 2.2)))
    arrow(pen, (868, 330), (934, 330), px=2.6, color=NAVY, head=11)

    # ---- the simulator, with room for a few things only -----------------------------------
    bx0, bx1 = 950, 1146
    cv.add_patch(FancyBboxPatch((bx0, 196), bx1 - bx0, 268,
                                boxstyle="round,pad=0,rounding_size=14", facecolor=PAPER,
                                edgecolor=NAVY, linewidth=lw(2.4), zorder=1))
    cv.text((bx0 + bx1) / 2, 488, "the simulator", ha="center", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold", color=NAVY)
    for y in (410, 330, 250):
        cv.add_patch(FancyBboxPatch((bx0 + 22, y - 28), bx1 - bx0 - 44, 56,
                                    boxstyle="round,pad=0,rounding_size=10", facecolor=SOFT,
                                    edgecolor=FAINT, linewidth=lw(1.2), linestyle=(0, (3, 2.4)),
                                    zorder=2))
        cv.text((bx0 + bx1) / 2, y, "?", ha="center", va="center", fontsize=fs(FS_TITLE),
                color=NAVY, fontweight="bold")
    cv.plot([40, 1140], [136, 136], color=LINE, lw=lw(1.2))
    cv.text(590, 92, "For interacting species, what does theory say we cannot leave out?",
            ha="center", va="center", fontsize=fs(FS_LABEL), fontweight="bold")
    status(cv, "schematic", "a question; the next slide answers it")
    save(fig, "s20b_generative")


# =================================================================== s21 ===
def lessons():
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    tiles = [("interface", "how the partners meet\nsets the outcome",
              lambda x, y: sys_pair(pen, x, y, 0.5), "s08_arms_race"),
             ("community", "strong competitors select strongly,\nin a direction competition does not give",
              lambda x, y: sys_bip(pen, x, y, 0.42), "s11c_circle"),
             ("space", "who is locally adapted\ndepends on the distance",
              lambda x, y: sys_space(pen, x + 4, y, 0.5), "s17_scale")]
    for k, (head, line, icon, thumb) in enumerate(tiles):
        x0 = 30 + k * 380
        icon(x0 + 44, 492)
        cv.text(x0 + 100, 484, head, ha="left", va="baseline", fontsize=fs(FS_TITLE),
                fontweight="bold", color=INDIGO)
        cv.text(x0, 440, line, ha="left", va="top", fontsize=fs(FS_SMALL),
                linespacing=1.25)
        f = FIG_DIR / f"{thumb}.png"
        w, h = 360, 360 * 560 / 1180
        if f.exists():
            cv.imshow(mpimg.imread(f), extent=[x0, x0 + w, 196, 196 + h], zorder=3,
                      aspect="auto", interpolation="lanczos")
        cv.add_patch(Rectangle((x0, 196), w, h, facecolor="none", edgecolor=LINE,
                               linewidth=lw(1.4), zorder=4))
    cv.set_xlim(0, fig._px[0])
    cv.set_ylim(0, fig._px[1])
    cv.plot([30, 1150], [150, 150], color=LINE, lw=lw(1.2))
    cv.text(590, 108, "Three things decide what coevolution leaves behind.",
            ha="center", va="center", fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(590, 74, "Next: a system where all three act at once.", ha="center",
            va="center", fontsize=fs(FS_SMALL), color=MUTED)
    status(cv, "schematic", "thumbnails: slides 8, 11 and 17")
    save(fig, "s21_lessons")


# =================================================================== s23 ===
def g_box(cv, x, y, s=1.0):
    cv.add_patch(FancyBboxPatch((x - 15 * s, y - 12 * s), 30 * s, 24 * s,
                                boxstyle=f"round,pad=0,rounding_size={6 * s}",
                                facecolor=PAPER, edgecolor=TEAL, linewidth=lw(2.0), zorder=8))
    t = cv.text(x, y, r"$g$", ha="center", va="center", fontsize=fs(max(FS_SMALL * s, FS_TINY)),
                color=TEAL, zorder=9)
    t._allow_overlap = True


def cartoon_sum(cv, pen, x, y, lead="host trait  =", k=1.0):
    """host trait = [g] + microbes + [g] x microbes + ...  (a cartoon, not a model)"""
    if lead:
        cv.text(x, y, lead, ha="left", va="center", fontsize=fs(FS_LABEL), fontweight="bold")
        x += 12.2 * len(lead) + 18
    else:
        x += 16 * k
    g_box(cv, x, y, k)
    plus = dict(ha="center", va="center", fontsize=fs(max(FS_LABEL * k, FS_SMALL)), color=MUTED)
    cv.text(x + 36 * k, y, "+", **plus)
    for i in range(3):
        microbe(pen, x + (70 + 26 * i) * k, y, 7.5 * k, i, z=8)
    cv.text(x + 158 * k, y, "+", **plus)
    g_box(cv, x + 194 * k, y, k)
    cv.text(x + 222 * k, y, "\u00d7", **plus)
    microbe(pen, x + 248 * k, y, 7.5 * k, 2, z=8)
    cv.text(x + 286 * k, y, "+ \u2026", **plus)
    return x + 310 * k


def qgmmt():
    """The trait can be built any way at all; what is measured is host genomes,
    microbiomes and the trait, host by host; fitting them together partitions
    the trait's variance into G_A, M_A and their covariance C_A."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(9)
    cartoon_sum(cv, pen, 56, 510)
    cv.text(640, 510, "however the trait is built,", ha="left", va="center",
            fontsize=fs(FS_SMALL), color=MUTED)

    # ---- what we measure, host by host ------------------------------------------
    cv.text(56, 440, "What we measure, host by host", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    from matplotlib.colors import LinearSegmentedColormap
    y0, hh = 170, 230
    for kx, (lab, tint, col, wcol) in enumerate([("host loci", TEAL_T, TEAL, 150),
                                                 ("microbial taxa", OCHRE_T, OCHRE, 150)]):
        ax = px_axes(fig, 92 + kx * 172, y0, wcol, hh)
        M = rng.integers(0, 3, (12, 6)) if kx == 0 else rng.gamma(0.8, 1.0, (12, 6))
        ax.imshow(M, cmap=LinearSegmentedColormap.from_list("m", [PAPER, tint, col]),
                  aspect="auto", interpolation="nearest")
        blank(ax)
        cv.text(92 + kx * 172 + wcol / 2, y0 - 20, lab, ha="center", va="center",
                fontsize=fs(FS_TINY), color=col, fontweight="bold")
    ax = px_axes(fig, 440, y0, 64, hh)
    ax.barh(np.arange(12), rng.uniform(0.25, 1.0, 12)[::-1], height=0.62, color=INK)
    ax.set_ylim(-0.5, 11.5)
    blank(ax)
    cv.text(472, y0 - 20, "trait", ha="center", va="center", fontsize=fs(FS_TINY),
            fontweight="bold")
    cv.text(72, y0 + hh / 2, "one row per host", ha="center", va="center", rotation=90,
            fontsize=fs(FS_TINY), color=MUTED)
    arrow(pen, (530, y0 + hh / 2), (606, y0 + hh / 2), px=2.6, color=NAVY, head=11)
    cv.text(568, y0 + hh / 2 + 20, "fit\ntogether", ha="center", va="bottom",
            fontsize=fs(FS_TINY), color=NAVY, linespacing=1.15)

    # ---- the partition ----------------------------------------------------------------
    X, W, by, bh = 640, 480, 300, 62
    cv.text(X, 440, "we can partition its variance", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    x = X
    parts = [(r"$G_A$", 0.40, "g", "host genes", TEAL),
             (r"$M_A$", 0.34, "m", "microbes", OCHRE),
             (r"$C_A$", 0.26, "c", "genes and microbes\noccurring together", INK)]
    for sym, f_, kind, words, col in parts:
        variance_block(cv, x, by, W * f_, bh, kind)
        cv.text(x + W * f_ / 2, by + bh + 22, sym, ha="center", va="center",
                fontsize=fs(FS_LABEL))
        cv.text(x + W * f_ / 2, by - 12, words, ha="center", va="top", fontsize=fs(FS_TINY),
                color=col if col != INK else MUTED, linespacing=1.2,
                fontweight="bold" if col != INK else "normal")
        x += W * f_
    cv.text(X + W * 0.87, 214, "a covariance,\nnot an interaction", ha="center", va="top",
            fontsize=fs(FS_TINY), color=INK, fontweight="bold", linespacing=1.2)
    cv.plot([56, 1124], [116, 116], color=LINE, lw=lw(1.2))
    cv.text(590, 78, "But explaining trait variation now does not tell us which of it is inherited.",
            ha="center", va="center", fontsize=fs(FS_LABEL), fontweight="bold")
    status(cv, "published", "Week et al. 2025, Evolution 79:2487 \u2014 the framework, drawn without values")
    save(fig, "s23_qgmmt")


# =================================================================== s28 ===
def social():
    """Social contact is microbial dispersal: hosts are the patches, contacts
    are the dispersal events, and one contact can move several taxa at once.
    Framework in development with Aura Raulo; schematic, no results."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(6)
    H = [(96, 330), (212, 440), (236, 262), (366, 356), (488, 448), (504, 270)]
    E = [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (3, 5), (4, 5)]
    for i, j in E:
        link(pen, H[i], H[j], px=2.2, color=MUTED, z=2)
    for k, (hx, hy) in enumerate(H):
        host(pen, hx, hy, 44)
        ptsm = microbe_positions(hx, hy, 44, 7, seed=30 + k)
        tx = rng.choice(3, size=len(ptsm), p=[0.5, 0.25, 0.25] if k < 3 else [0.2, 0.3, 0.5])
        for (px_, py_), t in zip(ptsm, tx):
            microbe(pen, px_, py_, 4.6, int(t), z=7)
    # one contact, several taxa together
    p, q = np.array(H[2]), np.array(H[3])
    mid, u = (p + q) / 2, (q - p) / np.hypot(*(q - p))
    nrm = np.array([-u[1], u[0]])
    cv.add_patch(Ellipse(mid + nrm * 2, 74, 34, angle=np.degrees(np.arctan2(u[1], u[0])),
                         facecolor=TEAL_T, edgecolor=TEAL, linewidth=lw(1.8), zorder=5))
    for k in range(3):
        c_ = mid + u * (k - 1) * 20 + nrm * 2
        microbe(pen, c_[0], c_[1], 6.0, k, z=8)
    arrow(pen, tuple(mid + u * 44 + nrm * 2), tuple(mid + u * 70 + nrm * 2), px=2.4,
          color=TEAL, head=10, z=6)
    cv.text(292, 204, "one contact,\nseveral taxa at once", ha="left",
            va="center", fontsize=fs(FS_SMALL), color=TEAL, fontweight="bold",
            linespacing=1.15)
    cv.plot([mid[0] - 4, 330], [mid[1] - 20, 224], color=TEAL, lw=lw(1.2), zorder=4)
    cv.text(40, 520, "Hosts are the patches; contacts are the dispersal events", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(40, 496, "a social network is a microbial metacommunity", ha="left",
            va="baseline", fontsize=fs(FS_TINY), color=MUTED)

    # ---- the consequence: taxa do not disperse independently --------------------------
    X = 640
    cv.plot([X - 40, X - 40], [50, 500], color=LINE, lw=lw(1.2))
    cv.text(X - 14, 520, "So taxa do not disperse independently", ha="left",
            va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    n_tax, T, rate = 8, 60, 0.11
    ind = rng.random((n_tax, T)) < rate
    events = rng.random(T) < rate / 0.75
    clu = events[None, :] & (rng.random((n_tax, T)) < 0.75)
    for k, (M, col, head, line) in enumerate([
            (ind, INDIGO, "if each taxon moved on its own clock", "independent dispersal"),
            (clu, TEAL, "when one contact moves many at once", "the same average rate, clustered")]):
        y0 = 300 - k * 212
        ax = px_axes(fig, X + 6, y0, 500, 120)
        if k == 1:
            for t_ in np.where(events)[0]:
                ax.axvspan(t_ - 0.5, t_ + 0.5, color=TEAL_T, lw=0, zorder=1)
        yy, xx = np.where(M)
        ax.plot(xx, yy, ls="none", marker="s", ms=5.6, mfc=col, mec="none", zorder=3)
        ax.set_xlim(-1, T)
        ax.set_ylim(-0.8, n_tax - 0.2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("time", fontsize=fs(FS_TINY))
        ax.set_ylabel("taxa", fontsize=fs(FS_TINY))
        cv.text(X + 6, y0 + 150, head, ha="left", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", color=col)
        cv.text(X + 6, y0 + 130, line, ha="left", va="center", fontsize=fs(FS_TINY),
                color=MUTED)
    status(cv, "schematic", "the social-microbiome metacommunity framework, in development with Aura Raulo; illustration, no results")
    save(fig, "s28a_social")


def geography():
    """Microbial geography exists before host selection.

    The background model in its first, most neutral form (after the spatial
    summary): no host movement and no effect of the microbe on hosts.  In each
    place hosts gain a microbe from the local environment at a rate that
    differs among places, and lose it at one rate everywhere, so prevalence
    settles at gain / (gain + loss), place by place.

    Consequence: if host allele frequencies also differ among places, for any
    reason, then allele and microbe are associated in the pooled sample
    although they are independent within every place and no interaction was
    assumed.  With the QGMMT components of one locus (frequency p, effect 1)
    and one microbe (prevalence q, effect 1), and D = Cov(g, m):

        M_A = q(1 - q)  = mean_a[q_a (1 - q_a)] + Var_a(q_a)
        C_A = 2 D       = 2 Cov_a(p_a, q_a)       (independence within places)

    which are the laws of total variance and covariance.  The bars use the five
    places drawn on the left.  This is a sketch for planned work, not a result
    of the papers; values are illustrative.
    """
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    gain = np.array([0.2, 0.5, 1.0, 2.0, 4.0])
    loss = 1.0
    q = gain / (gain + loss)
    p = np.array([0.15, 0.30, 0.50, 0.70, 0.85])            # host allele frequency, for other reasons
    nA = len(q)
    n_host = 6

    # ---- the model, place by place ------------------------------------------------------
    cv.text(40, 520, "Gain and loss, place by place", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(40, 494, "no host movement, no effect on hosts: just the local environment",
            ha="left", va="baseline", fontsize=fs(FS_TINY), color=MUTED)
    quad = landscape(pen, 40, 96, 470, 300, skew=0.20, n_lines=3, seed=5, px=1.6)
    uv = [(0.13, 0.80), (0.30, 0.50), (0.50, 0.76), (0.68, 0.44), (0.87, 0.72)]
    P = [on_tile(quad, *c_) for c_ in uv]
    for a_, pt in enumerate(P):
        k = int(round(q[a_] * n_host))
        cv.add_patch(Ellipse((pt[0], pt[1] - 34), 96, 26, facecolor=PAPER, edgecolor=FAINT,
                             linewidth=lw(1.2), zorder=3))
        for j in range(int(3 + 2.2 * gain[a_])):                 # more of it in richer places
            ang, rr = 2.4 * j + a_, 0.25 + 0.7 * ((j * 37 + a_ * 11) % 10) / 10
            microbe(pen, pt[0] + 40 * rr * np.cos(ang), pt[1] - 34 + 9 * rr * np.sin(ang), 2.8,
                    1, state="muted", z=4)
        arrow(pen, (pt[0] - 26, pt[1] - 20), (pt[0] - 26, pt[1] + 2), px=1.2 + 1.1 * gain[a_] ** 0.6,
              color=TEAL, head=6 + 1.5 * gain[a_] ** 0.5, z=5)
        arrow(pen, (pt[0] + 26, pt[1] + 2), (pt[0] + 26, pt[1] - 20), px=1.6, color=INK,
              head=7, z=5)
        for j in range(n_host):
            hx, hy = pt[0] - 30 + 12 * j, pt[1] + 18 + (5 if j % 2 else -5)
            host(pen, hx, hy, 7.5, edge_px=1.2)
            if j < k:
                microbe(pen, hx, hy, 3.2, 1, z=8)
    lx, ly = 56, 70
    arrow(pen, (lx, ly - 9), (lx, ly + 9), px=2.4, color=TEAL, head=8)
    cv.text(lx + 14, ly, "gain from the local environment (differs by place)", ha="left",
            va="center", fontsize=fs(FS_TINY), color=TEAL)
    arrow(pen, (lx, ly - 16), (lx, ly - 34), px=1.6, color=INK, head=7)
    cv.text(lx + 14, ly - 25, "loss (the same everywhere)", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INK)
    cv.text(40, 448, "so the share of hosts carrying it settles at", ha="left",
            va="center", fontsize=fs(FS_SMALL))
    cv.text(40, 424, "gain \u00f7 (gain + loss),  a different value in each place", ha="left",
            va="center", fontsize=fs(FS_SMALL), fontweight="bold", color=TEAL)

    # ---- a spurious association ----------------------------------------------------------
    X = 640
    cv.plot([X - 44, X - 44], [40, 540], color=LINE, lw=lw(1.2))
    cv.text(X - 16, 520, "An association nobody selected for", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    ax = px_axes(fig, X + 30, 318, 200, 160)
    ax.plot(p, q, ls="none", marker="o", ms=11, mfc=TEAL, mec=PAPER, mew=lw(1.4), zorder=3)
    c1 = np.polyfit(p, q, 1)
    ax.plot([0.05, 0.95], np.polyval(c1, [0.05, 0.95]), color=MUTED, lw=lw(1.4),
            ls=(0, (3, 2)), zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xlabel("host allele frequency", fontsize=fs(FS_TINY), labelpad=-6)
    ax.set_ylabel("microbe prevalence", fontsize=fs(FS_TINY), labelpad=-6)
    ax.text(0.05, 0.92, "one point per place", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    D_within = 0.0
    D_pooled = float(np.mean(p * q) - p.mean() * q.mean())
    r_pooled = D_pooled / np.sqrt(p.mean() * (1 - p.mean()) * q.mean() * (1 - q.mean()))
    tx = X + 262
    cv.text(tx, 458, "if allele frequencies also differ", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(tx, 438, "among places, for any reason:", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(tx, 398, "within any one place", ha="left", va="center", fontsize=fs(FS_TINY))
    cv.text(tx, 376, "allele and microbe: independent", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold")
    cv.text(tx, 340, "all hosts pooled", ha="left", va="center", fontsize=fs(FS_TINY))
    cv.text(tx, 318, f"associated  (r = {r_pooled:.2f} here)", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold", color=TEAL)

    # ---- what that does to the variance components ------------------------------------------
    cv.plot([X - 16, 1150], [266, 266], color=LINE, lw=lw(1.2))
    cv.text(X - 16, 240, "Through QME and the variance components", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold")
    status_mark(cv, "planned", X + 378, 240, color=NAVY)
    cv.text(X + 392, 240, "a first sketch", ha="left", va="center", fontsize=fs(FS_TINY),
            color=NAVY)
    within = float(np.mean(q * (1 - q)))
    among = float(np.var(q))
    assert abs(within + among - q.mean() * (1 - q.mean())) < 1e-12      # total variance
    scale_px = 900.0                                   # pixels per unit of variance
    y1, y2, bh = 168, 96, 26
    cv.text(X - 16, y1 + bh / 2, r"$M_A$", ha="left", va="center", fontsize=fs(FS_LABEL))
    x0 = X + 44
    variance_block(cv, x0, y1, within * scale_px, bh, "m")
    cv.add_patch(Rectangle((x0 + within * scale_px, y1), among * scale_px, bh,
                           facecolor=OCHRE, edgecolor=INK, linewidth=lw(1.2), zorder=3))
    cv.text(x0 + within * scale_px / 2, y1 + bh + 12, "variation within places", ha="center",
            va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(x0 + within * scale_px + among * scale_px + 10, y1 + bh / 2,
            "+ differences among places", ha="left", va="center", fontsize=fs(FS_TINY),
            color=OCHRE, fontweight="bold")
    cv.text(X - 16, y2 + bh / 2, r"$C_A$", ha="left", va="center", fontsize=fs(FS_LABEL))
    variance_block(cv, x0, y2, 2 * D_pooled * scale_px, bh, "c")
    cv.text(x0 + 2 * D_pooled * scale_px + 10, y2 + bh / 2 + 9,
            "allele frequency and prevalence", ha="left", va="center", fontsize=fs(FS_TINY),
            color=INK, fontweight="bold")
    cv.text(x0 + 2 * D_pooled * scale_px + 10, y2 + bh / 2 - 10,
            "varying together among places", ha="left", va="center", fontsize=fs(FS_TINY),
            color=INK, fontweight="bold")
    cv.text(X - 16, 56, "background, present before any selection on hosts; next: host movement",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    print(f"  [check] geography: q = {np.round(q, 2).tolist()}; M_A within {within:.3f} + among "
          f"{among:.3f}; pooled D = {D_pooled:.3f}, r = {r_pooled:.2f}")
    status(cv, "planned", "the spatial background model, first step: gain and loss only; unit effects, illustrative values")
    save(fig, "s28b_geography")


# ================================================================== s28c ===
def _sheet(cv, z, x0, y0, w, h, skew, cmap, z_=3):
    """A small tilted sheet (see `_layer` in s14_17_space.py)."""
    from matplotlib.transforms import Affine2D
    tr = Affine2D.from_values(w, 0, skew * h, h, x0, y0) + cv.transData
    cv.imshow(z, extent=[0, 1, 0, 1], origin="lower", cmap=cmap, vmin=-2.4, vmax=2.4,
              transform=tr, interpolation="bilinear", aspect="auto", zorder=z_)
    cv.add_patch(Polygon([(x0, y0), (x0 + w, y0), (x0 + w + skew * h, y0 + h), (x0 + skew * h, y0 + h)],
                         closed=True, facecolor="none", edgecolor=INK, linewidth=lw(1.2),
                         zorder=z_ + 1))


def _smooth(n, ell, rng):
    k = 2 * np.pi * np.fft.fftfreq(n, d=1.0 / n)
    kx, ky = np.meshgrid(k, k)
    z = np.real(np.fft.ifft2(np.fft.fft2(rng.normal(size=(n, n))) / (2.0 / ell ** 2 + kx ** 2 + ky ** 2)))
    return (z - z.mean()) / z.std()


def directions():
    """A transition, in three rows, each with pictures.  Top: the ingredients
    gathered in Part II and the question of how to combine them to find the
    signature of host-microbiome interactions.  Middle: the background any such
    attempt has to account for: where hosts live (spurious associations from
    geography) and whom they meet (host social structure, with Aura Raulo and
    Guilhem Sommeria-Klein); together they shape how microbiomes are inherited.
    Bottom: the same question for coevolution, with its own two pieces of
    background: selection by the abiotic environment that varies across the
    landscape (a small copy of the optima slide), and shared dispersal history
    (two species that spread along the same route have allele frequencies that
    covary across places without any selection between them).
    Planned and ongoing work; nothing here is a result."""
    from matplotlib.colors import LinearSegmentedColormap
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(7)
    cv.text(40, 530, "How do we combine these ingredients to find the signature of host\u2013microbiome interactions?",
            ha="left", va="center", fontsize=fs(FS_LABEL), fontweight="bold", color=TEAL)
    # ---- the ingredients ---------------------------------------------------------------------
    y = 462
    x = 56
    for w_, kind in ((62, "g"), (52, "m"), (40, "c")):
        variance_block(cv, x, y - 14, w_, 28, kind)
        x += w_
    cv.text(133, y - 36, "a description of variation", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    gx, gy = 300, y - 22
    hp = [(1, 0), (1, 1), (2, 2), (2, 3)]
    mp = [(1, 0), (3, 1), (0, 2), (3, 3)]
    for r in range(4):
        for c in range(4):
            cv.add_patch(Circle((gx + c * 16, gy + r * 15), 3.8, facecolor=PAPER, edgecolor=FAINT,
                                linewidth=lw(0.9), zorder=2))
    cv.plot([gx + c * 16 for c, r in hp], [gy + r * 15 for c, r in hp], color=LINE, lw=lw(6.5),
            solid_capstyle="round", zorder=3)
    cv.plot([gx + c * 16 for c, r in mp], [gy + r * 15 for c, r in mp], color=TEAL, lw=lw(1.5),
            ls=(0, (2.2, 1.6)), marker="s", ms=3.6, mfc=TEAL, mec="none", zorder=5)
    cv.text(gx + 24, y - 36, "what is inherited", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    for k, (hx, has_g, has_m) in enumerate([(470, 1, 1), (522, 0, 1)]):
        host(pen, hx, y + 2, 21, edge_px=1.8)
        for dx, letter, col, on in ((-8, "g", TEAL, has_g), (9, "m", OCHRE, has_m)):
            cv.add_patch(FancyBboxPatch((hx + dx - 8, y - 5), 16, 14,
                                        boxstyle="round,pad=0,rounding_size=4",
                                        facecolor=PAPER if on else SOFT,
                                        edgecolor=col if on else FAINT, linewidth=lw(1.4),
                                        zorder=9))
            t = cv.text(hx + dx, y + 2, f"${letter}$" if on else "\u2013", ha="center", va="center",
                        fontsize=fs(13), color=col if on else FAINT, zorder=10)
            t._allow_overlap = True
    cv.text(496, y - 36, "a mechanistic model", ha="center", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    cv.add_patch(FancyBboxPatch((636, y - 18), 92, 40, boxstyle="round,pad=0,rounding_size=10",
                                facecolor=PAPER, edgecolor=TEAL, linewidth=lw(2.0), zorder=2))
    cv.text(682, y + 2, "QME", ha="center", va="center", fontsize=fs(FS_SMALL), fontweight="bold")
    cv.text(682, y - 36, "a reduction that scales", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    for xp in (236, 406, 590):
        cv.text(xp, y + 2, "+", ha="center", va="center", fontsize=fs(FS_TITLE), color=FAINT)
    arrow(pen, (770, y + 2), (850, y + 2), px=3.0, color=TEAL, head=13)
    for kx, (tint, col) in enumerate(((TEAL_T, TEAL), (OCHRE_T, OCHRE))):
        ax = px_axes(fig, 880 + kx * 92, y - 24, 82, 52)
        M = rng.integers(0, 3, (6, 8)) if kx == 0 else rng.gamma(0.8, 1.0, (6, 8))
        ax.imshow(M, cmap=LinearSegmentedColormap.from_list("m", [PAPER, tint, col]),
                  aspect="auto", interpolation="nearest")
        blank(ax)
    cv.text(1078, y + 2, "?", ha="center", va="center", fontsize=fs(40), color=TEAL,
            fontweight="bold")
    cv.text(967, y - 36, "their signature in the data", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)

    # ---- the background to account for ----------------------------------------------------------
    cv.plot([40, 1140], [404, 404], color=LINE, lw=lw(1.2))
    cv.text(40, 376, "Not without the background", ha="left", va="center",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(40, 350, "structure that creates associations by itself", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(40, 290, "Together they decide how", ha="left", va="center", fontsize=fs(FS_SMALL))
    cv.text(40, 266, "microbiomes are inherited.", ha="left", va="center",
            fontsize=fs(FS_SMALL))
    gain = np.array([0.25, 0.9, 4.0])
    q = gain / (gain + 1.0)
    quad = landscape(pen, 392, 240, 250, 100, skew=0.22, n_lines=2, seed=5, px=1.3)
    for a_, (u, v) in enumerate([(0.17, 0.68), (0.50, 0.32), (0.83, 0.70)]):
        pt = on_tile(quad, u, v)
        k = int(round(q[a_] * 5))
        for j in range(5):
            hx, hy = pt[0] - 24 + 12 * j, pt[1] + (4 if j % 2 else -4)
            host(pen, hx, hy, 7, edge_px=1.1)
            if j < k:
                microbe(pen, hx, hy, 3.0, 1, z=8)
    cv.text(392, 362, "where hosts live", ha="left", va="center", fontsize=fs(FS_SMALL),
            fontweight="bold", color=TEAL)
    cv.text(392, 224, "geography: spurious associations", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INK)
    H = [(800, 308), (868, 262), (930, 318), (1000, 264), (1064, 312), (1112, 258)]
    E = [(0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (4, 5), (3, 5)]
    for i, j in E:
        link(pen, H[i], H[j], px=1.8, color=MUTED, z=2)
    for k, (hx, hy) in enumerate(H):
        host(pen, hx, hy, 18)
        tx = rng.choice(3, size=4, p=[0.5, 0.25, 0.25] if k < 3 else [0.2, 0.3, 0.5])
        for (px_, py_), t in zip(microbe_positions(hx, hy, 18, 4, seed=50 + k), tx):
            microbe(pen, px_, py_, 3.0, int(t), z=7)
    cv.text(780, 362, "whom they meet: social structure", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold", color=TEAL)
    cv.text(780, 224, "with Aura Raulo and Guilhem Sommeria-Klein", ha="left", va="center",
            fontsize=fs(FS_TINY), color=INK)

    # ---- and the same question for coevolution, with its own two kinds of background --------------
    cv.plot([40, 1140], [200, 200], color=LINE, lw=lw(1.2))
    for k, line in enumerate(["Similarly, for coevolution:", "how do we account for", "background selection and",
                              "shared dispersal history?"]):
        cv.text(40, 164 - 30 * k, line, ha="left", va="center", fontsize=fs(FS_LABEL),
                fontweight="bold", color=NAVY)
    # background selection: a mean-trait sheet pulled toward an optimum sheet that varies in space
    env = LinearSegmentedColormap.from_list("env", ["#FFFFFF", "#C9D1D9", "#4A5560"])
    th = _smooth(80, 0.30, rng)
    zz = 0.7 * th + 0.7 * _smooth(80, 0.18, rng)
    x0 = 390
    _sheet(cv, th, x0, 72, 150, 38, 0.6, env, z_=3)
    _sheet(cv, zz, x0, 134, 150, 38, 0.6, "cividis", z_=6)
    for u in (0.25, 0.5, 0.75):
        xa = x0 + 150 * u + 12
        arrow(pen, (xa, 94), (xa, 150), px=1.3, color=MUTED, head=6, z=5)
    cv.set_xlim(0, fig._px[0])
    cv.set_ylim(0, fig._px[1])
    cv.text(x0 + 88, 52, "background selection", ha="center", va="center",
            fontsize=fs(FS_TINY), color=INK, fontweight="bold")
    # shared dispersal history: two species spread along the same route, so their allele
    # frequencies (shade) change together from place to place, with no selection between them
    cmA = LinearSegmentedColormap.from_list("a", ["#E3E0FA", INDIGO, "#1E1670"])
    cmB = LinearSegmentedColormap.from_list("b", ["#F6E3C2", OCHRE, "#5A3306"])
    xs_, yy = [636, 742, 848], 124
    freq = [0.12, 0.5, 0.9]
    for k, (x_, f_) in enumerate(zip(xs_, freq)):
        cv.add_patch(Ellipse((x_, yy), 62, 46, facecolor=PAPER, edgecolor=FAINT, linewidth=lw(1.2),
                             zorder=2))
        cv.add_patch(Circle((x_ - 12, yy), 9.5, facecolor=cmA(f_), edgecolor=INK,
                            linewidth=lw(1.2), zorder=5))
        cv.add_patch(Rectangle((x_ + 4, yy - 9), 18, 18, facecolor=cmB(f_), edgecolor=INK,
                               linewidth=lw(1.2), zorder=5))
        if k < 2:
            arrow(pen, (x_ + 36, yy + 9), (xs_[k + 1] - 36, yy + 9), px=2.2, color=INDIGO, head=9)
            arrow(pen, (x_ + 36, yy - 9), (xs_[k + 1] - 36, yy - 9), px=2.2, color=OCHRE, head=9)
    cv.text(742, 172, "two species spread along the same route", ha="center", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    cv.text(742, 52, "shared dispersal history", ha="center", va="center", fontsize=fs(FS_TINY),
            color=INK, fontweight="bold")
    arrow(pen, (892, yy), (930, yy), px=2.2, color=NAVY, head=9)
    ax = px_axes(fig, 968, 84, 92, 82)
    ax.plot([0, 1], [0, 1], color=FAINT, lw=lw(1.2), ls=(0, (2, 2)))
    for f_ in freq:
        ax.plot([f_], [f_ + 0.04 * (0.5 - f_)], marker="o", ms=9, mfc=INK, mec=PAPER, mew=lw(1.0))
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, 1.08)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("species 1", fontsize=fs(FS_TINY), color=INDIGO, labelpad=3)
    ax.set_ylabel("species 2", fontsize=fs(FS_TINY), color=OCHRE, labelpad=3)
    cv.text(1072, 150, "they\ncovary", ha="left", va="center", fontsize=fs(FS_TINY),
            color=INK, linespacing=1.15)
    cv.text(1072, 104, "without\ncoevolving", ha="left", va="center", fontsize=fs(FS_TINY),
            color=CRIMSON, fontweight="bold", linespacing=1.15)
    status(cv, "planned", "a question for ongoing and planned work; the panels are schematics, not results")
    save(fig, "s28_directions")


# =================================================================== s30 ===
def _genome(cv, x, y, w=54, col=INDIGO, seed=0, n=7):
    rng = np.random.default_rng(seed)
    for k in range(n):
        cv.add_patch(Rectangle((x + k * w / n, y), w / n - 1.5, 9,
                               facecolor=col if rng.uniform() < 0.5 else PAPER,
                               edgecolor=col, linewidth=lw(1.0), zorder=6))


def two_questions():
    """Two open questions, one per row, mostly as pictures: what is sampled
    across the landscape, the mechanism we cannot see, and a short question.
    No method is named."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(4)
    rows = [dict(y=306, col=INDIGO, q="Are two species shaping\none another?\nWhere in their genomes?",
                 sub="or is it shared environment,\ndispersal, history?"),
            dict(y=82, col=TEAL, q="Is the microbiome changing\nhost evolution?\nThrough what?",
                 sub="or is it microbial ecology,\ngeography, demography?")]
    for r_ in rows:
        y, col = r_["y"], r_["col"]
        quad = landscape(pen, 36, y, 380, 176, skew=0.22, n_lines=2, seed=5, px=1.5)
        for k, (u, v) in enumerate([(0.13, 0.70), (0.37, 0.28), (0.62, 0.72), (0.86, 0.30)]):
            p = on_tile(quad, u, v)
            if col == INDIGO:                       # two species, each with a genome
                species(pen, p[0] - 20, p[1] + 12, 8, "A", z=7, edge_px=1.2)
                species(pen, p[0] + 20, p[1] + 12, 8, "B", z=7, edge_px=1.2)
                _genome(cv, p[0] - 46, p[1] - 12, 44, INDIGO, seed=k)
                _genome(cv, p[0] + 2, p[1] - 12, 44, OCHRE, seed=10 + k)
            else:                                   # hosts with genome and microbiome
                for j in range(2):
                    hx = p[0] - 20 + 40 * j
                    host(pen, hx, p[1] + 4, 14, edge_px=1.4)
                    for (mx, my), t in zip(microbe_positions(hx, p[1] + 4, 14, 3, seed=k + j),
                                           rng.integers(0, 3, 3)):
                        microbe(pen, mx, my, 3.0, int(t), z=8)
                _genome(cv, p[0] - 22, p[1] + 24, 44, TEAL, seed=20 + k)
        cv.text(36, y + 196, "what we can sample, across places", ha="left", va="center",
                fontsize=fs(FS_TINY), color=MUTED)
        arrow(pen, (444, y + 88), (478, y + 88), px=2.6, color=NAVY, head=11)
        cv.add_patch(FancyBboxPatch((494, y + 14), 250, 150,
                                    boxstyle="round,pad=0,rounding_size=14", facecolor=PAPER,
                                    edgecolor=MUTED, linewidth=lw(1.6), linestyle=(0, (4, 3)),
                                    zorder=1))
        if col == INDIGO:                           # a hidden phenotypic interface
            ax = px_axes(fig, 522, y + 58, 194, 84)
            dd = np.linspace(-1.15, 1.15, 200)
            ax.plot(dd, np.exp(-7.5 * (dd - 0.5) ** 2), color=INDIGO, lw=lw(2.6))
            ax.plot(dd, np.exp(-2.5 * (dd + 0.5) ** 2), color=OCHRE, lw=lw(2.6),
                    ls=(0, (3.2, 1.6)))
            ax.set_xlim(-1.15, 1.15)
            ax.set_ylim(0, 1.12)
            blank(ax)
            ax.spines["bottom"].set_visible(True)
            cv.text(619, y + 36, "an interface we cannot see", ha="center", va="center",
                    fontsize=fs(FS_TINY), color=MUTED)
        else:
            cartoon_sum(cv, pen, 500, y + 104, lead="", k=0.74)
            cv.text(619, y + 60, "genes, microbes and their\ninteractions shape host fitness",
                    ha="center", va="center", fontsize=fs(FS_TINY), color=MUTED,
                    linespacing=1.2)
        arrow(pen, (760, y + 88), (794, y + 88), px=2.6, color=NAVY, head=11)
        cv.text(814, y + 110, r_["q"], ha="left", va="center", fontsize=fs(FS_LABEL),
                fontweight="bold", color=col, linespacing=1.2)
        cv.text(814, y + 36, r_["sub"], ha="left", va="center", fontsize=fs(FS_TINY),
                color=MUTED, linespacing=1.2)
    cv.plot([36, 1144], [288, 288], color=LINE, lw=lw(1.2))
    cv.text(590, 46, "What models, data and spatial replication make these mechanisms distinguishable?",
            ha="center", va="center", fontsize=fs(FS_LABEL), fontweight="bold", color=NAVY)
    status(cv, "planned", "two proposed directions; no results")
    save(fig, "s30_two_questions")


def main():
    question()
    generative()
    lessons()
    qgmmt()
    social()
    geography()
    directions()
    two_questions()


if __name__ == "__main__":
    main()
