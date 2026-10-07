"""
Slide 24 -- which microbial variation can contribute across generations?

Ancestral concordance, redrawn from Week et al. (2025) Evolution
79:2487-2502, Fig. 1, on the lineage grid of rev15/figure2_results.py (C).
Rows are host generations (the present at the bottom); the band is the focal
host's own ancestry; the dashed line is the ancestry of one of its microbes.

This slide only says where ancestry runs.  Whether selection on hosts reaches
that variation is the next slide (`source_pool`).
"""
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

import qgmmt_sim as qg
from glyphs import Pen, arrow, host, microbe
from variance_blocks import variance_block
from semstyle import *  # noqa: F403

HOST = [2, 2, 1, 1, 2, 3, 3]           # focal host's ancestors, present -> past
TRACKS = {
    "lineal": HOST,
    "non-lineal": [2, 4, 3, 0, 1, 4, 3],
    "novel": [2, 5, 5, 5, 5, 5, 5],    # column 5 = outside the host population
}
SAYS = {
    "lineal": ("follows the host's own ancestors", "yes"),
    "non-lineal": ("among hosts, but not this host's ancestors",
                   "only if acquired after selection"),
    "novel": ("from outside the host population", "no"),
}
BAND = "#C9D0D6"


def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    sx, sy, env_gap = 36, 43, 14
    y0 = 150
    ys = [y0 + g * sy for g in range(len(HOST))]

    # key
    ky = 514
    cv.plot([150, 206], [ky, ky], color=BAND, lw=lw(12), solid_capstyle="round")
    cv.add_patch(Circle((178, ky), 7, fc=MUTED, ec="white", lw=lw(1.2), zorder=3))
    cv.text(222, ky, "the focal host's ancestry", va="center", fontsize=fs(FS_SMALL))
    cv.plot([486, 542], [ky, ky], color=TEAL, lw=lw(2.6), ls=(0, (2.2, 1.2)))
    cv.add_patch(Rectangle((514 - 6.5, ky - 6.5), 13, 13, fc=TEAL, ec="none", zorder=3))
    cv.text(558, ky, "the ancestry of one of its microbes", va="center",
            fontsize=fs(FS_SMALL))
    cv.add_patch(Circle((922, ky), 13, fc="none", ec=INK, lw=lw(2.0)))
    cv.text(944, ky, "focal host, today", va="center", fontsize=fs(FS_SMALL))

    for gi, (name, track) in enumerate(TRACKS.items()):
        gx = 196 + gi * 340
        xs = [gx + i * sx for i in range(5)] + [gx + 5 * sx + env_gap]
        cv.add_patch(Rectangle((xs[5] - 15, ys[0] - 20), 30, ys[-1] - ys[0] + 40,
                               fc=SOFT, ec=LINE, lw=lw(1.2), zorder=0))
        cv.text(xs[5], ys[-1] + 28, "outside", ha="center", va="bottom",
                fontsize=fs(FS_TINY), color=MUTED)
        for y in ys:
            for x in xs[:5]:
                cv.add_patch(Circle((x, y), 7, fc="white", ec=FAINT, lw=lw(1.2),
                                    zorder=1))
        hx = [xs[i] for i in HOST]
        cv.plot(hx, ys, color=BAND, lw=lw(13), zorder=2, solid_capstyle="round",
                solid_joinstyle="round")
        for x, y in zip(hx, ys):
            cv.add_patch(Circle((x, y), 8, fc=MUTED, ec="white", lw=lw(1.2), zorder=3))
        mx = [xs[i] for i in track]
        cv.plot(mx, ys, color=TEAL, lw=lw(2.6), ls=(0, (2.2, 1.2)), zorder=4)
        for x, y in zip(mx, ys):
            cv.add_patch(Rectangle((x - 6, y - 6), 12, 12, fc=TEAL, ec="white",
                                   lw=lw(1.0), zorder=5))
        cv.add_patch(Circle((xs[HOST[0]], ys[0]), 18, fc="none", ec=INK, lw=lw(2.0),
                            zorder=6))
        xc = (xs[0] + xs[5]) / 2
        cv.text(xc, 92, name, ha="center", va="baseline", fontsize=fs(FS_TITLE),
                fontweight="bold")
        cv.text(xc, 64, SAYS[name][0], ha="center", va="baseline",
                fontsize=fs(FS_TINY), color=MUTED)
    arrow(pen, (118, ys[0]), (118, ys[-1]), px=1.4, color=MUTED, head=9)
    cv.text(96, (ys[0] + ys[-1]) / 2, "generations back", rotation=90, ha="center",
            va="center", fontsize=fs(FS_SMALL), color=MUTED)
    status(cv, "published", "Week et al. 2025, Evolution 79:2487 \u2014 Fig. 1 redrawn; three patterns of ancestry, not three modes of transmission")
    save(fig, "s24_ancestral_concordance")


# ============================================================ simulations ===
STYLE = {"G": dict(color=INK, ls=(0, (1, 1.6)), lw=LW_DATA, name="genes only"),
         "GL": dict(color=TEAL, ls="-", lw=LW_DATA + 0.6, name="+ lineal microbes"),
         "GLN": dict(color=OCHRE, ls=(0, (4, 1.8)), lw=LW_DATA + 0.6, name="+ non-lineal"),
         "GLNV": dict(color=MUTED, ls=(0, (6, 1.5, 1.2, 1.5)), lw=LW_DATA - 0.6, name="+ novel")}


def _fig3_data():
    out = {}
    for k, arch in enumerate(qg.ARCH):
        for ns in (True, False):
            out[f"{arch}_{'post' if ns else 'pre'}"] = qg.timeseries(
                arch, ns, T=10, reps=20, seed=100 + 2 * k + ns)
    return out


def response():
    """Fig. 3 of Week et al. (2025): the response to selection on hosts when
    the trait is built from more and more kinds of factor, with non-lineal
    microbes taken from hosts after or before selection."""
    D = cached("qgmmt_fig3", _fig3_data)
    fig = figure_px(FULL)
    cv = canvas(fig)
    gen = np.arange(1, 11)
    ymax = 1.08 * max(D[k].mean(0).max() for k in D)
    for x0, tag, title, line in [
            (112, "post", "Non-lineal microbes from hosts that passed selection",
             "they add to the response; novel microbes do not"),
            (690, "pre", "Non-lineal microbes from hosts before selection",
             "now they add nothing; only genes and lineal microbes respond")]:
        cv.text(x0, 500, title, ha="left", va="baseline", fontsize=fs(FS_LABEL),
                fontweight="bold")
        cv.text(x0, 474, line, ha="left", va="baseline", fontsize=fs(FS_TINY),
                color=MUTED)
        ax = px_axes(fig, x0, 104, 440, 340)
        ends = {}
        for arch in ("GLNV", "GLN", "GL", "G"):
            y = D[f"{arch}_{tag}"]
            m, sd = y.mean(0), y.std(0, ddof=1)
            st = STYLE[arch]
            ax.fill_between(gen, m - sd, m + sd, color=st["color"], alpha=0.10, lw=0)
            ax.plot(gen, m, color=st["color"], lw=lw(st["lw"]), ls=st["ls"],
                    zorder=3 if arch != "GLNV" else 4)
            ends[arch] = m[-1]
        # label the line ends, merging lines that finish together
        groups, used = [], set()
        for arch in ("G", "GL", "GLN", "GLNV"):
            if arch in used:
                continue
            g = [b for b in ("G", "GL", "GLN", "GLNV")
                 if b not in used and abs(ends[b] - ends[arch]) < 0.09 * ymax]
            used |= set(g)
            groups.append(g)
        for g in groups:
            yv = np.mean([ends[b] for b in g])
            lab = ", ".join(STYLE[b]["name"] for b in g)
            ax.text(9.75, yv + 0.035 * ymax, lab, ha="right", va="bottom",
                    fontsize=fs(FS_TINY), color=STYLE[g[0]]["color"] if len(g) == 1 else INK,
                    fontweight="bold")
        ax.set_xlim(1, 10)
        ax.set_ylim(0, ymax)
        ax.set_xticks([1, 4, 7, 10])
        ax.set_yticks([0, 50, 100])
        ax.set_xlabel("host generations")
        if tag == "post":
            ax.set_ylabel("change in mean trait")
    # the published figure's ordering, checked on every build
    e = {k: D[k].mean(0)[-1] for k in D}
    assert e["GLN_post"] > 1.5 * e["GL_post"] > 1.5 * e["G_post"]
    assert abs(e["GLN_pre"] - e["GL_pre"]) < 0.2 * e["GL_pre"]
    assert abs(e["GLNV_post"] - e["GLN_post"]) < 0.2 * e["GLN_post"]
    status(cv, "published",
           "Week et al. 2025, Evolution 79:2487 \u2014 Fig. 3, re-run from a Python port of the repository's simulation (20 replicates, mean \u00b1 s.d.)")
    save(fig, "s24x_response")


def _mini(pen, x, y, carrier, state="on", r=16, bold=False):
    host(pen, x, y, r, state=state, edge_px=2.8 if bold else 1.6)
    if carrier:
        microbe(pen, x, y, 5.6, taxon=1, z=8, state="on" if state == "on" else "muted")


def _mini_lineage(cv, pen, x0, y0, kind, dx=30, dy=30):
    """A small version of the ancestry grid: five host generations, the focal
    host's line of descent, and the ancestry of one of its microbes."""
    host_path = [(1, 0), (1, 1), (2, 2), (2, 3), (1, 4)]
    paths = {"lineal": host_path,
             "non-lineal": [(1, 0), (3, 1), (0, 2), (3, 3), (2, 4)],
             "novel": [(1, 0), (4.5, 1), (4.5, 2), (4.5, 3), (4.5, 4)]}
    pos = lambda c, r: (x0 + c * dx, y0 + r * dy)
    cv.add_patch(Rectangle((x0 + 4.5 * dx - 11, y0 - 12), 22, 4 * dy + 24, facecolor=SOFT,
                           edgecolor=LINE, linewidth=lw(1.0), zorder=1))
    for r in range(5):
        for c in range(4):
            cv.add_patch(Circle(pos(c, r), 5.5, facecolor=PAPER, edgecolor=FAINT,
                                linewidth=lw(1.0), zorder=2))
    hp = np.array([pos(*p) for p in host_path])
    cv.plot(hp[:, 0], hp[:, 1], color=LINE, lw=lw(9), solid_capstyle="round", zorder=3)
    cv.plot(hp[:, 0], hp[:, 1], ls="none", marker="o", ms=7, mfc=MUTED, mec="none", zorder=4)
    mp = np.array([pos(*p) for p in paths[kind]])
    cv.plot(mp[:, 0], mp[:, 1], color=TEAL, lw=lw(1.8), ls=(0, (2.4, 1.8)), zorder=5)
    cv.plot(mp[:, 0], mp[:, 1], ls="none", marker="s", ms=5.5, mfc=TEAL, mec="none", zorder=6)
    cv.add_patch(Circle(pos(1, 0), 11, facecolor="none", edgecolor=INK, linewidth=lw(1.6),
                        zorder=7))


def expectation():
    """What we expect before any model: selection acts on hosts, so microbial
    variation whose ancestry runs with the host lineage should respond, novel
    microbes cannot, and for non-lineal microbes it is not obvious."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    cv.text(590, 516, "Selection acts on hosts. Which microbial variation can respond?",
            ha="center", va="center", fontsize=fs(FS_LABEL), fontweight="bold")
    cards = [(206, "lineal", "follows the selected hosts'\nown line of descent", "yes",
              "expected to respond", TEAL),
             (590, "non-lineal", "among hosts, but not\nthese hosts' ancestors", "?",
              "it is not obvious", INDIGO),
             (974, "novel", "from outside the\nhost population", "no",
              "cannot, by definition", MUTED)]
    for xc, name, where, mark, verdict, col in cards:
        hot = mark == "?"
        cv.add_patch(FancyBboxPatch((xc - 170, 76), 340, 396,
                                    boxstyle="round,pad=0,rounding_size=14",
                                    facecolor=INDIGO_T if hot else PAPER,
                                    edgecolor=INDIGO if hot else LINE,
                                    linewidth=lw(2.6 if hot else 1.4), zorder=0))
        cv.text(xc, 436, name, ha="center", va="center", fontsize=fs(FS_TITLE),
                fontweight="bold")
        cv.text(xc, 398, where, ha="center", va="center", fontsize=fs(FS_TINY), color=MUTED,
                linespacing=1.2)
        _mini_lineage(cv, pen, xc - 76, 232, name)
        if mark == "?":
            cv.text(xc, 164, "?", ha="center", va="center", fontsize=fs(46), color=INDIGO,
                    fontweight="bold")
        else:
            cv.add_patch(Circle((xc, 164), 26, facecolor=col if mark == "yes" else PAPER,
                                edgecolor=col, linewidth=lw(2.6), zorder=3))
            if mark == "yes":
                cv.plot([xc - 11, xc - 3, xc + 12], [164, 155, 175], color=PAPER, lw=lw(4.2),
                        solid_capstyle="round", solid_joinstyle="round", zorder=4)
            else:
                for sgn in (1, -1):
                    cv.plot([xc - 10, xc + 10], [164 - 10 * sgn, 164 + 10 * sgn], color=col,
                            lw=lw(4.2), solid_capstyle="round", zorder=4)
        cv.text(xc, 112, verdict, ha="center", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", color=INK if mark != "?" else INDIGO)
    cv.text(590, 44, "grey: the hosts' line of descent; dotted: the ancestry of one of their microbes; shaded column: outside the host population",
            ha="center", va="center", fontsize=fs(FS_TINY), color=MUTED)
    status(cv, "published", "Week et al. 2025, Evolution 79:2487 \u2014 the three ancestry classes and what one would expect of each")
    save(fig, "s24a_expectation")


def source_pool():
    """For a non-lineal microbe the answer turns on one thing: is it passed to
    the next generation before selection acts on hosts, or after?  Before, its
    variation carries no signal of selection; after, it does."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    cv.text(40, 520, "A non-lineal microbe: is it passed on before selection on hosts, or after?",
            ha="left", va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    cols = {"parents": 250, "sel": 520, "off": 790}
    for key, head in (("parents", "parents"), ("sel", "after selection on hosts"),
                      ("off", "offspring")):
        cv.text(cols[key], 484, head, ha="center", va="center", fontsize=fs(FS_SMALL),
                color=MUTED)
    grid = [(-44, 22), (0, 22), (44, 22), (-44, -22), (0, -22), (44, -22)]
    carriers = {0, 3, 4}
    rows = [(338, "before", False, {0, 2, 5}, "no signal of selection", "no response"),
            (160, "after", True, {0, 1, 2, 3, 4, 5}, "carries the signal", "a response")]
    for y, when, filtered, kid_carriers, out1, out2 in rows:
        col = TEAL if filtered else INK
        cv.text(40, y + 14, "passed on", ha="left", va="center", fontsize=fs(FS_SMALL),
                color=MUTED)
        cv.text(40, y - 14, when, ha="left", va="center", fontsize=fs(FS_TITLE),
                fontweight="bold", color=col)
        for n, (dx, dy) in enumerate(grid):
            _mini(pen, cols["parents"] + dx, y + dy, n in carriers)
            keep = n in carriers                      # hosts with the microbe do better
            _mini(pen, cols["sel"] + dx, y + dy, n in carriers,
                  state="on" if keep else "ghost", bold=keep)
            _mini(pen, cols["off"] + dx, y + dy, n in kid_carriers)
        arrow(pen, (cols["parents"] + 76, y), (cols["sel"] - 76, y), px=2.2, color=MUTED,
              head=10)
        arrow(pen, (cols["sel"] + 76, y - 10), (cols["off"] - 76, y - 10), px=2.2,
              color=MUTED, head=10)
        cv.text((cols["sel"] + cols["off"]) / 2, y - 28, "hosts", ha="center", va="center",
                fontsize=fs(FS_TINY), color=MUTED)
        src = cols["sel"] if filtered else cols["parents"]
        arrow(pen, (src, y + 50), (cols["off"], y + 50), px=2.8, color=TEAL,
              rad=-0.24 if filtered else -0.15, head=11, ls=(0, (3.2, 2.2)), z=9)
        cv.text((src + cols["off"]) / 2, y + 50 + (50 if filtered else 62), "microbes",
                ha="center", va="center", fontsize=fs(FS_SMALL), color=TEAL,
                fontweight="bold")
        cv.text(900, y + 14, out1, ha="left", va="center", fontsize=fs(FS_SMALL), color=MUTED)
        cv.text(900, y - 14, out2, ha="left", va="center", fontsize=fs(FS_LABEL),
                fontweight="bold", color=col)
    cv.plot([40, 1140], [258, 258], color=LINE, lw=lw(1.2))
    cv.text(40, 62, "The same species of microbe can do either: what matters is whether its variation has been through selection.",
            ha="left", va="center", fontsize=fs(FS_SMALL))
    status(cv, "published", "Week et al. 2025, Evolution 79:2487 \u2014 schematic of the result of their Figs 3 and 5 (re-runs in the backup)")
    save(fig, "s24b_source_pool")


def bridge():
    """From description to mechanism.  Left: what the framework so far is, a
    statistical description of host variation in terms of genes and microbes,
    and a source of hypotheses.  Right: what is needed next, a mechanistic
    model simple enough to give clear first answers."""
    from variance_blocks import variance_block as vb
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(2)
    # ---- a description -----------------------------------------------------------------
    cv.add_patch(FancyBboxPatch((40, 96), 470, 404, boxstyle="round,pad=0,rounding_size=14",
                                facecolor=PAPER, edgecolor=LINE, linewidth=lw(1.6), zorder=0))
    cv.text(64, 464, "So far: a statistical description", ha="left", va="center",
            fontsize=fs(FS_LABEL), fontweight="bold")
    x = 64
    for sym, f_, kind in ((r"$G_A$", 0.40, "g"), (r"$M_A$", 0.34, "m"), (r"$C_A$", 0.26, "c")):
        vb(cv, x, 380, 420 * f_, 44, kind)
        cv.text(x + 420 * f_ / 2, 362, sym, ha="center", va="center", fontsize=fs(FS_SMALL))
        x += 420 * f_
    cv.text(64, 320, "host variation, in terms of genes and microbes", ha="left",
            va="center", fontsize=fs(FS_SMALL))
    for k, name in enumerate(("lineal", "non-lineal", "novel")):
        _mini_lineage(cv, pen, 84 + k * 146, 150, name, dx=19, dy=19)
    cv.text(64, 264, "and which of it could respond to selection", ha="left", va="center",
            fontsize=fs(FS_SMALL))
    cv.text(64, 118, "a way to describe, and a source of hypotheses", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    # ---- the step -------------------------------------------------------------------------
    arrow(pen, (530, 300), (650, 300), px=3.0, color=NAVY, head=13)
    cv.text(590, 330, "how does it\nactually happen?", ha="center", va="bottom",
            fontsize=fs(FS_SMALL), color=NAVY, fontweight="bold", linespacing=1.15)
    # ---- a mechanism ----------------------------------------------------------------------
    cv.add_patch(FancyBboxPatch((670, 96), 470, 404, boxstyle="round,pad=0,rounding_size=14",
                                facecolor=PAPER, edgecolor=NAVY, linewidth=lw(2.4), zorder=0))
    cv.text(694, 464, "Next: a simple mechanistic model", ha="left", va="center",
            fontsize=fs(FS_LABEL), fontweight="bold", color=NAVY)
    pts = [(760, 380), (880, 404), (1010, 372), (1080, 300), (800, 280), (930, 300),
           (1020, 226), (860, 196), (740, 196)]
    has_g = [1, 0, 1, 0, 0, 1, 1, 0, 1]
    has_m = [1, 1, 0, 0, 0, 1, 0, 0, 0]
    for (hx, hy), g_, m_ in zip(pts, has_g, has_m):
        host(pen, hx, hy, 27, edge_px=2.0)
        _chip(cv, hx - 10, hy + 2, "g", TEAL, on=bool(g_))
        if m_:
            _chip(cv, hx + 12, hy + 2, "m", OCHRE, on=True)
    for i, j in ((0, 4), (1, 2), (5, 3)):
        arrow(pen, pts[i], pts[j], px=2.2, color=OCHRE, head=9, shrink=(31, 31),
              ls=(0, (3.2, 2.2)), z=6)
    cv.text(694, 142, "one host allele, one microbe, a population in trouble", ha="left",
            va="center", fontsize=fs(FS_SMALL))
    cv.text(694, 118, "simple enough to give clear first answers", ha="left", va="center",
            fontsize=fs(FS_TINY), color=MUTED)
    status(cv, "schematic", "from Week et al. 2025 (Evolution) to the rescue model (in prep.)")
    save(fig, "s24d_bridge")


def _chip(cv, x, y, letter, col, on=True, s=1.0):
    """A small labelled chip: [g] the host allele, [m] the microbe."""
    cv.add_patch(FancyBboxPatch((x - 11 * s, y - 9 * s), 22 * s, 18 * s,
                                boxstyle=f"round,pad=0,rounding_size={4.5 * s}",
                                facecolor=PAPER if on else SOFT, edgecolor=col if on else FAINT,
                                linewidth=lw(1.8 if on else 1.1), zorder=9))
    t = cv.text(x, y, f"${letter}$" if on else "\u2013", ha="center", va="center",
                fontsize=fs(FS_TINY * s), color=col if on else FAINT, zorder=10)
    t._allow_overlap = True


def _fig5_data():
    out = {}
    for ns in (False, True):
        obs, pred = qg.prediction_runs(ns, nr=60, sr=6, orr=6, seed=7 + ns)
        tag = "post" if ns else "pre"
        out[f"obs_{tag}"] = obs
        for k in pred:
            out[f"pred_{k}_{tag}"] = pred[k]
    return out


def prediction():
    """Fig. 5 of Week et al. (2025): the one-generation response observed in
    the simulation against the response predicted from additive variance,
    counting more and more kinds of factor as transmissible."""
    D = cached("qgmmt_fig5", _fig5_data)
    fig = figure_px(FULL)
    cv = canvas(fig)
    lo = min(v.min() for v in D.values()) - 0.5
    hi = max(v.max() for v in D.values()) + 0.5
    heads = {"G": "genes", "GL": "+ lineal", "GLN": "+ non-lineal", "GLNV": "+ novel"}
    cv.text(250, 524, "Response predicted when these count as inherited:",
            ha="left", va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")
    for r, (tag, rowlab) in enumerate([("pre", "non-lineal microbes\nfrom hosts\nbefore selection"),
                                       ("post", "non-lineal microbes\nfrom hosts that\npassed selection")]):
        y0 = 296 - r * 216
        cv.text(24, y0 + 90, rowlab, ha="left", va="center", fontsize=fs(FS_TINY),
                color=MUTED, linespacing=1.2)
        err = {k: np.mean(np.abs(D[f"pred_{k}_{tag}"] - D[f"obs_{tag}"])) for k in heads}
        best = min(err, key=err.get)
        for c, k in enumerate(heads):
            x0 = 250 + c * 228
            ax = px_axes(fig, x0, y0, 190, 180)
            ax.plot([lo, hi], [lo, hi], color=FAINT, lw=lw(1.4), zorder=1)
            ax.plot(D[f"obs_{tag}"], D[f"pred_{k}_{tag}"], ls="none", marker="o", ms=5,
                    mfc=TEAL if k == best else MUTED, mec=PAPER, mew=lw(0.8), zorder=3)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_xticks([5, 10])
            ax.set_yticks([5, 10])
            if r == 0:
                ax.set_xticklabels([])
                cv.text(x0 + 95, 492, heads[k], ha="center", va="baseline",
                        fontsize=fs(FS_SMALL), fontweight="bold")
            else:
                ax.set_xlabel("observed", fontsize=fs(FS_SMALL))
            if c == 0:
                ax.set_ylabel("predicted", fontsize=fs(FS_SMALL))
            else:
                ax.set_yticklabels([])
            if k == best:
                for sp_ in ax.spines.values():
                    sp_.set_visible(True)
                    sp_.set_color(TEAL)
                    sp_.set_linewidth(lw(2.4))
                ax.text(0.95, 0.07, "matches", transform=ax.transAxes, ha="right",
                        va="bottom", fontsize=fs(FS_TINY), color=TEAL, fontweight="bold")
        want = "GL" if tag == "pre" else "GLN"
        assert best == want, (tag, best)
    status(cv, "published",
           "Week et al. 2025, Evolution 79:2487 \u2014 Fig. 5, re-run from a Python port (60 parameter draws; the paper used 100)")
    save(fig, "s24c_prediction")


def main_all():
    main()
    expectation()
    source_pool()
    bridge()
    response()
    prediction()


if __name__ == "__main__":
    main_all()
