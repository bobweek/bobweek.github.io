"""
Slides 25-27 -- microbial rescue and quasi-microbial equilibrium (QME).

Every number comes from scripts/qme_rescue_model.py, copied unchanged from
rev15 (the model code of the manuscript in preparation).

s25_rescue_setup   the model: growth of a host = -r + s_g g + s_m m + e g m,
                   where g, m in {0, 1} mark the host allele and the microbe;
                   the microbe spreads between hosts at rate beta and is lost
                   at rate ell, so its prevalence settles at 1 - ell/beta.
s26_rescue_result  the statement's Fig. 3A at slide size.  Scenarios of the
                   manuscript: beta = 1, ell = 0.4, r = 0.04, c = 0.05,
                   s_g = 0.05, p_g(0) = 0.02, N(0) = 0.8; total microbial
                   benefit 0.05, placed either in s_m (direct) or in e
                   (interaction); "microbe removed" is the host-only model.
                   QME-reduced trajectories.  The build asserts the numbers
                   quoted on the figure.
s27_qme            A: the QME manifold of the statement (berkeley_statement
                   figure3, panel A): prevalence against association, traced
                   as the allele goes from absent to fixed (beta = 10, ell = 4,
                   s_g = 0.03, s_m = 0.05, e = 0.10).  Added here: two paths
                   of the full model started off the manifold; the build
                   asserts each is back on it within two generations.
                   B: the analogy with quasi-linkage equilibrium (Fig. 3B).
"""
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from scipy.integrate import solve_ivp

from glyphs import Pen, arrow, host, microbe, tick_timeline
from qme_rescue_model import (D_qme, Params, frequency_rhs, pm_qme,
                              solve_host_only, solve_qme)
from semstyle import *  # noqa: F403
from variance_blocks import variance_block

T_END = 300.0
DASH = {"host": (0, (3.0, 2.2)), "direct": "-", "inter": (0, (7, 1.8))}
COL = {"host": MUTED, "direct": OCHRE, "inter": TEAL}


def rescue_params(sm, e):
    return Params(beta=1.0, ell0=0.4, phi=0.0, sg=0.05, sm=sm, e=e, r=0.04, c=0.05)


def rescue_trajectories():
    scen = {"host": rescue_params(0.05, 0.0), "direct": rescue_params(0.05, 0.0),
            "inter": rescue_params(0.0, 0.05)}
    out = {}
    for k, p in scen.items():
        solver = solve_host_only if k == "host" else solve_qme
        t, z = solver(p, pg0=0.02, N0=0.8, t_end=T_END, n=1500)
        out[k] = (t, z[0], z[1])
    return out


# =================================================================== s25 ===
def _chip(cv, x, y, letter, col, on=True, s=1.0):
    """A small labelled chip: [g] the host's rescue allele, [m] the microbe."""
    cv.add_patch(FancyBboxPatch((x - 13 * s, y - 10 * s), 26 * s, 20 * s,
                                boxstyle=f"round,pad=0,rounding_size={5 * s}",
                                facecolor=PAPER if on else SOFT, edgecolor=col if on else FAINT,
                                linewidth=lw(1.8 if on else 1.2), zorder=9))
    t = cv.text(x, y, f"${letter}$" if on else "\u2013", ha="center", va="center",
                fontsize=fs(max(FS_TINY * s, 13)), color=col if on else FAINT, zorder=10)
    t._allow_overlap = True


def g_chip(cv, x, y, on=True, s=1.0):
    _chip(cv, x, y, "g", TEAL, on, s)


def m_chip(cv, x, y, on=True, s=1.0):
    _chip(cv, x, y, "m", OCHRE, on, s)


def setup():
    """A simple model of microbially assisted rescue, said three ways.
    Left: a population in which the microbe [m] is caught by contact with
    carriers and lost again (an SIS model of carriage) while the host allele
    [g] is inherited (haploid population genetics).  The inset shows the part
    of loss that happens at birth: an offspring receives the microbe with
    probability v and misses it with probability 1 - v; the rest of loss is
    clearance.  Top right: the four kinds of host and the two transitions of a
    host, with their rates (after Fig. 1A of the manuscript).  Bottom right:
    host fitness in the same cartoon as the trait on the QGMMT slide."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    rng = np.random.default_rng(14)
    cv.text(40, 520, "A host allele that is inherited, a microbe that is caught and lost",
            ha="left", va="baseline", fontsize=fs(FS_LABEL), fontweight="bold")

    # ---- the population --------------------------------------------------------------
    cv.add_patch(FancyBboxPatch((40, 96), 560, 392, boxstyle="round,pad=0,rounding_size=14",
                                facecolor=PAPER, edgecolor=LINE, linewidth=lw(1.4), zorder=0))
    pts = []
    tries = 0
    while len(pts) < 9 and tries < 20000:
        tries += 1
        p = np.array([rng.uniform(86, 554), rng.uniform(290, 448)])
        if all(np.hypot(*(p - q)) > 86 for q in pts):
            pts.append(p)
    pts = np.array(pts)
    R = 27
    n = len(pts)
    order = np.argsort(pts[:, 0])
    allele = np.zeros(n, bool)
    allele[order[[0, 3, 4, 7]]] = True
    carrier = np.zeros(n, bool)
    carrier[order[[1, 3, 6, 8]]] = True
    dd = np.hypot(pts[:, None, 0] - pts[None, :, 0], pts[:, None, 1] - pts[None, :, 1])
    gains = []
    for j in np.argsort(np.where(~carrier, dd[:, carrier].min(axis=1), 1e9))[:3]:
        i = int(np.where(carrier)[0][np.argmin(dd[j, carrier])])
        gains.append((i, int(j)))
    loser = int(order[8])
    for k, p in enumerate(pts):
        host(pen, p[0], p[1], R, edge_px=2.0)
        g_chip(cv, p[0] - 11, p[1] + 2, on=bool(allele[k]), s=0.82)
        if carrier[k] and k != loser:
            m_chip(cv, p[0] + 12, p[1] + 2, s=0.82)
    for i, j in gains:
        arrow(pen, pts[i], pts[j], px=2.4, color=OCHRE, head=10, shrink=(R + 4, R + 4),
              ls=(0, (3.2, 2.2)), z=6)
    p = pts[loser]                                         # one host clears it
    arrow(pen, (p[0] + 16, p[1] + 10), (p[0] + 40, p[1] + 36), px=1.8, color=INK, head=8, z=9)
    _chip(cv, p[0] + 52, p[1] + 46, "m", FAINT, on=True, s=0.82)
    # the legend of processes
    cv.plot([62, 96], [214, 214], color=OCHRE, lw=lw(2.4), ls=(0, (3.2, 2.2)))
    cv.plot([96], [214], marker=">", ms=7, mfc=OCHRE, mec=OCHRE)
    cv.text(110, 214, "caught by contact with a carrier", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold", color=OCHRE)
    cv.text(110, 193, r"at rate $\beta\,p_m$  ($p_m$: share carrying it)",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.plot([62, 96], [160, 160], color=INK, lw=lw(2.4))
    cv.plot([96], [160], marker=">", ms=7, mfc=INK, mec=INK)
    cv.text(110, 160, "lost again", ha="left", va="center", fontsize=fs(FS_SMALL),
            fontweight="bold")
    cv.text(110, 139, "cleared, or not passed on at birth", ha="left",
            va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(62, 112, "microbe: an SIS model.   allele: haploid.", ha="left",
            va="center", fontsize=fs(FS_TINY), color=INK)
    # at birth: the microbe is passed on with probability v
    bx, by = 504, 222
    cv.add_patch(FancyBboxPatch((bx - 84, 104), 168, 150,
                                boxstyle="round,pad=0,rounding_size=10", facecolor=SOFT,
                                edgecolor="none", zorder=1))
    cv.text(bx - 72, 240, "at birth", ha="left", va="center", fontsize=fs(FS_TINY),
            color=MUTED)
    host(pen, bx, by - 10, 22, edge_px=1.8)
    g_chip(cv, bx - 9, by - 9, s=0.68)
    m_chip(cv, bx + 10, by - 9, s=0.68)
    for dx, has, lab in ((-44, True, r"$v$"), (44, False, r"$1 - v$")):
        kx, ky = bx + dx, 134
        arrow(pen, (bx + 0.35 * dx, by - 34), (kx, ky + 24), px=1.6, color=MUTED, head=7)
        host(pen, kx, ky, 20, edge_px=1.6)
        g_chip(cv, kx - 8, ky + 1, s=0.62)
        m_chip(cv, kx + 9, ky + 1, on=has, s=0.62)
        cv.text(kx + (-30 if dx < 0 else 34), ky + 34, lab, ha="center", va="center",
                fontsize=fs(FS_TINY))

    # ---- four kinds of host, and their transitions ---------------------------------------
    X = 650
    cv.text(X, 488, "Four kinds of host", ha="left", va="center", fontsize=fs(FS_SMALL),
            fontweight="bold")
    cols, rws = [X + 96, X + 386], [410, 300]
    labels = {(0, 0): r"$0$", (1, 0): r"$s_m$", (0, 1): r"$s_g$", (1, 1): r"$s_g + s_m + e$"}
    for ci, x in enumerate(cols):
        for ri, y in enumerate(rws):
            best = (ci, ri) == (1, 1)
            cv.add_patch(FancyBboxPatch((x - 96, y - 40), 192, 80,
                                        boxstyle="round,pad=0,rounding_size=10",
                                        facecolor=TEAL_T if best else PAPER,
                                        edgecolor=TEAL if best else LINE,
                                        linewidth=lw(2.0 if best else 1.3), zorder=1))
            host(pen, x - 56, y, 28, tint=PAPER if best else None, edge_px=2.0)
            g_chip(cv, x - 68, y + 2, on=ri == 1, s=0.82)
            m_chip(cv, x - 44, y + 2, on=ci == 1, s=0.82)
            cv.text(x + 34, y + 12, "adds", ha="center", va="center", fontsize=fs(FS_TINY),
                    color=MUTED)
            cv.text(x + 34, y - 12, labels[(ci, ri)], ha="center", va="center",
                    fontsize=fs(FS_SMALL))
    xa, xb = cols[0] + 102, cols[1] - 102
    for y in rws:
        arrow(pen, (xa, y + 12), (xb, y + 12), px=2.2, color=OCHRE, head=10,
              ls=(0, (3.2, 2.2)))
        cv.text((xa + xb) / 2, y + 28, r"$\beta\,p_m$", ha="center", va="center",
                fontsize=fs(FS_TINY), color=OCHRE)
        arrow(pen, (xb, y - 12), (xa, y - 12), px=2.2, color=INK, head=10)
        cv.text((xa + xb) / 2, y - 28, r"$\ell$", ha="center", va="center",
                fontsize=fs(FS_TINY), color=INK)
    cv.text(X, 240, "what each kind adds to a growth rate that is otherwise negative",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)

    # ---- fitness, in the grammar of the trait slide -----------------------------------------
    cv.plot([X, 1150], [216, 216], color=LINE, lw=lw(1.2))
    y = 160
    cv.text(X, y, "host fitness  =", ha="left", va="center", fontsize=fs(FS_SMALL),
            fontweight="bold")
    x = X + 168
    g_chip(cv, x, y)
    cv.text(x + 58, y, "+", ha="center", va="center", fontsize=fs(FS_LABEL), color=MUTED)
    m_chip(cv, x + 116, y)
    cv.text(x + 176, y, "+", ha="center", va="center", fontsize=fs(FS_LABEL), color=MUTED)
    g_chip(cv, x + 234, y)
    cv.text(x + 262, y, "\u00d7", ha="center", va="center", fontsize=fs(FS_SMALL), color=MUTED)
    m_chip(cv, x + 290, y)
    for xt, sym, word in [(x, r"$s_g$", "the allele"), (x + 116, r"$s_m$", "direct effect"),
                          (x + 262, r"$e$", "interaction")]:
        cv.text(xt, y - 34, sym, ha="center", va="center", fontsize=fs(FS_SMALL))
        cv.text(xt, y - 58, word, ha="center", va="center", fontsize=fs(FS_TINY), color=MUTED)
    cv.text(X, 66, "the microbe can help directly, or only together with the allele",
            ha="left", va="center", fontsize=fs(FS_TINY), color=INK)
    status(cv, "inprep", "after Fig. 1A of the rescue manuscript (scripts/qme_rescue_model.py); the population is a cartoon")
    save(fig, "s25_rescue_setup")


# =================================================================== s26 ===
def result():
    fig = figure_px(FULL)
    cv = canvas(fig)
    tr = rescue_trajectories()
    Nmin = {k: float(v[2].min()) for k, v in tr.items()}
    t50 = {k: float(v[0][np.argmax(v[1] >= 0.5)]) for k, v in tr.items()}
    assert abs(Nmin["direct"] - 0.21) < 0.01 and abs(Nmin["host"] - 0.03) < 0.01
    assert abs(t50["inter"] - 48) < 1.5 and abs(t50["host"] - 78) < 1.5
    assert abs(t50["direct"] - t50["host"]) < 1.5       # a direct effect does not speed it

    cv.text(112, 508, "Host abundance", ha="left", va="baseline",
            fontsize=fs(FS_TITLE), fontweight="bold")
    ax = px_axes(fig, 112, 100, 520, 370)
    for k in ("host", "direct", "inter"):
        t, pg, N = tr[k]
        ax.plot(t, N, color=COL[k], lw=lw(LW_DATA + (0.6 if k != "host" else 0)),
                ls=DASH[k], zorder=3 if k != "host" else 2)
    ax.set_xlim(0, T_END)
    ax.set_ylim(0, 0.92)
    ax.set_yticks([0, 0.4, 0.8])
    ax.set_xticks([0, 100, 200, 300])
    ax.set_xlabel("time (host generations)")
    ax.set_ylabel("host abundance")
    ax.text(292, 0.115, "microbe removed", ha="right", va="bottom",
            fontsize=fs(FS_SMALL), color=MUTED)
    ax.text(128, 0.16, r"direct effect $s_m$", ha="left", va="bottom",
            fontsize=fs(FS_SMALL), color=OCHRE, fontweight="bold")
    ax.text(104, 0.55, r"interaction $e$", ha="right", va="center",
            fontsize=fs(FS_SMALL), color=TEAL, fontweight="bold")
    kd = int(np.argmin(tr["direct"][2]))
    ax.annotate("a direct effect buffers the crash", xy=(tr["direct"][0][kd], Nmin["direct"]),
                xytext=(150, 0.44), ha="left", va="center", fontsize=fs(FS_SMALL),
                color=INK, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=INK, lw=lw(1.4), shrinkA=3,
                                shrinkB=5))

    cv.text(722, 508, "Frequency of the host allele", ha="left", va="baseline",
            fontsize=fs(FS_TITLE), fontweight="bold")
    bx = px_axes(fig, 722, 100, 430, 370)
    for k in ("host", "direct", "inter"):
        t, pg, N = tr[k]
        bx.plot(t, pg, color=COL[k], lw=lw(LW_DATA + (0.6 if k != "host" else 0)),
                ls=DASH[k], zorder=3 if k != "host" else 2)
    bx.axhline(0.5, color=LINE, lw=lw(1.2), zorder=1)
    bx.set_xlim(0, 200)
    bx.set_ylim(0, 1.04)
    bx.set_yticks([0, 0.5, 1])
    bx.set_xticks([0, 100, 200])
    bx.set_xlabel("time (host generations)")
    bx.set_ylabel("allele frequency")
    bx.annotate("an interaction speeds\nthe allele's spread", xy=(t50["inter"], 0.5),
                xytext=(112, 0.22), ha="left", va="center", fontsize=fs(FS_SMALL),
                color=INK, fontweight="bold", linespacing=1.15,
                arrowprops=dict(arrowstyle="-", color=INK, lw=lw(1.4), shrinkA=3,
                                shrinkB=5))
    bx.text(196, 0.90, "direct effect and\nmicrobe removed\ncoincide", ha="right",
            va="top", fontsize=fs(FS_TINY), color=MUTED, linespacing=1.2)
    status(cv, "inprep",
           "QME-reduced rescue model; same total microbial benefit (0.05) placed in $s_m$ or in $e$")
    save(fig, "s26_rescue_result")


# =================================================================== s27 ===
P_MANIFOLD = Params(beta=10.0, ell0=4.0, phi=0.0, sg=0.03, sm=0.05, e=0.10,
                    r=0.04, c=0.05)          # the statement's manifold panel


def off_manifold_paths(p, starts, t_end=120.0):
    """Full model (p_m, p_g, D) from points displaced off the QME manifold."""
    out = []
    for pg0, dpm, dD in starts:
        y0 = [float(pm_qme(pg0, p)) + dpm, pg0, float(D_qme(pg0, p)) + dD]
        sol = solve_ivp(lambda t, z: frequency_rhs(t, z, p), (0, t_end), y0,
                        method="LSODA", dense_output=True, rtol=1e-10, atol=1e-13)
        z2 = sol.sol(2.0)                       # collapsed within two generations
        assert abs(z2[0] - float(pm_qme(z2[1], p))) < 1e-4
        assert abs(z2[2] - float(D_qme(z2[1], p))) < 3e-5
        out.append(sol)
    return out


def qme():
    """Read left to right.  Left: the same move in two fields, and what it
    delivers in each.  Quasi-linkage equilibrium lets fast associations among
    loci be written in terms of slow allele frequencies, which is how
    multilocus population genetics yields the genetic variances of
    quantitative genetics.  Quasi-microbial equilibrium does the same with
    fast transmission and loss, and so can yield the components of
    microbiome-mediated quantitative genetics (work in progress, after
    qme_qgmmt_growth_rate_derivations.md).  Right: the QME curve of the
    rescue model (statement's manifold parameters), drawn without the fast
    relaxation onto it, and the point that the reduction is not tied to that
    model."""
    fig = figure_px(FULL)
    cv = canvas(fig)
    pen = Pen(cv)
    cv.text(40, 520, "The same move, in two fields", ha="left", va="baseline",
            fontsize=fs(FS_TITLE), fontweight="bold")
    rows = [dict(y=390, col=CRIMSON, name="Multilocus genetics", fast="recombination",
                 clos="quasi-linkage\nequilibrium", kinds=("g",),
                 gives="the genetic variances of\nquantitative genetics", wip=False),
            dict(y=190, col=TEAL, name="Host\u2013microbe dynamics", fast="transmission, loss",
                 clos="quasi-microbial\nequilibrium", kinds=("g", "m", "c"),
                 gives="the components of microbiome-\nmediated quantitative genetics",
                 wip=True)]
    for r in rows:
        y, col = r["y"], r["col"]
        cv.text(40, y + 72, r["name"], ha="left", va="center", fontsize=fs(FS_LABEL),
                fontweight="bold", color=col)
        for word, yy, n, c_ in [(r["fast"], y + 22, 14, col), ("selection", y - 14, 3, INK)]:
            cv.text(40, yy + 15, word, ha="left", va="center", fontsize=fs(FS_TINY))
            tick_timeline(pen, 40, yy, 120, n, color=c_, tick=9, px=1.8, arrow_head=7)
        cv.text(40, y - 44, "fast against slow", ha="left", va="center", fontsize=fs(FS_TINY),
                color=MUTED)
        arrow(pen, (186, y + 4), (222, y + 4), px=2.4, color=NAVY, head=10)
        cv.add_patch(FancyBboxPatch((236, y - 30), 150, 68,
                                    boxstyle="round,pad=0,rounding_size=10", facecolor=PAPER,
                                    edgecolor=col, linewidth=lw(2.0), zorder=1))
        cv.text(311, y + 4, r["clos"], ha="center", va="center", fontsize=fs(FS_SMALL),
                fontweight="bold", linespacing=1.15)
        arrow(pen, (400, y + 4), (436, y + 4), px=2.4, color=NAVY, head=10,
              ls=(0, (2.6, 2.0)) if r["wip"] else "-")
        x = 450
        w = {"g": 150} if len(r["kinds"]) == 1 else {"g": 66, "m": 56, "c": 42}
        for kind in r["kinds"]:
            variance_block(cv, x, y - 12, w[kind], 34, kind)
            cv.text(x + w[kind] / 2, y + 36, {"g": r"$G_A$", "m": r"$M_A$", "c": r"$C_A$"}[kind],
                    ha="center", va="center", fontsize=fs(FS_SMALL))
            x += w[kind]
        cv.text(450, y - 30, r["gives"], ha="left", va="top", fontsize=fs(FS_TINY),
                color=INK, linespacing=1.2)
        if r["wip"]:
            status_mark(cv, "inprep", 456, y + 72, color=NAVY)
            cv.text(470, y + 72, "work in progress", ha="left", va="center",
                    fontsize=fs(FS_TINY), color=NAVY)
    cv.plot([40, 640], [290, 290], color=LINE, lw=lw(1.2))
    cv.text(40, 74, "The same reduction applies wherever microbes turn over faster than hosts evolve.",
            ha="left", va="center", fontsize=fs(FS_TINY), color=MUTED)

    # ---- the curve, in the rescue model, and what it does to the variance components ------
    X = 710
    cv.plot([X - 34, X - 34], [60, 520], color=LINE, lw=lw(1.2))
    cv.text(X, 520, "In the rescue model", ha="left", va="baseline", fontsize=fs(FS_LABEL),
            fontweight="bold")
    p = P_MANIFOLD
    grid = np.linspace(0, 1, 501)
    pmM, DM = pm_qme(grid, p), D_qme(grid, p)
    ax = px_axes(fig, X + 84, 336, 300, 152)
    ax.plot(pmM, DM, color=TEAL, lw=lw(5.0), zorder=3, solid_capstyle="round")
    for k in (95, 330):
        ax.annotate("", xy=(pmM[k + 16], DM[k + 16]), xytext=(pmM[k], DM[k]),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.9,head_width=0.5",
                                    color=TEAL, lw=0, mutation_scale=16), zorder=4)
    for k, lab, ha, dx, dy in [(0, "allele absent", "left", 0.00016, -0.00012),
                               (250, "allele at 50%", "center", 0, 0.00015),
                               (500, "allele fixed", "right", -0.00016, -0.00012)]:
        ax.plot(pmM[k], DM[k], "o", ms=8, mfc="white", mec=INK, mew=lw(2.0), zorder=7,
                clip_on=False)
        ax.text(pmM[k] + dx, DM[k] + dy, lab, ha=ha, va="center", fontsize=fs(FS_TINY),
                color=INK)
    span = pmM[-1] - pmM[0]
    ax.set_xlim(pmM[0] - 0.22 * span, pmM[-1] + 0.22 * span)
    ax.set_ylim(-2.6e-4, 9.4e-4)
    ax.set_xticks([])
    ax.set_yticks([0.0])
    ax.set_xlabel(r"share of hosts carrying the microbe, $p_m$", fontsize=fs(FS_TINY))
    ax.set_ylabel("allele\u2013microbe\nassociation, $D$", fontsize=fs(FS_TINY), linespacing=1.1)

    # The QGMMT components of variance in host growth rate, for one allele and one
    # microbe (the derivation note, sections 8-9): the additive effects are the
    # partial regression coefficients of growth rate on allele and microbe.
    pg = grid[1:-1]
    pm_, D_ = pmM[1:-1], DM[1:-1]
    Vg, Vm = pg * (1 - pg), pm_ * (1 - pm_)
    pgm = pg * pm_ + D_
    cg = p.sg * Vg + p.sm * D_ + p.e * pgm * (1 - pg)        # Cov(g, r)
    cm = p.sg * D_ + p.sm * Vm + p.e * pgm * (1 - pm_)       # Cov(m, r)
    det = Vg * Vm - D_ ** 2
    gam, om = (cg * Vm - cm * D_) / det, (cm * Vg - cg * D_) / det
    GA, MA, CA = gam ** 2 * Vg, om ** 2 * Vm, 2 * gam * om * D_
    assert np.allclose(gam, p.sg + p.e * pm_, rtol=0.02) and np.allclose(om, p.sm + p.e * pg, rtol=0.02)
    cv.plot([X, 1150], [276, 276], color=LINE, lw=lw(1.2))
    cv.text(X, 250, "so the variance in host growth rate splits as", ha="left", va="center",
            fontsize=fs(FS_SMALL), fontweight="bold")
    status_mark(cv, "inprep", X + 392, 250, color=NAVY)
    # at QME the association is small, so the additive effects reduce to
    # s_g + e p_m (allele) and s_m + e p_g (microbe); the build checks this above
    for y, txt in [(214, r"$G_A = (s_g + e\,p_m)^2\; p_g(1 - p_g)$"),
                   (184, r"$M_A = (s_m + e\,p_g)^2\; p_m(1 - p_m)$"),
                   (154, r"$C_A = 2\,(s_g + e\,p_m)(s_m + e\,p_g)\, D$")]:
        cv.text(X, y, txt, ha="left", va="center", fontsize=fs(FS_SMALL))
    cv.text(X, 104, r"$p_m$ and $D$ follow the curve" + "\n" + r"above, so all three depend" + "\n" + r"on $p_g$ alone",
            ha="left", va="center", fontsize=fs(FS_TINY), color=TEAL, linespacing=1.2,
            fontweight="bold")
    bx = px_axes(fig, X + 262, 62, 172, 68)
    top = max(GA.max(), MA.max())
    kC = 10 ** np.floor(np.log10(0.6 * top / CA.max()))
    bx.plot(pg, GA, color=TEAL, lw=lw(2.6))
    bx.plot(pg, MA, color=OCHRE, lw=lw(2.6), ls=(0, (3.2, 1.6)))
    bx.plot(pg, CA * kC, color=INK, lw=lw(2.0), ls=(0, (1, 1.6)))
    bx.text(0.20, GA[100] + 0.12 * top, r"$G_A$", ha="center", va="bottom", fontsize=fs(FS_TINY),
            color=TEAL)
    bx.text(0.80, MA[400] + 0.08 * top, r"$M_A$", ha="right", va="bottom", fontsize=fs(FS_TINY),
            color=OCHRE)
    bx.text(0.62, CA[310] * kC - 0.04 * top, rf"$C_A \times {int(kC)}$", ha="left", va="top",
            fontsize=fs(FS_TINY), color=INK)
    bx.set_xlim(0, 1)
    bx.set_ylim(0, top * 1.12)
    bx.set_xticks([0, 1])
    bx.set_yticks([])
    bx.spines["left"].set_visible(False)
    bx.set_xlabel(r"allele frequency, $p_g$", fontsize=fs(FS_TINY), labelpad=-8)
    print(f"  [check] QME components at p_g = 0.5: G_A = {GA[249]:.2e}, M_A = {MA[249]:.2e}, C_A = {CA[249]:.2e} (C_A drawn x{int(kC)})")
    status(cv, "inprep", "the QME reduction (rescue manuscript, in prep.); components from the QGMMT definitions on its curve: work in progress")
    save(fig, "s27_qme")


def main():
    setup()
    result()
    qme()


if __name__ == "__main__":
    main()
