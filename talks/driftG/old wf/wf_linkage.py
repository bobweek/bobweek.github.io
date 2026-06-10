#!/usr/bin/env python3
"""
wf_linkage.py
=============
Wright-Fisher genetic drift builds linkage disequilibrium (LD) in a finite,
non-recombining population.  Produces an animated GIF with two panels that
echo the static "genomic linkage" figure:

    LEFT  : no recombination  -> whole haplotypes are copied each generation,
            so drift drags every locus along one shared genealogy and allelic
            associations build up (the population collapses onto a few
            correlated haplotypes -> high r^2).

    RIGHT : free recombination -> each locus is inherited from an independently
            chosen parent, so every locus has its own genealogy, associations
            are reshuffled every generation, and the population stays in
            linkage equilibrium (r^2 ~ 0) while loci drift & fix independently.

Both panels start from the SAME initial population (alleles assigned
independently across individuals = linkage equilibrium), so you literally
watch the left panel develop the correlations the anchor figure shows while
the right panel stays scrambled.

Model (haploid, discrete non-overlapping Wright-Fisher, N constant, NO
mutation):
    - genome = L loci; a handful are biallelic "segregating" loci (coloured),
      the rest are monomorphic (white) genome background.
    - each generation the next population is N independent draws.  For each
      offspring we walk the loci copying from a parent; between adjacent loci
      we switch to a freshly chosen random parent with probability c
      (c = 0  -> clonal copy of one haplotype; c = 1 -> every locus an
      independent parent = free recombination).  This is a sequential
      copying / Markov recombination model; for large N the per-gap
      probability that two adjacent loci come from different lineages is ~c,
      so adjacent-locus LD decays by ~(1-c) per generation from recombination.

Outputs: wf_linkage.gif
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.animation import FuncAnimation, PillowWriter

# --------------------------------------------------------------------------
# palette  (allele bands sampled from the anchor figure; chrome from norns)
# --------------------------------------------------------------------------
BG       = "#0E0F12"   # near-black ground
BAR      = "#F5F6F7"   # genome bar (white)
BAR_EDGE = "#0E0F12"   # bar outline (reads as a thin dark gutter on dark bg)
INK      = "#E8E9EC"   # primary text
INK_DIM  = "#7C828C"   # secondary text
RULE     = "#2A2D33"   # hairlines / header segments
ACCENT   = "#C7CCD4"   # active header segment / emphasis

# two-hue allele pairs, one pair per segregating locus (allele 0 / allele 1)
LOCUS_COLORS = [
    ("#3D7FE0", "#6A48A6"),   # locus A : blue   / purple   (anchor)
    ("#93C079", "#F0C530"),   # locus B : green  / gold     (anchor)
    ("#3FB0A6", "#E8794A"),   # locus C : teal   / coral
]

# --------------------------------------------------------------------------
# model
# --------------------------------------------------------------------------
def init_pop(N, L, seg_loci, rng):
    """allele matrix (N, L); non-seg loci = -1 (white); seg loci ~ Bernoulli(1/2)."""
    g = np.full((N, L), -1, dtype=np.int8)
    for loc in seg_loci:
        g[:, loc] = rng.integers(0, 2, N)
    return g


def wf_step(g, c, rng):
    """One Wright-Fisher generation with per-gap recombination prob c in [0,1]."""
    N, L = g.shape
    if c <= 0.0:                       # clonal: copy whole haplotypes (fast path)
        return g[rng.integers(0, N, N)].copy()
    if c >= 1.0:                       # free recombination: per-locus parent
        idx = rng.integers(0, N, size=(N, L))
        return np.take_along_axis(g, idx, axis=0).copy()
    # intermediate: sequential copying with switch prob c at each gap
    out = np.empty_like(g)
    parent = rng.integers(0, N, N)                       # parent for locus 0
    out[:, 0] = g[parent, 0]
    for loc in range(1, L):
        switch = rng.random(N) < c
        parent = np.where(switch, rng.integers(0, N, N), parent)
        out[:, loc] = g[parent, loc]
    return out


def r_squared(g, i, j):
    """biallelic r^2 between loci i, j; None if either locus is fixed."""
    a, b = g[:, i], g[:, j]
    pA, pB = a.mean(), b.mean()
    if pA in (0.0, 1.0) or pB in (0.0, 1.0):
        return None
    pAB = np.mean((a == 1) & (b == 1))
    D = pAB - pA * pB
    return D * D / (pA * (1 - pA) * pB * (1 - pB))


def n_haplotypes(g, seg_loci):
    return len(np.unique(g[:, seg_loci], axis=0))


# --------------------------------------------------------------------------
# simulate a full run (record every generation for both panels)
# --------------------------------------------------------------------------
def simulate(N, L, seg_loci, focal, n_gen, seed):
    rng_L = np.random.default_rng(seed)          # drift stream, left panel
    rng_R = np.random.default_rng(seed + 1)      # drift stream, right panel
    rng_0 = np.random.default_rng(seed + 7)      # shared initial population
    g0 = init_pop(N, L, seg_loci, rng_0)

    histL, histR = [g0.copy()], [g0.copy()]
    gL, gR = g0.copy(), g0.copy()
    for _ in range(n_gen):
        gL = wf_step(gL, 0.0, rng_L)             # no recombination
        gR = wf_step(gR, 1.0, rng_R)             # free recombination
        histL.append(gL.copy())
        histR.append(gR.copy())
    return histL, histR


def all_loci_polymorphic(g, seg_loci):
    return all(0 < g[:, loc].mean() < 1 for loc in seg_loci)


def evaluate(histL, histR, seg_loci, focal):
    """
    Find a stop generation g* where the NO-RECOMB panel has collapsed onto
    exactly two haplotypes that still differ at every segregating locus
    (so every locus pair is perfectly correlated -> the anchor's left panel),
    while the FREE-RECOMB panel is still scrambled (many haplotypes, low r^2).
    Returns (score, g_star) or (None, None) if the run never lands there.
    """
    i, j = focal
    started_low = histL[0] is not None and (r_squared(histL[0], i, j) or 0) < 0.25
    if not started_low:
        return None, None
    for g in range(1, len(histL)):
        gl = histL[g]
        # need 2 left haplotypes, both alleles present at every seg locus,
        # and balanced enough to read as "two groups"
        if n_haplotypes(gl, seg_loci) == 2 and all_loci_polymorphic(gl, seg_loci):
            uniq, counts = np.unique(gl[:, seg_loci], axis=0, return_counts=True)
            balance = counts.min() / counts.sum()        # 0..0.5
            gr = histR[g]
            right_hap = n_haplotypes(gr, seg_loci)
            right_r2 = r_squared(gr, i, j) or 0.0
            if right_hap >= 4 and right_r2 < 0.30 and balance >= 0.22:
                # prefer a gradual build (g* in a mid window), balanced split,
                # and a still-scrambled right panel
                window = -0.05 * max(0, 16 - g) - 0.04 * max(0, g - 30)
                score = (balance
                         + 0.04 * right_hap
                         - 0.5 * right_r2
                         + window)
                return score, g
    return None, None


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def draw_panel(ax, g, x0, w, y0, h, seg_loci, focal):
    """Draw a population as stacked genome bars within [x0,x0+w] x [y0,y0+h]."""
    N, L = g.shape
    seg_idx = {loc: k for k, loc in enumerate(seg_loci)}
    row_h = h / N
    bar_h = row_h * 0.74
    pad = row_h * 0.13
    cell_w = w / L
    for r in range(N):
        by = y0 + h - (r + 1) * row_h + pad
        # genome bar (white, rounded) -------------------------------------
        bar = FancyBboxPatch(
            (x0, by), w, bar_h,
            boxstyle="round,pad=0,rounding_size=" + str(bar_h * 0.45),
            linewidth=1.1, edgecolor=BAR_EDGE, facecolor=BAR, mutation_aspect=1.0,
            zorder=2,
        )
        ax.add_patch(bar)
        # coloured allele bands at segregating loci -----------------------
        for loc in seg_loci:
            allele = g[r, loc]
            col = LOCUS_COLORS[seg_idx[loc] % len(LOCUS_COLORS)][allele]
            bx = x0 + loc * cell_w + cell_w * 0.12
            bw = cell_w * 0.76
            inset = bar_h * 0.16
            ax.add_patch(Rectangle((bx, by + inset), bw, bar_h - 2 * inset,
                                   facecolor=col, edgecolor="none", zorder=3))
    # focal-pair carets above the panel ----------------------------------
    for loc in focal:
        cx = x0 + (loc + 0.5) * cell_w
        ax.plot([cx - cell_w * 0.16, cx, cx + cell_w * 0.16],
                [y0 + h + h * 0.018, y0 + h + h * 0.052, y0 + h + h * 0.018],
                color=INK_DIM, lw=1.0, zorder=4, solid_capstyle="round")


def render_gif(histL, histR, seg_loci, focal, out_path,
               fps=7, hold=14):
    N, L = histL[0].shape
    n_frames = len(histL)
    order = list(range(n_frames)) + [n_frames - 1] * hold   # hold final state

    fig = plt.figure(figsize=(11.0, 5.6), dpi=110)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # static panel geometry
    pxw, gap = 0.40, 0.10
    pLx, pRx = 0.055, 0.055 + pxw + gap
    py0, pyh = 0.10, 0.66

    i, j = focal

    def frame(fi):
        ax.clear()
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        ax.set_facecolor(BG)
        gen = order[fi]
        gL, gR = histL[gen], histR[gen]

        # title + segmented header rule (norns-style tabs) ----------------
        ax.text(0.055, 0.945, "drift builds linkage", color=INK,
                fontsize=18, fontweight="bold", family="DejaVu Sans", va="center")
        ax.text(0.055, 0.895,
                "same start, Wright-Fisher resampling, no mutation  ·  "
                "carets mark the focal locus pair",
                color=INK_DIM, fontsize=9.5, va="center")
        segs = 8
        for s in range(segs):
            xs = 0.62 + s * (0.33 / segs)
            ax.plot([xs, xs + 0.33 / segs - 0.006], [0.95, 0.95],
                    color=ACCENT if s < (gen * segs)//max(1, n_frames-1) else RULE,
                    lw=2.0, solid_capstyle="butt")
        ax.text(0.955, 0.905, f"gen {gen:>2d}", color=INK_DIM, fontsize=10,
                ha="right", va="center", family="DejaVu Sans Mono")

        # panels ----------------------------------------------------------
        draw_panel(ax, gL, pLx, pxw, py0, pyh, seg_loci, focal)
        draw_panel(ax, gR, pRx, pxw, py0, pyh, seg_loci, focal)

        # panel labels ----------------------------------------------------
        ax.text(pLx, py0 + pyh + 0.075, "NO RECOMBINATION", color=INK,
                fontsize=11, fontweight="bold", va="center")
        ax.text(pRx, py0 + pyh + 0.075, "FREE RECOMBINATION", color=INK,
                fontsize=11, fontweight="bold", va="center")

        # readouts --------------------------------------------------------
        def readout(x, g):
            r2 = r_squared(g, i, j)
            nh = n_haplotypes(g, seg_loci)
            r2s = "fixed" if r2 is None else f"{r2:0.2f}"
            ax.text(x, py0 - 0.045,
                    f"r\u00b2 = {r2s}", color=INK, fontsize=12,
                    family="DejaVu Sans Mono", va="center")
            ax.text(x + 0.165, py0 - 0.045,
                    f"{nh} haplotype" + ("s" if nh != 1 else ""),
                    color=INK_DIM, fontsize=10.5, family="DejaVu Sans Mono",
                    va="center")
            # r^2 meter bar
            mlen = pxw * 0.42
            mx = x
            my = py0 - 0.085
            ax.add_patch(Rectangle((mx, my), mlen, 0.012, facecolor=RULE,
                                   edgecolor="none"))
            frac = 0.0 if r2 is None else r2
            ax.add_patch(Rectangle((mx, my), mlen * frac, 0.012,
                                   facecolor=ACCENT, edgecolor="none"))
        readout(pLx, gL)
        readout(pRx, gR)
        return []

    anim = FuncAnimation(fig, frame, frames=len(order), interval=1000 / fps,
                         blit=False)
    anim.save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    N        = 24
    L        = 9
    seg_loci = [1, 4, 7]
    focal    = (1, 4)
    n_gen    = 110

    best = None  # (score, seed, g_star, histL, histR)
    for seed in range(600):
        hL, hR = simulate(N, L, seg_loci, focal, n_gen, seed)
        sc, gstar = evaluate(hL, hR, seg_loci, focal)
        if sc is None:
            continue
        if best is None or sc > best[0]:
            best = (sc, seed, gstar, hL, hR)
    if best is None:
        raise SystemExit("no clean realization found; widen the search.")
    sc, seed, gstar, hL, hR = best
    print(f"chosen seed={seed}  g*={gstar}  score={sc:.3f}")

    # truncate to the money-shot window and hold there
    hL, hR = hL[: gstar + 1], hR[: gstar + 1]
    render_gif(hL, hR, seg_loci, focal, "/home/claude/wf_linkage.gif",
               fps=6, hold=18)
    print("wrote wf_linkage.gif")


if __name__ == "__main__":
    main()
