#!/usr/bin/env python3
"""
wf_correlations.py
==================
Random genetic drift builds correlations among loci in a finite, NON-recombining
population.  Produces an animated GIF with two regions:

    LEFT   the population itself -- N haploid sequences, every locus biallelic
           and polymorphic, each allele drawn as a bright / dark shade of that
           locus' hue.  Under drift the sequences coalesce onto a few shared
           haplotypes, so the columns organise into correlated (and anti-
           correlated) blocks.

    RIGHT  the dynamics of the genetic correlations themselves.  One line per
           locus pair (there are L-choose-2 of them) tracking the signed allelic
           correlation r_ij over generations.  All start scattered near zero and
           fan out toward +/-1 as drift commits the population to a shared
           genealogy.  A line stops when one of its loci fixes (r undefined).

Model (haploid, discrete non-overlapping Wright-Fisher, N constant, NO mutation,
NO recombination):
    - genome = L biallelic loci, all polymorphic at the start.  Initial state is
      linkage equilibrium: every locus is an independent Bernoulli(1/2) across
      individuals, so all pairwise correlations start ~0.
    - each generation the next population is N whole-haplotype copies of randomly
      chosen parents (clonal Wright-Fisher resampling).  With no recombination,
      the entire genome rides one shared genealogy, so as lineages coalesce the
      loci become correlated purely from drift -- this is the finite-population
      build-up of linkage disequilibrium (Hill-Robertson), here in the fully
      linked limit.

The signed correlation of allelic state between loci i and j is
    r_ij = D_ij / sqrt(p_i(1-p_i) p_j(1-p_j)),   D_ij = p_ij - p_i p_j,
the usual LD correlation coefficient; r^2 is its square.  It is undefined once
either locus fixes.

Note this is transient: with no mutation the only true equilibrium is fixation
(a single haplotype, all loci monomorphic).  The interesting picture is the
intermediate regime, where the population sits on ~2 complementary haplotypes
and the surviving polymorphic loci are perfectly correlated.

Output: wf_correlations.gif
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, PillowWriter
from itertools import combinations

# --------------------------------------------------------------------------
# palette : Rosé Pine Dawn (light)
# --------------------------------------------------------------------------
BASE    = "#faf4ed"   # ground
SURFACE = "#fffaf3"   # panel
OVERLAY = "#f2e9e1"   # cell gutters
MUTED   = "#9893a5"   # faint
SUBTLE  = "#797593"   # secondary text
TEXT    = "#575279"   # primary text / ink
HL_MED  = "#dfdad9"   # hairlines

# accent hues (love, gold, rose, pine, foam, iris)
ACCENTS = ["#b4637a", "#ea9d34", "#d7827e", "#286983", "#56949f", "#907aa9"]


def _hex2rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[k:k + 2], 16) for k in (0, 2, 4)], dtype=float)


def _rgb2hex(c):
    c = np.clip(c, 0, 255).astype(int)
    return "#%02x%02x%02x" % (c[0], c[1], c[2])


def _mix(a, b, t):
    return _hex2rgb(a) * (1 - t) + _hex2rgb(b) * t


# per-locus allele shades: allele 0 = deep, allele 1 = bright tint of the hue.
def locus_shades(hue):
    deep   = _rgb2hex(_mix(hue, TEXT, 0.18))      # saturated, slightly deepened
    bright = _rgb2hex(_mix(hue, "#ffffff", 0.62))  # pale tint, still clearly the hue
    return _hex2rgb(deep), _hex2rgb(bright)


# --------------------------------------------------------------------------
# model
# --------------------------------------------------------------------------
def init_pop(N, L, rng):
    """(N, L) binary matrix; every locus an independent Bernoulli(1/2) = LE."""
    return rng.integers(0, 2, size=(N, L), dtype=np.int8)


def wf_step(g, rng):
    """One clonal Wright-Fisher generation: copy whole haplotypes (no recomb)."""
    N = g.shape[0]
    return g[rng.integers(0, N, N)].copy()


def corr(g, i, j):
    """Signed allelic correlation r_ij; np.nan if either locus is fixed."""
    a, b = g[:, i].astype(float), g[:, j].astype(float)
    pa, pb = a.mean(), b.mean()
    if pa in (0.0, 1.0) or pb in (0.0, 1.0):
        return np.nan
    D = (a * b).mean() - pa * pb
    denom = np.sqrt(pa * (1 - pa) * pb * (1 - pb))
    if abs(denom)>0:
        D /= denom
    else:
        D /= abs(D) 
    return D


def n_haplotypes(g):
    return len(np.unique(g, axis=0))


def n_polymorphic(g):
    p = g.mean(axis=0)
    return int(np.sum((p > 0) & (p < 1)))


# --------------------------------------------------------------------------
# simulate + score
# --------------------------------------------------------------------------
def simulate(N, L, n_gen, seed):
    rng = np.random.default_rng(seed)
    g = init_pop(N, L, rng)
    hist = [g.copy()]
    for _ in range(n_gen):
        g = wf_step(g, rng)
        hist.append(g.copy())
    return hist


def find_money_shot(hist):
    """
    First generation where the population has collapsed to <=2 haplotypes while
    keeping as many polymorphic (hence perfectly correlated) loci as possible.
    Returns (score, g_star) or (None, None).
    """
    best = None
    for g in range(1, len(hist)):
        pop = hist[g]
        nh = n_haplotypes(pop)
        if nh > 2:
            continue
        m = n_polymorphic(pop)
        if m < 4:
            continue
        # balance of the two surviving haplotypes (0..0.5)
        _, counts = np.unique(pop, axis=0, return_counts=True)
        balance = counts.min() / counts.sum() if len(counts) == 2 else 0.5
        # prefer many polymorphic loci and a balanced split; keep g* in a window
        # that leaves a legible build-up (not instant, not dragged out)
        window = -0.03 * max(0, 6 - g) - 0.006 * max(0, g - 40)
        score = m + 2.0 * balance + window
        if best is None or score > best[0]:
            best = (score, g)
        break  # first qualifying generation is the cleanest "just collapsed"
    return best if best is not None else (None, None)


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def build_corr_history(hist, pairs):
    """(n_frames, n_pairs) signed correlations; np.nan where undefined."""
    C = np.full((len(hist), len(pairs)), np.nan)
    for t, pop in enumerate(hist):
        for k, (i, j) in enumerate(pairs):
            C[t, k] = corr(pop, i, j)
    return C


def render(hist, out_path, fps=14, hold_ms=1700, max_frames=70):
    from PIL import Image

    N, L = hist[0].shape
    pairs = list(combinations(range(L), 2))
    Cfull = build_corr_history(hist, pairs)
    n_gen = len(hist) - 1

    # sample generations down to <= max_frames (keep first & last)
    if n_gen + 1 > max_frames:
        idx = np.unique(np.linspace(0, n_gen, max_frames).round().astype(int))
    else:
        idx = np.arange(n_gen + 1)

    # locus hues + their two shades
    hues   = [ACCENTS[k % len(ACCENTS)] for k in range(L)]
    shades = [locus_shades(h) for h in hues]                 # (deep, bright) RGB
    pair_colors = [hues[i] for (i, j) in pairs]              # colour line by 1st locus

    fig = plt.figure(figsize=(10.6, 5.7), dpi=100)
    fig.patch.set_facecolor(BASE)
    axP = fig.add_axes([0.045, 0.085, 0.40, 0.74])           # population
    axC = fig.add_axes([0.560, 0.135, 0.405, 0.62])          # correlations

    def draw_pop(pop):
        axP.clear()
        axP.set_facecolor(SURFACE)
        order = np.lexsort(pop.T[::-1])                      # group haplotypes
        sp = pop[order]
        img = np.empty((N, L, 3))
        for loc in range(L):
            deep, bright = shades[loc]
            img[:, loc, :] = np.where(sp[:, loc][:, None] == 1, bright, deep) / 255.0
        axP.imshow(img, aspect="auto", interpolation="nearest",
                   extent=[0, L, 0, N], zorder=1)
        for x in range(L + 1):                               # norns-style gutters
            axP.plot([x, x], [0, N], color=BASE, lw=1.4, zorder=2)
        axP.set_xlim(0, L); axP.set_ylim(0, N)
        axP.set_xticks([]); axP.set_yticks([])
        for s in axP.spines.values():
            s.set_visible(False)

    def draw_corr(gen):
        axC.clear()
        axC.set_facecolor(SURFACE)
        axC.set_xlim(0, n_gen); axC.set_ylim(-1.05, 1.05)
        axC.axhline(0, color=MUTED, lw=1.0, zorder=1)
        for yv in (-1, 1):
            axC.axhline(yv, color=HL_MED, lw=1.0, ls=(0, (2, 3)), zorder=1)
        axC.set_yticks([-1, -0.5, 0, 0.5, 1])
        axC.tick_params(colors=SUBTLE, labelsize=8.5, length=0)
        axC.set_xlabel("generation", color=SUBTLE, fontsize=9.5)
        axC.set_ylabel("Genetic correlation  ρ", color=SUBTLE, fontsize=9.5)
        for side, sp in axC.spines.items():
            sp.set_visible(side in ("left", "bottom")); sp.set_color(HL_MED)

        xs = np.arange(gen + 1)
        segs, cols = [], []                                  # break lines at NaN
        for k, base_c in enumerate(pair_colors):
            good = np.isfinite(Cfull[: gen + 1, k])
            run = []
            for t in range(gen + 1):
                if good[t]:
                    run.append((xs[t], Cfull[t, k]))
                elif len(run) >= 2:
                    segs.append(run); cols.append(base_c); run = []
                else:
                    run = []
            if len(run) >= 2:
                segs.append(run); cols.append(base_c)
        if segs:
            axC.add_collection(LineCollection(segs, colors=cols, linewidths=1.3,
                                              alpha=0.65, capstyle="round", zorder=3))

    # static titles
    # fig.text(0.045, 0.955, "drift builds genetic correlations", color=TEXT,
    #          fontsize=18, fontweight="bold", va="center")
    # fig.text(0.045, 0.915,
    #          "Wright-Fisher resampling of whole haplotypes  \u00b7  no mutation, "
    #          "no recombination  \u00b7  every locus polymorphic",
    #          color=SUBTLE, fontsize=9.5, va="center")

    base_texts = len(fig.texts)

    def draw_frame(gen):
        draw_pop(hist[gen])
        draw_corr(gen)
        del fig.texts[base_texts:]                           # refresh dynamic labels
        fig.text(0.045, 0.862, "POPULATION", color=TEXT, fontsize=11,
                 fontweight="bold", va="center")
        fig.text(0.445, 0.862, f"gen {gen}", color=SUBTLE, fontsize=10,
                 family="DejaVu Sans Mono", va="center", ha="right")
        fig.text(0.560, 0.862, "CORRELATIONS", color=TEXT, fontsize=11,
                 fontweight="bold", va="center")
        # fig.text(0.690, 0.862, "one line per locus pair", color=SUBTLE,
        #          fontsize=9.5, va="center")
        rs = Cfull[gen][np.isfinite(Cfull[gen])]
        mean_abs = float(np.mean(np.abs(rs))) if rs.size else 0.0
        committed = int(np.sum(np.abs(rs) >= 0.9)) if rs.size else 0
        fig.text(0.560, 0.815, f"mean |ρ| = {mean_abs:0.2f}", color=TEXT,
                 fontsize=9.5, family="DejaVu Sans Mono", va="center")
        # fig.text(0.965, 0.815, f"{committed} pairs at |ρ| \u2265 0.9", color=SUBTLE,
        #          fontsize=9.5, family="DejaVu Sans Mono", va="center", ha="right")

    # render each sampled generation to a PIL frame
    frames, durations = [], []
    for gen in idx:
        draw_frame(int(gen))
        fig.canvas.draw()
        w, h = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
        frames.append(Image.fromarray(buf[..., :3].copy()))
        durations.append(int(1000 / fps))
    plt.close(fig)

    durations[-1] = hold_ms                                  # hold the money shot
    frames[0].save(out_path, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, disposal=2, optimize=False)


# --------------------------------------------------------------------------
def main():
    N      = 1000          # population size (much larger than before)
    L      = 10          # loci, all polymorphic -> 45 pairs
    n_gen  = 1000         # search horizon

    best = None          # (score, seed, g_star, hist)
    for seed in range(400):
        hist = simulate(N, L, n_gen, seed)
        sc, gstar = find_money_shot(hist)
        if sc is None:
            continue
        if best is None or sc > best[0]:
            best = (sc, seed, gstar, hist)
    if best is None:
        raise SystemExit("no clean realization found; widen the search.")
    sc, seed, gstar, hist = best
    print(f"chosen seed={seed}  g*={gstar}  score={sc:.3f}  "
          f"polymorphic@g*={n_polymorphic(hist[gstar])}  "
          f"haplotypes@g*={n_haplotypes(hist[gstar])}")

    hist = hist[: gstar + 1]                      # stop at the money shot
    render(hist, "wf_correlations.gif", fps=14, hold_ms=1700, max_frames=70)
    print("wrote wf_correlations.gif")


if __name__ == "__main__":
    main()
