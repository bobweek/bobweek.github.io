#!/usr/bin/env python3
"""
wf_correlations.py
==================
Random genetic drift builds correlations among loci in a finite, NON-recombining
population.  Animated GIF with two regions:

    LEFT   the population -- N haploid sequences, every locus biallelic and
           polymorphic, each allele a bright / dark shade of that locus' hue.
    RIGHT  the dynamics of the genetic correlations -- one line per locus pair
           (L-choose-2 of them), the signed allelic correlation rho_ij over
           generations, fanning from ~0 toward +/-1 as drift commits the
           population to a shared genealogy.

Model: haploid Wright-Fisher, N constant, NO mutation, NO recombination.  Start
in linkage equilibrium (each locus an independent Bernoulli(1/2)).  Each
generation the next population is N whole-haplotype copies of random parents.

    rho_ij = D_ij / sqrt(p_i(1-p_i) p_j(1-p_j)),   D_ij = p_ij - p_i p_j

is the usual LD correlation coefficient; undefined once either locus fixes.

--------------------------------------------------------------------------------
WHAT CHANGED (and why the old script "didn't run to final time"):

  * It was NOT the correlation denominator.  corr() returns np.nan as soon as a
    locus is fixed (pa or pb in {0,1}), so the denominator is never divided when
    it is zero -- the guard added for that case was unreachable (and if reached,
    D/=abs(D) is 0/0 when D==0, and would fabricate a +/-1 value for a fixed
    locus whose correlation is genuinely undefined).

  * The real reason: main() did `hist = hist[: gstar + 1]`, truncating the run
    at the FIRST generation the population hit <=2 haplotypes (the "money shot").
    So the animation always ended at that collapse, never at fixation, no matter
    how large n_gen was.

  * Timescale: the time to collapse to <=2 haplotypes (then to a single one)
    scales like ~2N generations.  At N=1000 that is ~2000-3000 generations, so
    n_gen=1000 was also far too short, and a 400-seed search over full histories
    with np.unique per generation is very slow.

  Fixes below: RUN_TO controls the stopping rule ("fixation" = run to final
  time, the default; "money_shot" = the old truncated behaviour).  Each run
  early-stops the instant the population fixes (cheap all-rows-equal test), so
  even N=1000 is fast.  The seed search is light and only chooses which
  realization to show; it no longer caps how far the chosen run plays.

Output: wf_correlations.gif
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from itertools import combinations

# --------------------------------------------------------------------------
# palette : Rosé Pine Dawn (light)
# --------------------------------------------------------------------------
BASE    = "#faf4ed"
SURFACE = "#fffaf3"
MUTED   = "#9893a5"
SUBTLE  = "#797593"
TEXT    = "#575279"
HL_MED  = "#dfdad9"
ACCENTS = ["#b4637a", "#ea9d34", "#d7827e", "#286983", "#56949f", "#907aa9"]


def _hex2rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[k:k + 2], 16) for k in (0, 2, 4)], dtype=float)


def _rgb2hex(c):
    c = np.clip(c, 0, 255).astype(int)
    return "#%02x%02x%02x" % (c[0], c[1], c[2])


def _mix(a, b, t):
    return _hex2rgb(a) * (1 - t) + _hex2rgb(b) * t


def locus_shades(hue):
    deep   = _rgb2hex(_mix(hue, TEXT, 0.18))
    bright = _rgb2hex(_mix(hue, "#ffffff", 0.62))
    return _hex2rgb(deep), _hex2rgb(bright)


# --------------------------------------------------------------------------
# model
# --------------------------------------------------------------------------
def init_pop(N, L, rng):
    return rng.integers(0, 2, size=(N, L), dtype=np.int8)


def wf_step(g, rng):
    N = g.shape[0]
    return g[rng.integers(0, N, N)].copy()


def corr(g, i, j):
    """Signed allelic correlation; np.nan if either locus is fixed.

    No denominator guard is needed: a fixed locus has pa or pb in {0,1}, which
    returns nan here, so the denominator is never divided when it would be zero.
    """
    a, b = g[:, i].astype(float), g[:, j].astype(float)
    pa, pb = a.mean(), b.mean()
    if pa in (0.0, 1.0) or pb in (0.0, 1.0):
        return np.nan
    D = (a * b).mean() - pa * pb
    return D / np.sqrt(pa * (1 - pa) * pb * (1 - pb))


def dprime(g, i, j):
    """Lewontin's D' = D / D_max for a biallelic pair; np.nan if either fixed.

    D = p_11 - p_i p_j (gametic covariance / linkage disequilibrium).  D_max is
    the largest |D| compatible with the marginal allele frequencies:
        D > 0 :  D_max = min(p_i(1-p_j), (1-p_i)p_j)
        D < 0 :  D_max = min(p_i p_j,   (1-p_i)(1-p_j))
    so D' in [-1, 1], reaching +/-1 exactly when one gamete type is absent
    (complete LD) -- regardless of how skewed the allele frequencies are.  This
    is why a pair heading to fixation reads as +/-1 here even though Pearson r
    (which divides by the frequency SDs) gets squeezed toward 0.
    """
    a, b = g[:, i].astype(float), g[:, j].astype(float)
    p, q = a.mean(), b.mean()
    if p in (0.0, 1.0) or q in (0.0, 1.0):
        return np.nan
    D = (a * b).mean() - p * q
    Dmax = min(p * (1 - q), (1 - p) * q) if D >= 0 else min(p * q, (1 - p) * (1 - q))
    return np.nan if Dmax == 0 else D / Dmax


def is_fixed(g):
    """True once the whole population is a single haplotype (cheap, O(N*L))."""
    return bool((g == g[0]).all())


def n_haplotypes(g):
    return len(np.unique(g, axis=0))


def n_polymorphic(g):
    p = g.mean(axis=0)
    return int(np.sum((p > 0) & (p < 1)))


# --------------------------------------------------------------------------
# simulate (early-stops at fixation, so "final time" is reached cheaply)
# --------------------------------------------------------------------------
def simulate(N, L, seed, n_gen_cap):
    rng = np.random.default_rng(seed)
    g = init_pop(N, L, rng)
    hist = [g.copy()]
    for _ in range(n_gen_cap):
        g = wf_step(g, rng)
        hist.append(g.copy())
        if is_fixed(g):
            break
    return hist


def first_two_haplotype_gen(hist):
    """First generation with <=2 haplotypes, plus the polymorphic-loci count
    there (used to pick a seed and to set the emphasis hold).  (None, 0) if it
    never happens within the simulated history."""
    for g in range(1, len(hist)):
        if n_haplotypes(hist[g]) <= 2:
            return g, n_polymorphic(hist[g])
    return None, 0


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def build_corr_history(hist, pairs, stat=corr):
    C = np.full((len(hist), len(pairs)), np.nan)
    for t, pop in enumerate(hist):
        for k, (i, j) in enumerate(pairs):
            C[t, k] = stat(pop, i, j)
    return C


def render(hist, out_path, emphasis_gen=None, fps=14, hold_ms=1700,
           emphasis_ms=1200, max_frames=80, fold=False, display_rows=60,
           ref_locus=None, measure="dprime"):
    from PIL import Image

    stat = dprime if measure == "dprime" else corr
    sym  = "D\u2032" if measure == "dprime" else "\u03c1"     # D' or rho
    N, L = hist[0].shape
    pairs = list(combinations(range(L), 2))
    Cfull = build_corr_history(hist, pairs, stat)
    n_gen = len(hist) - 1

    # Honest sampling: uniform in GENERATIONS, so time on screen is proportional
    # to real time.  (The earlier activity-weighted sampler distorted this -- it
    # starved the long two-haplotype drift phase of frames, making the haplotype
    # frequencies look frozen, and it skipped generations where a single line was
    # still moving, manufacturing apparent jumps.)  With enough frames each step
    # is a true single-or-few-generation increment; raise max_frames toward n_gen
    # for exact per-generation fidelity, or lower N for shorter runs.
    if n_gen + 1 > max_frames:
        idx = set(np.linspace(0, n_gen, max_frames).round().astype(int).tolist())
    else:
        idx = set(range(n_gen + 1))
    idx.update({0, n_gen})
    if emphasis_gen is not None:
        idx.add(int(emphasis_gen))
    idx = sorted(idx)

    hues   = [ACCENTS[k % len(ACCENTS)] for k in range(L)]
    shades = [locus_shades(h) for h in hues]
    # which pair-lines to draw, and their colors
    if ref_locus is None:
        draw_list = [(k, hues[pairs[k][0]]) for k in range(len(pairs))]
    else:                                            # only ref vs each other locus
        draw_list = []
        for k, (i, j) in enumerate(pairs):
            if i == ref_locus or j == ref_locus:
                partner = j if i == ref_locus else i
                draw_list.append((k, hues[partner]))

    fig = plt.figure(figsize=(10.6, 5.7), dpi=100)
    fig.patch.set_facecolor(BASE)
    axP = fig.add_axes([0.045, 0.085, 0.40, 0.74])
    axC = fig.add_axes([0.560, 0.135, 0.405, 0.62])

    def draw_pop(pop):
        axP.clear()
        axP.set_facecolor(SURFACE)
        order = np.lexsort(pop.T[::-1])
        sp = pop[order]
        # Show a legible number of whole sequences.  The correlations use the
        # full population; here we display a representative, haplotype-sorted
        # SAMPLE of rows so individual sequences are resolvable even at large N
        # (1000 rows can't be drawn distinctly in ~400 px).
        Nfull = sp.shape[0]
        R = min(Nfull, display_rows)
        if Nfull > R:
            sp = sp[np.linspace(0, Nfull - 1, R).round().astype(int)]
        img = np.empty((R, L, 3))
        for loc in range(L):
            deep, bright = shades[loc]
            img[:, loc, :] = np.where(sp[:, loc][:, None] == 1, bright, deep) / 255.0
        axP.imshow(img, aspect="auto", interpolation="nearest",
                   extent=[0, L, 0, R], zorder=1)
        # columns flush (no vertical gutters); tiny gaps BETWEEN rows so the eye
        # reads each sequence left-to-right
        for y in range(1, R):
            axP.plot([0, L], [y, y], color=BASE, lw=1.1, zorder=2)
        axP.set_xlim(0, L); axP.set_ylim(0, R)
        axP.set_xticks([]); axP.set_yticks([])
        for s in axP.spines.values():
            s.set_visible(False)

    def draw_corr(gen):
        axC.clear()
        axC.set_facecolor(SURFACE)
        axC.set_xlim(0, n_gen)
        lname = "Lewontin's D\u2032" if measure == "dprime" else "Genetic correlation"
        if fold:
            axC.set_ylim(-0.03, 1.05)
            axC.axhline(1, color=HL_MED, lw=1.0, ls=(0, (2, 3)), zorder=1)
            axC.axhline(0, color=MUTED, lw=1.0, zorder=1)
            axC.set_yticks([0, 0.5, 1])
            axC.set_ylabel(f"{lname}  |{sym}|", color=SUBTLE, fontsize=9.5)
        else:
            axC.set_ylim(-1.05, 1.05)
            axC.axhline(0, color=MUTED, lw=1.0, zorder=1)
            for yv in (-1, 1):
                axC.axhline(yv, color=HL_MED, lw=1.0, ls=(0, (2, 3)), zorder=1)
            axC.set_yticks([-1, -0.5, 0, 0.5, 1])
            ylab = (f"{lname}  {sym}" if ref_locus is None
                    else f"{sym} with locus {ref_locus + 1}")
            axC.set_ylabel(ylab, color=SUBTLE, fontsize=9.5)
        axC.tick_params(colors=SUBTLE, labelsize=8.5, length=0)
        axC.set_xlabel("generation", color=SUBTLE, fontsize=9.5)
        for side, sp in axC.spines.items():
            sp.set_visible(side in ("left", "bottom")); sp.set_color(HL_MED)

        xs = np.arange(gen + 1)
        segs, cols = [], []
        for k, base_c in draw_list:
            good = np.isfinite(Cfull[: gen + 1, k])
            run = []
            for t in range(gen + 1):
                if good[t]:
                    v = abs(Cfull[t, k]) if fold else Cfull[t, k]
                    run.append((xs[t], v))
                elif len(run) >= 2:
                    segs.append(run); cols.append(base_c); run = []
                else:
                    run = []
            if len(run) >= 2:
                segs.append(run); cols.append(base_c)
        if segs:
            axC.add_collection(LineCollection(segs, colors=cols, linewidths=1.3,
                                              alpha=0.65, capstyle="round", zorder=3))

    base_texts = len(fig.texts)

    def draw_frame(gen):
        draw_pop(hist[gen])
        draw_corr(gen)
        del fig.texts[base_texts:]
        fig.text(0.045, 0.862, "POPULATION", color=TEXT, fontsize=11,
                 fontweight="bold", va="center")
        shown = min(N, display_rows)
        ptag = f"gen {gen}" if shown >= N else f"{shown} of {N} seqs  ·  gen {gen}"
        fig.text(0.445, 0.862, ptag, color=SUBTLE, fontsize=10,
                 family="DejaVu Sans Mono", va="center", ha="right")
        base_hdr = "LINKAGE D\u2032" if measure == "dprime" else "CORRELATIONS"
        hdr = base_hdr if ref_locus is None else f"{sym} vs locus {ref_locus + 1}"
        fig.text(0.560, 0.862, hdr, color=TEXT, fontsize=11,
                 fontweight="bold", va="center")
        drawn = np.array([Cfull[gen, k] for k, _ in draw_list])
        rs = drawn[np.isfinite(drawn)]
        mean_abs = float(np.mean(np.abs(rs))) if rs.size else 0.0
        nh = n_haplotypes(hist[gen])
        tag = "fixed \u2014 one haplotype" if nh == 1 else f"mean |{sym}| = {mean_abs:0.2f}"
        fig.text(0.560, 0.815, tag, color=TEXT, fontsize=9.5,
                 family="DejaVu Sans Mono", va="center")

    frames, durations = [], []
    for gen in idx:
        draw_frame(int(gen))
        fig.canvas.draw()
        w, h = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
        frames.append(Image.fromarray(buf[..., :3].copy()))
        durations.append(int(1000 / fps))
    plt.close(fig)

    durations[-1] = hold_ms                              # brief hold so the loop reads
    frames[0].save(out_path, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, disposal=2, optimize=False)


# --------------------------------------------------------------------------
def main():
    N         = 1000          # population size
    L         = 10            # loci, all polymorphic -> 45 pairs
    RUN_TO    = "fixation"    # "fixation" (run to final time) or "money_shot"
    MEASURE   = "dprime"      # "dprime" (Lewontin's D'=D/D_max) or "r" (Pearson)
    FOLD_ABS  = False         # True -> plot |.| on a 0..1 half-axis (no +/- mirror)
    PLOT_MODE = "reference"   # "all_pairs" (every L-choose-2 line) or
                              # "reference" (rho of each locus vs one reference locus;
                              #  the non-redundant view -- in the rank-1 collapse these
                              #  L-1 lines determine every other pair, no mirror band)
    N_SEEDS   = 24            # light search: only chooses which realization to show
    n_gen_cap = 8 * N         # safety cap; runs early-stop at fixation anyway
    ref_locus = 0 if PLOT_MODE == "reference" else None

    best = None               # (score, seed, hist, g_money)
    for seed in range(N_SEEDS):
        hist = simulate(N, L, seed, n_gen_cap)
        g_money, m = first_two_haplotype_gen(hist)
        score = m                                        # rich +/-1 fan at collapse
        if best is None or score > best[0]:
            best = (score, seed, hist, g_money)
    _, seed, hist, g_money = best

    fixed_gen = (len(hist) - 1) if is_fixed(hist[-1]) else None
    print(f"chosen seed={seed}  generations simulated={len(hist) - 1}  "
          f"fixation gen={fixed_gen}  2-haplotype gen={g_money}  "
          f"polymorphic@2hap={0 if g_money is None else n_polymorphic(hist[g_money])}")

    if RUN_TO == "money_shot" and g_money is not None:
        render(hist[: g_money + 1], "wf_correlations.gif", emphasis_gen=None,
               fold=FOLD_ABS, ref_locus=ref_locus, max_frames=160, measure=MEASURE)
    else:                                                # run to final time
        render(hist, "wf_correlations.gif", emphasis_gen=None,
               fold=FOLD_ABS, ref_locus=ref_locus, max_frames=160, measure=MEASURE)
    print("wrote wf_correlations.gif")


if __name__ == "__main__":
    main()
