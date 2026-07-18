"""
Two figures for the opening slide "Demographic stochasticity is super important".

  TOP  extinction.png      — population abundance under demographic stochasticity.
        Stochastic logistic  dn = (r - c n) n dt + sqrt(v n) dB  (the model's own
        abundance noise, sqrt(v n) dB). Deterministically the size is constant at
        carrying capacity, yet demographic noise drives some replicates to n = 0.

  BOTTOM fixation.png      — neutral allele-frequency drift.
        Wright-Fisher diffusion  dp = sqrt(p(1-p)/Ne) dB  from p0 = 1/2.
        Deterministically the frequency is constant, yet drift drives replicates to
        fixation (p = 1) or loss (p = 0): the "fixation/loss" faces of extinction.

Design: deck palette (white slides, ink #1b1f24, cyan #007c91 = persists,
pink #c2185b = absorbed), MULTIPLE replicate paths (the point is the spread of
fates), minimal annotation, transparent background, wide-and-short to sit in the
slide where the placeholder iframe was.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# ---- palette (from daytime.css) ----
INK  = '#1b1f24'
CYAN = '#007c91'   # persists  (survives / fixes)
PINK = '#c2185b'   # absorbed  (extinct / lost)
GRAY = '#9aa5b1'   # still segregating
AXIS = '#b7c0c9'   # spines / boundaries
TICK = '#5f6b76'   # tick labels
MEANL= '#cbd4dc'   # faint deterministic reference

mpl.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 15,
    'text.color': INK, 'axes.labelcolor': INK,
    'xtick.color': TICK, 'ytick.color': TICK,
    'axes.edgecolor': AXIS, 'axes.linewidth': 1.1, 'svg.fonttype': 'none',
})

FIGSIZE = (9.0, 2.7)     # wide & short; ~1800x540 px at dpi 200

# ---------------- simulations ----------------
def logistic_paths(n0, K, r, v, T, dt, npaths, seed):
    c = r / K; rng = np.random.default_rng(seed); steps = int(T / dt)
    t = np.linspace(0, T, steps + 1); N = np.empty((npaths, steps + 1)); N[:, 0] = n0
    for k in range(steps):
        n = N[:, k]
        nxt = n + (r - c * n) * n * dt + np.sqrt(np.maximum(n, 0) * v) * np.sqrt(dt) * rng.standard_normal(npaths)
        nxt[nxt < 0] = 0; nxt[n <= 0] = 0
        N[:, k + 1] = nxt
    return t, N

def wf_paths(p0, Ne, T, dt, npaths, seed):
    rng = np.random.default_rng(seed); steps = int(T / dt)
    t = np.linspace(0, T, steps + 1); P = np.empty((npaths, steps + 1)); P[:, 0] = p0
    for k in range(steps):
        p = P[:, k]
        nxt = np.clip(p + np.sqrt(np.maximum(p * (1 - p), 0) / Ne) * np.sqrt(dt) * rng.standard_normal(npaths), 0, 1)
        done = (p <= 0) | (p >= 1); nxt[done] = p[done]
        P[:, k + 1] = nxt
    return t, P

def first_hit(traj, lo=None, hi=None):
    for i, val in enumerate(traj):
        if lo is not None and val <= lo: return i
        if hi is not None and val >= hi: return i
    return None

def tidy(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.tick_params(length=3)

# ---------------- figure 1: extinction ----------------
def make_extinction(fname, seed=6):
    n0 = K = 30.0; r, v, T, dt, npaths = 0.6, 3.5, 110.0, 0.03, 11
    t, N = logistic_paths(n0, K, r, v, T, dt, npaths, seed)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axhline(K, color=MEANL, lw=1.2, ls=(0, (5, 4)), zorder=1)      # deterministic (constant)
    ax.axhline(0, color=AXIS, lw=1.1, zorder=2)                        # extinction boundary
    for j in range(npaths):
        h = first_hit(N[j], lo=0); extinct = h is not None
        ax.plot(t, N[j], color=(PINK if extinct else CYAN), lw=1.6, alpha=0.9,
                solid_capstyle='round', zorder=3)
        if extinct:
            ax.plot(t[h], 0, 'o', color=PINK, ms=5.5, zorder=4,
                    markeredgecolor='white', markeredgewidth=0.7)
    ax.set_xlim(0, T); ax.set_ylim(-2.5, N.max() * 1.06)
    ax.set_xlabel('time  \u2192', labelpad=2); ax.set_ylabel('population size,  $n$')
    ax.set_xticks([]); ax.set_yticks([0, K]); ax.set_yticklabels(['$0$', '$n_0$'])
    tidy(ax)
    ax.text(0.985, 0.93, 'surviving', transform=ax.transAxes, ha='right', va='top',
            color=CYAN, fontsize=13.5, fontweight='bold')
    ax.text(0.015, 0.05, 'extinctions', transform=ax.transAxes, ha='left', va='bottom',
            color=PINK, fontsize=13.5, fontweight='bold')
    fig.tight_layout(pad=0.4)
    fig.savefig(fname, dpi=200, transparent=True, bbox_inches='tight')
    plt.close(fig)
    return sum(first_hit(N[j], lo=0) is not None for j in range(npaths)), npaths

# ---------------- figure 2: fixation / loss ----------------
def make_fixation(fname, seed=0):
    p0, Ne, T, dt, npaths = 0.5, 0.6, 1.0, 0.0004, 11
    t, P = wf_paths(p0, Ne, T, dt, npaths, seed)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axhline(0.5, color=MEANL, lw=1.2, ls=(0, (5, 4)), zorder=1)     # deterministic (constant)
    ax.axhline(1.0, color=AXIS, lw=1.1, zorder=2)
    ax.axhline(0.0, color=AXIS, lw=1.1, zorder=2)
    nfix = nloss = 0
    for j in range(npaths):
        hf, hl = first_hit(P[j], hi=1), first_hit(P[j], lo=0)
        if hf is not None:
            col = CYAN; nfix += 1
            ax.plot(t[hf], 1, 'o', color=CYAN, ms=5.5, zorder=4, markeredgecolor='white', markeredgewidth=0.7)
        elif hl is not None:
            col = PINK; nloss += 1
            ax.plot(t[hl], 0, 'o', color=PINK, ms=5.5, zorder=4, markeredgecolor='white', markeredgewidth=0.7)
        else:
            col = GRAY
        ax.plot(t, P[j], color=col, lw=1.6, alpha=0.9, solid_capstyle='round', zorder=3)
    ax.set_xlim(0, T); ax.set_ylim(-0.06, 1.06)
    ax.set_xlabel('time  \u2192', labelpad=2); ax.set_ylabel('allele frequency,  $p$')
    ax.set_xticks([]); ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(['$0$', '$1/2$', '$1$'])
    tidy(ax)
    ax.text(0.985, 0.95, 'fixation', transform=ax.transAxes, ha='right', va='top',
            color=CYAN, fontsize=13.5, fontweight='bold')
    ax.text(0.985, 0.05, 'loss', transform=ax.transAxes, ha='right', va='bottom',
            color=PINK, fontsize=13.5, fontweight='bold')
    fig.tight_layout(pad=0.4)
    fig.savefig(fname, dpi=200, transparent=True, bbox_inches='tight')
    plt.close(fig)
    return nfix, nloss, npaths

if __name__ == '__main__':
    ne, npp = make_extinction('extinction.png')
    nf, nl, _ = make_fixation('fixation.png')
    print(f'extinction.png : {ne}/{npp} populations extinct, {npp-ne} surviving')
    print(f'fixation.png   : {nf} fixed, {nl} lost, {npp-nf-nl} segregating (of {npp})')
