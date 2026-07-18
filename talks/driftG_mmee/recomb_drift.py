"""
Drift builds cross-trait genetic correlation via linkage disequilibrium alone
(NO pleiotropy), and recombination slows -- but does not stop -- the collapse
rho_G -> +/-1, until recombination swamps drift (N_e * c >> 1).

Haploid Wright-Fisher, N gametes, L biallelic loci with equal additive effects.
Loci INTERLEAVED: even loci -> trait A, odd loci -> trait B (disjoint sets,
so every cross-trait covariance is LD; there is zero pleiotropy by construction).
Uniform inter-locus recombination fraction c; two-parent recombination.
Initialised at frequency 1/2 in linkage equilibrium => rho_G(0) ~ 0.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

rng = np.random.default_rng(7)

# ---- parameters -------------------------------------------------------------
N     = 200        # gametes (haploid pop size ~ N_e)
L     = 20         # loci total (10 per trait, interleaved)
R     = 120        # replicate populations
GEN   = 500        # generations
C_VALUES = [0.0, 0.005, 0.05, 0.5]   # recombination fractions -> N*c = 0, 1, 10, 100
VAR_EPS  = 1e-9

evens = np.arange(0, L, 2)   # trait A loci
odds  = np.arange(1, L, 2)   # trait B loci
Rar   = np.arange(R)[:, None]

def rho_series(c):
    """Return array (R, GEN+1) of signed rho_G trajectories for recombination c."""
    pop = (rng.random((R, N, L)) < 0.5).astype(np.int8)   # LE, freq 1/2
    out = np.full((R, GEN + 1), np.nan)
    out[:, 0] = _rho(pop)
    for g in range(1, GEN + 1):
        idxA = rng.integers(0, N, size=(R, N))
        idxB = rng.integers(0, N, size=(R, N))
        parA = pop[Rar, idxA, :]                 # (R, N, L)
        parB = pop[Rar, idxB, :]
        # crossover parity mask along the chromosome
        switches = rng.random((R, N, L - 1)) < c
        parity = np.zeros((R, N, L), dtype=np.int8)
        parity[:, :, 1:] = np.cumsum(switches, axis=2) % 2   # 0 -> parent A, 1 -> parent B
        pop = np.where(parity == 1, parB, parA)
        out[:, g] = _rho(pop)
    return out

def _rho(pop):
    gA = pop[:, :, evens].sum(axis=2).astype(float)   # (R, N)
    gB = pop[:, :, odds ].sum(axis=2).astype(float)
    gA -= gA.mean(axis=1, keepdims=True)
    gB -= gB.mean(axis=1, keepdims=True)
    vA = (gA * gA).mean(axis=1)
    vB = (gB * gB).mean(axis=1)
    cov = (gA * gB).mean(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = cov / np.sqrt(vA * vB)
    r[(vA < VAR_EPS) | (vB < VAR_EPS)] = np.nan
    return np.clip(r, -1, 1)

results = {c: rho_series(c) for c in C_VALUES}

# ---- figure -----------------------------------------------------------------
BG   = "#faf4ed"; TXT = "#575279"; GRID = "#dfdad1"
PINK = "#c2185b"; VIOLET = "#5a4cc7"; ACID = "#1b8c98"; PINE = "#797593"
panel_col = {0.0: PINK, 0.005: VIOLET, 0.05: ACID, 0.5: PINE}

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "text.color": TXT, "axes.labelcolor": TXT, "xtick.color": TXT,
    "ytick.color": TXT, "axes.edgecolor": GRID, "font.size": 13,
})

fig, axes = plt.subplots(1, len(C_VALUES), figsize=(14.5, 4.1), sharey=True)
gens = np.arange(GEN + 1)
for ax, c in zip(axes, C_VALUES):
    col = panel_col[c]
    r = results[c]
    for i in range(R):                                   # signed trajectories (the spread)
        ax.plot(gens, r[i], color=col, alpha=0.14, lw=0.8)
    rms = np.sqrt(np.nanmean(r ** 2, axis=0))            # typical |rho| = sqrt(E[rho^2])
    ax.plot(gens,  rms, color=col, lw=2.6)
    ax.plot(gens, -rms, color=col, lw=2.6)
    ax.axhline( 1, color=GRID, lw=1, ls=(0, (4, 4)))
    ax.axhline(-1, color=GRID, lw=1, ls=(0, (4, 4)))
    ax.axhline( 0, color=GRID, lw=0.8)
    final = float(np.sqrt(np.nanmean(r[:, -1] ** 2)))
    ax.text(0.96, 0.06, f"$\\sqrt{{\\mathbb{{E}}[\\rho_G^2]}}={final:.2f}$",
            transform=ax.transAxes, ha="right", va="bottom", color=col, fontsize=12)
    ax.set_ylim(-1.08, 1.08)
    ax.set_xlim(0, GEN)
    ax.set_xlabel("generation")
    tag = "drift only" if c == 0 else ("free recombination" if c == 0.5 else "")
    ttl = f"$c = {c:g}$   ($N_e c = {N*c:g}$)"
    ax.set_title(ttl + (f"\n{tag}" if tag else "\n "), color=col, fontsize=14, pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylabel(r"genetic correlation  $\rho_G$")
fig.suptitle("Drift drives $\\rho_G\\!\\to\\!\\pm1$ through linkage disequilibrium alone — "
             "recombination slows the collapse, but only free recombination ($N_e c\\gg1$) suppresses it",
             color=TXT, fontsize=15, y=1.03)
fig.text(0.5, -0.06,
         "Haploid Wright–Fisher, $N=200$, 20 non-pleiotropic loci (10 per trait), "
         "80 replicates, $\\rho_G(0)=0$.  Thin lines: replicates.  Bold: $\\pm\\sqrt{\\mathbb{E}[\\rho_G^2]}$.",
         ha="center", color=PINE, fontsize=11)
fig.tight_layout()
fig.savefig("/home/claude/recomb_drift.png", dpi=170, bbox_inches="tight")
print("saved; final sqrt(E[rho^2]) by c:",
      {c: round(float(np.sqrt(np.nanmean(results[c][:, -1] ** 2))), 3) for c in C_VALUES})
