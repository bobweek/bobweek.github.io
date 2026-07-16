#
# G_drift_2trait.py
#
# Drift-only G-matrix dynamics for two traits (paper eq. 21 with mutation and
# selection switched off), simulated with Euler-Maruyama, then animated as
# covariance ellipses (top: G itself; bottom: the correlation matrix).
#
# The stochastic matrix equation (pure drift, MVN traits):
#
#     dG = -(v/n) G dt  +  √(v/n) · √Γ : dB_G ,
#
#   • √Γ is a 4th-order tensor with  √Γ_ijkl = (√G_ik √G_jl + √G_il √G_jk)/√2,
#     where √G is the symmetric matrix square root (√G √Gᵀ = G); Γ = √Γ:√Γ.
#   • dB_G is a SYMMETRIC matrix-valued Brownian motion with covariance heuristic
#     (dB_G)_ij (dB_G)_kl = (δ_ik δ_jl + δ_il δ_jk)/2 dt
#     i.e. diagonal entries have variance dt, off-diagonals variance dt/2.
#   • The double contraction (√Γ : dB_G)_ij = Σ_kl √Γ_ijkl (dB_G)_kl is the only
#     subtle part. Because dB_G is symmetric and √G is symmetric, the two tensor
#     products (otimesu, otimesl in the repo's Tensors.jl version) collapse to
#     LEFT and RIGHT matrix multiplication by √G:
#
#         √Γ : dB_G  =  √2 · √G · dB_G · √G .
#
#     (verified: this reproduces Cov(dG_ij,dG_kl)/dt = (v/n)(G_ik G_jl + G_il G_jk).)
#
# So one Euler-Maruyama step is
#     G ← G - (v/n) G Δt + √(v/n)·√(2Δt) · √G X √G ,
# with X symmetric, X_ii ~ N(0,1), X_ij ~ N(0,1/2). Symmetrize and project to the
# PD cone each step for numerical safety. This matches drift𝐆 in solvers/diy.
#
# Drift erodes overall genetic variance (G shrinks, the deterministic -(v/n)G
# term) while simultaneously driving the genetic correlation to ±1 (the noise),
# so each replicate's ellipse shrinks AND collapses onto a line — the orientation
# change the normalized (correlation) ellipses make plain.
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.animation import FuncAnimation, PillowWriter

def sqrtm_sym(G):
    w, V = np.linalg.eigh(G); w = np.clip(w, 1e-12, None)
    return (V * np.sqrt(w)) @ V.T

def simulate(G0, v, n0, T, N, R, seed=0):
    rng = np.random.default_rng(seed); dt = T / N; delta = v / n0
    hist = np.zeros((R, N + 1, 2, 2))
    for r in range(R):
        G = G0.copy(); hist[r, 0] = G
        for k in range(1, N + 1):
            sG = sqrtm_sym(G)
            a, b = rng.standard_normal(), rng.standard_normal()
            c = rng.standard_normal() / np.sqrt(2)        # off-diagonal: variance 1/2
            X = np.array([[a, c], [c, b]])                # symmetric noise matrix
            stoch = np.sqrt(2) * (sG @ X @ sG)            # = √Γ : dB_G  (per √dt)
            G = G - delta * G * dt + np.sqrt(delta) * np.sqrt(dt) * stoch
            G = 0.5 * (G + G.T)                            # symmetrize
            w, V = np.linalg.eigh(G); G = (V * np.clip(w, 1e-9, None)) @ V.T  # keep PD
            hist[r, k] = G
    return hist

def ellipse_axes(M):
    w, V = np.linalg.eigh(M); w = np.clip(w, 0, None)
    angle = np.degrees(np.arctan2(V[1, 1], V[0, 1]))
    return 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]), angle      # full width, height, angle

def correlation(G):
    d = np.sqrt(np.diag(G)); return G / np.outer(d, d)

# ---- run: 3 replicates, drift only, start at G = I ----
hist = simulate(np.eye(2), v=1.0, n0=80.0, T=150.0, N=1500, R=3, seed=7)
N = hist.shape[1] - 1
frames = np.linspace(0, N, 100).astype(int)

# ---- styling ----
BG, TXT, SUBTLE = "#faf4ed", "#575279", "#dfdad9"
REPCOL = ["#286983", "#b4637a", "#907aa9"]
plt.rcParams.update({"axes.edgecolor": SUBTLE, "text.color": TXT,
                     "figure.facecolor": BG, "axes.facecolor": BG,
                     "axes.labelcolor": TXT})
fig, axes = plt.subplots(2, 3, figsize=(9, 6))
LIM = 1.6
arts = []
for col in range(3):
    for row in range(2):
        ax = axes[row, col]
        ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.add_patch(Ellipse((0, 0), 2, 2, fill=False, ec=SUBTLE, lw=1, ls=(0, (4, 4))))
        el = Ellipse((0, 0), 1, 1, fc=REPCOL[col], ec=REPCOL[col], alpha=0.35, lw=2)
        ax.add_patch(el); arts.append(el)
axes[0, 0].set_ylabel("G-matrix", fontsize=12)
axes[1, 0].set_ylabel("correlation matrix", fontsize=12)
fig.tight_layout()

def update(fi):
    k = frames[fi]; idx = 0
    for col in range(3):
        G = hist[col, k]
        for M in (G, correlation(G)):
            w, h, a = ellipse_axes(M)
            arts[idx].width, arts[idx].height, arts[idx].angle = w, h, a
            idx += 1
    return arts

anim = FuncAnimation(fig, update, frames=len(frames), interval=80, blit=False)
anim.save("G_drift_dynamics.gif", writer=PillowWriter(fps=14))
print("saved G_drift_dynamics.gif")
