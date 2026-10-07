"""
Slide 05 -- a single-species aside: drift reorients G.

Week (2026) J. Theor. Biol. 625:112428.  Sample paths of the paper's
Brownian-motion-gradient G equation (Eq. 21 / 23 with M = 0, no selection,
fixed n):
    dG = -(v/n) G dt + sqrt(v/n) * sqrt(2) * A dB A^T,   A A^T = G,
dB symmetric with Var(dB_ii) = dt, Var(dB_ij) = dt/2 (Eq. 22).
Time is in units of n/v (= N_e when v = 1).  The build checks the simulated
correlations against Eq. 42.  The deterministic expectation G_0 exp(-vt/n)
keeps the correlation fixed (Section 3.2).

Simulation code and settings are those of rev15/figure2_results.py, panel A.
Ellipse *shape* is the simulated G of the highlighted path; ellipse *size* is
schematic (SHRINK): true size would be 3 % by t = 7 and unreadable.
"""
import numpy as np
from matplotlib.patches import Ellipse

from semstyle import *  # noqa: F403

SHRINK = 0.114


def simulate_G(n_paths=40, rho0=0.30, vn=1e-3, dt=0.5, tau_max=8.0, seed=21):
    rng = np.random.default_rng(seed)
    n = int(round(tau_max / (vn * dt)))
    G = np.tile(np.array([[1.0, rho0], [rho0, 1.0]]), (n_paths, 1, 1))
    Gs = np.empty((n + 1, n_paths, 2, 2))
    Gs[0] = G
    s = np.sqrt(vn * dt)
    for k in range(n):
        A = np.linalg.cholesky(G)
        dB = rng.normal(size=(n_paths, 2, 2))
        off = dB[:, 0, 1] / np.sqrt(2.0)
        dB[:, 0, 1] = off
        dB[:, 1, 0] = off
        G = G - vn * G * dt + s * np.sqrt(2.0) * A @ dB @ np.transpose(A, (0, 2, 1))
        G = 0.5 * (G + np.transpose(G, (0, 2, 1)))
        Gs[k + 1] = G
    tau = np.arange(n + 1) * vn * dt
    rho = Gs[:, :, 0, 1] / np.sqrt(Gs[:, :, 0, 0] * Gs[:, :, 1, 1])
    return tau, Gs, rho


def check_against_eq42(rho0=0.30, vn=1e-3, dt=0.5, tau_max=8.0):
    rng = np.random.default_rng(11)
    u = np.full(20000, np.arctanh(rho0))
    for _ in range(int(round(tau_max / (vn * dt)))):
        u += 0.5 * vn * np.tanh(u) * dt + np.sqrt(vn * dt) * rng.normal(size=u.size)
    r42 = np.tanh(u)
    _, _, rho_big = simulate_G(n_paths=2000, seed=99)
    a, b = np.mean(abs(rho_big[-1]) > .9), np.mean(abs(r42) > .9)
    print(f"  [check] P(|rho|>0.9) at tau={tau_max}: G-eq {a:.3f} vs Eq.42 {b:.3f}")
    assert abs(a - b) < 0.03, "simulated G disagrees with Eq. 42"


def draw_cov(ax, G, scale=1.0, **kw):
    G = G / np.trace(G)
    vals, vecs = np.linalg.eigh(G)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    w, h = 2 * 0.95 * scale * np.sqrt(np.maximum(vals, 1e-4))
    ax.add_patch(Ellipse((0, 0), w, h, angle=ang, clip_on=False, **kw))


def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    tau, Gs, rho = simulate_G()
    check_against_eq42()
    dtau = tau[1] - tau[0]
    k15, k3 = int(round(1.5 / dtau)), int(round(3.0 / dtau))
    hi = next(j for j in range(rho.shape[1])
              if rho[k15, j] > -0.1 and -0.75 < rho[k3, j] < -0.3 and rho[-1, j] < -0.9)

    X0, W, Y0, Hh = 124, 1000, 88, 248
    bx = px_axes(fig, X0, Y0, W, Hh)
    th = slice(None, None, 4)
    for j in range(rho.shape[1]):
        if j != hi:
            bx.plot(tau[th], rho[th, j], color=FAINT, lw=lw(1.0), alpha=0.8, zorder=2)
    bx.axhline(rho[0, 0], color=INK, lw=lw(2.4), ls=(0, (3.2, 1.8)), zorder=4)
    bx.plot(tau[th], rho[th, hi], color=CRIMSON, lw=lw(3.4), zorder=5)
    bx.set_xlim(0, tau[-1])
    bx.set_ylim(-1.06, 1.06)
    bx.set_yticks([-1, 0, 1])
    bx.set_xticks([0, 2, 4, 6, 8])
    bx.set_xlabel(r"time (units of $N_e$)")
    bx.set_ylabel("genetic correlation")
    box = dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.88)
    t = bx.text(7.9, rho[0, 0] + 0.07, "the deterministic expectation", ha="right",
                va="bottom", fontsize=fs(FS_SMALL), color=INK, bbox=box)
    t.set_zorder(6)
    t = bx.text(5.3, -0.62, "one population", ha="left", va="center",
                fontsize=fs(FS_LABEL), color=CRIMSON, fontweight="bold", bbox=box)
    t.set_zorder(6)
    t = bx.text(7.9, 0.70, "39 others, same start", ha="right", va="center",
                fontsize=fs(FS_SMALL), color=MUTED, bbox=box)
    t.set_zorder(6)

    # G of the highlighted population at three times, above their moments
    cy, size_px = 448, 132
    for ts in (0.0, 3.0, 7.0):
        k = int(round(ts / dtau))
        xc = X0 + W * ts / tau[-1]
        ea = px_axes(fig, xc - size_px / 2, cy - size_px / 2, size_px, size_px)
        ea.set_xlim(-1, 1)
        ea.set_ylim(-1, 1)
        ea.set_aspect("equal")
        blank(ea)
        size = np.exp(-SHRINK * ts)
        ea.plot([-0.95, 0.95], [0, 0], color=LINE, lw=lw(1.0), zorder=0)
        ea.plot([0, 0], [-0.95, 0.95], color=LINE, lw=lw(1.0), zorder=0)
        if ts > 0:
            draw_cov(ea, Gs[0, hi], scale=size, facecolor="none", edgecolor=INK,
                     lw=lw(2.0), ls=(0, (3.0, 1.6)))
        draw_cov(ea, Gs[k, hi], scale=size, facecolor=CRIMSON_T, edgecolor=CRIMSON,
                 lw=lw(3.0), alpha=0.95)
        cv.plot([xc, xc], [Y0 + Hh + 4, cy - size_px / 2 + 6], color=FAINT,
                lw=lw(1.2), ls=(0, (1, 2.2)))
    cv.text(246, 470, r"$\mathbf{G}$ of that population", ha="left", va="center",
            fontsize=fs(FS_LABEL), color=CRIMSON, fontweight="bold")
    cv.text(246, 440, "trait 1 against trait 2", ha="left", va="center",
            fontsize=fs(FS_SMALL), color=MUTED)
    cv.text(600, 470, "dashed: the expectation", ha="left", va="center",
            fontsize=fs(FS_LABEL), color=INK)
    cv.text(600, 440, "smaller, but the same shape", ha="left", va="center",
            fontsize=fs(FS_SMALL), color=MUTED)
    status(cv, "published",
           "Week 2026, J. Theor. Biol. 625:112428 \u2014 Eq. 21, drift alone; ellipse shape simulated, ellipse size schematic")
    save(fig, "s05_drift_G")


if __name__ == "__main__":
    main()
