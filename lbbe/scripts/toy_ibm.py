"""
toy_ibm.py -- a toy individual-based simulation of a predator and its prey,
used only to illustrate "events among individuals -> a stochastic process for
the populations" on slides 4 and 18.  It is not a model from any of the papers.

Individuals live on a plane (a unit square with wrap-around edges, so the
frame is a window, not a boundary).  Each has a position and a one-dimensional
trait.

  environment  a smooth field of habitat quality that keeps changing: a sum of
               a few long waves whose phases drift at different speeds.  Where
               quality is low prey die faster, so sparse areas and the edges
               of dense ones move about through time.
  prey         born at rate b; die at rate  c * (prey nearby)
                                          + a * (predators nearby) * match
                                          + poor habitat + s_h * trait^2
  predator     born at rate eps * a * (prey nearby) * match;
               die at rate m + s_p * trait^2
  match        the chance that a predator with trait x catches a prey with
               trait y falls off as exp(-g (x - y)^2).  Each individual's rate
               uses that kernel averaged over the other species' current trait
               distribution (taken as normal), so prey deaths and predator
               births balance.  Prey are selected away from the predators'
               mean and predators toward the prey's; with weak stabilising
               selection holding both near zero the chase cannot run away: the
               prey mean swings from side to side and the predator mean
               follows (Red Queen cycles), kept going by drift and by the
               cycling abundances.
  inheritance  an offspring's trait is the average of its parent's and a
               random mate's, plus a normal deviation, which keeps each
               species' trait distribution single-peaked.

"Nearby" is the count in the individual's cell of a grid.  Offspring settle close to
the parent; everyone also moves a little each step.
"""
import numpy as np


class Environment:
    """Habitat quality in [0, 1] at position xy and time t."""

    def __init__(self, seed=5, kmax=2, speed=1.0, cover=0.62, sharp=0.22):
        rng = np.random.default_rng(seed)
        ks = [(i, j) for i in range(-kmax, kmax + 1) for j in range(0, kmax + 1)
              if (i, j) != (0, 0) and not (j == 0 and i < 0)]
        self.k = np.array(ks, float)
        self.amp = rng.normal(size=len(ks)) / np.hypot(self.k[:, 0], self.k[:, 1]) ** 1.3
        self.phi = rng.uniform(0, 2 * np.pi, len(ks))
        self.om = speed * rng.choice([-1, 1], len(ks)) * rng.uniform(0.05, 0.16, len(ks))
        g = (np.arange(64) + 0.5) / 64
        xx, yy = np.meshgrid(g, g, indexing="ij")
        z = self.raw(np.c_[xx.ravel(), yy.ravel()], 0.0)
        self.cut, self.sd, self.sharp = np.quantile(z, 1 - cover), z.std(), sharp

    def raw(self, xy, t):
        ph = 2 * np.pi * (xy @ self.k.T) + self.phi + self.om * t
        return np.cos(ph) @ self.amp

    def quality(self, xy, t):
        return 1.0 / (1.0 + np.exp(-(self.raw(xy, t) - self.cut) / (self.sharp * self.sd)))


def simulate(T=150.0, dt=0.04, seed=3, G=14, n_h0=6000, n_p0=1800,
             g=20.0, mu=0.20, s_h=0.44, s_p=0.06, verbose=False):
    rng = np.random.default_rng(seed)
    b, c, a, eps, m = 1.0, 0.0012, 0.060, 0.33, 0.6
    move, settle, bad = 0.022, 0.022, 2.6
    env = Environment()

    def scatter(n):
        out = np.empty((0, 2))
        while len(out) < n:
            xy = rng.uniform(0, 1, (2 * n, 2))
            out = np.vstack([out, xy[rng.random(2 * n) < env.quality(xy, 0.0)]])
        return out[:n]

    H = dict(xy=scatter(n_h0), z=rng.normal(0.30, 0.20, n_h0))
    P = dict(xy=scatter(n_p0), z=rng.normal(-0.10, 0.20, n_p0))
    cell = lambda xy: ((xy[:, 0] * G).astype(int) % G) * G + ((xy[:, 1] * G).astype(int) % G)
    out = dict(t=[], NH=[], NP=[], mH=[], mP=[], sH=[], sP=[], frames=[], env=env)
    n_steps = int(round(T / dt))
    keep_every = max(1, n_steps // 480)
    frame_every = keep_every * 8                     # snapshots are heavy: keep 60
    for k in range(n_steps + 1):
        t = k * dt
        if k % keep_every == 0:
            out["t"].append(t)
            for key, S in (("H", H), ("P", P)):
                out["N" + key].append(len(S["z"]))
                out["m" + key].append(S["z"].mean() if len(S["z"]) else np.nan)
                out["s" + key].append(S["z"].std() if len(S["z"]) else np.nan)
        if k % frame_every == 0:
            out["frames"].append((t, H["xy"].astype(np.float32).copy(), H["z"].astype(np.float32).copy(),
                                  P["xy"].astype(np.float32).copy(), P["z"].astype(np.float32).copy()))
        if len(H["z"]) == 0 or len(P["z"]) < 20 or len(H["z"]) > 90000:
            break
        cH, cP = cell(H["xy"]), cell(P["xy"])
        nH = np.bincount(cH, minlength=G * G)
        nP = np.bincount(cP, minlength=G * G)
        zH, zP, vH, vP = H["z"].mean(), P["z"].mean(), H["z"].var(), P["z"].var()
        kP, kH = 1 + 2 * g * vP, 1 + 2 * g * vH          # the kernel, averaged over the other species
        dH = (c * nH[cH] + a * nP[cH] * np.exp(-g * (H["z"] - zP) ** 2 / kP) / np.sqrt(kP)
              + bad * (1 - env.quality(H["xy"], t)) + s_h * H["z"] ** 2)
        bP = eps * a * nH[cP] * np.exp(-g * (P["z"] - zH) ** 2 / kH) / np.sqrt(kH)
        dP = m + s_p * P["z"] ** 2
        new = {}
        for S, birth, death in ((H, b, dH), (P, bP, dP)):
            n = len(S["z"])
            born = rng.random(n) < birth * dt
            dies = rng.random(n) < death * dt
            kid_xy = S["xy"][born] + rng.normal(0, settle, (born.sum(), 2))
            mate = S["z"][rng.integers(0, n, born.sum())]
            kid_z = 0.5 * (S["z"][born] + mate) + rng.normal(0, mu, born.sum())
            xy = np.vstack([S["xy"][~dies], kid_xy])
            xy = (xy + rng.normal(0, move * np.sqrt(dt), xy.shape)) % 1.0
            new[id(S)] = (xy, np.r_[S["z"][~dies], kid_z])
        H["xy"], H["z"] = new[id(H)]
        P["xy"], P["z"] = new[id(P)]
    for key in ("t", "NH", "NP", "mH", "mP", "sH", "sP"):
        out[key] = np.array(out[key])
    return out


def _report(o, tag=""):
    t, zh, zp = o["t"], o["mH"], o["mP"]
    cross = int(np.sum(np.diff(np.sign(zh - zp)) != 0))
    flips = int(np.sum(np.diff(np.sign(zh)) != 0))
    print(f"{tag} T={t[-1]:.0f} prey {o['NH'].min()}..{o['NH'].max()} pred {o['NP'].min()}..{o['NP'].max()} "
          f"zH [{zh.min():+.2f},{zh.max():+.2f}] sd {np.nanmean(o['sH']):.2f}; zP [{zp.min():+.2f},{zp.max():+.2f}]; "
          f"sign flips of zH: {flips}; crossings zH-zP: {cross}")


if __name__ == "__main__":
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    grid = [dict(), dict(seed=8)]
    fig, axs = plt.subplots(len(grid), 2, figsize=(13, 2.4 * len(grid)))
    for r, kw in enumerate(grid):
        t0 = time.time()
        o = simulate(**kw)
        _report(o, f"{kw} ({time.time() - t0:.1f}s)")
        axs[r, 0].plot(o["t"], o["NH"], color="#93B200"); axs[r, 0].plot(o["t"], o["NP"], color="#5B36C4")
        axs[r, 1].plot(o["t"], o["mH"], color="#93B200"); axs[r, 1].plot(o["t"], o["mP"], color="#5B36C4")
        axs[r, 1].axhline(0, color="k", lw=0.5); axs[r, 0].set_title(str(kw), fontsize=8)
    fig.savefig("/home/claude/work/sheets/ibm_scan.png", dpi=70, bbox_inches="tight")
    fr = o["frames"]
    ks = [int(len(fr) * f) for f in (0.02, 0.2, 0.4, 0.6, 0.8, 0.99)]
    fig, axs = plt.subplots(1, 6, figsize=(18, 3.2))
    for ax, k in zip(axs, ks):
        t_, hxy, hz, pxy, pz = fr[k]
        ax.scatter(hxy[:, 0], hxy[:, 1], s=0.5, c="#B5C400", lw=0)
        ax.scatter(pxy[:, 0], pxy[:, 1], s=0.5, c="#5B4BC4", lw=0)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        ax.set_title(f"t={t_:.0f}  H={len(hz)} P={len(pz)}", fontsize=9)
    fig.savefig("/home/claude/work/sheets/ibm_test.png", dpi=80, bbox_inches="tight")
