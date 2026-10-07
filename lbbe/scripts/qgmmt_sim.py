"""
qgmmt_sim.py -- the simulation model of Week et al. (2025) Evolution
79:2487-2502, ported from fcts_strcts.jl, figs34.jl and fig5.jl of the qgmmt
repository.

Each host carries L diploid loci and three sets of microbial taxa:
  lineal      offspring abundance ~ Poisson(mid-parent abundance)
  non-lineal  offspring abundance ~ Poisson(abundance in one donor host); the
              donor is drawn from the selected parents (Nsel = True, "post-
              selected") or from the whole parental population (Nsel = False,
              "pre-selected")
  novel       drawn afresh every generation, ~ Poisson(K)
The trait is additive in allele counts and abundances plus noise; fitness is
exp(s z); parents are sampled in proportion to fitness.

Differences from the Julia code, none of which change the model:
  * vectorised over hosts;
  * the additive variances A used for prediction are computed as the variance
    of the corresponding additive value.  The Julia code obtains them by
    regressing that value on allele counts and abundances, which returns the
    same number because the value is an exact linear function of them;
  * fewer replicates by default (see the figure script), cached in data/.
The selection strength is the repository's s = 0.1, which reproduces the scale
of the published Fig. 3 (the paper's text quotes 1e-3).
"""
from dataclasses import dataclass, replace

import numpy as np


@dataclass
class Par:
    EP: float = 100.0      # expected total phenotypic variance
    h2: float = 0.2        # shares of it: host genes,
    l2: float = 0.2        #   lineal microbes,
    nu2: float = 0.2       #   non-lineal microbes,
    e2: float = 0.2        #   novel (environmental) microbes
    p: float = 0.5
    L: int = 100
    SL: int = 100
    SN: int = 100
    SE: int = 100
    n: int = 1000
    s: float = 0.1
    K: float = 100.0
    Nsel: bool = True


def init(par, rng):
    EE = max(0.0, (1 - par.h2 - par.l2 - par.nu2 - par.e2) * par.EP)
    vG = par.h2 * par.EP / (2 * par.p * (1 - par.p))
    pp = dict(
        gamma=rng.normal(0, np.sqrt(vG / par.L), par.L),
        wL=rng.normal(0, np.sqrt(par.l2 * par.EP / par.K / par.SL), par.SL),
        wN=rng.normal(0, np.sqrt(par.nu2 * par.EP / par.K / par.SN), par.SN),
        wE=rng.normal(0, np.sqrt(par.e2 * par.EP / par.K / par.SE), par.SE),
        sdE=np.sqrt(EE))
    pd = dict(
        g=rng.random((par.n, par.L, 2)) < par.p,
        mL=rng.poisson(par.K, (par.n, par.SL)),
        mN=rng.poisson(par.K, (par.n, par.SN)),
        mE=rng.poisson(par.K, (par.n, par.SE)))
    pd["z"] = value(pp, pd, "GLNV") + rng.normal(0, pp["sdE"], par.n)
    return pp, pd


def value(pp, pd, factors):
    """Additive value from the named factors: G, L, N, V."""
    v = (pd["g"].sum(2)) @ pp["gamma"]
    if "L" in factors:
        v = v + pd["mL"] @ pp["wL"]
    if "N" in factors:
        v = v + pd["mN"] @ pp["wN"]
    if "V" in factors:
        v = v + pd["mE"] @ pp["wE"]
    return v


def selection(par, pd, rng):
    w = np.exp(par.s * (pd["z"] - pd["z"].max()))
    w /= w.sum()
    pairs = rng.choice(par.n, size=(par.n, 2), p=w)
    W = np.bincount(pairs.ravel(), minlength=par.n) / 2.0     # realised fitness
    return pairs, W


def offspring(par, pp, pd, pairs, rng):
    n, L = par.n, par.L
    donors = rng.choice(pairs.ravel() if par.Nsel else np.arange(n), size=n)
    mN = rng.poisson(pd["mN"][donors])
    fr = rng.random((n, L)) < 0.5
    c1 = rng.integers(0, 2, n)[:, None]
    c2 = rng.integers(0, 2, n)[:, None]
    g1, g2 = pd["g"][pairs[:, 0]], pd["g"][pairs[:, 1]]
    take = lambda g, c: np.take_along_axis(g, np.broadcast_to(c[:, :, None], (n, L, 1)), 2)[:, :, 0]
    ch1 = np.where(fr, take(g1, c1), take(g1, c2))
    ch2 = np.where(fr, take(g2, c2), take(g2, c1))
    new = dict(g=np.stack([ch1, ch2], 2),
               mL=rng.poisson(0.5 * (pd["mL"][pairs[:, 0]] + pd["mL"][pairs[:, 1]])),
               mN=mN, mE=rng.poisson(par.K, (n, par.SE)))
    new["z"] = value(pp, new, "GLNV") + rng.normal(0, pp["sdE"], n)
    return new


# the four nested trait architectures of figs34.jl: each active factor has an
# expected variance of 25
ARCH = {
    "G": dict(EP=25.0, h2=0.99, l2=0.001, nu2=0.001, e2=0.001, SL=1, SN=1, SE=1),
    "GL": dict(EP=50.0, h2=0.49, l2=0.49, nu2=0.001, e2=0.001, SN=1, SE=1),
    "GLN": dict(EP=75.0, h2=0.33, l2=0.33, nu2=0.33, e2=0.001, SE=1),
    "GLNV": dict(EP=100.0, h2=0.25, l2=0.25, nu2=0.25, e2=0.25),
}


def timeseries(arch, Nsel, T=10, reps=20, seed=0):
    """Mean trait over T generations (relative to the first), reps x T."""
    rng = np.random.default_rng(seed)
    par = replace(Par(), Nsel=Nsel, **ARCH[arch])
    out = np.empty((reps, T))
    for r in range(reps):
        pp, pd = init(par, rng)
        zb = [pd["z"].mean()]
        for _ in range(T - 1):
            pairs, _ = selection(par, pd, rng)
            pd = offspring(par, pp, pd, pairs, rng)
            zb.append(pd["z"].mean())
        out[r] = np.array(zb) - zb[0]
    return out


def prediction_runs(Nsel, nr=60, sr=6, orr=6, seed=0):
    """Observed one-generation change against A*beta for the four nested
    additive variances (fig5.jl), full GLNV architecture."""
    rng = np.random.default_rng(seed)
    par = replace(Par(), Nsel=Nsel, **ARCH["GLNV"])
    obs = np.empty(nr)
    pred = {k: np.empty(nr) for k in ARCH}
    for i in range(nr):
        pp, pd = init(par, rng)
        A = {k: value(pp, pd, k).var(ddof=1) for k in ARCH}
        zbar, P = pd["z"].mean(), pd["z"].var(ddof=1)
        dz, beta = [], []
        for _ in range(sr):
            pairs, W = selection(par, pd, rng)
            dz.append(np.mean([offspring(par, pp, pd, pairs, rng)["z"].mean()
                               for _ in range(orr)]) - zbar)
            beta.append(np.cov(W, pd["z"])[0, 1] / P)
        obs[i] = np.mean(dz)
        for k in ARCH:
            pred[k][i] = A[k] * np.mean(beta)
    return obs, pred
