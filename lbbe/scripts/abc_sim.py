"""
abc_sim.py -- a small, vectorised port of the simulation in
ABC_Coevo_3_1_19.cpp (Nuismer & Week 2019, PLoS Comput. Biol. 15:e1006988),
used here only to illustrate the method.

Life cycle, as in the C++ source, on the mean traits X, Y of two species in
each of `npops` populations:

  selection   W_x(x, y) = exp(-g_x (x - theta_x)^2) / (1 + exp(-A_x (x - y)))
              W_y(y, x) = exp(-g_y (y - theta_y)^2) / (1 + exp(-A_y (y - x)))
              averaged over Gaussian phenotype distributions (variances V_x,
              V_y) by numerical integration over +-5 s.d.
  drift       + Normal(0, V / N)
  movement    island model, rates m_x, m_y
  mating      X <- h X' + (1 - h) X

What this allows that the likelihood of Week & Nuismer (2019) does not:
selection of any strength (the logistic interaction term is not linearised),
gene flow among populations (m > 0), and abiotic optima theta that differ
among populations.

Summary statistics, as in the source: the mean of each species' population
means, the standard deviation of each, and their correlation.  A draw is kept
if all five fall within the source's thresholds of the data:
0.05 + 0.1 |mean|, 0.05 + 0.2 sd, 0.05 + 0.2 |rho|.

Simplifications for the illustration: a fixed number of generations instead
of the source's equilibrium test, a coarser integration grid, and only the
two strengths of coevolutionary selection (A_x, A_y) are inferred; the
background parameters are fixed at the modes of the source's priors.
"""
import numpy as np

BG = dict(Vx=0.71, Vy=3.72, hx=0.412, hy=0.7375, gx=0.05, gy=0.05,
          Nx=30870.0, Ny=1818.0, mx=3.07e-6, my=1.0e-6,
          ETx=5.911818, ETy=6.248235, VTx=0.15, VTy=1.27, RoT=0.0)


def simulate(Ax, Ay, rng, npops=13, gens=200, ngrid=30, bg=BG):
    """Ax, Ay: arrays of length B.  Returns the five summary statistics, B x 5,
    and the final population means X, Y (B x npops)."""
    Ax = np.atleast_1d(np.asarray(Ax, float))[:, None]
    Ay = np.atleast_1d(np.asarray(Ay, float))[:, None]
    B = Ax.shape[0]
    zx, zy = rng.normal(size=(B, npops)), rng.normal(size=(B, npops))
    Tx = bg["ETx"] + np.sqrt(bg["VTx"]) * zx
    Ty = bg["ETy"] + np.sqrt(bg["VTy"]) * (zx * bg["RoT"] + zy * np.sqrt(1 - bg["RoT"] ** 2))
    X = Tx + rng.random((B, npops))
    Y = Ty + rng.random((B, npops))
    u = (np.arange(ngrid) + 0.5) / ngrid * 10 - 5            # midpoints over +-5 s.d.
    f = np.exp(-u ** 2 / 2)
    f /= f.sum()
    sx, sy = np.sqrt(bg["Vx"]), np.sqrt(bg["Vy"])
    for _ in range(gens):
        x = X[..., None] + sx * u                             # B, P, n
        y = Y[..., None] + sy * u
        d = x[..., :, None] - y[..., None, :]                 # B, P, n, n
        sel_x = np.exp(-bg["gx"] * (x - Tx[..., None]) ** 2)
        sel_y = np.exp(-bg["gy"] * (y - Ty[..., None]) ** 2)
        EWx = sel_x * ((1.0 / (1.0 + np.exp(-Ax[..., None, None] * d))) @ f)
        EWy = sel_y * np.einsum("bpij,i->bpj", 1.0 / (1.0 + np.exp(Ay[..., None, None] * d)), f)
        Xp = (x * f * EWx).sum(-1) / (f * EWx).sum(-1)
        Yp = (y * f * EWy).sum(-1) / (f * EWy).sum(-1)
        Xp = Xp + rng.normal(0, np.sqrt(bg["Vx"] / bg["Nx"]), Xp.shape)
        Yp = Yp + rng.normal(0, np.sqrt(bg["Vy"] / bg["Ny"]), Yp.shape)
        Xp = (1 - bg["mx"]) * Xp + bg["mx"] * (Xp.sum(1, keepdims=True) - Xp) / (npops - 1)
        Yp = (1 - bg["my"]) * Yp + bg["my"] * (Yp.sum(1, keepdims=True) - Yp) / (npops - 1)
        X = bg["hx"] * Xp + (1 - bg["hx"]) * X
        Y = bg["hy"] * Yp + (1 - bg["hy"]) * Y
    return summaries(X, Y), X, Y


def summaries(X, Y):
    mx_, my_ = X.mean(1), Y.mean(1)
    sx_, sy_ = X.std(1, ddof=1), Y.std(1, ddof=1)
    ro = ((X - mx_[:, None]) * (Y - my_[:, None])).sum(1) / (X.shape[1] - 1) / (sx_ * sy_)
    return np.c_[mx_, my_, sx_, sy_, ro]


def accept(S, data):
    thr = np.array([0.05 + 0.1 * abs(data[0]), 0.05 + 0.1 * abs(data[1]),
                    0.05 + 0.2 * data[2], 0.05 + 0.2 * data[3],
                    0.05 + 0.2 * abs(data[4])])
    return (np.abs(S - data) <= thr).all(1)


def run_abc(true=(1.5, 0.5), n_sims=4000, chunk=250, seed=3, **kw):
    rng = np.random.default_rng(seed)
    data, X, Y = simulate([true[0]], [true[1]], rng, **kw)
    data = data[0]
    Ax, Ay = rng.uniform(0, 3, n_sims), rng.uniform(0, 3, n_sims)   # the source's priors
    S = np.vstack([simulate(Ax[i:i + chunk], Ay[i:i + chunk], rng, **kw)[0]
                   for i in range(0, n_sims, chunk)])
    keep = accept(S, data)
    return dict(true=np.array(true), data=data, X=X[0], Y=Y[0], Ax=Ax, Ay=Ay, S=S,
                keep=keep)
