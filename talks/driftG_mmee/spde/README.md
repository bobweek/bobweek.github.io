# SPDE trait-density applet

Interactive companion to the branching-diffusion applet. It integrates, in real
time, the **stochastic PDE** that the measure-valued martingale characterization
reduces to in **one trait dimension** (Supplementary Material, Eq. 63):

```
∂ν/∂t = m(ν,z) ν + ½ M ∂²ν/∂z²  +  √(v ν) · ξ(z,t)
m(ν,z) = r − ½ ψ (θ − z)²  − c n        (b = 0: no directional selection)
n = ∫ ν(z) dz
```

`ξ(z,t)` is space–time white noise (the supplement's Ẇₜ(z)); the `√(v ν)`
prefactor is the demographic-stochasticity / drift term. The main panel plots the
current density `ν(z,t)` — abundance density on y, trait value on x — as it settles
under the balance of stabilizing selection ↔ mutation (variance) and intrinsic
growth ↔ competition ↔ the cost of trait variance (abundance). No density history
is kept; only the instantaneous field is shown.

## The Hilbert-space idea

The three check boxes overlay a test function `f(z)` on the density (dashed) and
open a window tracing the associated **inner product** `⟨f, ν⟩ = ∫ f(z) ν(z) dz` —
i.e. the population moment obtained by projecting the measure onto `f`:

| overlay      | window statistic (main-text formula, with `p(z) → ν(z)/n`) |
|--------------|-------------------------------------------------------------|
| `f(z) = 1`   | **n**  = ∫ ν(z) dz                                          |
| `f(z) = z`   | **z̄** = ∫ z · ν(z)/n dz                                     |
| `f(z) = z²`  | **P**  = ∫ (z − z̄)² · ν(z)/n dz                             |

Symbols are set in bold to keep the multivariate reading in view (in `d` dimensions
`z̄` is a vector and `P` a matrix; `n` is a scalar). Each window also draws the
**deterministic equilibrium** as a dashed reference line — `n̂ = 40`, `θ = 0`,
`P̂ = 1` — so you can watch the moments converge to the mutation–selection–
competition balance. (`P` hovers slightly *above* `P̂`: that gap is the genuine
demographic-noise/drift contribution to the variance, not a numerical artefact.)

## Fixed parameters

`r = 1.25   ψ = 0.5   M = 0.5   c = 0.025   θ = 0   v = 0.6`

Deterministic equilibrium: variance `P̂ = √(M/ψ) = 1`, mean `ẑ̄ = θ = 0`,
abundance `n̂ = (r − ½√(ψM))/c = 40`. These are baked in (not sliders), as
requested.

## Controls

- **check boxes** — toggle each overlay + its moment window
- **▶ / ⏸** — play / pause (starts paused; no animation runs in the deck until asked)
- **↺ reset** — restart from the initial condition
- keyboard: `space` play/pause · `r` reset · `1` `2` `3` toggle the three overlays

## Embedding in the Quarto/reveal.js deck

Drop this `spde/` folder next to `branching/` and embed exactly as the branching
applet, e.g. on the **“A Hilbert Space Perspective”** slide:

```markdown
## A Hilbert Space Perspective

<iframe src="spde/index.html" width="980" height="560"
        style="border:0; display:block; margin:0 auto;"
        title="SPDE trait density"></iframe>
```

## Notes on the visualization

- **Overlays are drawn on the shape (test-function) scale, not the density's
  y-axis** — `f(z)` lives in a different space than `ν(z)`. The meaningful features
  are the *flatness* of `f = 1`, the *zero-crossing at θ* of `f = z`, and the
  *minimum at θ* of `f = z²`; the inner product `∫ f ν` is what each window plots.
- **Numerics.** Euler–Maruyama on a 200-cell grid over `z ∈ [−4.5, 4.5]`,
  sub-step `Δt = 1.5e-3`, 35 sub-steps/frame, Neumann (zero-flux) boundaries.
  Space–time white noise discretizes to a per-cell noise s.d. of
  `√(v ν_i / Δz)·√(Δt)`, which reproduces the abundance SDE
  `dn = m̄ n dt + √(v n) dB` (Eq. 18). Validated: the deterministic limit hits
  `P̂, ẑ̄, n̂` to <0.1%, and the stochastic `sd(n)` matches the Ornstein–Uhlenbeck
  prediction `√(v/2c)`.
- Self-contained; only dependency is p5.js 1.x from a CDN (as in the branching applet).
```
