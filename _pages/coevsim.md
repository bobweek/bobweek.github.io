---
layout: page
permalink: /coevsim/
title: CoevSim
description:
nav: true
nav_order: 3
---

<script src="https://cdn.jsdelivr.net/npm/p5@latest/lib/p5.min.js"></script>

<div style="display:flex;justify-content:center;">
  <iframe
    src="/sims/spatial_coev/"
    width="400"
    height="700"
    style="border:0;display:block;"
    loading="lazy"
    allowfullscreen>
  </iframe>
</div>

---

A visualizer of spatially distributed coevolution among three species. The
three engage in a rock-paper-scissors dynamic: pairwise interaction outcomes
are determined by phenotypic matching, so each species tracks the next around
the cycle while fleeing the previous one, and none can win. The trait value of
each species is the red, green, and blue content of a pixel, and the community
at each spatial location is summarized by the resulting RGB color.

Adjust relative selection strengths to see emergent color palettes, and
relative dispersal abilities to see emergent spatial patterns. Two further
forces are off by default: stabilizing selection pulls each species toward its
own spatial optimum (switch the view to Optima to see where those sit), and
drift perturbs each population by sqrt(G/Ne) per generation, which can be made
to strengthen toward the range margin where Ne is small.

Controls:

- Sliders for dispersal, coevolution, stabilizing selection, and drift
- Tap/click the lattice, or the play button, to pause and resume
- **Step** advances a single generation
- **-/+** coarsen or refine the lattice
- **Random**, **At optima**, and **Uniform** reseed the field
- Press **h** to hide or show the sliders

Keyboard shortcuts (space, s, r, h, [ , ]) need the frame focused — click it once first.