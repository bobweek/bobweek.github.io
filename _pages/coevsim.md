---
layout: page
permalink: /coevsim/
title: CoevSim
description:
nav: true
nav_order: 3
---

<center>
<iframe src="https://openprocessing.org/sketch/396443/embed/" width="850" height="700"></iframe>
</center>

This is an old visualizer I wrote of spatially distributed coevolution among three species. Very thankful the [openprocessing.org](https://openprocessing.org/) website is still supported and maintained! In the simulation, the three species engange in a rock-paper-scissors dynamic. Pairwise interaction outcomes are determined by phenotypic matching. This results in cyclical coevolutionary chases where the trait of one species attempts to escape another, and that other is attempting to escape the third, which finally the third attempts to escape the first. The trait value for each species is associated with the red, green, and blue content in each pixel (spatial location) and the community is summarized by the resulting RGB/HEX encoded color. I recommend adjusting relative selection strengths to observe emergent color palettes, and relative dispersal abilities to observe emergent spatial patterns.