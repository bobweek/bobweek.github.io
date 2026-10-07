"""
variance_blocks.py -- one way of drawing the three additive components of
Week et al. (2025) wherever they appear:

  G_A   host genes                     solid teal
  M_A   microbes                       ochre, hatched
  C_A   their covariance               teal and ochre stripes interleaved:
        host-genetic and microbial effects occurring together non-randomly
        (a covariance, not a gene-by-microbe interaction)
"""
from matplotlib.patches import Rectangle

from semstyle import INK, OCHRE, OCHRE_T, PAPER, TEAL, lw


def variance_block(cv, x, y, w, h, kind, z=2):
    if kind == "g":
        cv.add_patch(Rectangle((x, y), w, h, facecolor=TEAL, edgecolor="none", zorder=z))
    elif kind == "m":
        cv.add_patch(Rectangle((x, y), w, h, facecolor=PAPER, edgecolor=OCHRE,
                               hatch="////", linewidth=lw(1.2), zorder=z))
    else:
        n = max(4, int(round(w / 9)))
        for i in range(n):
            cv.add_patch(Rectangle((x + i * w / n, y), w / n, h,
                                   facecolor=TEAL if i % 2 == 0 else OCHRE_T,
                                   edgecolor="none", zorder=z))
            if i % 2:
                cv.add_patch(Rectangle((x + i * w / n, y), w / n, h, facecolor="none",
                                       edgecolor=OCHRE, hatch="////", linewidth=0,
                                       zorder=z))
    cv.add_patch(Rectangle((x, y), w, h, facecolor="none", edgecolor=INK,
                           linewidth=lw(1.2), zorder=z + 1))
