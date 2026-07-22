"""
Opening-slide visual: individuals vary in offspring number (reproductive variance).

A one-generation reproduction diagram: a row of identical parent individuals, each
connected to the offspring it produced. The number varies -- some leave none, some
leave several -- with no axes, numbers, or distributions. This is the individual-level
picture whose population-level consequence is demographic stochasticity.

Deck palette (daytime.css): white slides, ink #1b1f24, on-brand teal #1b8c98.
Single colour for every individual: the variation is in NUMBER, not type.
Transparent background so it drops onto the white slide.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

INK   = '#1b1f24'
TEAL  = '#1b8c98'   # every individual (parents + offspring): same kind
LINE  = '#b2bcc6'   # parent -> offspring links
LABEL = '#5f6b76'

mpl.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'})

# one clear, punchy spread: two leave none, two leave several
counts = [2, 0, 4, 1, 0, 3, 1]        # 7 parents -> 11 offspring
GAP = 1.15                             # space between families

# ---- lay families out left-to-right, each parent centred over its brood ----
families = []
x = 0.0
for c in counts:
    if c == 0:
        families.append((x + 0.5, [])); x += 1 + GAP
    else:
        block = [x + i + 0.5 for i in range(c)]
        families.append((x + c / 2.0, block)); x += c + GAP
x_right = x - GAP

parent_x   = [f[0] for f in families]
offspring_x = [o for f in families for o in f[1]]
Y_P, Y_O = 1.7, 0.0

fig, ax = plt.subplots(figsize=(10.5, 3.5))

# links first (behind the circles)
for px, block in families:
    for ox in block:
        ax.plot([px, ox], [Y_P, Y_O], color=LINE, lw=1.4, alpha=0.9,
                solid_capstyle='round', zorder=1)

# individuals (scatter keeps them perfectly round at any aspect); parents a touch larger
ax.scatter(offspring_x, [Y_O]*len(offspring_x), s=470, color=TEAL,
           edgecolors='white', linewidths=1.6, zorder=3)
ax.scatter(parent_x, [Y_P]*len(parent_x), s=830, color=TEAL,
           edgecolors='white', linewidths=1.8, zorder=3)

# soft row labels on the left (the only text)
xL = -1.4
ax.text(xL, Y_P, 'parents',   ha='right', va='center', color=LABEL, fontsize=15)
ax.text(xL, Y_O, 'offspring', ha='right', va='center', color=LABEL, fontsize=15)

ax.set_xlim(xL - 3.0, x_right + 0.6)
ax.set_ylim(-0.55, 2.25)
ax.axis('off')
fig.tight_layout(pad=0.3)
fig.savefig('repro_variance.png', dpi=200, transparent=True, bbox_inches='tight')
plt.close(fig)
print(f'saved repro_variance.png : {len(counts)} parents, offspring counts {counts} '
      f'(total {sum(counts)}), mean {np.mean(counts):.2f}, var {np.var(counts):.2f}')

# ---- trim transparent border to a tight crop (+ small margin) ----
from PIL import Image
im = Image.open('repro_variance.png')
a = np.array(im)[..., 3]
ys, xs = np.where(a > 0)
m = 24
box = (max(xs.min()-m,0), max(ys.min()-m,0), min(xs.max()+m, im.width), min(ys.max()+m, im.height))
im.crop(box).save('repro_variance.png')
print('trimmed to', Image.open('repro_variance.png').size)
