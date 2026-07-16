from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import matplotlib.pyplot as plt

plt.ioff()

def style_axes(ax, axis_color="black"):
    # no tick marks or labels
    ax.set_xticks([])
    ax.set_yticks([])

    # show only left and bottom axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    ax.spines["left"].set_linewidth(1.6)
    ax.spines["bottom"].set_linewidth(1.6)
    ax.spines["left"].set_color(axis_color)
    ax.spines["bottom"].set_color(axis_color)

    ax.set_facecolor("none")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

out = Path("scatter-icons-side-by-side.png")

fig, axes = plt.subplots(1, 2, figsize=(4.8, 1.9))
fig.patch.set_alpha(0)

# --------------------------
# left: positive relationship
# --------------------------
np.random.seed(4)

n = 34
x1 = np.linspace(0.05, 0.95, n)
y1 = 0.18 + 0.68 * x1 + np.random.normal(scale=0.15, size=n)
y1 = np.clip(y1, 0.03, 0.97)

m1, b1 = np.polyfit(x1, y1, 1)
xx1 = np.linspace(0.03, 0.97, 200)

axes[0].scatter(x1, y1, s=20, color="black", alpha=0.85)
axes[0].plot(xx1, m1 * xx1 + b1, color="black", linewidth=2.0)
style_axes(axes[0], axis_color="black")

# --------------------------
# right: no relationship, but still with fitted line
# --------------------------
pink = "#c2185b"

np.random.seed(7)

n = 34
x2 = np.random.uniform(0.05, 0.95, n)
y2 = np.random.normal(loc=0.5, scale=0.20, size=n)
y2 = np.clip(y2, 0.05, 0.95)

m2, b2 = np.polyfit(x2, y2, 1)
xx2 = np.linspace(0.03, 0.97, 200)

axes[1].scatter(x2, y2, s=20, color=pink, alpha=0.85)
axes[1].plot(xx2, m2 * xx2 + b2, color=pink, linewidth=2.0)
style_axes(axes[1], axis_color=pink)

plt.tight_layout(pad=0.15, w_pad=5.4)

fig.savefig(
    out,
    dpi=220,
    transparent=True,
    bbox_inches="tight",
    pad_inches=0.02,
)

plt.close(fig)

print(f"saved {out.resolve()}")