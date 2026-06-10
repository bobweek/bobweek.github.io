from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import matplotlib.pyplot as plt

plt.ioff()

# ------------------------------------------------------------
# Simple host–microbiome rescue model
#
# U(t): hosts without beneficial microbe
# I(t): hosts with beneficial microbe
# N(t) = U(t) + I(t)
#
# Uninfected hosts decline after environmental change.
# Microbe-carrying hosts have improved growth.
# The microbe is introduced at t_intro and also spreads
# horizontally from I to U.
# ------------------------------------------------------------

# output file
out = Path("microbial_rescue_curve.png")

# --------------------------
# parameters
# --------------------------
T = 120.0         # total time
dt = 0.05         # time step
times = np.arange(0, T + dt, dt)

K = 1.0           # carrying capacity (scaled)

# growth rates under environmental stress
r_U = -0.20
r_I = 0.08

# horizontal spread of beneficial microbe
beta = 1.0

# introduction event
t_intro = 8.0
I_intro = 0.005   # small microbe-carrying introduction pulse

# initial conditions
U = np.zeros_like(times)
I = np.zeros_like(times)
U[0] = 0.82
I[0] = 0.0

# --------------------------
# simulate by Euler method
# --------------------------
introduced = False

for k in range(len(times) - 1):
    t = times[k]

    # one-time introduction of beneficial microbe
    if (not introduced) and (t >= t_intro):
        I[k] += I_intro
        introduced = True

    N = U[k] + I[k]
    N_safe = max(N, 1e-12)

    # density-regulated growth
    grow_U = r_U * U[k] * (1 - N / K)
    grow_I = r_I * I[k] * (1 - N / K)

    # horizontal spread of beneficial microbe
    trans = beta * U[k] * I[k] / N_safe

    dU = grow_U - trans
    dI = grow_I + trans

    U[k + 1] = max(U[k] + dt * dU, 0.0)
    I[k + 1] = max(I[k] + dt * dI, 0.0)

# handle case where intro lands exactly on final step
if (not introduced) and (times[-1] >= t_intro):
    I[-1] += I_intro

N = U + I

# --------------------------
# identify "rescue" point
# --------------------------
# approximate rescue turnaround as the minimum of total population
imin = np.argmin(N)
t_min = times[imin]
N_min = N[imin]

# first time after the minimum that total population clearly rises
# and exceeds the minimum by a small amount
post = np.where((times > t_min) & (N > N_min + 0.03))[0]
irescue = post[0] if len(post) > 0 else imin
t_rescue = times[irescue]
N_rescue = N[irescue]

# --------------------------
# plot
# --------------------------
fig, ax = plt.subplots(figsize=(8.2, 4.8))

# total population = main rescue curve
ax.plot(times, N, linewidth=3.4, color="#31748f", label="Total host population", zorder=4)

# optional component curves
ax.plot(times, U, linewidth=2.0, color="#6e6a86", alpha=0.95, linestyle="--",
        label="Without beneficial microbe", zorder=2)
ax.plot(times, I, linewidth=2.0, color="#c2185b", alpha=0.95,
        label="With beneficial microbe", zorder=3)

# mark introduction time
N_intro = np.interp(t_intro, times, N)
ax.scatter([t_intro], [N_intro], s=60, color="#c2185b", zorder=5)

# ax.annotate(
#     "",
#     xy=(t_intro, N_intro),
#     xytext=(t_intro - 22, N_intro + 0.18),
#     arrowprops=dict(arrowstyle="->", lw=1.8, color="#c2185b"),
#     fontsize=11,
#     color="#8e1244"
# )

# mark rescue turnaround
# ax.scatter([t_min], [N_min], s=55, color="#f6c177", edgecolor="black", zorder=6)
# ax.annotate(
#     "population reaches minimum\nthen begins recovery",
#     xy=(t_min, N_min),
#     xytext=(t_min + 8, N_min - 0.10),
#     arrowprops=dict(arrowstyle="->", lw=1.5, color="#f6c177"),
#     fontsize=10,
#     color="#7a5a00"
# )

# subtle marker for "rescue underway"
# ax.scatter([t_rescue], [N_rescue], s=40, color="#ebbcba", edgecolor="black", zorder=6)
# ax.annotate(
#     "rescue underway",
#     xy=(t_rescue, N_rescue),
#     xytext=(t_rescue + 8, N_rescue + 0.08),
#     arrowprops=dict(arrowstyle="->", lw=1.4, color="#ebbcba"),
#     fontsize=10,
#     color="#8b5e5c"
# )


# little microbe icon near the recovery phase
def add_microbe(ax, x=0.82, y=0.80, size=0.04,
                body_color="#c2185b", accent_color="#8e1244"):
    tr = ax.transAxes  # draw in axes coordinates, not data coordinates

    body = Circle((x, y), size, facecolor=body_color, edgecolor=accent_color,
                  linewidth=1.8, zorder=7, transform=tr, clip_on=False)
    ax.add_patch(body)

    # cilia / spikes
    angles = np.linspace(0, 2*np.pi, 12, endpoint=False)
    for ang in angles:
        x0 = x + size * np.cos(ang)
        y0 = y + size * np.sin(ang)
        x1 = x + 1.45 * size * np.cos(ang)
        y1 = y + 1.45 * size * np.sin(ang)
        ax.plot([x0, x1], [y0, y1],
                color=accent_color, lw=1.0, zorder=6,
                transform=tr, clip_on=False)

    # eyes
    eye1 = Circle((x - 0.35*size, y + 0.22*size), 0.18*size,
                  color="white", zorder=8, transform=tr, clip_on=False)
    eye2 = Circle((x + 0.35*size, y + 0.22*size), 0.18*size,
                  color="white", zorder=8, transform=tr, clip_on=False)
    pupil1 = Circle((x - 0.30*size, y + 0.18*size), 0.08*size,
                    color="black", zorder=9, transform=tr, clip_on=False)
    pupil2 = Circle((x + 0.30*size, y + 0.18*size), 0.08*size,
                    color="black", zorder=9, transform=tr, clip_on=False)

    ax.add_patch(eye1)
    ax.add_patch(eye2)
    ax.add_patch(pupil1)
    ax.add_patch(pupil2)

    # smile
    tt = np.linspace(-0.9, 0.9, 60)
    xs = x + 0.42*size*np.sin(tt)
    ys = y - 0.28*size - 0.18*size*np.cos(tt)
    ax.plot(xs, ys, color=accent_color, lw=1.1, zorder=9,
            transform=tr, clip_on=False)

add_microbe(ax, x=0.84, y=0.80, size=0.075)

# aesthetics
ax.set_xlim(times.min(), times.max())
ax.set_ylim(0, 1.05)
# ax.set_xlabel("Time", fontsize=12)
# ax.set_ylabel("Population size", fontsize=12)
# ax.set_title("Microbiome-mediated evolutionary rescue", fontsize=14, pad=10)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", labelsize=10)
ax.legend(frameon=False, loc="lower right", fontsize=10)

plt.tight_layout()

fig.patch.set_alpha(0)
ax.set_facecolor("none")

fig.savefig(
    out,
    dpi=220,
    bbox_inches="tight",
    transparent=True,
    pad_inches=0.02
)

# fig.savefig(
#     out,
#     dpi=220,
#     bbox_inches="tight",
#     transparent=False
# )

plt.close(fig)

print(f"saved {out.resolve()}")