"""
Slide 06 -- what is the interface between interacting species?

Three panels of Fig. 1 of Pauw, Stofberg & Waterman (2009) Evolution
63:268-279, cropped from assets/pauw2009_fig1.png (a published figure at low
resolution: check reuse terms before circulating the slides).

  A, B  Moegistorhynchus longirostris at Lapeirousia anceps
  C     nectar consumed (open circles, solid quadratic fit) and pollen
        deposited (filled circles, dashed fit) against proboscis length minus
        tube length

Panels D and E of the original are left out so the slide stays on the
interaction and on outcomes as a function of the trait difference.

Added here: the headings, the two legend lines, and one question.  The solid
fit rises, levels off near a lead of about 18 mm and turns down; the callout
asks: a cost, or diminishing returns?
It is a question put to the audience, not a claim of the paper.
"""
import matplotlib.image as mpimg
import numpy as np
from matplotlib.patches import Circle, Rectangle

from semstyle import *  # noqa: F403

PHOTO = HERE.parent / "assets" / "pauw2009_fig1.png"
# crop boxes in source pixels (x0, x1, y0, y1), found from the white gutters
BOX = {"A": (0, 150, 0, 275), "B": (154, 346, 0, 275), "C": (0, 346, 303, 559)}
# where each goes on the slide: x, y (lower left) and width, in slide pixels
PLACE = {"A": (30, 150, 184), "B": (226, 150, 236), "C": (556, 92, 520)}
# the top of the solid fit in panel C, in that panel's source pixels (x, y from its top)
PEAK = (250, 15)


def main():
    fig = figure_px(FULL)
    cv = canvas(fig)
    peak = None
    if PHOTO.exists():
        img = mpimg.imread(PHOTO)
        for key, (x0, x1, y0, y1) in BOX.items():
            sub = img[y0:y1, x0:x1]
            px_, py_, w = PLACE[key]
            h = w * sub.shape[0] / sub.shape[1]
            cv.imshow(sub, extent=[px_, px_ + w, py_, py_ + h], zorder=3, aspect="auto",
                      interpolation="lanczos")
            if key in "AB":
                cv.add_patch(Rectangle((px_, py_), w, h, facecolor="none", edgecolor=INK,
                                       linewidth=lw(1.2), zorder=4))
            else:
                k = w / sub.shape[1]
                peak = (px_ + PEAK[0] * k, py_ + h - PEAK[1] * k)
        cv.set_xlim(0, fig._px[0])
        cv.set_ylim(0, fig._px[1])
    cv.text(30, 118, "Moegistorhynchus longirostris", ha="left", va="center",
            fontsize=fs(FS_SMALL), style="italic", color=INDIGO)
    cv.text(30, 92, "at Lapeirousia anceps", ha="left", va="center",
            fontsize=fs(FS_SMALL), style="italic", color=OCHRE)
    cv.text(30, 520, "A fly and a flower, in the field", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(556, 520, "Who gains depends on the difference", ha="left", va="baseline",
            fontsize=fs(FS_LABEL), fontweight="bold")
    cv.text(556, 66, "open circles, solid line: nectar the fly takes", ha="left",
            va="center", fontsize=fs(FS_SMALL), color=INDIGO)
    cv.text(556, 42, "filled circles, dashed line: pollen the flower receives", ha="left",
            va="center", fontsize=fs(FS_SMALL), color=OCHRE)
    if peak is not None:                      # the question the quadratic fit raises
        cv.add_patch(Circle(peak, 15, facecolor="none", edgecolor=CRIMSON, linewidth=lw(2.4),
                            zorder=6))
        tx = peak[0] + 58                                  # to the right of the ring
        cv.text(tx, 514, "a cost, or", ha="left", va="center", fontsize=fs(FS_SMALL),
                color=CRIMSON, fontweight="bold")
        cv.text(tx, 492, "diminishing returns?", ha="left", va="center",
                fontsize=fs(FS_SMALL), color=CRIMSON, fontweight="bold")
        end = (tx - 8, 500)                                # stops short of the words
        d = np.hypot(end[0] - peak[0], end[1] - peak[1])
        start = (peak[0] + 15 * (end[0] - peak[0]) / d, peak[1] + 15 * (end[1] - peak[1]) / d)
        cv.plot([start[0], end[0]], [start[1], end[1]], color=CRIMSON, lw=lw(1.6), zorder=6)
    status(cv, "published",
           "Pauw, Stofberg & Waterman 2009, Evolution 63:268 \u2014 their Fig. 1, panels A\u2013C; the question in red is added")
    save(fig, "s06_interface")


if __name__ == "__main__":
    main()
