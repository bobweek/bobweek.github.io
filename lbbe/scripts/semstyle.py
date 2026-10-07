"""
semstyle.py -- shared style for the seminar figures.

Design contract (inherited from rev15/statement_style.py, rescaled to a slide)
-----------------------------------------------------------------------------
* The deck is 1280 x 720.  Every figure is drawn on a canvas measured in
  *slide pixels* and saved without tight cropping, so one unit in the script is
  one pixel on the slide and type sizes are literal: FS_LABEL = 20 means
  20 px on the projected slide.  PNGs are written at 2x for sharp projection.
* One typeface for labels and maths (Inter -> TeX Gyre Heros -> Liberation).
* Colour is additive, never the only channel: every distinction is also
  carried by marker shape, line style, fill or hatching.  SEMINAR_MONO=1
  renders pure greyscale; a greyscale proof is written on every build.
* Every build runs a collision check (overlapping text, text off canvas).
  SEMINAR_STRICT=1 turns a problem into a failed build.
* Every figure carries a status line: published / in prep. / planned /
  schematic.  Nothing unpublished is drawn as a finding.
"""
from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE.parent / "figures"

PX = 1 / 100.0                    # inches per slide pixel (figure dpi = 100)
PT = 0.72                         # points per slide pixel
OUT_DPI = 200                     # PNGs at 2x

MONO = os.environ.get("SEMINAR_MONO", "0") == "1"
STRICT = os.environ.get("SEMINAR_STRICT", "0") == "1"

# canvases (slide px).  FULL sits under a slide title; BLEED is the whole slide.
FULL = (1180, 560)
BLEED = (1280, 720)
HALF = (580, 540)

# ---------------------------------------------------------------- palette --
INK = "#151A1F"
MUTED = "#5F6A73"
FAINT = "#A4ADB5"
LINE = "#D5DCE2"
SOFT = "#F3F5F7"
PAPER = "#FFFFFF"
NAVY = "#243A5E"

_ACCENTS = {
    "CRIMSON": ("#A8123F", "#2B2B2B"),   # single-species quantitative genetics
    "INDIGO": ("#5A4CC7", "#2B2B2B"),    # coevolution; partner X; pollinators
    "OCHRE": ("#B0741A", "#7A7A7A"),     # second partner Y; plants; contrast
    "TEAL": ("#0B6F9C", "#2B2B2B"),      # host-microbiome
}
_TINTS = {
    "CRIMSON_T": ("#F8DBE3", "#E6E6E6"),
    "INDIGO_T": ("#E9E6FA", "#E6E6E6"),
    "OCHRE_T": ("#F6E6CE", "#EFEFEF"),
    "TEAL_T": ("#DCEEF6", "#E6E6E6"),
}
for _name, (_col, _grey) in {**_ACCENTS, **_TINTS}.items():
    globals()[_name] = _grey if MONO else _col

# ------------------------------------------------------- type sizes (px) --
FS_TITLE = 26       # in-figure headings
FS_LABEL = 20       # axis labels, node names
FS_SMALL = 17       # annotations
FS_TINY = 15        # status line, citations; nothing smaller is used
LW_AXIS = 1.4       # px
LW_DATA = 3.4       # px

FONT_STACK = ["Inter", "TeX Gyre Heros", "Liberation Sans", "DejaVu Sans"]
SERIF_STACK = ["TeX Gyre Pagella", "Palatino Linotype", "Palatino", "DejaVu Serif"]


def fs(px: float) -> float:
    """Font size in points for a height of `px` slide pixels."""
    return px * PT


def lw(px: float) -> float:
    """Line width in points for a width of `px` slide pixels."""
    return px * PT


def apply_style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 100,
        "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "axes.facecolor": "none",
        "axes.edgecolor": FAINT,
        "axes.labelcolor": INK,
        "axes.linewidth": lw(LW_AXIS),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": fs(FS_LABEL),
        "axes.labelpad": 4,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.serif": SERIF_STACK,
        "font.size": fs(FS_LABEL),
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": fs(FS_SMALL),
        "ytick.labelsize": fs(FS_SMALL),
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.major.width": lw(LW_AXIS),
        "ytick.major.width": lw(LW_AXIS),
        "xtick.major.pad": 4,
        "ytick.major.pad": 4,
        "lines.linewidth": lw(LW_DATA),
        "lines.solid_capstyle": "round",
        "lines.dash_capstyle": "butt",
        "hatch.linewidth": 0.8,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Inter",
        "mathtext.it": "Inter:italic",
        "mathtext.bf": "Inter:bold",
        "mathtext.sf": "Inter",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


# ------------------------------------------------------------------- cache --
DATA_DIR = HERE.parent / "data"


def cached(name: str, fn):
    """Run `fn()` (which returns a dict of arrays) once and keep the result in
    data/<name>.npz.  Set SEMINAR_RECOMPUTE=1 to run it again."""
    import numpy as np
    DATA_DIR.mkdir(exist_ok=True)
    f = DATA_DIR / f"{name}.npz"
    if f.exists() and os.environ.get("SEMINAR_RECOMPUTE", "0") != "1":
        return dict(np.load(f, allow_pickle=False))
    out = fn()
    np.savez_compressed(f, **out)
    return out


# ------------------------------------------------------------------ layout --
def figure_px(size=FULL):
    """A figure of `size` slide pixels."""
    apply_style()
    w, h = size
    fig = plt.figure(figsize=(w * PX, h * PX))
    fig._px = (w, h)
    return fig


def px_axes(fig, x, y, w, h, **kw):
    """Axes placed by lower-left corner and size, in slide pixels."""
    W, H = fig._px
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def canvas(fig):
    """Full-figure axes whose data coordinates are slide pixels."""
    W, H = fig._px
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    ax.set_zorder(-1)
    return ax


def camera_axes(fig, cx, cy, width):
    """Full-figure axes looking at a scene: centred on (cx, cy), showing
    `width` scene units across.  Returns (ax, k) where k = slide px per scene
    unit, to be passed to glyphs so strokes scale with the zoom."""
    W, H = fig._px
    k = W / width
    height = H / k
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(cx - width / 2, cx + width / 2)
    ax.set_ylim(cy - height / 2, cy + height / 2)
    ax.axis("off")
    ax.set_zorder(-1)
    return ax, k


def blank(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ------------------------------------------------------------------ status --
_STATUS = {"published": "published", "inprep": "in prep.",
           "planned": "planned", "schematic": "schematic"}


def status_mark(cv, kind: str, x: float, y: float, r: float = 6.0,
                color=None, z=20):
    """The four status marks, drawn as patches (no font dependence):
    filled = published, half = in prep., dashed open = planned,
    open square = schematic."""
    from matplotlib.patches import Circle, Rectangle, Wedge
    col = color or MUTED
    if kind == "published":
        cv.add_patch(Circle((x, y), r, facecolor=col, edgecolor=col,
                            linewidth=lw(1.3), zorder=z))
    elif kind == "inprep":
        cv.add_patch(Circle((x, y), r, facecolor=PAPER, edgecolor=col,
                            linewidth=lw(1.3), zorder=z))
        cv.add_patch(Wedge((x, y), r, 90, 270, facecolor=col, edgecolor="none",
                           zorder=z + 0.1))
    elif kind == "planned":
        cv.add_patch(Circle((x, y), r, facecolor=PAPER, edgecolor=col,
                            linewidth=lw(1.3), linestyle=(0, (1.6, 1.4)),
                            zorder=z))
    elif kind == "schematic":
        cv.add_patch(Rectangle((x - r * 0.9, y - r * 0.9), 1.8 * r, 1.8 * r,
                               facecolor=PAPER, edgecolor=col,
                               linewidth=lw(1.3), zorder=z))
    else:
        raise KeyError(kind)


def status(cv, kind: str, text: str = "", x: float = 10, y: float = 8):
    """Status line in the lower-left corner: mark, kind, and source.

    published  a result from a published paper (give the citation in `text`)
    inprep     a real result from a manuscript in preparation
    planned    planned work: only ever drawn as a schematic
    schematic  a picture of an idea; no data or simulation is claimed
    """
    status_mark(cv, kind, x + 6, y + 7)
    s = _STATUS[kind] + (f"  \u00B7  {text}" if text else "")
    return cv.text(x + 19, y, s, ha="left", va="bottom", fontsize=fs(FS_TINY),
                   color=MUTED)


# ---------------------------------------------------------------- quality --
def collision_report(fig, pad_px: float = 1.0, ignore=()):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    items = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        if t in ignore or getattr(t, "_allow_overlap", False):
            continue
        bb = t.get_window_extent(r)
        if bb.width < 1 or bb.height < 1:
            continue
        items.append((t.get_text().replace("\n", " / ")[:40], bb))
    problems = []
    for name, bb in items:
        if (bb.x0 < fb.x0 - 0.5 or bb.y0 < fb.y0 - 0.5 or
                bb.x1 > fb.x1 + 0.5 or bb.y1 > fb.y1 + 0.5):
            problems.append(f"off-canvas: '{name}'")
    for (n1, b1), (n2, b2) in combinations(items, 2):
        ox = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
        oy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
        if ox > pad_px and oy > pad_px:
            problems.append(f"overlap: '{n1}' x '{n2}'")
    return problems


def save(fig, stem: str, outdir: Path | None = None):
    outdir = Path(outdir or FIG_DIR)
    outdir.mkdir(parents=True, exist_ok=True)
    probs = collision_report(fig)
    for p in probs:
        print(f"  [{stem}] {p}")
    if probs and STRICT:
        raise SystemExit(f"{stem}: {len(probs)} layout problem(s)")
    fig.savefig(outdir / f"{stem}.png", dpi=OUT_DPI)
    fig.savefig(outdir / f"{stem}.pdf")
    try:
        from PIL import Image
        proof = outdir / "proofs"
        proof.mkdir(exist_ok=True)
        Image.open(outdir / f"{stem}.png").convert("L").save(
            proof / f"{stem}_grey.png")
    except Exception as exc:  # pragma: no cover
        print(f"  [{stem}] greyscale proof skipped: {exc}")
    plt.close(fig)
    print(f"{stem}: {fig._px[0]} x {fig._px[1]} px, {len(probs)} layout problem(s)")
