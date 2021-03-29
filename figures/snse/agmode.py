import matplotlib.pyplot as plt
from plotutils import LARGE_FIGURE_WIDTH, named_arrow
from pathlib import Path
from skimage.io import imread

IMAGES = Path("images")

fig, ax = plt.subplots(1, 1, figsize=(2, 4))
ax.axis("off")

im = imread(IMAGES / "Ag_mode.png")
im = im[800:6730, 2000:5000, :]
ax.imshow(im)

# lattice vectors
arrow_kwds = dict(
    x=-0.05,
    y=-0.05,
    length_includes_head=True,
    width=0.001,
    head_width=0.007,
    fc="k",
    transform=ax.transAxes,
    clip_on=False,
)
arrow_length = 0.2
named_arrow(
    ax,
    dx=0,
    dy=arrow_length / 2,  # aspect ratio of 0.5
    text=r"$\mathbf{a}$",
    toffset=(-0.01, 0),
    tkwds=dict(va="center", ha="right", transform=ax.transAxes),
    **arrow_kwds,
)
named_arrow(
    ax,
    dx=arrow_length,
    dy=0,
    text=r"$\mathbf{c}$",
    toffset=(0, -0.01),
    tkwds=dict(ha="center", va="top", transform=ax.transAxes),
    **arrow_kwds,
)
plt.subplots_adjust(top=1.0, bottom=0.0, left=0.1, right=0.9, hspace=0.2, wspace=0.2)
