import matplotlib.pyplot as plt
from dissutils import MEDIUM_FIGURE_WIDTH, tag_axis, named_arrow
from pathlib import Path
from skimage.io import imread

IMAGES = Path("images")

fig, (ax_pnma, ax_cmcm) = plt.subplots(
    1, 2, figsize=(MEDIUM_FIGURE_WIDTH, MEDIUM_FIGURE_WIDTH)
)

for ax in [ax_pnma, ax_cmcm]:
    ax.axis("off")

for ax, filename in zip([ax_pnma, ax_cmcm], [IMAGES / "pnma.png", IMAGES / "cmcm.png"]):
    im = imread(filename)
    im = im[800:6730, 2000:5000, :]
    ax.imshow(im)

# lattice vectors
arrow_kwds = dict(
    x=1,
    y=0,
    length_includes_head=True,
    width=0.001,
    head_width=0.007,
    fc="k",
    transform=ax_pnma.transAxes,
    clip_on=False,
)
arrow_length = 0.2
named_arrow(
    ax_pnma,
    dx=0,
    dy=arrow_length / 2,  # aspect ratio of 0.5
    text=r"$\mathbf{a}$",
    toffset=(-0.01, 0),
    tkwds=dict(va="center", ha="right", transform=ax_pnma.transAxes),
    **arrow_kwds,
)
named_arrow(
    ax_pnma,
    dx=arrow_length,
    dy=0,
    text=r"$\mathbf{c}$",
    toffset=(0, -0.01),
    tkwds=dict(ha="center", va="top", transform=ax_pnma.transAxes),
    **arrow_kwds,
)

tag_axis(ax_pnma, "a)", x=-0.05, y=1.05)
tag_axis(ax_cmcm, "b)", x=-0.05, y=1.05)
plt.subplots_adjust(
    top=0.95, bottom=0.05, left=0.05, right=0.95, hspace=0.2, wspace=0.2
)
