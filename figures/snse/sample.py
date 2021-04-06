from skimage.io import imread
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

from dissutils import LARGE_FIGURE_WIDTH, ImageGrid, tag_axis, discrete_colors

DATADIR = Path("images")
DOWNSAMPLING = 5
COLOR = discrete_colors(1)[0]

sample_tape_lowmag = DATADIR / "snse_surface_tape_lowmag.jpg"
sample_cut_lowmag = DATADIR / "snse_trimmed_lowmag.jpg"
sample_on_grid = DATADIR / "snse_section_grid.jpg"

images = [sample_tape_lowmag, sample_cut_lowmag, sample_on_grid]
scale_bar_left = 2 * [2017] + [1813]

fig = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 0.28 * LARGE_FIGURE_WIDTH))
axes = ImageGrid(fig, 111, nrows_ncols=(1, 3), cbar_mode="none")

for (ax, fname, left, tag) in zip(axes, images, scale_bar_left, "abc"):
    im = imread(fname)[::DOWNSAMPLING, ::DOWNSAMPLING, :]

    scale_bar = mpatches.Rectangle(
        xy=(left / DOWNSAMPLING, 1848 / DOWNSAMPLING),
        width=(2432 - left) / DOWNSAMPLING,
        height=35 / DOWNSAMPLING,
        ec="none",
        fc=COLOR,
    )

    ax.annotate(
        text=f"100 μm",
        xy=(scale_bar.xy[0] + scale_bar.get_width() / 2, scale_bar.xy[1]),
        transform=ax.transData,
        ha="center",
        va="bottom",
        color=COLOR,
    )
    ax.imshow(im)
    ax.add_patch(scale_bar)
    tag_axis(ax, tag + ")")
    ax.axis("off")

plt.tight_layout()
