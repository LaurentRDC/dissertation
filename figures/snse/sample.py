from skimage.io import imread
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from skued import diffshow

from plotutils import FIGURE_WIDTH, ImageGrid, tag_axis

DATADIR = Path("images")
DOWNSAMPLING = 5

sample_tape_lowmag = DATADIR / "snse_surface_tape_lowmag.jpg"
sample_tape_highmag = DATADIR / "snse_surface_tape_highmag.jpg"
sample_cut_lowmag = DATADIR / "snse_trimmed_lowmag.jpg"
sample_on_grid = DATADIR / "snse_section_grid.jpg"

images = [sample_tape_lowmag, sample_tape_highmag, sample_cut_lowmag, sample_on_grid]
scale_bar_left = 3 * [2017] + [1813]
scales_um = [100, 25, 100, 100]

fig = plt.figure(figsize=(FIGURE_WIDTH, 0.6 * FIGURE_WIDTH))
axes = ImageGrid(fig, 111, nrows_ncols=(2, 2), cbar_mode="none")

for (ax, fname, left, scale, tag) in zip(
    axes, images, scale_bar_left, scales_um, "abcd"
):
    im = imread(fname)[::DOWNSAMPLING, ::DOWNSAMPLING, :]

    scale_bar = mpatches.Rectangle(
        xy=(left / DOWNSAMPLING, 1848 / DOWNSAMPLING),
        width=(2432 - left) / DOWNSAMPLING,
        height=35 / DOWNSAMPLING,
        ec="none",
        fc="w",
    )

    ax.annotate(
        text=f"{scale} μm",
        xy=(scale_bar.xy[0] + scale_bar.get_width() / 2, scale_bar.xy[1]),
        transform=ax.transData,
        ha="center",
        va="bottom",
        color="w",
    )
    ax.imshow(im)
    ax.add_patch(scale_bar)
    tag_axis(ax, tag + ")")
    ax.axis("off")


plt.subplots_adjust(bottom=0, top=1)
