from math import floor
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from iris import DiffractionDataset
from skimage.filters import gaussian
from matplotlib.ticker import FixedFormatter, FixedLocator
from plotutils import FIGURE_WIDTH, ImageGrid, tag_axis
from skued import autocenter

INPUT = Path("data") / "snse"
DOWNSAMPLING = 4

with DiffractionDataset(INPUT / "static.hdf5") as source:
    mask_exf = source.valid_mask[::DOWNSAMPLING, ::DOWNSAMPLING]
    im_exf = source.diff_data(0)[::DOWNSAMPLING, ::DOWNSAMPLING]

with DiffractionDataset(INPUT / "overnight4.hdf5") as source:
    mask_umt = source.valid_mask[::DOWNSAMPLING, ::DOWNSAMPLING]
    im_umt = source.diff_eq()[::DOWNSAMPLING, ::DOWNSAMPLING]

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 2))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 2), cbar_location="top")

# Determine the smallest image center - side distance
# This way, we ensure that both diffraction patterns will appear to be
# the same size.
side_lengths = list()
centers = list()
for im, mask in zip([im_umt, im_exf], [mask_umt, mask_exf]):
    r, c = autocenter(im, mask).astype(np.int)
    centers.append((r, c))
    side_lengths.append(floor(min([c, abs(c - im.shape[1]), r, abs(r - im.shape[0])])))
side_length = min(side_lengths)

for ax, (r, c), im, mask, vmax, crop, label in zip(
    grid,
    centers,
    [im_umt, im_exf],
    [mask_umt, mask_exf],
    [4000, 150],
    [
        0,
        150,
    ],
    ["a)", "b)"],
):
    xs, ys = (
        slice(r - side_length, r + side_length),
        slice(c - side_length, c + side_length),
    )
    im = im[xs, ys]
    mask = mask[xs, ys]
    im = np.maximum(im, 0)
    im = gaussian(im, sigma=4 / DOWNSAMPLING)

    m = ax.imshow(
        im, vmin=0, vmax=vmax, cmap="magma", extent=[0, im.shape[1], im.shape[0], 0]
    )
    tag_axis(ax, text=label)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    width = 170 / DOWNSAMPLING
    height = 1.2 * im.shape[0] / 2
    x, y = im.shape[1] // 2 - width / 2, 0

    beamblock_patch = mpatches.Rectangle(
        xy=(x, y),
        width=width,
        height=height,
        edgecolor="k",
        facecolor="w",
    )
    move_in = 1  # needs to be pixel-perfect. Adjusted for 600 DPI
    crossover_patch = mpatches.Rectangle(
        xy=(x + move_in, y - 10),
        width=width - 2 * move_in,
        height=height + 10 - move_in,
        fill=True,
        color="w",
        edgecolor="none",
        zorder=10,
        clip_on=False,
    )
    ax.add_patch(beamblock_patch)
    ax.add_patch(crossover_patch)

cbar = grid[0].cax.colorbar(
    m,
    ticks=FixedLocator(locs=[0, 150]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")

plt.subplots_adjust(bottom=0.01)
