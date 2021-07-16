from math import floor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedFormatter, FixedLocator
from skued import autocenter

from dissutils import LARGE_FIGURE_WIDTH, ImageGrid, tag_axis

DOWNSAMPLING = 4

wse2 = np.load(Path("data") / "conclusion" / "wse2_monolayer.npy")
wse2_mask = np.load(Path("data") / "conclusion" / "wse2_monolayer_mask.npy")

mos2 = np.load(Path("data") / "conclusion" / "mos2_monolayer.npy")
mos2_mask = np.load(Path("data") / "conclusion" / "mos2_monolayer_mask.npy")

figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2))
grid = ImageGrid(figure, 111, nrows_ncols=(1, 2), cbar_location="top")

for ax, im, mask, label in zip(
    grid, [wse2, mos2], [wse2_mask, mos2_mask], ["a)", "b)"]
):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    r, c = autocenter(im, mask=mask).astype(int)

    im[:, 1024::] += 6  # correct for CCD bias

    im = im[::DOWNSAMPLING, ::DOWNSAMPLING]
    mask = mask[::DOWNSAMPLING, ::DOWNSAMPLING]
    r = r // DOWNSAMPLING
    c = c // DOWNSAMPLING

    # Determine the smallest center -> side distance, and crop around that
    side_length = floor(min([c, abs(c - im.shape[1]), r, abs(r - im.shape[0])]))
    xs, ys = (
        slice(r - side_length, r + side_length),
        slice(c - side_length, c + side_length),
    )
    im = im[xs, ys]

    m = ax.imshow(im, vmin=0, vmax=2000, cmap="CMRmap_r")

    width = 170 / DOWNSAMPLING
    height = 1.12 * im.shape[0] / 2
    x, y = im.shape[1] // 2 - width / 2, 0

    beamblock_patch = Rectangle(
        xy=(x, y),
        width=width,
        height=height,
        edgecolor="k",
        facecolor="w",
    )
    move_in = 2  # needs to be pixel-perfect. Adjusted for 600 DPI
    crossover_patch = Rectangle(
        xy=(x + move_in, y - 10),
        width=width - 2 * move_in,
        height=height + 10 - move_in,
        fill=True,
        color="w",
        zorder=10,
        clip_on=False,
    )
    ax.add_patch(beamblock_patch)
    ax.add_patch(crossover_patch)

    tag_axis(ax, label)

cbar = grid[0].cax.colorbar(
    m,
    ticks=FixedLocator(locs=[0, 2000]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")
