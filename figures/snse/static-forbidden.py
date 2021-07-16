import itertools as it
from math import floor
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from skimage.filters import gaussian

from dissutils import LARGE_FIGURE_WIDTH, ImageGrid
from dissutils.snse import overnight4

INPUT = Path("data") / "snse"
DOWNSAMPLING = 4


def pnma_valid(h, k, l):
    """
    True for Miller indices associated with allowed reflections in Pnma space group.
    False otherwise.
    """
    even = lambda n: n % 2 == 0

    return (h == 0) and (
        even(k + l) or (even(k) and (l == 0)) or ((k == 0) and even(l))
    )


with DiffractionDataset(INPUT / "overnight4.hdf5") as source:
    mask = source.valid_mask[::DOWNSAMPLING, ::DOWNSAMPLING]
    im = source.diff_eq()[::DOWNSAMPLING, ::DOWNSAMPLING]
    c, r = (np.asarray(source.center) / DOWNSAMPLING).astype(int)
    side_length = floor(min([c, abs(c - im.shape[1]), r, abs(r - im.shape[0])]))


fig = plt.figure(figsize=(LARGE_FIGURE_WIDTH / 2, LARGE_FIGURE_WIDTH / 2))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 1), cbar_location="top")
ax = grid[0]

xs, ys = (
    slice(r - side_length, r + side_length),
    slice(c - side_length, c + side_length),
)
im = im[xs, ys]
mask = mask[xs, ys]
im = np.maximum(im, 0)
im = gaussian(im, sigma=4 / DOWNSAMPLING)

m = ax.imshow(im, vmin=0, vmax=4000, cmap="magma")
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

off_x = im
for k, l in it.product(range(-3, 4), repeat=2):
    if pnma_valid(0, k, l):
        continue

    yi, xi = overnight4.miller_to_arrindex(0, k, l)
    xi /= DOWNSAMPLING
    yi /= DOWNSAMPLING

    ax.add_patch(
        mpatches.Circle(
            xy=(yi - (c - side_length), xi - (r - side_length)),
            radius=10,
            fc="none",
            ec="r",
        )
    )


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
    ticks=FixedLocator(locs=[0, 4000]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")

plt.subplots_adjust(bottom=0.01)
