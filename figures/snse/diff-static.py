from math import floor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from plotutils import FIGURE_WIDTH, ImageGrid, tag_axis
from skued import autocenter

INPUT = Path("data") / "snse"

with DiffractionDataset(INPUT / "static.hdf5") as source:
    mask_exf = source.valid_mask
    im_exf = source.diff_data(0)

with DiffractionDataset(INPUT / "overnight4.hdf5") as source:
    mask_umt = source.valid_mask
    im_umt = source.diff_eq()

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 2))
grid = ImageGrid(
    fig, 111, nrows_ncols=(1, 2), cbar_location="top", aspect=False, share_all=False
)

for ax, im, mask, vmax, crop, label in zip(
    grid, [im_exf, im_umt], [mask_exf, mask_umt], [150, 4000], [300, 0], ["a)", "b)"]
):

    r, c = autocenter(im, mask).astype(np.int)

    # Determine the smallest center -> side distance, and crop around that
    side_length = floor(min([c, abs(c - 2048), r, abs(r - 2048)])) - crop
    xs, ys = (
        slice(r - side_length, r + side_length),
        slice(c - side_length, c + side_length),
    )
    im = im[xs, ys]
    im = np.maximum(im, 0)

    m = ax.imshow(
        im,
        vmin=0,
        vmax=vmax,
        cmap="magma",
        extent=[
            0,
            2048,
            2048,
            0,
        ],  # Important to be the same for images which are not the same size
    )
    tag_axis(ax, text=label)
    ax.axis("off")

cbar = grid[0].cax.colorbar(
    m,
    ticks=FixedLocator(locs=[0, 4000]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")

plt.subplots_adjust(bottom=0.01)
