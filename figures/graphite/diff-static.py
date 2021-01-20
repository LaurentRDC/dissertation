import itertools as it
from math import floor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from plotutils import (
    FIGURE_WIDTH,
    GRAPHITE_ANGLE,
    GRAPHITE_CAMERA_LENGTH,
    ImageGrid,
    draw_hexagon_field,
    tag_axis,
)
from skimage.transform import rotate
from skued import detector_scattvectors, nfold, autocenter

with DiffractionDataset(
    Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"
) as source:
    mask = source.valid_mask
    b4t0 = source.diff_eq()

r, c = autocenter(im=b4t0, mask=mask).astype(np.int)

xx, yy = np.meshgrid(np.arange(0, 2048), np.arange(0, 2048))
rr = np.sqrt(np.square(xx - c) + np.square(yy - r))


b4t0_symmetrized = np.array(b4t0, copy=True)
b4t0_symmetrized = nfold(b4t0, mod=6, center=(c, r), mask=mask)
b4t0_symmetrized[rr < 125] = 0

b4t0[:] = rotate(b4t0, angle=GRAPHITE_ANGLE, center=(c, r), mode="reflect")
b4t0_symmetrized[:] = rotate(
    b4t0_symmetrized, angle=GRAPHITE_ANGLE, center=(c, r), mode="reflect"
)

qx, qy, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6,
    center=(r, c),
)

# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([c, abs(c - 2048), r, abs(r - 2048)]))
xs, ys = (
    slice(r - side_length, r + side_length),
    slice(c - side_length, c + side_length),
)
b4t0 = b4t0[xs, ys]
b4t0_symmetrized = b4t0_symmetrized[xs, ys]
qx = qx[ys, xs]
qy = qy[ys, xs]

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 2))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 2), cbar_location="top")

for ax, im, label in zip(grid, [b4t0, b4t0_symmetrized], ["a)", "b)"]):
    m = ax.imshow(
        im,
        vmin=0,
        vmax=200,
        cmap="magma",
        extent=[qx.min(), qx.max(), qy.min(), qy.max()],
    )
    draw_hexagon_field(
        ax,
        radius=1.7,
        crystal=Crystal.from_pwscf(Path("data") / "graphite" / "output.out"),
        reflections=it.product(range(-4, 4), range(-4, 4), [0]),
        color="w",
        linewidth=0.5,
        linestyle=":",
    )
    tag_axis(ax, text=label)

for ax in grid:
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

cbar = grid[0].cax.colorbar(
    m,
    ticks=FixedLocator(locs=[0, 200]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")

plt.subplots_adjust(bottom=0.01)
