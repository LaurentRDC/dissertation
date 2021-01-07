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
    GRAPHITE_CENTER,
    ImageGrid,
    draw_hexagon_field,
    tag_axis,
)
from skimage.transform import rotate
from skued import combine_masks, detector_scattvectors, nfold

xc, yc = GRAPHITE_CENTER

xx, yy = np.meshgrid(np.arange(0, 2048), np.arange(0, 2048))
rr = np.sqrt(np.square(xx - xc) + np.square(yy - yc))

beamblock = np.ones((2048, 2048), dtype=np.bool)
beamblock[0:1260, 900:1130] = False

artifact_mask = np.ones((2048, 2048), dtype=np.bool)
artifact_mask[1084::, 437:482] = False
artifact_mask[0:932, 1296:1324] = False

mask = combine_masks(beamblock, artifact_mask)

with DiffractionDataset(
    Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"
) as source:
    b4t0 = source.diff_eq()

b4t0_symmetrized = np.array(b4t0, copy=True)
b4t0_symmetrized = nfold(b4t0, mod=6, center=GRAPHITE_CENTER, mask=mask)
b4t0_symmetrized[rr < 125] = 0

b4t0[:] = rotate(b4t0, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect")
b4t0_symmetrized[:] = rotate(
    b4t0_symmetrized, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect"
)

qx, qy, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6,
    center=(yc, xc),
)

# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([xc, abs(xc - 2048), yc, abs(yc - 2048)]))
xs, ys = (
    slice(yc - side_length, yc + side_length),
    slice(xc - side_length, xc + side_length),
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
