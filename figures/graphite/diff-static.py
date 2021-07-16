import itertools as it
from math import floor
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from skimage.transform import rotate
from skued import detector_scattvectors, nfold

from dissutils import (
    GRAPHITE_ANGLE,
    GRAPHITE_CAMERA_LENGTH,
    LARGE_FIGURE_WIDTH,
    ImageGrid,
    draw_hexagon,
    draw_hexagon_field,
    tag_axis,
)

DOWNSAMPLING = 4

with DiffractionDataset(
    Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"
) as source:
    mask = source.valid_mask[::DOWNSAMPLING, ::DOWNSAMPLING]
    b4t0 = source.diff_data(source.time_points[0])[::DOWNSAMPLING, ::DOWNSAMPLING]
    c, r = (np.asarray(source.center) / DOWNSAMPLING).astype(int)

xx, yy = np.meshgrid(np.arange(0, b4t0.shape[1]), np.arange(0, b4t0.shape[0]))
rr = np.sqrt(np.square(xx - c) + np.square(yy - r))


b4t0_symmetrized = np.array(b4t0, copy=True)
b4t0_symmetrized = nfold(b4t0, mod=6, center=(c, r), mask=mask)
b4t0_symmetrized[rr < 125 / DOWNSAMPLING] = 0

b4t0[:] = rotate(b4t0, angle=GRAPHITE_ANGLE, center=(c, r), mode="reflect")
b4t0_symmetrized[:] = rotate(
    b4t0_symmetrized, angle=GRAPHITE_ANGLE, center=(c, r), mode="reflect"
)

qx, qy, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=b4t0.shape,
    pixel_size=DOWNSAMPLING * 14e-6,
    center=(r, c),
)

# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([c, abs(c - b4t0.shape[1]), r, abs(r - b4t0.shape[0])]))
xs, ys = (
    slice(r - side_length, r + side_length),
    slice(c - side_length, c + side_length),
)
b4t0 = b4t0[xs, ys]
b4t0_symmetrized = b4t0_symmetrized[xs, ys]
qx = qx[ys, xs]
qy = qy[ys, xs]

fig = plt.figure(figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 2), cbar_location="top")

for ax, im, label in zip(grid, [b4t0, b4t0_symmetrized], ["a)", "b)"]):
    m = ax.imshow(
        im,
        vmin=0,
        vmax=200,
        cmap="inferno",
        extent=[qx.min(), qx.max(), qy.min(), qy.max()],
    )
    tag_axis(ax, text=label)

draw_hexagon_field(
    grid[1],
    radius=1.7,
    crystal=Crystal.from_pwscf(Path("data") / "graphite" / "graphite.out"),
    reflections=it.product(range(-5, 5), range(-5, 5), [0]),
    color="w",
    linewidth=0.3,
    alpha=0.5,
)

# Beam block patch
# assuming that the image is centered
dk = abs(qx[1, 1] - qx[0, 0])
width = dk * 250 / DOWNSAMPLING
height = 1.5 * dk * b4t0.shape[0] / 2
x, y = -width / 2, -1.2

beamblock_patch = mpatches.Rectangle(
    xy=(x, y),
    width=width,
    height=height,
    edgecolor="k",
    facecolor="w",
)
move_in = 2 * dk  # needs to be pixel-perfect. Adjusted for 600 DPI
crossover_patch = mpatches.Rectangle(
    xy=(x + move_in, y + move_in),
    width=width - 2 * move_in,
    height=height + 10 * dk - move_in,
    fill=True,
    color="w",
    edgecolor="none",
    zorder=10,
    clip_on=False,
)
grid[0].add_patch(beamblock_patch)
grid[0].add_patch(crossover_patch)


draw_hexagon(
    grid[-1],
    radius=1.7,
    center=(0, 0),
    color="w",
    facecolor="w",
)

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
