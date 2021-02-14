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
    draw_hexagon,
    tag_axis,
)
from skimage.transform import rotate
from skimage.filters import gaussian
from skued import detector_scattvectors, nfold

DATASET = Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"
DOWNSAMPLING = 4

with DiffractionDataset(DATASET, mode="r") as source:
    b4t0 = source.diff_eq()[::DOWNSAMPLING, ::DOWNSAMPLING]
    mask = source.valid_mask[::DOWNSAMPLING, ::DOWNSAMPLING]
    c, r = (np.asarray(source.center) / DOWNSAMPLING).astype(int)

qx, qy, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=b4t0.shape,
    pixel_size=DOWNSAMPLING * 14e-6,
    center=(r, c),
)


# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([c, abs(c - b4t0.shape[1]), r, abs(r - b4t0.shape[0])]))
xs = slice(r - side_length, r + side_length)
ys = slice(c - side_length, c + side_length)

qx = qx[ys, xs]
qy = qy[ys, xs]
qq = np.sqrt(qx ** 2 + qy ** 2)

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH))
grid = ImageGrid(fig, 111, nrows_ncols=(2, 2), cbar_location="top")

with DiffractionDataset(DATASET) as dset:
    for time, ax, letter in zip([0.5, 1.5, 5, 100], grid, "abcd"):

        image = nfold(
            dset.diff_data(time)[::DOWNSAMPLING, ::DOWNSAMPLING] - b4t0,
            mod=6,
            center=(c, r),
            mask=mask,
            fill_value=np.nan,
        )
        gaussian(image, sigma=4 / DOWNSAMPLING, output=image)
        image[:] = rotate(image, angle=GRAPHITE_ANGLE, center=(c, r), mode="reflect")
        image = image[xs, ys]
        image[qq < 1.5] = 0
        m = ax.imshow(
            image,
            cmap="seismic",
            vmin=-0.4,
            vmax=0.4,
            extent=[qx.min(), qx.max(), qy.min(), qy.max()],
        )

        draw_hexagon(
            ax,
            radius=1.7,
            center=(0, 0),
            color="w",
            facecolor="w",
        )

        draw_hexagon_field(
            ax,
            radius=1.7,
            crystal=Crystal.from_pwscf(Path("data") / "graphite" / "output.out"),
            reflections=it.product(
                [-5, -4, -3, -2, -1, 0], [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4], [0]
            ),
            color="k",
            alpha=0.5,
        )

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        tag = f"{letter}) {1e3 * time:0.0f} fs" if time < 1 else f"{letter}) {time} ps"
        tag_axis(ax, text=tag)

cbar = ax.cax.colorbar(
    m,
    ticks=FixedLocator(locs=[-0.4, 0, 0.4]),
    format=FixedFormatter(["-1", "0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity change [a.u.]")

plt.subplots_adjust(bottom=0.01)
