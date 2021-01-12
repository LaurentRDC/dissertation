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
from skimage.filters import gaussian
from skued import detector_scattvectors, nfold

DATASET = Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"

xc, yc = GRAPHITE_CENTER
qx, qy, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6,
    center=(yc, xc),
)

with DiffractionDataset(DATASET, mode='r') as source:
    b4t0 = source.diff_eq()
    mask = source.valid_mask
    

# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([xc, abs(xc - 2048), yc, abs(yc - 2048)]))
xs = slice(yc - side_length, yc + side_length)
ys = slice(xc - side_length, xc + side_length)

qx = qx[ys, xs]
qy = qy[ys, xs]
qq = np.sqrt(qx ** 2 + qy ** 2)

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH))
grid = ImageGrid(fig, 111, nrows_ncols=(2, 2), cbar_location="top")

with DiffractionDataset(DATASET) as dset:
    for time, ax, letter in zip([0.5, 1.5, 5, 100], grid, "abcd"):

        image = nfold(
            dset.diff_data(time) - b4t0,
            mod=6,
            center=GRAPHITE_CENTER,
            mask=mask,
            fill_value=np.nan,
        )
        gaussian(image, sigma=4, output=image)
        image[:] = rotate(
            image, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect"
        )
        image = image[xs, ys]
        image[qq < 1.5] = 0
        m = ax.imshow(
            image,
            cmap="seismic",
            vmin=-0.4,
            vmax=0.4,
            extent=[qx.min(), qx.max(), qy.min(), qy.max()],
        )

        draw_hexagon_field(
            ax,
            radius=1.7,
            crystal=Crystal.from_pwscf(Path("data") / "graphite" / "output.out"),
            reflections=it.product(
                [-4, -3, -2, -1, 0], [-4, -3, -2, -1, 0, 1, 2, 3], [0]
            ),
            color="k",
            linestyle=":",
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