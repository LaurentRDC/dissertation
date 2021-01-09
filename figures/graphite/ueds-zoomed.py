import itertools as it
from math import floor, sqrt
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
from skued import detector_scattvectors, indices_to_text, nfold, combine_masks

DATADIR = Path("data") / "graphite"
DATASET = DATADIR / "graphite_time_corrected_iris5.hdf5"

qx, qy, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6,
    center=GRAPHITE_CENTER,
)
qq = np.sqrt(qx ** 2 + qy ** 2)


beamblock = np.ones((2048, 2048), dtype=np.bool)
beamblock[0:1260, 900:1130] = False

artifact_mask = np.ones((2048, 2048), dtype=np.bool)
artifact_mask[1084::, 437:482] = False
artifact_mask[0:932, 1296:1324] = False

mask = combine_masks(beamblock, artifact_mask)

with DiffractionDataset(DATASET) as source:
    b4t0 = source.diff_eq()

# Find location of (100) and (200) BZs
graphite = Crystal.from_pwscf(DATADIR / "output.out")
q010 = graphite.scattering_vector((0, 1, 0))
q020 = graphite.scattering_vector((0, 2, 0))
hex_radius = sqrt(3) * np.linalg.norm(q020 - q010) / 2.7

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 2))

with DiffractionDataset(DATASET) as dset:
    image = dset.diff_data(100) - b4t0

gaussian(image, sigma=4, output=image)
image[:] = nfold(image, mod=6, center=GRAPHITE_CENTER, mask=mask, fill_value=np.nan)
image[:] = rotate(image, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect")
image[qq < 1.5] = 0

for ax, center, indices in zip([ax1, ax2], [q010, q020], [(0, 1, 0), (0, 2, 0)]):
    m = ax.imshow(
        image,
        cmap="seismic",
        vmin=-0.4,
        vmax=0.4,
        extent=[qx.min(), qx.max(), qy.max(), qy.min()],
    )

    draw_hexagon_field(
        ax,
        radius=1.7,
        center=(0, 0),
        crystal=graphite,
        reflections=it.product([-3, -2, -1, 0, 1, 2, 3], [-3, -2, -1, 0, 1, 2, 3], [0]),
        color="k",
        linestyle=":",
    )

    cx, cy, _ = center

    ax.scatter(cx, cy, c="k", s=10)
    ax.set_xlim([cx - hex_radius, cx + hex_radius])
    ax.set_ylim([cy - hex_radius, cy + hex_radius])

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    tag_axis(
        ax,
        indices_to_text(*indices),
        x=0.5,
        y=0.95,
        verticalalignment="bottom",
        horizontalalignment="center",
    )
