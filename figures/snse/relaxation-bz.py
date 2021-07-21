import itertools as it
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from crystals import Crystal
from iris import DiffractionDataset
from matplotlib.patches import Circle, Ellipse
from matplotlib.ticker import PercentFormatter

from dissutils import MEDIUM_FIGURE_WIDTH, ImageGrid
from dissutils.snse import overnight4

# To perform analysis on another dataset, simply change the following
# line:
CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")
_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters

INNER_RADIUS = 15
OUTER_RADIUS = 30

# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [
    (0, -1, 3),
    (0, 2, 0),
    (0, -2, 0),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -3, 5),
]

EXTENT = np.linspace(start=-1 / 2, stop=1 / 2, num=64, endpoint=True)

with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points
    eq = dset.diff_eq()
    t1 = np.argmin(np.abs(timedelays - 4))
    t2 = np.argmin(np.abs(timedelays - 6))
    IMAGE = np.mean(dset.diffraction_group["intensity"][:, :, t1:t2], axis=2)
IMAGE -= eq
IMAGE /= eq


def diffuse_amplitude(h, k, l):

    result = np.zeros(shape=(len(EXTENT), len(EXTENT)), dtype=float)

    # Reflections are not perfectly aligned
    xoffset = -5
    yoffset = -5

    for y, z in it.product(EXTENT, repeat=2):

        yi, xi = overnight4.miller_to_arrindex(h, k + y, l + z)

        result[np.argmin(np.abs(EXTENT - y)), np.argmin(np.abs(EXTENT - z))] = np.mean(
            IMAGE[
                xi - INNER_RADIUS + xoffset : xi + INNER_RADIUS + xoffset,
                yi - INNER_RADIUS + yoffset : yi + INNER_RADIUS + yoffset,
            ]
        )

    return result


result = sum(diffuse_amplitude(h, k, l) for (h, k, l) in INDICES_DIFFUSE) / len(
    INDICES_DIFFUSE
)

# The image shows a deeply negative value where the Debye-Waller has occured
# Therefore, better to scale according to the max, and not abs().max()
vmax = result.max()
fig = plt.figure(figsize=(MEDIUM_FIGURE_WIDTH, MEDIUM_FIGURE_WIDTH))
(ax,) = ImageGrid(fig, 111, nrows_ncols=(1, 1), cbar_location="top", cbar_mode="each")
m = ax.imshow(
    result,
    vmin=-vmax,
    vmax=vmax,
    EXTENT=[EXTENT.min(), EXTENT.max(), EXTENT.max(), EXTENT.min()],
    cmap="RdBu_r",
)
ax.cax.colorbar(mappable=m, format=PercentFormatter(xmax=1))

ax.cax.set_xlabel("$\Delta I / I_0$ [a.u.]")
ax.set_xlabel("$k_y / b^{*}$")
ax.set_ylabel("$k_z / c^{*}$")

ax.xaxis.set_major_locator(ticker.FixedLocator([-0.5, 0, 0.5]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(["-½", "0", "½"]))
ax.yaxis.set_major_locator(ticker.FixedLocator([-0.5, 0, 0.5]))
ax.yaxis.set_major_formatter(ticker.FixedFormatter(["-½", "0", "½"]))

for label, pos in [
    ("T", (1 / 2, 1 / 2)),
    ("Y", (1 / 2, 0)),
    ("Z", (0, 1 / 2)),
]:
    ax.add_patch(
        Circle(
            xy=pos,
            radius=0.015,
            edgecolor="None",
            facecolor="k",
            zorder=np.inf,
            clip_on=False,
        )
    )
    ax.annotate(
        label,
        xy=pos,
        xytext=np.asarray(pos) + np.asarray([-0.05, -0.05]),
    )

ax.annotate(
    "$\mathbf{\Gamma}$",
    xy=(0, 0),
    xytext=(0.4, 0.65),
    textcoords="axes fraction",
    arrowprops=dict(arrowstyle="->"),
)

# Show regions 1 and 2 from the main text
for r in [0.114, 0.228]:  # Two radii in A^-1
    ax.add_patch(
        Ellipse(
            xy=(0, 0),
            width=r / bstar,
            height=r / cstar,
            edgecolor="k",
            linewidth=1,
            fill=False,
        )
    )

ax.set_aspect(cstar / bstar)

plt.subplots_adjust(
    top=0.884, bottom=0.133, left=0.218, right=0.934, hspace=0.2, wspace=0.2
)
