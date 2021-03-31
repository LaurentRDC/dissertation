import matplotlib.pyplot as plt
import numpy as np
import itertools as it
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Ellipse
from pathlib import Path
from iris import DiffractionDataset
from plotutils import MEDIUM_FIGURE_WIDTH, ImageGrid
from plotutils.snse_datasets import overnight4
from crystals import Crystal

# To perform analysis on another dataset, simply change the following
# line:
CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")
_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters

INNER_RADIUS = 15
OUTER_RADIUS = 30

# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [
    (0, -1, 3),
    (0, 0, 2),
    (0, 2, 0),
    (0, 0, -2),
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
    t1 = np.argmin(np.abs(timedelays - 5))
    t2 = np.argmin(np.abs(timedelays - 15))
    IMAGE = np.mean(dset.diffraction_group["intensity"][:, :, t1:t2], axis=2)
IMAGE -= eq
IMAGE /= eq


def diffuse_amplitude(h, k, l):

    result = np.zeros(shape=(len(EXTENT), len(EXTENT)), dtype=float)

    for y, z in it.product(EXTENT, repeat=2):

        yi, xi = overnight4.miller_to_arrindex(h, k + y, l + z)

        result[np.argmin(np.abs(EXTENT - y)), np.argmin(np.abs(EXTENT - z))] = np.mean(
            IMAGE[
                xi - INNER_RADIUS : xi + INNER_RADIUS,
                yi - INNER_RADIUS : yi + INNER_RADIUS,
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
ax.set_xlabel("$\mathbf{k}_y / \mathbf{a_2}$ [a.u.]")
ax.set_ylabel("$\mathbf{k}_z / \mathbf{a_3}$ [a.u.]")

for label, pos in [
    ("T", (1 / 2, 1 / 2)),
    ("Y", (1 / 2, 0)),
    ("Z", (0, 1 / 2)),
]:
    ax.scatter(*pos, color="k", clip_on=False)
    ax.annotate(
        label,
        xy=pos,
        xytext=np.asarray(pos) + np.asarray([-0.05, -0.05]),
    )

ax.annotate(
    "$\Gamma$",
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
