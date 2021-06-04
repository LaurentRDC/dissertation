import matplotlib.pyplot as plt
import numpy as np
import itertools as it
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Ellipse, Circle
import matplotlib.ticker as ticker
from pathlib import Path
from iris import DiffractionDataset
from skued import azimuthal_average
import scipy.optimize as opt
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors
from dissutils.snse import overnight4
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
    (0, 2, 0),
    (0, -2, 0),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -3, 5),
]

EXTENT = np.linspace(start=-1 / 2, stop=1 / 2, num=64, endpoint=True)


def polaron(q, A, rp, C):
    return A * q * rp ** 3 / (1 + (q * rp) ** 2) ** 2 + C


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

_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters
kx, ky = np.meshgrid(
    np.linspace(-bstar / 2, bstar / 2, num=len(EXTENT)),
    np.linspace(-cstar / 2, cstar / 2, num=len(EXTENT)),
)
kk = np.hypot(kx, ky)
_, profile = azimuthal_average(result, center=np.asarray(result.shape) // 2)
_, kr = azimuthal_average(kk, center=np.asarray(result.shape) // 2)


profile /= profile.max()
params, pcov = opt.curve_fit(
    polaron,
    kr,
    profile,
    bounds=([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf]),
)

figure, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))

ax.plot(
    kr,
    profile,
    marker="o",
    color=discrete_colors(1)[0],
    markersize=2,
    linestyle="None",
)
kr_ = np.linspace(kr.min(), kr.max(), num=1024)
ax.plot(kr_, polaron(kr_, *params), linestyle="-", color="k")

Aerr, rerr, Cerr = np.diag(np.sqrt(pcov))
A, r, C = params
ax.fill_between(
    kr_,
    y1=polaron(kr_, A - Aerr, r - rerr, C - Cerr),
    y2=polaron(kr_, A + Aerr, r + rerr, C + Cerr),
    facecolor="gray",
    edgecolor="k",
    linestyle="dashed",
    linewidth=1,
    alpha=0.2,
)

ax.set_xlim([kr.min(), kr.max()])
ax.set_ylabel("Diffuse intensity change [a.u.]")
ax.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")
plt.tight_layout()