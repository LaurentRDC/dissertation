from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import skued
from crystals import Crystal
from iris import DiffractionDataset
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis
from dissutils.snse import overnight4
from skimage.filters import gaussian

DATADIR = Path("data") / "snse"
CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")

INDICES_DIFFUSE = [
    (0, -1, 3),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -1, 7),
    (0, -3, 5),
    (0, -5, 3),
]

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, 0, *args, **kwargs)


def polaron(q, A, rp):
    return A * q * rp ** 3 / (1 + (q * rp) ** 2) ** 2


ks, amplitudes, amplitudes_err = np.loadtxt(
    DATADIR / "fast-diffuse.csv", delimiter=",", unpack=True
)

params, pcov = opt.curve_fit(
    polaron,
    ks,
    amplitudes,
    sigma=amplitudes_err,
    bounds=([-np.inf, 0], [np.inf, np.inf]),
)

figure, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))

ax.plot(
    ks,
    amplitudes,
    marker="o",
    color=discrete_colors(1)[0],
    markersize=2,
    linestyle="None",
)
ax.fill_between(
    ks,
    y1=amplitudes + amplitudes_err,
    y2=amplitudes - amplitudes_err,
    facecolor="gray",
    edgecolor="k",
    linestyle="dashed",
    linewidth=1,
    alpha=0.2,
)

ks_ = np.linspace(ks.min(), ks.max(), num=1024)
ax.plot(ks_, polaron(ks_, *params), linestyle="-", color="k")

ax.set_xlim([ks.min(), ks.max()])
ax.set_ylabel("Amplitude [a.u.]")
ax.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")

plt.tight_layout()