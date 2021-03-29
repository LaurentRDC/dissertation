from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset
from plotutils.snse_datasets import overnight4
from plotutils import MEDIUM_FIGURE_WIDTH, discrete_colors

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

RADIUS = 30

# Determine the peak position and the time-series integration bounding box

INDICES_BRAGG_C = [
    (0, 0, 2),
    (0, 0, 1),
    (0, 0, 3),
    (0, 0, -3),
    (0, 0, -2),
    (0, 0, 4),
    (0, 0, -4),
    (0, 0, 6),
]

INDICES_BRAGG_B = [(0, 2, 0), (0, -2, 0), (0, -4, 0), (0, -6, 0)]


# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, *args, **kwargs)


timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["debye-waller-b"] = np.zeros_like(timedelays)
    for h, k, l in INDICES_BRAGG_B:
        yj, xj = overnight4.miller_to_arrindex(h, k, l)
        q2 = np.linalg.norm(CRYSTAL.scattering_vector((h, k, l))) ** 2
        timeseries["debye-waller-b"] += (
            dset.time_series_selection(
                skued.DiskSelection(
                    shape=dset.resolution,
                    center=(xj, yj),
                    radius=RADIUS,
                )
            )
            / q2
        )

    timeseries["debye-waller-c"] = np.zeros_like(timedelays)
    for h, k, l in INDICES_BRAGG_C:
        yj, xj = overnight4.miller_to_arrindex(h, k, l)
        q2 = np.linalg.norm(CRYSTAL.scattering_vector((h, k, l))) ** 2
        timeseries["debye-waller-c"] += (
            dset.time_series_selection(
                skued.DiskSelection(
                    shape=dset.resolution,
                    center=(xj, yj),
                    radius=RADIUS,
                )
            )
            / q2
        )

timeseries["debye-waller-b"] /= np.mean(timeseries["debye-waller-b"][timedelays < 0])
timeseries["debye-waller-c"] /= np.mean(timeseries["debye-waller-c"][timedelays < 0])

figure, diffuse_ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))
diffuse_ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
diffuse_ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

plot_params = dict(
    markersize=3,
    linestyle="None",
    elinewidth=0.5,
)

for k, label, color, marker in zip(
    ["debye-waller-b", "debye-waller-c"],
    [r"$\mathbf{q} ~ || ~ \mathbf{b}^*$", r"$\mathbf{q} ~ || ~ \mathbf{c}^*$"],
    discrete_colors(2),
    ["x", "^"],
):
    ts = timeseries[k]

    params, pcov = opt.curve_fit(
        biexponential,
        timedelays,
        ts,
        p0=(0, 0.01, 0.001, 0.2, 5, 1),
        bounds=([-1, -1, -1, 0, 3.6, 0.9], [1, 1, 1, 0.4, 20, 1.1]),
    )

    fit_curve = biexponential(timedelays, *params)

    diffuse_ax.plot(timedelays, fit_curve, linewidth=1, color=color)

    # Rolling average for display purposes ONLY
    diffuse_ax.errorbar(
        x=timedelays,
        y=ts,
        yerr=scipy.stats.sem(ts[timedelays < 0]),
        color=color,
        label=label,
        marker=marker,
        markersize=3,
        linestyle="None",
        elinewidth=0.5,
    )

diffuse_ax.legend(ncol=2, loc="upper right", edgecolor="none")


diffuse_ax.set_xlim([-1.6, 30])
diffuse_ax.set_yticks([0.99, 1.00])

diffuse_ax.set_xlabel("Time-delay [ps]")
diffuse_ax.set_ylabel("$\Delta I/I_0$ [a.u.]")

plt.tight_layout()
