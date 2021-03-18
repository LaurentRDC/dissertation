from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset

from plotutils import FIGURE_WIDTH, discrete_colors
from plotutils.snse_datasets import overnight4

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

INNER_RADIUS = 15
OUTER_RADIUS = 30

INDICES_DIFFUSE_C = [
    (0, -1, 3),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -1, 7),
]

INDICES_DIFFUSE_B = [
    (0, -4, 0),
    (0, -3, 1),
    (0, -3, -1),
    (0, -5, 1),
]

timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["gamma-c"] = np.zeros_like(timedelays)
    for indices in INDICES_DIFFUSE_C:
        q2 = np.linalg.norm(CRYSTAL.scattering_vector(indices)) ** 2

        yi, xi = overnight4.miller_to_arrindex(*indices)
        diffuse_selection = skued.RingSelection(
            shape=dset.resolution,
            center=(xi, yi),
            inner_radius=INNER_RADIUS,
            outer_radius=OUTER_RADIUS,
        )
        timeseries["gamma-c"] += dset.time_series_selection(diffuse_selection) / q2

    timeseries["gamma-b"] = np.zeros_like(timedelays)
    for indices in INDICES_DIFFUSE_B:
        q2 = np.linalg.norm(CRYSTAL.scattering_vector(indices)) ** 2

        yi, xi = overnight4.miller_to_arrindex(*indices)
        diffuse_selection = skued.RingSelection(
            shape=dset.resolution,
            center=(xi, yi),
            inner_radius=INNER_RADIUS,
            outer_radius=OUTER_RADIUS,
        )
        timeseries["gamma-b"] += dset.time_series_selection(diffuse_selection) / q2


# Normalize all time-series to pre-time-zero
for k, ts in timeseries.items():
    timeseries[k] /= np.mean(ts[timedelays < 0])

figure, ax = plt.subplots(1, 1, figsize=(4, 3))

for ts, color, label in zip(
    ["gamma-b", "gamma-c"],
    discrete_colors(2),
    [
        "$\mathbf{q} ~ || ~ \mathbf{b}^{\star}$",
        "$\mathbf{q} ~ || ~ \mathbf{c}^{\star}$",
    ],
):
    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
    ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)
    ax.errorbar(
        x=timedelays,
        y=timeseries[ts],
        yerr=scipy.stats.sem(timeseries[ts][timedelays < 0]),
        marker="o",
        markersize=2,
        color=color,
        label=label,
        linestyle="None",
        elinewidth=0.5,
    )


ax.legend(ncol=2, loc="center", edgecolor="none", bbox_to_anchor=(0.5, 1.1))
ax.set_xlim([-1, 10])
ax.set_ylim([0.985, 1.01])

ax.set_xlabel("Time-delay [ps]")
ax.set_ylabel("$\Delta I/I_0$ [a.u.]")

plt.tight_layout()
