from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset

from plotutils import FIGURE_WIDTH, discrete_colors, tag_axis
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

figure, (ax_b, ax_c) = plt.subplots(2, 1, sharex=True, figsize=(4.25, 4))

for ax, ts, color, label in zip(
    [ax_b, ax_c],
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
        linestyle="None",
        elinewidth=0.5,
    )

    tag_axis(
        ax, text=label, y=0.95, x=0.975, horizontalalignment="right", edgecolor="w"
    )
    ax.set_ylabel("$\Delta I/I_0$ [a.u.]")

ax_b.set_xlim([-1.6, 12])

ax_b.xaxis.set_visible(False)
ax_c.set_xlabel("Time-delay [ps]")

plt.tight_layout()
