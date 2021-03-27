"""
Observation of time-series
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset
from plotutils.snse_datasets import overnight4

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

INNER_RADIUS = 15
OUTER_RADIUS = 30


# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [
    (0, -1, 3),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -1, 7),
    (0, -3, 5),
    (0, -5, 3),
]
INDICES_BRAGG = [(0, 1, 2), (0, -1, -2), (0, 0, 1), (0, 0, -1), (0, 1, 0), (0, -1, 0)]

inferno = plt.get_cmap("inferno")
DWCOLOR, GAMMA_COLOR, AVERAGE_COLOR = inferno(30), inferno(100), inferno(200)

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def exponential(time, *args, **kwargs):
    return skued.exponential(time, *args, **kwargs)


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, *args, **kwargs)


figure, diffuse_ax = plt.subplots(1, 1, figsize=(4.25, 3))

diffuse_ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
diffuse_ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["debye-waller"] = np.zeros_like(timedelays)
    for h, k, l in INDICES_BRAGG:
        q2 = np.linalg.norm(CRYSTAL.scattering_vector((h, k, l))) ** 2
        yj, xj = overnight4.miller_to_arrindex(h, k, l)
        timeseries["debye-waller"] += (
            dset.time_series_selection(
                skued.DiskSelection(
                    shape=dset.resolution,
                    center=(xj, yj),
                    radius=OUTER_RADIUS,
                )
            )
            / q2
        )

    timeseries["average"] = np.zeros_like(timedelays)
    timeseries["gamma"] = np.zeros_like(timedelays)
    for indices in INDICES_DIFFUSE:
        q2 = np.linalg.norm(CRYSTAL.scattering_vector(indices)) ** 2

        yi, xi = overnight4.miller_to_arrindex(*indices)
        diffuse_selection = skued.RingSelection(
            shape=dset.resolution,
            center=(xi, yi),
            inner_radius=INNER_RADIUS,
            outer_radius=OUTER_RADIUS,
        )
        timeseries["gamma"] += dset.time_series_selection(diffuse_selection) / q2

        bg_selection = skued.RingSelection(
            shape=dset.resolution,
            center=(xi, yi),
            inner_radius=int(2.5 * OUTER_RADIUS),
            outer_radius=int(4 * OUTER_RADIUS),
        )
        timeseries["average"] += dset.time_series_selection(bg_selection) / q2

# Normalize all time-series to pre-time-zero
for k, ts in timeseries.items():
    timeseries[k] /= np.mean(ts[timedelays < 0])

# Diffuse rises still contain nearby Bragg dynamics
timeseries["gamma"] -= timeseries["debye-waller"]
timeseries["gamma"] += 1

# -----------------------------------------------------------------------------
# fit and plot Debye-Waller dynamics
# -----------------------------------------------------------------------------

dwparams, _ = opt.curve_fit(
    biexponential,
    timedelays[timedelays < 20],
    timeseries["debye-waller"][timedelays < 20],
    p0=(0, 0.01, 0.02, 0.2, 3.7, timeseries["debye-waller"].min()),
    bounds=([-0.2, -1, -1, 0, 1, 0.9], [0.2, 1, 1, 1, 5, 1.1]),
)

dwfit_curve = biexponential(timedelays, *dwparams)
plot_params = dict(
    marker="x",
    markersize=3,
    color=DWCOLOR,
    label="(1) Bragg",
    linestyle="None",
    elinewidth=0.5,
)
diffuse_ax.plot(timedelays, dwfit_curve, linewidth=1, color=DWCOLOR)

diffuse_ax.errorbar(
    x=timedelays,
    y=timeseries["debye-waller"],
    yerr=scipy.stats.sem(timeseries["debye-waller"][timedelays < 0]),
    **plot_params,
)


# -----------------------------------------------------------------------------
# Gamma points
# -----------------------------------------------------------------------------

# Fit the renormalization to rise + fall
params, _ = opt.curve_fit(
    biexponential,
    timedelays[timedelays < 15],
    timeseries["gamma"][timedelays < 15],
    p0=(0, -0.01, 0.01, 0.1, 3.7, 1),
)

diffuse_ax.plot(
    timedelays, biexponential(timedelays, *params), linewidth=1, color=GAMMA_COLOR
)

diffuse_ax.errorbar(
    x=timedelays,
    y=timeseries["gamma"],
    yerr=scipy.stats.sem(timeseries["gamma"][timedelays < 0]),
    marker="o",
    markersize=2,
    color=GAMMA_COLOR,
    linestyle="None",
    label="(2) $\Gamma$",
    elinewidth=0.5,
)

# -----------------------------------------------------------------------------
# Determine the average diffuse increase
# -----------------------------------------------------------------------------
plot_params = dict(
    marker="D",
    markersize=2,
    color=AVERAGE_COLOR,
    label="(3) background",
    linestyle="None",
    elinewidth=0.5,
)

diffuse_ax.errorbar(
    x=timedelays,
    y=timeseries["average"],
    yerr=scipy.stats.sem(timeseries["average"][timedelays < 0]),
    **plot_params,
)

# Fit average diffuse increase
avparams, _ = opt.curve_fit(
    exponential,
    timedelays[timedelays < 20],
    timeseries["average"][timedelays < 20],
    p0=(0, -0.02, 3, 1),
)
diffuse_ax.plot(
    timedelays, exponential(timedelays, *avparams), linewidth=1, color=AVERAGE_COLOR
)


diffuse_ax.set_xlim([-1.6, 15])
diffuse_ax.set_ylim([0.98, 1.022])

diffuse_ax.legend(
    loc="center",
    ncol=1,
    bbox_to_anchor=(0.75, 0.3),
    bbox_transform=diffuse_ax.transAxes,
    edgecolor="none",
)

diffuse_ax.set_xlabel("Time-delay [ps]")
diffuse_ax.set_ylabel("Scattered intensity change [a.u.]")
diffuse_ax.yaxis.set_ticks([0.98, 0.99, 1.00, 1.01, 1.02])

plt.tight_layout()
