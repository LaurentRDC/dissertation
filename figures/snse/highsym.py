import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
from skimage.filters import gaussian
from pathlib import Path
from plotutils import FIGURE_WIDTH, discrete_colors, tag_axis
from plotutils.snse_datasets import overnight4
import skued
from iris import DiffractionDataset

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
    (0, -4, 0),
    (0, -3, 1),
    (0, -3, -1),
    (0, -5, 1),
]

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def exponential(time, *args, **kwargs):
    return skued.exponential(time, *args, **kwargs)


timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["bg"] = np.zeros_like(timedelays)
    timeseries["Y"] = np.zeros_like(timedelays)
    timeseries["Z"] = np.zeros_like(timedelays)
    timeseries["T"] = np.zeros_like(timedelays)
    for indices in INDICES_DIFFUSE:
        h, k, l = indices

        yi, xi = overnight4.miller_to_arrindex(h, k, l)
        timeseries["bg"] += dset.time_series_selection(
            skued.RingSelection(
                shape=dset.resolution,
                center=(xi, yi),
                inner_radius=int(2.5 * OUTER_RADIUS),
                outer_radius=int(4 * OUTER_RADIUS),
            )
        )

        for sign in [1, -1]:
            yi_Y, xi_Y = overnight4.miller_to_arrindex(h, k + sign * (1 / 2), l)
            yi_Z, xi_Z = overnight4.miller_to_arrindex(h, k, l + sign * (1 / 2))
            yi_T, xi_T = overnight4.miller_to_arrindex(
                h, k + sign * (1 / 2), l + sign * (1 / 2)
            )
            timeseries["Y"] += dset.time_series_selection(
                skued.DiskSelection(
                    dset.resolution, center=(xi_Y, yi_Y), radius=OUTER_RADIUS
                )
            )
            timeseries["Z"] += dset.time_series_selection(
                skued.DiskSelection(
                    dset.resolution, center=(xi_Z, yi_Z), radius=OUTER_RADIUS
                )
            )
            timeseries["T"] += dset.time_series_selection(
                skued.DiskSelection(
                    dset.resolution, center=(xi_T, yi_T), radius=OUTER_RADIUS
                )
            )

# Normalize all time-series to pre-time-zero
for k, ts in timeseries.items():
    timeseries[k] /= np.mean(ts[timedelays < 0])


figure, axes = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(4.25, 5))

# background fit
bgparams, _ = opt.curve_fit(
    exponential,
    timedelays[timedelays < 20],
    timeseries["bg"][timedelays < 20],
    p0=(0, -0.02, 3, 1),
)
bg_curve = exponential(timedelays, *bgparams)


colors = discrete_colors(len(timeseries) - 1)
for ax, ts_name, color in zip(
    axes.flat,
    ["Y", "Z", "T"],
    colors,
):
    plot_params = dict(
        marker="o",
        markersize=2,
        color=color,
        label=ts_name,
        linestyle="None",
        elinewidth=0.5,
    )

    ax.errorbar(
        x=timedelays,
        y=timeseries[ts_name],
        yerr=scipy.stats.sem(timeseries[ts_name][timedelays < 0]),
        **plot_params,
    )

    # Fit average diffuse increase
    avparams, avpcov = opt.curve_fit(
        exponential,
        timedelays[timedelays < 20],
        timeseries[ts_name][timedelays < 20],
        p0=(0, -0.02, 3, 1),
    )
    avperr = np.sqrt(np.diag(avpcov))
    ax.plot(
        timedelays,
        exponential(timedelays, *avparams),
        linewidth=1,
        color=color,
    )

    ax.plot(
        timedelays, bg_curve, linewidth=1, linestyle="dashed", color="k", zorder=np.inf
    )

    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
    ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)
    ax.set_ylabel("$\Delta I/I_0$ [a.u.]")
    ax.set_yticks([1.00, 1.01])
    tag_axis(ax, ts_name, y=0.9, x=0.05, edgecolor="w")

axes[-1].set_xlim([-1.6, 12])
axes[-1].set_xlabel("Time-delay [ps]")

plt.tight_layout()
