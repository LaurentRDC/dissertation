import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
from skimage.filters import gaussian
from pathlib import Path
from plotutils import FIGURE_WIDTH, discrete_colors
from plotutils.snse_datasets import overnight4
import skued
from iris import DiffractionDataset

OUTER_RADIUS = 30

# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [(0, -1, 3), (0, 0, 4), (0, -2, 4), (0, -1, 5)]

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def exponential(time, *args, **kwargs):
    return skued.exponential(time, *args, **kwargs)


figure, ax_ts = plt.subplots(1, 1, figsize=(4.25, 4))

ax_ts.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["background"] = dset.time_series_selection(
        skued.RectSelection(dset.resolution, 680, 739, 1085, 1153)
    )

    timeseries["Y"] = np.zeros_like(timedelays)
    timeseries["Z"] = np.zeros_like(timedelays)
    timeseries["T"] = np.zeros_like(timedelays)
    for indices in INDICES_DIFFUSE:
        h, k, l = indices
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

colors = discrete_colors(len(timeseries))
for index, (ts_name, color, marker) in enumerate(
    zip(["background", "Y", "Z", "T"], colors, ["D", "*", "o", "^"])
):
    plot_params = dict(
        marker=marker,
        markersize=3,
        color=color,
        label=ts_name,
        linestyle="None",
        elinewidth=0.5,
    )

    vertical_offset = index * 0.015
    ax_ts.axhline(y=1 + vertical_offset, linestyle="dashed", color="k", linewidth=0.5)

    ax_ts.errorbar(
        x=timedelays,
        y=vertical_offset + timeseries[ts_name],
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
    ax_ts.plot(
        timedelays,
        vertical_offset + exponential(timedelays, *avparams),
        linewidth=1,
        color=color,
    )

ax_ts.legend(
    loc="center",
    ncol=len(timeseries),
    bbox_to_anchor=(0.5, 1.05),
    bbox_transform=ax_ts.transAxes,
    framealpha=1,
    edgecolor="none",
)

ax_ts.set_yticks([])
ax_ts.set_xlim([-1.6, 12])
ax_ts.set_xlabel("Time-delay [ps]")
ax_ts.set_ylabel("Scattering intensity change [a.u.]")

plt.tight_layout()
