from math import sqrt
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from skimage.filters import gaussian

from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis
from dissutils.snse import overnight4, static

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

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
    return skued.biexponential(time, *args, **kwargs)


with DiffractionDataset(static.path, mode="r") as static_dset:
    static_pattern = static_dset.diff_data(0)
    static_mask = static_dset.invalid_mask
static_pattern[:, 1024::] += 6  # Difference in bias between the two halves of the CCD
static_pattern[static_mask] = 0
static_pattern[:] = gaussian(static_pattern, 4)
# Removing some hot pixels
static_pattern = np.maximum(static_pattern, 0)

# Scaling between 0 and 1
static_pattern -= static_pattern.min()
static_pattern /= 200

figure, axes = plt.subplots(
    4,
    2,
    sharex="col",
    sharey="col",
    figsize=(LARGE_FIGURE_WIDTH, 6),
    gridspec_kw=dict(width_ratios=[3.5, 1]),
)

axs_trace = axes[:, 0]
axs_im = axes[:, 1]

timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    for r in [10, 15, 20, 25]:
        inner_radius = r
        outer_radius = 2 * r
        timeseries[r] = np.zeros_like(timedelays)
        for indices in INDICES_DIFFUSE:
            q2 = np.linalg.norm(CRYSTAL.scattering_vector(indices)) ** 2

            yi, xi = overnight4.miller_to_arrindex(*indices)
            diffuse_selection = skued.RingSelection(
                shape=dset.resolution,
                center=(xi, yi),
                inner_radius=inner_radius,
                outer_radius=outer_radius,
            )
            timeseries[r] += dset.time_series_selection(diffuse_selection) / q2

# Normalize all time-series to pre-time-zero
for k, ts in timeseries.items():
    timeseries[k] /= np.mean(ts[timedelays < 0])


kx, _ = overnight4.kgrid()
dk = kx[0, 1] - kx[0, 0]

for ax, ax_im, (r, ts), color in zip(
    axs_trace, axs_im, timeseries.items(), discrete_colors(len(timeseries))
):

    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
    ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

    params, pcov = opt.curve_fit(
        biexponential,
        timedelays[timedelays < 15],
        ts[timedelays < 15],
        p0=(0, -0.01, 0.01, 0.1, 3.7, 1),
    )
    ax.plot(timedelays, biexponential(timedelays, *params), linewidth=1, color=color)

    ax.errorbar(
        x=timedelays,
        y=ts,
        yerr=scipy.stats.sem(ts[timedelays < 0]),
        marker="o",
        markersize=2,
        color=color,
        linestyle="None",
        elinewidth=0.5,
    )
    ax.set_ylabel("$\Delta I/I_0$ [a.u.]")
    ax.xaxis.set_visible(False)
    tag_axis(
        ax,
        f"$r= {r * dk:.3f} ~ \AA^{{-1}}$",
        y=0.9,
        x=0.95,
        horizontalalignment="right",
        edgecolor="w",
    )

    kx, _ = static.kgrid()
    dk = kx[0, 1] - kx[0, 0]
    yc, xc = static.center
    yi_, xi_ = static.miller_to_arrindex(0, 0, 2)

    ax_im.xaxis.set_visible(False)
    ax_im.yaxis.set_visible(False)
    ax_im.imshow(static_pattern, vmin=0, vmax=1, cmap="CMRmap_r")
    ax_im.set_xlim([yi_ - 75, yi_ + 75])
    ax_im.set_ylim([xi_ - 75, xi_ + 75])

    inner_circ, outer_circ = skued.RingSelection(
        shape=static_pattern.shape,
        center=(xi_, yi_),
        inner_radius=r,
        outer_radius=2 * r,
    ).mpatch(edgecolor="k", linewidth=1, fill=False)
    ax_im.add_patch(inner_circ)
    ax_im.add_patch(outer_circ)


ax_im = axs_im[-1]
ax_im.xaxis.set_major_locator(
    FixedLocator(locs=[yi_ + (-1 / 3) / dk, yi_, yi_ + (1 / 3) / dk])
)
ax_im.xaxis.set_major_formatter(FixedFormatter(["-⅓", "0", "⅓"]))
ax_im.xaxis.set_visible(True)
ax_im.set_xlabel("$|\mathbf{k}|$ [$\AA^{-1}$]")
for x in [yi_ + (-1 / 3) / dk, yi_, yi_ + (1 / 3) / dk]:
    ax_im.axvline(
        x=x, linewidth=0.5, linestyle="--", color="k", ymax=4.57, clip_on=False
    )

axs_trace[0].set_xlim([-1.6, 15])
axs_trace[-1].xaxis.set_visible(True)
axs_trace[-1].set_xlabel("Time-delay [ps]")

plt.tight_layout()
