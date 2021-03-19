"""
Observation of time-series
"""
from math import sqrt
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from plotutils import discrete_colors, tag_axis
from plotutils.snse_datasets import overnight4, static
from skimage.filters import gaussian

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


figure, axes = plt.subplots(4, 1, sharex=True, sharey=True, figsize=(4.25, 5))


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
                inner_radius=r,
                outer_radius=2 * r,
            )
            timeseries[r] += dset.time_series_selection(diffuse_selection) / q2

# Normalize all time-series to pre-time-zero
for k, ts in timeseries.items():
    timeseries[k] /= np.mean(ts[timedelays < 0])


kx, _ = static.kgrid()
dk = kx[0, 1] - kx[0, 0]

for ax, (r, ts), color in zip(
    axes, timeseries.items(), discrete_colors(len(timeseries))
):

    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
    ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

    params, pcov = opt.curve_fit(
        biexponential,
        timedelays[timedelays < 15],
        ts[timedelays < 15],
        p0=(0, -0.01, 0.01, 0.1, 3.7, 1),
    )
    # perr = np.sqrt(np.diag(pcov))
    # tau1 = (params[3], perr[3])
    # tau2 = (params[4], perr[4])

    # ax.text(x=1.02, y=2/3, s=f'fast: ${tau1[0]:.3f} \pm {tau1[1]:.3f}$ ps', ha='left', va='center', transform=ax.transAxes)
    # ax.text(x=1.02, y=1/3, s=f'slow: ${tau2[0]:.3f} \pm {tau2[1]:.3f}$ ps', ha='left', va='center', transform=ax.transAxes)

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


axes[0].set_xlim([-1.6, 15])

axes[-1].xaxis.set_visible(True)
axes[-1].set_xlabel("Time-delay [ps]")


plt.tight_layout()
