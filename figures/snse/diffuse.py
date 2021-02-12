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
from plotutils import FIGURE_WIDTH, FONTSIZE, SNSE_CAMERA_LENGTH, ImageGrid, tag_axis
from plotutils.snse_datasets import overnight4, static
from skimage.filters import gaussian

INPUT = Path("data") / "snse"

CRYSTAL = Crystal.from_cif(INPUT / "snse_pnma.cif")

# Static dataset stays the same, illustration purposes only
CROP = 300  # pixels to crop out of the static picture, on all sides.
COLORMAP = "CMRmap_r"

INNER_RADIUS = 15
OUTER_RADIUS = 30


# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [(0, -1, 3), (0, 0, 4), (0, -2, 4), (0, -1, 5)]
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


figure = plt.figure(figsize=(FIGURE_WIDTH, 0.42 * FIGURE_WIDTH))

# The following shenanigans are based on
#   https://matplotlib.org/3.1.3/gallery/subplots_axes_and_figures/gridspec_nested.html#sphx-glr-gallery-subplots-axes-and-figures-gridspec-nested-py
gs0 = gridspec.GridSpec(1, 3, width_ratios=[2.2, 1, 2.2], figure=figure)

ax_im = figure.add_subplot(gs0[0])
diffuse_ax = figure.add_subplot(gs0[-1])

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1])

diff_context_ax = figure.add_subplot(gs00[0])
ax_linecut = figure.add_subplot(gs00[1])

divider = make_axes_locatable(ax_im)
cbar_ax = divider.append_axes("bottom", size=0.07, pad=0.05)

# Diffraction image is also used to show where time-series are selected.
for ax in (ax_im, diff_context_ax):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

diffuse_ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
diffuse_ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

# Static diffraction pattern
with DiffractionDataset(static.path, mode="r") as dset:
    static_pattern = dset.diff_data(0)
    static_mask = dset.invalid_mask
    static_pattern[static_mask] = 0
static_pattern[:] = gaussian(static_pattern, 4)
# Removing some hot pixels
static_pattern = np.maximum(static_pattern, 0)

# Scaling between 0 and 1
static_pattern -= static_pattern.min()
static_pattern /= 200

# -----------------------------------------------------------------------------
#   ACQUIRE TIME-SERIES
# -----------------------------------------------------------------------------
timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["average"] = dset.time_series_selection(
        skued.RectSelection(dset.resolution, 680, 739, 1085, 1153)
    )

    timeseries["debye-waller"] = np.zeros_like(timedelays)
    for h, k, l in INDICES_BRAGG:
        yj, xj = overnight4.miller_to_arrindex(h, k, l)
        timeseries["debye-waller"] += dset.time_series_selection(
            skued.RingSelection(
                shape=dset.resolution,
                center=(xj, yj),
                inner_radius=INNER_RADIUS,
                outer_radius=OUTER_RADIUS,
            )
        )

    for indices in INDICES_DIFFUSE:
        yi, xi = overnight4.miller_to_arrindex(*indices)
        diffuse_selection = skued.RingSelection(
            shape=dset.resolution,
            center=(xi, yi),
            inner_radius=INNER_RADIUS,
            outer_radius=OUTER_RADIUS,
        )
        timeseries[indices] = dset.time_series_selection(diffuse_selection)

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

# Diffuse rises still contain nearby Bragg dynamics
for indices in INDICES_DIFFUSE:
    timeseries[indices] -= timeseries["debye-waller"]
    timeseries[indices] += 1

m = ax_im.imshow(static_pattern, vmin=0, vmax=1, cmap=COLORMAP)

# Index pattern
for (h, k, l), offset in zip([(0, -2, 0), (0, 1, -1), (0, 0, 2)], [25, 25, 75]):
    rj, cj = static.miller_to_arrindex(h, k, l)
    ax_im.text(
        s=skued.indices_to_text(h, k, l),
        x=rj,
        y=cj + offset,
        fontsize=FONTSIZE - 2,
        transform=ax_im.transData,
        verticalalignment="top",
        horizontalalignment="center",
    )

plt.colorbar(mappable=m, cax=cbar_ax, orientation="horizontal")
cbar_ax.set_xlabel("Scattered intensity [a.u.]")
cbar_ax.xaxis.set_label_position("bottom")
cbar_ax.xaxis.set_ticks([0, 1])
cbar_ax.xaxis.tick_bottom()

# -----------------------------------------------------------------------------
# fit and plot Debye-Waller dynamics
# -----------------------------------------------------------------------------

dwparams, dwpcov = opt.curve_fit(
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
    label="Bragg (1)",
    linestyle="None",
    elinewidth=0.5,
)
diffuse_ax.plot(timedelays, dwfit_curve, linewidth=1, color=DWCOLOR)

# dwerr = np.std(timeseries["debye-waller"][timedelays < 0])
# diffuse_ax.fill_between(
#     timedelays,
#     timeseries["debye-waller"] - dwerr,
#     timeseries["debye-waller"] + dwerr,
#     where=timedelays > 0,
#     color=DWCOLOR,
#     alpha=0.2,
# )

diffuse_ax.errorbar(
    x=timedelays,
    y=timeseries["debye-waller"],
    yerr=scipy.stats.sem(timeseries["debye-waller"][timedelays < 0]),
    **plot_params,
)


# -----------------------------------------------------------------------------
# Combining all Gamma points
# -----------------------------------------------------------------------------

gamma_timeseries = np.zeros_like(timedelays)
for indices in INDICES_DIFFUSE:
    # The amplitude of increase should scale like q^2
    q2 = np.linalg.norm(CRYSTAL.scattering_vector(indices)) ** 2
    gamma_timeseries += timeseries[indices] / q2
gamma_timeseries /= np.mean(gamma_timeseries[timedelays < 0])

# Fit the renormalization to rise + fall
params, pcov = opt.curve_fit(
    biexponential,
    timedelays[timedelays < 15],
    gamma_timeseries[timedelays < 15],
    p0=(0, -0.01, 0.01, 0.1, 3.7, 1),
)
perr = np.sqrt(np.diag(pcov))

diffuse_ax.plot(
    timedelays, biexponential(timedelays, *params), linewidth=1, color=GAMMA_COLOR
)
# diffuse_ax.fill_between(
#     timedelays,
#     biexponential(timedelays, *(params - perr)),
#     biexponential(timedelays, *(params + perr)),
#     where=timedelays > 0,
#     color=GAMMA_COLOR,
#     alpha=0.2,
# )
diffuse_ax.errorbar(
    x=timedelays,
    y=gamma_timeseries,
    yerr=scipy.stats.sem(gamma_timeseries[timedelays < 0]),
    marker="o",
    markersize=2,
    color=GAMMA_COLOR,
    linestyle="None",
    label="$\Gamma$ (2)",
    elinewidth=0.5,
)

# -----------------------------------------------------------------------------
# Determine the average diffuse increase
# -----------------------------------------------------------------------------
plot_params = dict(
    marker="D",
    markersize=2,
    color=AVERAGE_COLOR,
    label="bg (3)",
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
avparams, avpcov = opt.curve_fit(
    exponential,
    timedelays[timedelays < 20],
    timeseries["average"][timedelays < 20],
    p0=(0, -0.02, 3, 1),
)
avperr = np.sqrt(np.diag(avpcov))
diffuse_ax.plot(
    timedelays, exponential(timedelays, *avparams), linewidth=1, color=AVERAGE_COLOR
)

# -----------------------------------------------------------------------------

# Show the shape of the selections use to contruct time-series
yi_, xi_ = static.miller_to_arrindex(0, 0, 2)
width = int(1.4 * OUTER_RADIUS)
diffuse_selection = skued.RingSelection(
    shape=(2048, 2048),
    center=(xi_, yi_),
    inner_radius=INNER_RADIUS,
    outer_radius=OUTER_RADIUS,
)
inner_circ, outer_circ = diffuse_selection.mpatch(
    edgecolor="k", linewidth=1, fill=False
)
diff_context_ax.imshow(static_pattern, vmin=0, vmax=1, cmap=COLORMAP)
diff_context_ax.set_ylim([xi_ - width, xi_ + width])
diff_context_ax.set_xlim([yi_ - width, yi_ + width])
diff_context_ax.add_patch(inner_circ)
diff_context_ax.add_patch(outer_circ)
diff_context_ax.axhline(y=xi_, linewidth=0.5, color="k")

step = 1
for index, label in enumerate(["(1)", "(2)", "(3)"]):
    diff_context_ax.text(
        s=label,
        color="w" if not index else "k",
        x=yi_ + index * step * INNER_RADIUS,
        y=xi_ + index * step * INNER_RADIUS,
        fontsize=FONTSIZE - 2,
        horizontalalignment="center",
        verticalalignment="center",
    )


# Draw connection between main figure and inset
mark_inset(
    ax_im,
    diff_context_ax,
    loc1=2,
    loc2=3,
    fc="none",
    ec="k",
    linewidth=0.5,
    linestyle="--",
    zorder=np.inf,
)

# -----------------------------------------------------------------------------
# Line cut to better describe the integration region
yi_, xi_ = static.miller_to_arrindex(0, 0, 2)
width = int(1.4 * OUTER_RADIUS)
kx, _ = static.kgrid()
dk = kx[0, 1] - kx[0, 0]
k = dk * (np.arange(start=yi_ - width, stop=yi_ + width) - yi_)  # inverse angstroms
linecut = static_pattern[xi_, yi_ - width : yi_ + width]

# Fit to extract width
def peak(f, A, c, fwhm_g, fwhm_l, o):
    return A * skued.pseudo_voigt(f, center=c, fwhm_g=fwhm_g, fwhm_l=fwhm_l) + o


params, pcov = opt.curve_fit(
    peak,
    k,
    linecut,
    p0=(1, 0, 0.05, 0.05, 0),
)

ax_linecut.plot(k, peak(k, *params), color="k", linewidth=0.5)
ax_linecut.scatter(x=k, y=linecut, s=5, c=plt.get_cmap(COLORMAP)(linecut))

ax_linecut.set_xlim([-dk * width, dk * width])
for r in [dk * INNER_RADIUS, -dk * INNER_RADIUS, dk * OUTER_RADIUS, -dk * OUTER_RADIUS]:
    ax_linecut.axvline(
        x=r, linewidth=0.5, linestyle="--", color="k", ymax=1.65, clip_on=False
    )

ax_linecut.yaxis.set_visible(False)
ax_linecut.set_xlabel("$|\mathbf{k}|$ [$\AA^{-1}$]")
ax_linecut.xaxis.set_major_locator(FixedLocator(locs=[-1 / 4, 0, 1 / 4]))
ax_linecut.xaxis.set_major_formatter(FixedFormatter(["-¼", "0", "¼"]))

# -----------------------------------------------------------------------------


diffuse_ax.set_xlim([-1, 12])
diffuse_ax.set_ylim([0.984, 1.02])

diffuse_ax.legend(
    loc="center",
    ncol=1,
    bbox_to_anchor=(0.75, 0.3),
    bbox_transform=diffuse_ax.transAxes,
    fontsize=FONTSIZE - 2,
    edgecolor="none",
)

diffuse_ax.set_xlabel("Time-delay [ps]")
diffuse_ax.set_ylabel("Scattered intensity change [a.u.]")
diffuse_ax.yaxis.set_label_position("right")
diffuse_ax.yaxis.tick_right()
diffuse_ax.yaxis.set_ticks([0.99, 1.00, 1.01, 1.02])

ax_im.set_xlim([CROP, 2048 - CROP])
ax_im.set_ylim([2048 - CROP, CROP])

# -----------------------------------------------------------------------------
# Cover the beam block with a tasteful patch
# This makes it look like the beam block is actually part of the page
# NOTE: the following dimensions are for the static dataset ONLY, are
#       need to be pixel perfect with respect to the DPI
beamblock_indices = np.argwhere(static_mask > 0)
x, y = (903, 0)
width = 173
height = 1135

beamblock_patch = mpatches.Rectangle(
    xy=(x, y),
    width=width,
    height=height,
    edgecolor="k",
    facecolor="w",
)
move_in = 7
crossover_patch = mpatches.Rectangle(
    xy=(x + move_in, y - 10),
    width=width - 2 * move_in,
    height=height + 10 - move_in,
    fill=True,
    color="w",
    zorder=10,
    clip_on=False,
)
ax_im.add_patch(beamblock_patch)
ax_im.add_patch(crossover_patch)

tag_axis(ax_im, "a)")
tag_axis(diff_context_ax, "b)", x=0.1, y=0.9)
tag_axis(ax_linecut, "c)", x=0.1, y=0.9)
tag_axis(diffuse_ax, "d)")

plt.subplots_adjust(
    top=0.98, bottom=0.19, left=0.02, right=0.89, hspace=0.18, wspace=0.07
)
