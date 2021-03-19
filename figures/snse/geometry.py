from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import skued
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from plotutils import FONTSIZE, tag_axis
from plotutils.snse_datasets import static
from skimage.filters import gaussian

INPUT = Path("data") / "snse"

CROP = 300  # pixels to crop out of the static picture, on all sides.
COLORMAP = "CMRmap_r"

INNER_RADIUS = 15
OUTER_RADIUS = 30


figure = plt.figure(figsize=(5, 3.5))

# The following shenanigans are based on
#   https://matplotlib.org/3.1.3/gallery/subplots_axes_and_figures/gridspec_nested.html#sphx-glr-gallery-subplots-axes-and-figures-gridspec-nested-py
gs0 = gridspec.GridSpec(1, 2, figure=figure, width_ratios=[2.2, 1])
gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0, 1])

ax_im = figure.add_subplot(gs0[0])
diff_context_ax = figure.add_subplot(gs00[0])
ax_linecut = figure.add_subplot(gs00[1], sharex=diff_context_ax)

divider = make_axes_locatable(ax_im)
cbar_ax = divider.append_axes("top", size=0.07, pad=0.05)

ax_im.xaxis.set_visible(False)
ax_im.yaxis.set_visible(False)

# Static diffraction pattern
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

# -----------------------------------------------------------------------------
# Static diffraction pattern
# -----------------------------------------------------------------------------

kx, ky = static.kgrid()
m = ax_im.imshow(
    static_pattern,
    vmin=0,
    vmax=1,
    cmap=COLORMAP,
    extent=[kx.min(), kx.max(), ky.max(), ky.min()],
)

plt.colorbar(mappable=m, cax=cbar_ax, orientation="horizontal")
cbar_ax.set_xlabel("Scattered intensity [a.u.]")
cbar_ax.xaxis.set_label_position("top")
cbar_ax.xaxis.set_ticks([0, 1])
cbar_ax.xaxis.tick_top()

# -----------------------------------------------------------------------------
# Inset
# -----------------------------------------------------------------------------
width = int(1.4 * OUTER_RADIUS)
kx, _ = static.kgrid()
dk = kx[0, 1] - kx[0, 0]

# Show the shape of the selections use to contruct time-series
yc, xc = static.center
yi_, xi_ = static.miller_to_arrindex(0, 0, 2)
kx_, ky_ = dk * (yi_ - yc), dk * (xi_ - xc)

diffuse_selection = skued.RingSelection(
    shape=static_pattern.shape,
    center=(ky_, kx_),
    inner_radius=dk * INNER_RADIUS,
    outer_radius=dk * OUTER_RADIUS,
)
inner_circ, outer_circ = diffuse_selection.mpatch(
    edgecolor="k", linewidth=1, fill=False
)
diff_context_ax.imshow(
    static_pattern,
    vmin=0,
    vmax=1,
    cmap=COLORMAP,
    extent=[kx.min(), kx.max(), ky.max(), ky.min()],
)
diff_context_ax.set_xlim([kx_ - dk * width, kx_ + dk * width])
diff_context_ax.set_ylim([ky_ - dk * width, ky_ + dk * width])
diff_context_ax.add_patch(inner_circ)
diff_context_ax.add_patch(outer_circ)
diff_context_ax.axhline(y=ky_, linewidth=0.5, color="k")

diff_context_ax.yaxis.set_visible(False)
diff_context_ax.set_xlabel("$|\mathbf{k}|$ [$\AA^{-1}$]")
diff_context_ax.xaxis.tick_top()
diff_context_ax.xaxis.set_label_position("top")
diff_context_ax.xaxis.set_major_locator(
    FixedLocator(locs=[kx_ + -1 / 4, kx_, kx_ + 1 / 4])
)
diff_context_ax.xaxis.set_major_formatter(FixedFormatter(["-¼", "0", "¼"]))

step = 1
for index, label in enumerate(["(1)", "(2)", "(3)"]):
    diff_context_ax.text(
        s=label,
        color="w" if not index else "k",
        x=kx_ + index * step * dk * INNER_RADIUS,
        y=ky_ + index * step * dk * INNER_RADIUS,
        fontsize=FONTSIZE - 1,
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
# Line cut
# -----------------------------------------------------------------------------

yi_, xi_ = static.miller_to_arrindex(0, 0, 2)
width = int(1.4 * OUTER_RADIUS)
kx, _ = static.kgrid()
dk = kx[0, 1] - kx[0, 0]
kc = dk * (yi_ - yc)

k = kc + dk * (
    np.arange(start=xi_ - width, stop=xi_ + width) - xi_
)  # inverse angstroms
linecut = static_pattern[xi_ - width : xi_ + width, yi_]

# Fit to extract width
def peak(f, A, c, fwhm_g, fwhm_l, o):
    return A * skued.pseudo_voigt(f, center=c, fwhm_g=fwhm_g, fwhm_l=fwhm_l) + o


params, pcov = opt.curve_fit(
    peak,
    k,
    linecut,
    p0=(1, kc, 0.05, 0.05, 0),
)

ax_linecut.plot(k, peak(k, *params), color="k", linewidth=0.5)
ax_linecut.scatter(x=k, y=linecut, s=5, c=plt.get_cmap(COLORMAP)(linecut))

ax_linecut.set_xlim([k.min(), k.max()])
for r in [dk * INNER_RADIUS, -dk * INNER_RADIUS, dk * OUTER_RADIUS, -dk * OUTER_RADIUS]:
    ax_linecut.axvline(
        x=kc + r, linewidth=0.5, linestyle="--", color="k", ymax=1.6, clip_on=False
    )

ax_linecut.xaxis.set_visible(False)
ax_linecut.yaxis.set_visible(False)
# -----------------------------------------------------------------------------

yc, xc = static.center
width = dk * (yc - CROP)
ax_im.set_xlim([-width, width])
ax_im.set_ylim([width, -width])

# -----------------------------------------------------------------------------
# Beam block patch
# -----------------------------------------------------------------------------

beamblock_indices = np.argwhere(static_mask > 0)
width = 0.15
height = 0.6
x, y = 0.5 - width / 2, 1


beamblock_patch = mpatches.Rectangle(
    xy=(x, y),
    width=width,
    height=-height,
    edgecolor="k",
    facecolor="w",
    transform=ax_im.transAxes,
)
move_in = 0.005
crossover_patch = mpatches.Rectangle(
    xy=(x + move_in, y + move_in),
    width=width - 2 * move_in,
    height=-height + 2 * move_in,
    fill=True,
    color="w",
    zorder=np.inf,
    clip_on=False,
    transform=ax_im.transAxes,
)
ax_im.add_patch(beamblock_patch)
ax_im.add_patch(crossover_patch)

tag_axis(ax_im, "a)")
tag_axis(diff_context_ax, "b)", x=0.1, y=0.9)
tag_axis(ax_linecut, "c)", x=0.1, y=0.9)

plt.subplots_adjust(
    top=0.86, bottom=0.07, left=0.06, right=0.9, hspace=0.115, wspace=0.065
)
