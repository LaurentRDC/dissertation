"""
Observation of time-series
"""
from math import sqrt, cos, sin, degrees
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import skued
from crystals import Crystal
from iris import DiffractionDataset
from matplotlib.ticker import FixedFormatter, FixedLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plotutils import FIGURE_WIDTH, FONTSIZE, SNSE_CAMERA_LENGTH, ImageGrid, tag_axis
from plotutils.snse_datasets import overnight4, static
from skimage.filters import gaussian

INPUT = Path("data") / "snse"

CRYSTAL = Crystal.from_cif(INPUT / "snse_pnma.cif")

# Static dataset stays the same, illustration purposes only
CROP = 500  # pixels to crop out of the static picture, on all sides.
COLORMAP = "CMRmap_r"

INNER_RADIUS = 15
OUTER_RADIUS = 30


# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [(0, -1, 3), (0, 0, 4), (0, -2, 4), (0, -1, 5)]
INDICES_BRAGG = [(0, 1, 2), (0, -1, -2), (0, 0, 1), (0, 0, -1), (0, 1, 0), (0, -1, 0)]
inferno = plt.get_cmap("inferno")
DWCOLOR, GAMMA_COLOR, AVERAGE_COLOR = inferno(30), inferno(100), inferno(200)


figure, ax_im = plt.subplots(1, 1, figsize=(4, 4))

divider = make_axes_locatable(ax_im)
cbar_ax = divider.append_axes("top", size=0.07, pad=0.05)

ax_im.xaxis.set_visible(False)
ax_im.yaxis.set_visible(False)

# Static diffraction pattern
with DiffractionDataset(static.path, mode="r") as dset:
    static_pattern = dset.diff_data(0)
    static_mask = dset.invalid_mask

static_pattern[:, 1024::] += 6  # Difference in bias between the two halves of the CCD
static_pattern[static_mask] = 0
static_pattern[:] = gaussian(static_pattern, 4)
# Removing some hot pixels
static_pattern = np.maximum(static_pattern, 0)

# Scaling between 0 and 1
static_pattern -= static_pattern.min()
static_pattern /= 200

kx, ky = static.kgrid()
dk = kx[0, 1] - kx[0, 0]
m = ax_im.imshow(
    static_pattern,
    vmin=0,
    vmax=1,
    cmap=COLORMAP,
)

# Index pattern
for (h, k, l), offset in zip([(0, -2, 0), (0, 0, -2)], [25, 25, 25]):
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


# -----------------------------------------------------------------------------
# Brillouin zones
# -----------------------------------------------------------------------------
r010, c010 = static.miller_to_arrindex(0, 1, 0)
r001, c001 = static.miller_to_arrindex(0, 0, 1)
r000, c000 = static.miller_to_arrindex(0, 0, 0)

bz_height = np.linalg.norm(np.array([r010 - r000, c010 - c000]))
bz_width = np.linalg.norm(np.array([r001 - r000, c001 - c000]))

v010 = np.array([r010, c010]) - np.array(static.center)
v010 /= np.linalg.norm(v010)
angle = np.arccos(np.dot(v010, np.array([0, 1])))
rotmat = np.array([[cos(angle), -sin(angle)], [sin(angle), cos(angle)]])

rj, cj = static.miller_to_arrindex(0, 0, -2)

xy = np.array([rj, cj]) - np.array([bz_width, bz_height]) / 2

ax_im.add_patch(
    mpatches.Rectangle(
        xy=np.array([rj, cj]) - rotmat @ np.array([bz_width / 2, bz_height / 2]),
        width=bz_width,
        height=bz_height,
        angle=degrees(angle),
        fc="none",
        ec="k",
        linestyle="solid",
        linewidth=0.1,
    )
)


# -----------------------------------------------------------------------------
# Show scattering vectors k, q, and G
# -----------------------------------------------------------------------------

rc, cc = static.center
rj, cj = static.miller_to_arrindex(0, 0, -2)

arrow_kwds = dict(
    arrowstyle="-|>", shrinkA=0.5, shrinkB=0.5, mutation_scale=6, fc="k", ec="k"
)

shadow_arrow_kwds = dict(
    arrowstyle="-",
    shrinkA=0,
    shrinkB=2,
    fc="k",
    ec="k",
    linestyle="dotted",
    linewidth=0.5,
    zorder=np.inf,
)
ax_im.add_patch(mpatches.FancyArrowPatch(posA=(rc, cc), posB=(rj, cj), **arrow_kwds))
ax_im.add_patch(
    mpatches.FancyArrowPatch(posA=(rc, cc), posB=(rj, cj), **shadow_arrow_kwds)
)
ax_im.text(
    x=(rc + rj) / 2, y=(cc + cj) / 2 + 20, s="$\mathbf{H}$", ha="center", va="top"
)

ax_im.add_patch(
    mpatches.FancyArrowPatch(posA=(rj, cj), posB=(rj + 20, cj - 100), **arrow_kwds)
)
ax_im.text(
    x=(2 * rj + 20) / 2 - 20,
    y=(2 * cj - 100) / 2,
    s="$\mathbf{k}_0$",
    ha="right",
    va="center",
)

ax_im.add_patch(
    mpatches.FancyArrowPatch(posA=(rc, cc), posB=(rj + 20, cj - 100), **arrow_kwds)
)
ax_im.add_patch(
    mpatches.FancyArrowPatch(
        posA=(rc, cc), posB=(rj + 20, cj - 100), **shadow_arrow_kwds
    )
)
ax_im.text(
    x=(rc + rj + 20) / 2,
    y=(cc + cj - 100) / 2,
    s="$\mathbf{q}$",
    ha="center",
    va="bottom",
)


# -----------------------------------------------------------------------------

plt.colorbar(mappable=m, cax=cbar_ax, orientation="horizontal")
cbar_ax.set_xlabel("Scattered intensity [a.u.]")
cbar_ax.xaxis.set_label_position("top")
cbar_ax.xaxis.set_ticks([0, 1])
cbar_ax.xaxis.tick_top()

ax_im.set_xlim([CROP, 2048 - CROP])
ax_im.set_ylim([2048 - CROP, CROP])

# -----------------------------------------------------------------------------
# Cover the beam block with a tasteful patch
# This makes it look like the beam block is actually part of the page
# NOTE: the following dimensions are for the static dataset ONLY, are
#       need to be pixel perfect with respect to the DPI
beamblock_indices = np.argwhere(static_mask > 0)
x, y = (903, 0)
width = 180
height = 1135

beamblock_patch = mpatches.Rectangle(
    xy=(x, y),
    width=width,
    height=height,
    edgecolor="k",
    facecolor="w",
)
move_in = 5
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

plt.tight_layout()
