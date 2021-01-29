# -*- coding: utf-8 -*-
import itertools as it
import sys
from pathlib import Path

import matplotlib.colors as cm
import matplotlib.pyplot as plt
import npstreams as ns
import numpy as np
import scipy.interpolate as interpolate
from crystals import Crystal
from matplotlib.ticker import FixedFormatter, FixedLocator
from mpl_toolkits.axes_grid1 import inset_locator, make_axes_locatable
from skimage.filters import gaussian

from plotutils import (
    FIGURE_WIDTH,
    ImageGrid,
    draw_hexagon_field,
    tag_axis,
)

MODE_ORDERING = {
    "LA": 0,
    "TA": 1,
    "ZA": 2,
    "LO1": 3,
    "LO2": 4,
    "LO3": 5,
    "TO1": 6,
    "TO2": 7,
    "TO3": 8,
    "ZO1": 9,
    "ZO2": 10,
    "ZO3": 11,
}
MODES = sorted(MODE_ORDERING.keys())
IN_PLANE_MODES = sorted(set(MODE_ORDERING.keys()) - {"ZA", "ZO1", "ZO2", "ZO3"})

INPUT = Path("data") / "graphite"


def oneph_weighted(mode):
    """ |F_{1j}|^2 / \omega_{j, k} """
    F1j = np.load(INPUT / "oneph" / f"{mode}_oneph.npy")
    omega = np.load(INPUT / "oneph" / f"{mode}_freq.npy")
    return F1j / omega


def threshold_oneph(threshold):
    """
    Show locations in reciprocal space where one one-ph structure factor dominates all others

    Returns
    -------
    majority : ndarray, shape (N, M)
        Mode number
    modes : list[str]
        Name of mode with majority index
    """
    weighted = np.stack(
        [oneph_weighted(mode) for mode in sorted(IN_PLANE_MODES)],
        axis=2,
    )
    weighted /= np.sum(weighted, axis=2, keepdims=True)
    weighted = np.nan_to_num(weighted)

    above_threshold = np.greater_equal(weighted, threshold)
    mask = np.any(above_threshold, axis=2)
    majority = np.argmax(
        above_threshold, axis=2
    )  # This operation implies that thresdhold > 50%

    # There are locations where no modes are above threshold.
    # In this case, argmax will pick the first argument (axis 0), which will be LA
    # Therefore, all mode indices are incremented (LA -> 1, TA -> 2, ...)
    # and a value of 0 means no majority. Since `mask` is False (0) where there is no majority,
    # a simple multiplication suffices
    majority += 1
    majority *= mask

    # Modes with majority might have indices 0, 3, 4, 7, 11
    # We want to map this to 0, 1, 2, 3, 4, ...
    values = sorted(np.unique(majority.ravel()))
    replacement = dict(enumerate(values))

    compact = np.zeros_like(majority, dtype=np.int)
    for i, v in replacement.items():
        compact[np.equal(majority, v)] = i

    in_plane_modes_with_nothing = [r"$\varnothing$"] + list(sorted(IN_PLANE_MODES))
    modes = [in_plane_modes_with_nothing[v] for v in values]

    return compact, modes


def draw_bragg_peaks(ax, reflections):
    """ Draw bragg peak locations with black scatter points. """
    cryst = Crystal.from_pwscf(INPUT / "output.out")
    astar, bstar, cstar = cryst.reciprocal_vectors
    bragg_peaks = np.vstack(
        [h * astar + k * bstar + l * cstar for (h, k, l) in reflections]
    )
    # ax.scatter(x=bragg_peaks[:, 0], y=bragg_peaks[:, 1], s=1, c="k")

    # Draw subtle hexagon field
    draw_hexagon_field(
        ax=ax,
        radius=1.7015921874416833,
        crystal=cryst,
        reflections=reflections,
        color="k",
        linestyle=":",
    )


reflections = list(
    filter(lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12))
)

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 2))
grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(1, 2),
    cbar_location="top",
)

ax = grid[0]
ax75 = grid[1]

majority, labels = threshold_oneph(threshold=0.5)
above_75, _ = threshold_oneph(threshold=0.75)

# print(
#     "Area with no majority: ",
#     100 * np.sum(np.less(majority, 1)) / np.size(majority),
# )
# print(
#     "Area with no 75% majority: ",
#     100 * np.sum(np.less(above_75, 1)) / np.size(above_75),
# )

bounds = np.arange(0, stop=len(labels) + 1)

# Base color for no majority should be white
# Also, cmaplist[1] and cmaplist[2] are very similar
# so we swap the unused cmaplist[0]
cmaplist = ["w", "red", "blue", "goldenrod"]
assert len(cmaplist) >= len(labels)

cmap = cm.ListedColormap(name="Modes", colors=cmaplist)
norm = cm.BoundaryNorm(bounds, cmap.N)

qx, qy = np.load(INPUT / "oneph" / "qx.npy"), np.load(INPUT / "oneph" / "qy.npy")
m = ax.imshow(
    majority, cmap=cmap, norm=norm, extent=[qx.min(), qx.max(), qy.min(), qy.max()]
)
ax75.imshow(
    above_75, cmap=cmap, norm=norm, extent=[qx.min(), qx.max(), qy.min(), qy.max()]
)

for ax_, threshold, letter in zip(grid, (50, 75), "ab"):
    ax_.xaxis.set_visible(False)
    ax_.yaxis.set_visible(False)
    tag_axis(ax_, text=f"{letter}) $>${str(threshold)}%")

    # Bragg peaks might extend beyond the rest of the image
    ax_.set_xlim([qx.min(), qx.max()])
    ax_.set_ylim([qy.min(), qy.max()])

    draw_bragg_peaks(ax_, reflections=reflections)

cbar = ax.cax.colorbar(
    m,
    ticks=FixedLocator(locs=np.arange(0, 10) + 0.5),
    format=FixedFormatter(labels),
)

plt.subplots_adjust(bottom=0.01)
