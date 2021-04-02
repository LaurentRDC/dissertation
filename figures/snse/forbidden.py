import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from dissutils import LARGE_FIGURE_WIDTH, tag_axis, discrete_colors
from dissutils.snse import overnight4
import skued
from iris import DiffractionDataset


INNER_RADIUS = 15
OUTER_RADIUS = 30

ALLOWED_PEAKS = [
    (0, -1, 3),
    (0, 1, 3),
    (0, -1, 5),
    (0, -2, 4),
]

FORBIDDEN_PEAKS = [
    (0, 0, 1),
    (0, -2, 5),
    (0, 1, 2),
    (0, -1, 6),
]

# To create a clearer distinction between allowed and forbidden reflections,
# let's remove the middle colors
COLORS = discrete_colors(10)
COLORS.pop(4)
COLORS.pop(4)

figure, axes = plt.subplots(
    4,
    2,
    sharex=True,
    sharey=True,
    figsize=(LARGE_FIGURE_WIDTH, 5.5),
    gridspec_kw=dict(hspace=0.05, wspace=0.025),
)

for indices, color, ax in zip(
    ALLOWED_PEAKS + FORBIDDEN_PEAKS,
    COLORS,
    axes[:, 0].tolist() + axes[:, 1].tolist(),
):
    h, k, l = indices
    yi, xi = overnight4.miller_to_arrindex(*indices)

    with DiffractionDataset(overnight4.path, mode="r") as dset:
        selection = skued.RingSelection(
            shape=dset.resolution,
            center=(xi, yi),
            inner_radius=INNER_RADIUS,
            outer_radius=OUTER_RADIUS,
        )
        timedelays = dset.time_points
        timeseries = dset.time_series_selection(selection)

    # Normalize time-series to the mean of pre-time-zero
    pret0 = np.mean(timeseries[timedelays < 0])
    timeseries /= pret0
    error = scipy.stats.sem(timeseries)
    ax.errorbar(
        x=timedelays,
        y=timeseries,
        yerr=error,
        linestyle="None",
        color=color,
        marker="o",
        markersize=2,
        label=skued.indices_to_text(*indices),
    )

    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=1)
    ax.axvline(x=0, linestyle="dashed", color="k", linewidth=1)
    tag_axis(
        ax,
        skued.indices_to_text(*indices),
        x=0.95,
        y=0.9,
        horizontalalignment="right",
        edgecolor="w",
    )

for ax in axes[-1, :]:
    ax.set_xlim([-1.5, 15])
    ax.set_xlabel("Time-delay [ps]")

for ax in axes[:, 0]:
    ax.set_ylabel("$\Delta I/I_0$ [a.u.]")

axes[0, 0].set_title("Allowed reflections")
axes[0, 1].set_title("Forbidden reflections")

plt.subplots_adjust(
    top=0.96, bottom=0.09, left=0.11, right=0.95, hspace=0.2, wspace=0.2
)
