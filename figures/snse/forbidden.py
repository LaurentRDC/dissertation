import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from plotutils import FIGURE_WIDTH, FONTSIZE, tag_axis
from plotutils.snse_datasets import overnight4
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
    (0, 5, 2),
    (0, 1, 2),
    (0, -1, 6),
    # (0,2,-1)
]


figure, (allowed_ax, forbidden_ax) = plt.subplots(
    1, 2, sharey=True, figsize=(FIGURE_WIDTH, 3.5)
)

for ax in (allowed_ax, forbidden_ax):
    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=1)
    ax.axvline(x=0, linestyle="dashed", color="k", linewidth=1)


for indices, color, marker, ax in zip(
    ALLOWED_PEAKS + FORBIDDEN_PEAKS,
    skued.spectrum_colors(len(ALLOWED_PEAKS) + len(FORBIDDEN_PEAKS)),
    Line2D.filled_markers,
    [allowed_ax] * len(ALLOWED_PEAKS) + [forbidden_ax] * len(FORBIDDEN_PEAKS),
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
        marker=marker,
        color=color,
        markersize=3,
        label=skued.indices_to_text(*indices),
    )

for ax in (allowed_ax, forbidden_ax):
    ax.set_xlim([-1.5, 15])
    ax.set_xlabel("Time-delay [ps]")

tag_axis(allowed_ax, "a)")
tag_axis(forbidden_ax, "b)")

allowed_ax.legend(
    ncol=2,
    loc="lower right",
    framealpha=1,
    edgecolor="none",
)
forbidden_ax.legend(
    ncol=2,
    loc="upper right",
    framealpha=1,
    edgecolor="none",
)

allowed_ax.set_ylabel("$\Delta I/I_0$ [a.u.]")

plt.tight_layout()
