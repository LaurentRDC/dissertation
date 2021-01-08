from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import numpy as np
from plotutils import FIGURE_WIDTH, draw_hexagon
from scipy.constants import physical_constants
from plotutils import CBAR_SIZE

INPUT = Path("data") / "graphite" / "weighted-dispersion"

HBAR = physical_constants["Planck constant over 2 pi in eV s"][0]
HZ_TO_EV = HBAR * 2 * np.pi

ACOUSTIC_CMAP = "Reds"
OPTICAL_CMAP = "Blues"

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, 0.35 * FIGURE_WIDTH))

LABELS = [
    r"$\mathbf{\Gamma}_{(010)}$",
    r"$\mathbf{M}$",
    r"$\mathbf{K}$",
    r"$\mathbf{\Gamma}_{(010)} ~ \leftrightarrow ~ \mathbf{\Gamma}_{(\bar{1}10)}$",
    r"$\mathbf{K}$",
    r"$\mathbf{M}$",
    r"$\mathbf{\Gamma}_{(\bar{1}10)}$",
]
NSTEPS = 1000

for mode in ["LA", "TA"]:
    frequencies = np.load(INPUT / f"{mode}_frequencies.npy")
    weights = np.load(INPUT / f"{mode}_oneph.npy")

    am = ax.scatter(
        x=range(np.size(frequencies)),
        y=frequencies * HZ_TO_EV * 1000,
        c=weights / weights.max(),
        s=3,
        cmap=ACOUSTIC_CMAP,
    )


for mode in ["LO2", "TO2"]:
    frequencies = np.load(INPUT / f"{mode}_frequencies.npy")
    weights = np.load(INPUT / f"{mode}_oneph.npy")

    om = ax.scatter(
        x=range(np.size(frequencies)),
        y=frequencies * HZ_TO_EV * 1000,
        c=weights / weights.max(),
        s=3,
        cmap=OPTICAL_CMAP,
    )

divider = make_axes_locatable(ax)
acoustic_cax = divider.append_axes("right", size=CBAR_SIZE, pad=0.03)
optical_cax = divider.append_axes("right", size=CBAR_SIZE, pad=0.03)

for cax in (acoustic_cax, optical_cax):
    cax.yaxis.set_visible(False)
    cax.xaxis.set_visible(False)

plt.colorbar(mappable=am, cax=acoustic_cax, orientation="vertical")
plt.colorbar(mappable=om, cax=optical_cax, orientation="vertical")
optical_cax.yaxis.tick_right()
optical_cax.yaxis.set_label_position("right")
optical_cax.yaxis.set_visible(True)
optical_cax.set_ylabel(r"$|F_{1j}(\mathbf{q}, \tau<0)|^2$ [a.u.]")

for i in range(1, len(LABELS)):
    ax.axvline(x=i * NSTEPS, color="k", linestyle="--", linewidth=1)

# Middle line traces out the boundary between Brillouin zones
ax.axvline(x=3 * NSTEPS, color="k", linewidth=3)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(
    ticker.FixedLocator([i * NSTEPS for i in range(len(LABELS))])
)
ax.xaxis.set_major_formatter(plt.FixedFormatter(LABELS))
ax.set_xlim([0, (len(LABELS) - 1) * NSTEPS])
ax.set_ylabel("Energy [meV]")
ax.set_ylim([0, 210])
ax.set_facecolor("lightgray")


# Custom legend
# Help here:
#   https://matplotlib.org/3.1.0/gallery/text_labels_and_annotations/custom_legends.html
legend_kwds = dict(marker="s", color="none", markeredgecolor="None")
legend_elements = [
    Line2D(
        [0],
        [0],
        markerfacecolor=plt.get_cmap(OPTICAL_CMAP)(256),
        label="Optical modes",
        **legend_kwds,
    ),
    Line2D(
        [0],
        [0],
        markerfacecolor=plt.get_cmap(ACOUSTIC_CMAP)(256),
        label="Acoustic modes",
        **legend_kwds,
    ),
]
ax.legend(
    handles=legend_elements,
    ncol=2,
    bbox_to_anchor=(0.5, 1.07),
    loc="center",
    fontsize=8,
)

# Annotate the modes
for label, pt in zip(
    ["TA", "LA", "TO2", "LO2"],
    [(0.08, 0.15), (0.025, 0.45), (0.42, 0.77), (0.35, 0.94)],
):
    ax.annotate(label, pt, color="k", textcoords=ax.transAxes, fontsize=8)
    ax.annotate(
        label,
        np.asarray(pt) + np.array([0.5, 0]),
        color="k",
        textcoords=ax.transAxes,
        fontsize=8,
    )


# Draw path in inset
path_ax = fig.add_axes([0.25, 0.11, 0.12, 0.33])
path_ax.axis("off")

# Draw path and hexagon
vertices1 = np.load(INPUT / "vertices1.npy")
vertices2 = np.load(INPUT / "vertices2.npy")

# Note that the hexagon "radius" is dependent on the fact that path goes from gamma (center) to K (furthest
# point from gamma)
# If the path ever changes, the definition of this radius will have to change
hex_radius = 1.7015921874416833
draw_hexagon(path_ax, center=vertices1[0][0:2], radius=hex_radius, color="k")
draw_hexagon(path_ax, center=vertices2[0][0:2], radius=hex_radius, color="k")
draw_hexagon(path_ax, center=(0, 0), radius=hex_radius, color="k")

# Last vertex is same as first
# Note that for display reasons, we remove the e2 offset
# from every point
for vertex, label, v_alignment, h_alignment in zip(
    vertices1[0:3].tolist() + vertices2[1::].tolist(),
    LABELS[0:3] + LABELS[4::],
    ["top", "bottom", "bottom", "bottom", "bottom", "top"],
    ["center", "right", "left", "right", "right", "center"],
):
    path_ax.scatter(
        x=vertex[0],
        y=vertex[1],
        s=2,
        c="k",
    )
    path_ax.annotate(
        label,
        (vertex[0], vertex[1]),
        color="k",
        verticalalignment=v_alignment,
        horizontalalignment=h_alignment,
        fontsize=8,
    )

# Show center of (000) reflections
path_ax.scatter(0, 0, s=2, c="k")
path_ax.annotate(
    "$\mathbf{\Gamma}_{(000)}$",
    (0, 0),
    color="k",
    verticalalignment="top",
    horizontalalignment="center",
    fontsize=8,
)
