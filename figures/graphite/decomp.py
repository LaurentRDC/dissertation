import itertools as it
from pathlib import Path
import matplotlib.path as mpath

import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import sqrt
import numpy as np
import scipy.interpolate as interpolate
from matplotlib.ticker import FixedFormatter, FixedLocator
from plotutils import FIGURE_WIDTH, ImageGrid, draw_hexagon, tag_axis
from skimage.filters import gaussian
from skued import nfold
from crystals import Crystal
from crystals.affine import change_of_basis

GAMMA_RADIUS = 0.45  # inverse angstroms
INPUT = Path("data") / "graphite" / "populations"
TIMES = [0.5, 1.5, 5, 100]  # ps
MODES = ["TO2", "TA", "LA"]
GRAPHITE = Crystal.from_pwscf(Path("data") / "graphite" / "output.out")


def draw_high_sym_points(ax):
    """Draw the graphite in-plane high-symmetry point Gamma, M, and K
    on an axis."""
    from_frac = change_of_basis(np.array(GRAPHITE.reciprocal_vectors), np.eye(3))

    pts = [(0, 0, 0), (0, 1 / 2, 0), (1 / 3, 1 / 3, 0)]  # gamma  # M  # K
    pts = [from_frac @ pt for pt in pts]

    labels = ["$\Gamma$", "M", "K"]
    offsets = [(-0.1, 0), (-0.1, -0.1), (-0.1, -0.1)]
    valignments = ["center", "top", "top"]

    for index, vertex in enumerate(pts):
        ax.scatter(vertex[0], vertex[1], s=5, c="w")
        ax.annotate(
            labels[index],
            (vertex[0], vertex[1]),
            xytext=vertex[0:2] + np.array(offsets[index]),
            color="w",
            verticalalignment=valignments[index],
            horizontalalignment="right",
        )


def draw_a1prime_integration(ax):
    """Draw the zone where the population is
    integrated to determine hte a1prime population"""
    from_frac = change_of_basis(np.array(GRAPHITE.reciprocal_vectors), np.eye(3))
    pt = from_frac @ (1 / 3, 1 / 3, 0)  # K-point

    circ = mpatches.Wedge(
        center=pt[0:2],
        r=0.3,
        theta1=180,
        theta2=300,
        facecolor="None",
        fill=False,
        edgecolor="w",
        linewidth=0.5,
    )
    ax.add_patch(circ)

    ax.annotate(
        "$A^{\prime}_1$",
        (pt[0], pt[1]),
        xytext=pt[0:2] + np.array([0.1, 0]),
        color="w",
        verticalalignment="center",
        horizontalalignment="left",
    )


populations = dict()
with h5py.File(INPUT / "population_timeseries.hdf5", mode="r") as f:
    times = f.attrs["times"]
    kx = np.array(f["kx"])
    ky = np.array(f["ky"])

    for mode_name in MODES:
        populations[mode_name] = np.array(f[mode_name])


def time_index(time):
    return np.argmin(np.abs(time - times))


images = dict()
for time in TIMES:

    images[time] = dict()
    for mode_name in MODES:
        images[time][mode_name] = populations[mode_name][:, :, time_index(time)]

fig = plt.figure(figsize=(FIGURE_WIDTH, 0.8 * FIGURE_WIDTH))
grid = ImageGrid(
    fig, rect=111, nrows_ncols=(len(MODES), len(TIMES)), cbar_location="top"
)


for ax, (mode_name, time) in zip(grid, it.product(MODES, TIMES)):
    image = images[time][mode_name] / max(images[t][mode_name].max() for t in TIMES)

    m = ax.imshow(
        image,
        extent=[kx.min(), kx.max(), ky.min(), ky.max()],
        cmap="hot",
        vmin=0,
        vmax=1,
    )
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    draw_hexagon(ax, center=(0, 0), radius=1.7, orientation=np.deg2rad(30))
    draw_hexagon(
        ax,
        center=(0, 0),
        radius=1.7 / 2.0,
        orientation=np.deg2rad(30),
        linestyle="dashed",
    )
    ax.add_patch(
        mpatches.Circle(
            xy=(0, 0),
            radius=GAMMA_RADIUS,
            facecolor="none",
            fill=False,
            edgecolor="w",
            linewidth=0.5,
        )
    )

    if time == 0.5 and mode_name == "TA":
        draw_high_sym_points(ax)

    if mode_name == "TO2":
        draw_a1prime_integration(ax)


for index, tag in enumerate(MODES):
    ax = grid[4 * index]
    ax.text(
        x=-0.1,
        y=0.5,
        s=tag,
        color="k",
        horizontalalignment="right",
        verticalalignment="center",
        transform=ax.transAxes,
    )

for index, time in enumerate(TIMES):
    ax = grid[8 + index]
    ax.text(
        x=0.5,
        y=-0.1,
        s=f"{1e3 * time:0.0f} fs" if time < 1 else f"{time} ps",
        color="k",
        horizontalalignment="center",
        verticalalignment="top",
        transform=ax.transAxes,
    )

cbar = grid[0].cax.colorbar(
    m, ticks=FixedLocator(locs=[0, 1]), format=FixedFormatter(["0", "1"])
)
cbar.ax.set_xlabel(r"Change in population $\Delta n_{j}(\mathbf{k}, \tau)$ [a.u.]")

plt.subplots_adjust(bottom=0.01)
