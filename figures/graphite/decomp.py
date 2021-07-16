import itertools as it
from math import sqrt
from pathlib import Path

import h5py
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interpolate
from crystals import Crystal
from crystals.affine import change_of_basis
from matplotlib.ticker import FixedFormatter, FixedLocator
from skued import nfold

from dissutils import LARGE_FIGURE_WIDTH, ImageGrid, draw_hexagon

GAMMA_RADIUS = 0.45  # inverse angstroms
INPUT = Path("data") / "graphite" / "populations"
TIMES = [0.5, 1.5, 5, 100]  # ps
MODES = ["TO2", "TA", "LA"]
GRAPHITE = Crystal.from_pwscf(Path("data") / "graphite" / "graphite.out")
COLORMAP = "hot"


def draw_high_sym_points(ax):
    """Draw the graphite in-plane high-symmetry point Gamma, M, and K
    on an axis."""
    from_frac = change_of_basis(np.array(GRAPHITE.reciprocal_vectors), np.eye(3))

    pts = [(0, 0, 0), (0, 1 / 2, 0), (1 / 3, 1 / 3, 0)]  # gamma  # M  # K
    pts = [from_frac @ pt for pt in pts]

    labels = ["$\Gamma$", "M", "K"]
    offsets = [(-0.1, 0), (-0.1, -0.1), (-0.1, -0.1)]
    valignments = ["center", "top", "top"]
    colors = ["k", "w", "w"]

    for index, vertex in enumerate(pts):
        ax.scatter(vertex[0], vertex[1], s=5, c=colors[index], zorder=np.inf)
        ax.annotate(
            labels[index],
            (vertex[0], vertex[1]),
            xytext=vertex[0:2] + np.array(offsets[index]),
            color=colors[index],
            verticalalignment=valignments[index],
            horizontalalignment="right",
        )


def draw_a1prime_integration(ax, with_label=True):
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

    if with_label:
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
    kx_ = np.array(f["kx"])
    ky_ = np.array(f["ky"])

    kmax = np.linalg.norm(kx_, axis=1).max() + 0.25
    kx, ky = np.meshgrid(np.linspace(-kmax, kmax, 256), np.linspace(-kmax, kmax, 256))
    kk = np.sqrt(kx ** 2 + ky ** 2)

    for mode_name in MODES:
        image = interpolate.griddata(
            points=np.hstack((kx_, ky_)),
            values=np.array(f[mode_name]),
            xi=(kx, ky),
            method="linear",  # fill_value has no effect if method = 'nearest'
            fill_value=0.0,
        )
        image = nfold(image, 6)
        image[kk < 0.45] = 0
        populations[mode_name] = image


def time_index(time):
    return np.argmin(np.abs(time - times))


images = dict()
for time in TIMES:

    images[time] = dict()
    for mode_name in MODES:
        images[time][mode_name] = populations[mode_name][:, :, time_index(time)]

fig = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 0.8 * LARGE_FIGURE_WIDTH))
grid = ImageGrid(
    fig, rect=111, nrows_ncols=(len(MODES), len(TIMES)), cbar_location="top"
)


for ax, (mode_name, time) in zip(grid, it.product(MODES, TIMES)):
    image = images[time][mode_name] / max(images[t][mode_name].max() for t in TIMES)

    m = ax.imshow(
        image,
        extent=[kx.min(), kx.max(), ky.min(), ky.max()],
        cmap=COLORMAP,
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
        linestyle=(
            0,
            (3, 3),
        ),  # loosely dashed: https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html
    )
    ax.add_patch(
        mpatches.Circle(
            xy=(0, 0),
            radius=GAMMA_RADIUS,
            facecolor="w",
            fill=True,
            edgecolor="k",
            linewidth=0.5,
        )
    )

    if time == 0.5 and mode_name == "TA":
        draw_high_sym_points(ax)

    if mode_name == "TO2":
        draw_a1prime_integration(ax, with_label=(time == 0.5))


for index, tag in enumerate(MODES):
    ax = grid[len(TIMES) * index]
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
cbar.ax.set_xlabel(
    r"Change in population $\Delta n_{\lambda}(\mathbf{k}, \tau)$ [a.u.]"
)

plt.subplots_adjust(bottom=0.01)
