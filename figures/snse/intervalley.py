import matplotlib.pyplot as plt
from pathlib import Path
import itertools as it
from math import sqrt, radians
from functools import reduce
from collections import namedtuple
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import numpy as np
import matplotlib.colors as cl
from dissutils import LARGE_FIGURE_WIDTH, tag_axis
import skued
from scipy.signal import correlate2d

DATADIR = Path("data") / "snse" / "estructure"


def signum(x):
    return (x > 0) - (x < 0)


def elec_cmap():
    cmap = plt.get_cmap("YlOrBr")
    colors = [cmap(i) for i in np.linspace(0, 1, num=256)]
    # The point of the plots in this script are to show regions where the images are 0
    # which is easiest to see if these regions are white
    colors[0] = (1, 1, 1, 1)
    return cl.ListedColormap(colors)


def hole_cmap():
    cmap = plt.get_cmap("Purples")
    colors = [cmap(i) for i in np.linspace(0, 1, num=256)]
    # The point of the plots in this script are to show regions where the images are 0
    # which is easiest to see if these regions are white
    colors[0] = (1, 1, 1, 1)
    return cl.ListedColormap(colors)


def zncc(arr, brr):
    """Zero-normalized cross-correlation"""
    return (
        2
        / (arr.size + brr.size)
        * (1 / (arr.std() * brr.std()))
        * correlate2d(
            arr - arr.mean(), brr - brr.mean(), mode="same", boundary="symmetric"
        )
    )


def pathways(ky, kz, energy, nlevels=32):
    """
    Decay pathways allowed by scattering DOWN in energy.
    """

    energy = np.array(energy, copy=True)
    energy -= energy.min()
    phonon_pathways = np.zeros_like(energy)

    dE = (energy.max() - energy.min()) / nlevels
    energy_levels = [
        (energy.min() + i * dE, energy.min() + (i + 1) * dE) for i in range(nlevels)
    ]
    for low, high in energy_levels:
        source = np.zeros_like(energy)
        drain = np.zeros_like(energy)
        # Only consider the dispersion between low and high energy
        source[
            np.logical_and(np.less_equal(energy, high), np.greater_equal(energy, low))
        ] = 1

        drain[np.less_equal(energy, low)] = 1
        phonon_pathways += zncc(source, drain)

    return np.abs(phonon_pathways)


ky = np.load(DATADIR / "ky.npy")
kz = np.load(DATADIR / "kz.npy")
elec_energy = np.load(DATADIR / "inplane-conduction.npy")
hole_energy = np.load(DATADIR / "inplane-valence.npy")

pathways_e = pathways(ky, kz, elec_energy)
pathways_h = pathways(ky, kz, -hole_energy)

# autocorrelation is always very bright in the center
center = (ky / 0.1) ** 2 + (kz / 0.1) ** 2 < 1
pathways_e[center] = 0
pathways_h[center] = 0

pathways_e -= pathways_e.min()
pathways_e /= pathways_e.max()

pathways_h -= pathways_h.min()
pathways_h /= pathways_h.max()


figure, (ax_e, ax_h) = plt.subplots(
    1,
    2,
    sharex=True,
    sharey=True,
    figsize=(LARGE_FIGURE_WIDTH, 0.6 * LARGE_FIGURE_WIDTH),
)

elec_cax = ax_e.inset_axes(bounds=[0.05, 1.04, 0.9, 0.05])
hole_cax = ax_h.inset_axes(bounds=[0.05, 1.04, 0.9, 0.05])

for ax in [ax_e, ax_h]:
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect(1)

me = ax_e.imshow(
    pathways_e,
    vmin=0,
    vmax=1,
    extent=[ky.min(), ky.max(), kz.min(), kz.max()],
    cmap=elec_cmap(),
)
mh = ax_h.imshow(
    pathways_h,
    vmin=0,
    vmax=1,
    extent=[ky.min(), ky.max(), kz.min(), kz.max()],
    cmap=hole_cmap(),
)


# Annotate high-symmetry points -----------------------------------------------
for i, ax in enumerate([ax_e, ax_h]):

    # Hide ugly pixelated center
    ax.add_patch(Circle(xy=(0, 0), radius=0.12, fc="w", ec="none", zorder=np.inf))

    ax.add_patch(
        Rectangle(xy=(-1, -1), width=2, height=2, fc="none", ec="k", clip_on=False)
    )

    for (x, y), label in zip(
        [(1, 0), (0, -1), (1, -1), (-1, 0), (0, 1), (-1, -1), (-1, 1), (1, 1)],
        ["$Y$", "$Z$", "$T$", "", "", "", "", ""],
    ):
        ax.add_patch(
            Circle(
                xy=(x, y),
                radius=0.03,
                edgecolor="None",
                facecolor="k",
                zorder=np.inf,
                clip_on=False,
            )
        )

        # only annotate the first subplot
        if i == 0:
            offset = np.array([signum(x) * 0.06, signum(y) * 0.06])
            ax.annotate(
                xy=(x, y),
                ha={-1: "right", 0: "center", 1: "left"}[signum(x)],
                va={-1: "top", 0: "center", 1: "bottom"}[signum(y)],
                text=label,
                xytext=np.array((x, y)) + offset,
            )

# Colorbars -------------------------------------------------------------------

for cax, mappable, label in zip(
    [elec_cax, hole_cax], [me, mh], ["$P_{e^{-}}$", "$P_{h^{+}}$"]
):
    plt.colorbar(mappable=mappable, cax=cax, orientation="horizontal")
    cax.xaxis.set_major_locator(ticker.FixedLocator([0, 1]))
    cax.set_xlabel(label)
    cax.xaxis.tick_top()
    cax.xaxis.set_label_position("top")

tag_axis(ax_e, "a)")
tag_axis(ax_h, "b)")
plt.subplots_adjust(
    top=0.875, bottom=0.037, left=0.029, right=0.972, hspace=0.2, wspace=0.128
)
