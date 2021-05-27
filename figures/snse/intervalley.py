import matplotlib.pyplot as plt
import itertools as it
from math import sqrt, radians
from functools import reduce
from collections import namedtuple
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from skimage.filters import gaussian
import numpy as np
import matplotlib.colors as cl
from dissutils import LARGE_FIGURE_WIDTH, tag_axis
import skued
from scipy.signal import correlate2d


def signum(x):
    return (x > 0) - (x < 0)


def elec_cmap():
    cmap = plt.get_cmap("YlOrBr")
    colors = [cmap(i) for i in np.linspace(0, 0.9, num=256)]
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


def electron_dispersion(ky, kz):
    elec_energy = np.zeros_like(ky)

    for (cy, cz), height, width in zip(
        [(2 / 3, 0), (-2 / 3, 0), (0, 0), (0, 2 / 3), (0, -2 / 3)],
        [1.0, 1.0, 0.8, 0.7, 0.7],
        [0.4, 0.4, 0.6, 0.4, 0.4],
    ):
        gauss = skued.gaussian(coordinates=[ky, kz], center=[cy, cz], fwhm=width)
        gauss /= gauss.max()
        elec_energy += height * gauss

    elec_energy /= elec_energy.max()
    elec_energy *= -1
    elec_energy -= elec_energy.min()
    return elec_energy


def hole_dispersion(ky, kz):
    hole_energy = np.zeros_like(ky)

    for (cy, cz), height, width in zip(
        [
            (-2 / 3, 0),
            (2 / 3, 0),
            (0, 0),
            (0, -3 / 4),
            (0, -3 / 4),
            (0, 3 / 4),
            (0, 3 / 4),
        ],
        # Pudding mold band represented by the difference of a tall wide gaussian
        # and a narrow, small gaussian
        [0.7, 0.7, 0.6, 1.0, -0.4, 1.0, -0.4],
        [0.4, 0.4, 0.7, 0.4, 0.08, 0.4, 0.08],
    ):
        gauss = skued.gaussian(coordinates=[ky, kz], center=[cy, cz], fwhm=width)
        gauss /= np.abs(gauss).max()
        hole_energy += height * gauss
    hole_energy /= hole_energy.max()

    hole_energy -= hole_energy.max()
    return hole_energy


def pathways(ky, kz, energy, nlevels=32):
    """
    Decay pathways allowed by scattering DOWN in energy.
    """
    energy = np.array(energy, copy=True)
    phonon_pathways = np.zeros_like(energy)

    dE = (energy.max() - energy.min()) / nlevels
    energy_levels = [
        (energy.min() + i * dE, energy.min() + (i + 1) * dE) for i in range(nlevels)
    ]

    for low, high in energy_levels:
        source = np.array(energy, copy=True)
        drain = np.array(energy, copy=True)
        source_mask = np.logical_or(
            np.less_equal(energy, low), np.greater_equal(energy, high)
        )
        source[source_mask] = 0

        drain[np.greater_equal(energy, low)] = 0
        phonon_pathways += correlate2d(source, drain, mode="same", boundary="symmetric")

    return np.abs(phonon_pathways)


ky, kz = np.meshgrid(*[np.linspace(-1, 1, num=92)] * 2)

# Only consider scattering from electrons/holes close to the local extrema
elec_energy = electron_dispersion(ky, kz)
elec_energy[elec_energy > elec_energy.min() + 0.4] = 0
hole_energy = hole_dispersion(ky, kz)
hole_energy[hole_energy < hole_energy.max() - 0.4] = 0

pathways_e = pathways(ky, kz, elec_energy)
pathways_e -= pathways_e.min()
pathways_e /= pathways_e.max()

pathways_h = pathways(ky, kz, -hole_energy)
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
