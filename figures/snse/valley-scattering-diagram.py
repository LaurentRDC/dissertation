import itertools as it
from collections import namedtuple
from functools import reduce
from math import radians, sqrt

import matplotlib.cm as cm
import matplotlib.colors as cl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import skued
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from skimage.filters import gaussian

from dissutils import LARGE_FIGURE_WIDTH, tag_axis

ARROWSTYLE = dict(arrowstyle="-|>", zorder=np.inf, mutation_scale=5, fc="k", ec="k")


def signum(x):
    return (x > 0) - (x < 0)


def rotmat(angle):
    r = radians(angle)
    return np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])


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


def elec_cmap():
    cmap = plt.get_cmap("YlOrBr")
    colors = [cmap(i) for i in np.linspace(0, 0.9, num=256)]
    colors[0] = (0, 0, 0, 0)
    return cl.ListedColormap(colors)


def hole_cmap():
    cmap = plt.get_cmap("Purples")
    colors = [cmap(i) for i in np.linspace(0, 1, num=256)]
    colors[0] = (0, 0, 0, 0)
    return cl.ListedColormap(colors)


ky, kz = np.meshgrid(*[np.linspace(-1, 1, num=512)] * 2)
hole_region = ky < -kz
elec_region = np.logical_not(hole_region)

elec_energy = electron_dispersion(ky, kz)
hole_energy = hole_dispersion(ky, kz)


elec_pop1 = np.zeros_like(elec_energy)
elec_pop1[elec_energy < 0.6 * elec_energy.max()] = 1
elec_pop1[:] = gaussian(elec_pop1, 20)

elec_pop2 = np.zeros_like(elec_energy)
elec_pop2[elec_energy < 0.4 * elec_energy.max()] = 1
elec_pop2[:] = gaussian(elec_pop2, 10)

elec_pop3 = np.zeros_like(elec_energy)
elec_pop3[elec_energy < 0.03 * elec_energy.max()] = 1
elec_pop3[:] = gaussian(elec_pop2, 10)

hole_pop1 = np.zeros_like(hole_energy)
hole_pop1[hole_energy > -0.6] = 1
hole_pop1[:] = gaussian(hole_pop1, 20)

hole_pop2 = np.zeros_like(hole_energy)
hole_pop2[hole_energy > -0.5] = 1
hole_pop2[:] = gaussian(hole_pop2, 10)

hole_pop3 = np.zeros_like(hole_energy)
hole_pop3[hole_energy > -0.35] = 1
hole_pop3[:] = gaussian(hole_pop3, 10)

elec_energy[hole_region] = elec_energy.max()
hole_energy[elec_region] = hole_energy.min()

figure, (ax_elec1, ax_elec2, ax_elec3) = plt.subplots(
    1, 3, figsize=(LARGE_FIGURE_WIDTH, 0.4 * LARGE_FIGURE_WIDTH)
)

elec_cax = ax_elec2.inset_axes(bounds=[0.05, 1.04, 0.9, 0.05])
hole_cax = ax_elec2.inset_axes(bounds=[0.05, -0.09, 0.9, 0.05])

for i, (ax, epop, hpop) in enumerate(
    zip(
        [ax_elec1, ax_elec2, ax_elec3],
        [elec_pop1, elec_pop2, elec_pop3],
        [hole_pop1, hole_pop2, hole_pop3],
    )
):

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect(1)

    # By adding a small constant, the "transparent" color
    # only appears in the opposite region
    epop += 0.001
    hpop += 0.001
    epop[hole_region] = 0
    hpop[elec_region] = 0

    me = ax.imshow(
        epop[::-1, :],
        vmin=0,
        vmax=1,
        cmap=elec_cmap(),
        extent=[ky.min(), ky.max(), kz.min(), kz.max()],
    )
    mh = ax.imshow(
        hpop[::-1, :],
        vmin=0,
        vmax=1,
        cmap=hole_cmap(),
        extent=[ky.min(), ky.max(), kz.min(), kz.max()],
    )

    ax.contour(
        ky,
        kz,
        elec_energy,
        colors="k",
        linestyles="solid",
        linewidths=0.5,
        levels=7,
        extent=[ky.min(), ky.max(), kz.min(), kz.max()],
    )
    ax.contour(
        ky,
        kz,
        hole_energy,
        colors="k",
        linestyles="solid",
        linewidths=0.5,
        levels=7,
        extent=[ky.min(), ky.max(), kz.min(), kz.max()],
    )

    # Split hole/electron regions
    y = np.linspace(-1, 1, num=1024)
    ax.plot(y, -y, "k", clip_on=True, linewidth=1)

    ax.text(x=0.9, y=0.9, s="e$^{-}$", va="top", ha="right")
    ax.text(x=-0.9, y=-0.9, s="h$^{+}$", va="bottom", ha="left")

# intravalley scattering ------------------------------------------------------
R60deg = rotmat(60)

for xy, v in zip(
    map(np.asarray, [(0, 0), (2 / 3, 0), (-2 / 3, 0), (0, 2 / 3)]),
    [(1 / 2, 1 / 2)] + [(0, 3 / 4)] * 6,
):
    for n in range(0, 6):
        start = 0.3 * reduce(lambda m1, m2: m1 @ m2, n * [R60deg], R60deg) @ np.array(v)
        ax_elec2.add_patch(FancyArrowPatch(posA=start + xy, posB=xy, **ARROWSTYLE))


# Pudding mold hole band is special
for n in range(0, 6):
    xy = (0, -3 / 4)
    start_out = (
        0.3
        * reduce(lambda m1, m2: m1 @ m2, n * [R60deg], R60deg)
        @ np.array([3 / 4, 0])
    )
    start_in = (
        0.2
        * reduce(lambda m1, m2: m1 @ m2, n * [R60deg], R60deg)
        @ np.array([0, 3 / 4])
    )
    ax_elec2.add_patch(
        FancyArrowPatch(
            posA=3 / 2 * start_out + xy, posB=1 / 2 * start_out + xy, **ARROWSTYLE
        )
    )
    ax_elec2.add_patch(
        FancyArrowPatch(posA=xy, posB=start_in + xy, shrinkA=0, **ARROWSTYLE)
    )

# Intervalley scattering ------------------------------------------------------
extrema = namedtuple("Extrema", ["yz", "E"])
elec_endpoints = [
    extrema(yz=(0, 2 / 3), E=2),
    extrema(yz=(0, 0), E=1),
    extrema(yz=(2 / 3, 0), E=0),
]
hole_endpoints = [
    extrema(yz=(0, -3 / 4), E=2),
    extrema(yz=(0, 0), E=0),
    extrema(yz=(-2 / 3, 0), E=1),
]

valid_elec_pathway = lambda t: t[0].E > t[1].E  # Electrons scatter down
valid_hole_pathway = lambda t: t[0].E < t[1].E  # Holes scatter up

# Filter to remove pathways that go up in energy
elec_pathways = list(filter(valid_elec_pathway, it.product(elec_endpoints, repeat=2)))
hole_pathways = list(filter(valid_hole_pathway, it.product(hole_endpoints, repeat=2)))

for (start, _), (end, _) in elec_pathways + hole_pathways:
    start = np.asarray(start)
    end = np.asarray(end)
    midpoint = (start + end) / 2

    for angle, alpha in zip(
        [-15, -7.5, 0, 7.5, 15], [0.3, 0.7, 1, 0.7, 0.3]
    ):  # degrees
        # Rotate arrows about the midpoint
        ax_elec3.add_patch(
            FancyArrowPatch(
                posA=midpoint + rotmat(angle) @ (start - midpoint),
                posB=midpoint + rotmat(angle) @ (end - midpoint),
                shrinkA=5,
                shrinkB=5,
                alpha=alpha,
                **ARROWSTYLE
            )
        )

# Annotate high-symmetry points -----------------------------------------------

for i, ax in enumerate([ax_elec1, ax_elec2, ax_elec3]):

    ax.add_patch(
        Rectangle(xy=(-1, -1), width=2, height=2, fc="none", ec="k", clip_on=False)
    )

    for (x, y), label in zip(
        [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1), (-1, 1), (1, -1)],
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

for cax, mappable, label, pos in zip(
    [elec_cax, hole_cax], [me, mh], ["$N_{e^{-}}$", "$N_{h^{+}}$"], ["top", "bottom"]
):
    plt.colorbar(mappable=mappable, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks([])
    cax.set_xlabel(label)
    cax.xaxis.set_label_position(pos)

ax_elec1.xaxis.set_visible(True)
ax_elec1.yaxis.set_visible(True)
ax_elec1.xaxis.set_major_locator(ticker.FixedLocator([-1, 0, 1]))
ax_elec1.xaxis.set_major_formatter(ticker.FixedFormatter(["-½", "0", "½"]))
ax_elec1.yaxis.set_major_locator(ticker.FixedLocator([-1, 0, 1]))
ax_elec1.yaxis.set_major_formatter(ticker.FixedFormatter(["-½", "0", "½"]))
ax_elec1.set_xlabel("$k_y/b^{*}$")
ax_elec1.set_ylabel("$k_z/c^{*}$")

tag_kwds = dict(x=0.12, y=0.88)
tag_axis(ax_elec1, "a)", **tag_kwds)
tag_axis(ax_elec2, "b)", **tag_kwds)
tag_axis(ax_elec3, "c)", **tag_kwds)

plt.tight_layout()
