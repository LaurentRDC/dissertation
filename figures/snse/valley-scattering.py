import matplotlib.pyplot as plt
import itertools as it
from math import sqrt
from functools import reduce
from collections import namedtuple
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from skimage.filters import gaussian
import numpy as np
import matplotlib.colors as cl
from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis

ARROWSTYLE = dict(arrowstyle="-|>", zorder=np.inf, mutation_scale=5, fc="k", ec="k")

_, PURPLE, _, _, _, _, YELLOW = discrete_colors(7)
ELEC_CMAP = cl.LinearSegmentedColormap.from_list("", [YELLOW, "w"])
HOLE_CMAP = cl.LinearSegmentedColormap.from_list("", ["w", PURPLE])

R60deg = np.array([[1 / 2, -sqrt(3) / 2], [sqrt(3) / 2, 1 / 2]])

ky, kz = np.meshgrid(*[np.linspace(-1, 1, num=512)] * 2)
hole_region = ky < -kz
elec_region = np.logical_not(hole_region)

elec_energy = np.zeros_like(ky)
hole_energy = np.zeros_like(ky)


for (cy, cz), height in zip(
    [(2 / 3, 0), (0, 0), (0, 2 / 3)],
    [1, 0.8, 0.7],
):
    elec_energy[
        np.argmin(np.abs(kz[:, 0] - cz)), np.argmin(np.abs(ky[0, :] - cy))
    ] = height

elec_energy[:] = gaussian(elec_energy, sigma=40)
elec_energy *= -1
elec_energy -= elec_energy.min()

elec_energy /= elec_energy.max()
elec_energy[hole_region] = elec_energy.max()

for (cy, cz), height in zip(
    [(-2 / 3, 0), (0, 0)],
    [0.7, 0.6],
):
    hole_energy[
        np.argmin(np.abs(kz[:, 0] - cz)), np.argmin(np.abs(ky[0, :] - cy))
    ] = height
hole_energy[:] = gaussian(hole_energy, sigma=40)
hole_energy /= hole_energy.max()

# Pudding mold band near 3/4Z
near_Z = np.sqrt(ky ** 2 + (kz + 3 / 4) ** 2)
pudding_mold_energy = np.zeros_like(hole_energy)
pudding_mold_energy[near_Z < 0.02] = -1.5
pudding_mold_energy[np.isclose(near_Z, 0.13, rtol=5e-2)] = 1
# pudding_mold_energy[near_Z < 0.19] = 0.15
pudding_mold_energy[:] = gaussian(pudding_mold_energy, sigma=20)
pudding_mold_energy /= pudding_mold_energy.max()

hole_energy += pudding_mold_energy
hole_energy -= hole_energy.max()
hole_energy[elec_region] = hole_energy.min()

figure, (ax_elec1, ax_elec2, ax_elec3) = plt.subplots(
    1, 3, figsize=(LARGE_FIGURE_WIDTH, 0.4 * LARGE_FIGURE_WIDTH)
)

cax = ax_elec1.inset_axes(bounds=[-0.03, 0.1, 0.05, 0.8])

for i, ax in enumerate([ax_elec1, ax_elec2, ax_elec3]):

    ax.axis("off")
    ax.set_aspect(1)

    elec_colors = [ELEC_CMAP(bc) for bc in np.linspace(0, 1, num=7)]
    hole_colors = [HOLE_CMAP(bc) for bc in np.linspace(0, 1, num=len(elec_colors))]
    hole_colors[0] = elec_colors[-1] = (0, 0, 0, 0)

    me = ax.contourf(ky, kz, elec_energy, colors=elec_colors, levels=len(elec_colors))
    mh = ax.contourf(ky, kz, hole_energy, colors=hole_colors, levels=len(hole_colors))

    # Split hole/electron regions
    y = np.linspace(-1, 1, num=1024)
    ax.plot(y, -y, "k", clip_on=True, linewidth=1)

    ax.add_patch(
        Rectangle(
            xy=(-1, -1),
            width=2,
            height=2,
            fc="none",
            ec="k",
            clip_on=False,
            zorder=np.inf,
        )
    )

    for (x, y), label in zip(
        [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1), (-1, 1), (1, -1)],
        ["$\Gamma$", "$Y$", "$Z$", "$T$"] + [""] * 5,
    ):
        ax.add_patch(
            Circle(
                xy=(x, y), radius=0.03, edgecolor="None", facecolor="k", zorder=np.inf
            )
        )

        # only annotate the first subplot
        if i == 0:
            ax.annotate(
                xy=(x, y),
                ha="left",
                va="bottom",
                text=label,
                xytext=np.array((x, y)) + np.array([0.06, 0.06]),
            )

    ax.text(x=0.9, y=0.9, s="e$^{-}$", va="top", ha="right")
    ax.text(x=-0.9, y=-0.9, s="h$^{+}$", va="bottom", ha="left")

cmap_colors = hole_colors + [(0.7, 0.7, 0.7, 0.7)] + elec_colors
mappable = cm.ScalarMappable(
    norm=cl.Normalize(vmin=-1, vmax=1),
    cmap=cl.ListedColormap(name="", colors=cmap_colors),
)
plt.colorbar(mappable=mappable, cax=cax)
cax.yaxis.set_major_locator(ticker.FixedLocator([0]))
cax.yaxis.set_major_formatter(ticker.FixedFormatter(["$E_f$"]))
cax.yaxis.tick_left()
cax.yaxis.set_label_position("left")


# intravalley scattering
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

# Intervalley scattering
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

# Electrons: positive quadrant
for (start, _), (end, _) in elec_pathways + hole_pathways:
    ax_elec3.add_patch(
        FancyArrowPatch(posA=start, posB=end, shrinkA=5, shrinkB=5, **ARROWSTYLE)
    )

tag_kwds = dict(x=0.12, y=0.88)
tag_axis(ax_elec1, "a)", **tag_kwds)
tag_axis(ax_elec2, "b)", **tag_kwds)
tag_axis(ax_elec3, "c)", **tag_kwds)

plt.subplots_adjust(
    top=1.0, bottom=0.006, left=0.056, right=1.0, hspace=0.135, wspace=0.052
)
