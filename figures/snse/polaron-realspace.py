from pathlib import Path
from math import pi
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle, FancyArrowPatch
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib
from crystals import Crystal, Atom, distance_cartesian
from dissutils import tag_axis, named_arrow, LARGE_FIGURE_WIDTH

NXE, NZE = 2, 5  # Supercell factor in a and c direction for electron polaron
NYH = NZH = 3  # Supercell factor in b and c directions for hole polaron
CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")
SNCOLOR, SECOLOR = plt.get_cmap("inferno")([0.2, 0.8])
RADII = {"Sn": 0.4, "Se": 0.3}
ATOM_PATCH_PARAMS = {
    "Sn": dict(fc=SNCOLOR, ec="k", linewidth=0.5),
    "Se": dict(fc=SECOLOR, ec="k", linewidth=0.5),
}
POLARON_SIZE_PARAMS = dict(ec="none", alpha=0.4, zorder=-np.inf)
DISP_ARROWS_PARAMS = dict(arrowstyle="<|-|>", mutation_scale=8, clip_on=False, fc="k")
LINE_PARAMS = dict(color="k", linewidth=1, zorder=-np.inf)
BORN_CHARGES = {
    "Sn": -1,
    "Se": 1,
}  # Sign of the Born effective charges from Caruso PRB 2019

ELECTRON_POLARON_DIMENSIONS = (1, 18)
HOLE_POLARON_DIMENSIONS = (4.5, 4.5)

HCOLOR, ECOLOR = plt.get_cmap("Purples")(0.7), plt.get_cmap("YlOrBr")(0.7)


def signum(x):
    return (x > 0) - (x < 0)


def gaussian1d(x, sx):
    """2D Gaussian centered at (0,0)"""
    return np.exp(-(0.5 * (x / sx) ** 2))


def gaussian2d(x, y, sx, sy):
    """2D Gaussian centered at (0,0)"""
    return np.exp(-(0.5 * (x / sx) ** 2 + 0.5 * (y / sy) ** 2))


def radial(x, y):
    direction = np.array([x, y])
    return direction / np.linalg.norm(direction)


def linear(x, y):
    return np.array([0.1 * signum(float(x)), signum(float(y))])


@dataclass
class AtomicBond:
    atm1: Atom
    atm2: Atom

    def matches(self, a1, a2):
        if (
            (self.atm1 == a1)
            and (self.atm2 == a2)
            or (self.atm1 == a2)
            and (self.atm2 == a1)
        ):
            return True
        return False


def bonds(cell):
    """Yield a list of bonds which connects atoms in `cell`."""
    atomic_bonds = list()

    def distance(atm1, atm2):
        """Only consider the distance between two distinct Sn and Se atoms."""
        # Prevent duplicate bonds
        if any(b.matches(atm1, atm2) for b in atomic_bonds):
            return float("inf")
        # Prevent bonds between same elements
        if atm1.element == atm2.element:
            return float("inf")
        d = distance_cartesian(atm1, atm2)
        if d < 1:
            return float("inf")
        return d

    atoms = list(sorted(cell))
    for i, atm in enumerate(atoms):
        closest_atm = min(atoms, key=lambda a: distance(atm, a))
        if distance(atm, closest_atm) < 3:
            atomic_bonds.append(AtomicBond(atm, closest_atm))

    yield from (
        (b.atm1.coords_cartesian, b.atm2.coords_cartesian) for b in atomic_bonds
    )


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 0.6 * LARGE_FIGURE_WIDTH))
gs = gridspec.GridSpec(2, 2, figure=figure, height_ratios=(1, 0.05))

ax_1d = figure.add_subplot(gs[0, 0])
ax_3d = figure.add_subplot(gs[0, 1])
ax_elements = figure.add_subplot(gs[1, :])

for ax in [ax_1d, ax_3d]:
    ax.axis("off")
    ax.set_aspect(1)

# 1D (electron) polaron -------------------------------------------------------
cell_elec = CRYSTAL.supercell(NXE, 1, NZE)
a, b, c, *_ = cell_elec.lattice_parameters
sx, sz = ELECTRON_POLARON_DIMENSIONS
ax_1d.add_patch(
    Rectangle(
        xy=(-sz / 2, -NXE * a / 2),
        width=sz,
        height=NXE * a,
        fc=ECOLOR,
        **POLARON_SIZE_PARAMS,
    )
)
for atm in cell_elec:
    x, y, z = atm.coords_cartesian
    # Center supercell about (0,0)
    x -= NXE * a / 2
    z -= NZE * c / 2

    born_eff_charge = BORN_CHARGES[atm.element]
    # Electron polaron attracts atoms with negative born effective charge
    uz = -0.5 * born_eff_charge * gaussian1d(z, sz / 2.335) * signum(float(z))

    ax_1d.add_patch(
        Circle(
            xy=(z + uz, x),
            **ATOM_PATCH_PARAMS[atm.element],
            radius=RADII[atm.element],
            clip_on=False,
            zorder=y,
        )
    )
    ax_1d.add_patch(
        Circle(
            xy=(z, x),
            **ATOM_PATCH_PARAMS[atm.element],
            radius=RADII[atm.element],
            clip_on=False,
            zorder=y - 1,
            alpha=0.3,
        )
    )

# Atomic bonds
for (x1, y1, z1), (x2, y2, z2) in bonds(cell_elec):
    ax_1d.add_patch(
        FancyArrowPatch(
            posA=(z1 - NZE * c / 2, x1 - NXE * a / 2),
            posB=(z2 - NZE * c / 2, x2 - NXE * a / 2),
            alpha=0.3,
            linewidth=0.5,
        )
    )

for i in range(-NZE // 2, (NZE + 1) // 2):
    ax_1d.axvline(x=i * c, **LINE_PARAMS)
assert NXE == 2  # The following line only makes sense if NXE == 2
ax_1d.axhline(y=0, **LINE_PARAMS)

ax_1d.set_xlim([-NZE * c / 2, NZE * c / 2])
ax_1d.set_ylim([-NXE * a / 2, NXE * a / 2])

for y in [-a / 2, a / 2]:
    ax_1d.add_patch(FancyArrowPatch(posA=(-3, y), posB=(3, y), **DISP_ARROWS_PARAMS))

# 3D (hole) polaron -----------------------------------------------------------
cell_hole = CRYSTAL.supercell(1, NYH, NZH)
a, b, c, *_ = cell_hole.lattice_parameters
sy, sz = HOLE_POLARON_DIMENSIONS

ax_3d.add_patch(
    Ellipse(xy=(0, 0), width=sy, height=sz, fc=HCOLOR, **POLARON_SIZE_PARAMS)
)
for atm in cell_hole:
    x, y, z = atm.coords_cartesian
    # Center supercell about (0,0)
    y -= NYH * b / 2
    z -= NZH * c / 2

    born_eff_charge = BORN_CHARGES[atm.element]
    # Hole polaron attracts atoms with negative born effective charge
    uy, uz = 0.3 * born_eff_charge * gaussian2d(y, z, sy, sz) * radial(y, z)

    ax_3d.add_patch(
        Circle(
            xy=(y + uy, z + uz),
            radius=RADII[atm.element] / 1.7,
            clip_on=False,
            zorder=x,
            **ATOM_PATCH_PARAMS[atm.element],
        )
    )
    ax_3d.add_patch(
        Circle(
            xy=(y, z),
            radius=RADII[atm.element] / 1.7,
            clip_on=False,
            zorder=x - a,
            alpha=0.3,
            **ATOM_PATCH_PARAMS[atm.element],
        )
    )

for i in range(-NYH // 2, (NYH + 1) // 2):
    ax_3d.axvline(x=i * b, **LINE_PARAMS)
for j in range(-NZH // 2, (NZH + 1) // 2):
    ax_3d.axhline(y=j * c, **LINE_PARAMS)

ax_3d.set_xlim([-NYH * b / 2, NYH * b / 2])
ax_3d.set_ylim([-NZH * c / 2, NZH * c / 2])

radius = 0.5 * (sy + sz) / 2
for angle in [i * pi / 2 + pi / 4 for i in range(4)]:
    direction = np.array([np.cos(angle), np.sin(angle)])

    ax_3d.add_patch(
        FancyArrowPatch(
            posA=(radius - 1) * direction,
            posB=(radius + 1) * direction,
            **DISP_ARROWS_PARAMS,
        )
    )

for ax, label in zip([ax_1d, ax_3d], ["a)", "b)"]):
    tag_axis(
        ax,
        label,
        x=0,
        y=1.03,
        verticalalignment="bottom",
        horizontalalignment="right",
    )

# Crystal axes
elec_arrow_kwds = dict(
    x=(NZE + 1) * c / 2,
    y=(NXE + 0.5) * a / 2,
    clip_on=False,
    width=0.001,
    head_width=0.5,
    head_length=0.5,
    fc="k",
)
named_arrow(
    ax_1d,
    dx=0,
    dy=-a,
    toffset=(0.5, 0),
    text="$\mathbf{a}$",
    tkwds=dict(
        ha="left",
        va="center",
    ),
    **elec_arrow_kwds,
)
named_arrow(
    ax_1d,
    dx=-c,
    dy=0,
    toffset=(0, 0.5),
    text="$\mathbf{c}$",
    tkwds=dict(
        ha="center",
        va="bottom",
    ),
    **elec_arrow_kwds,
)
hole_arrow_kwds = dict(
    x=-(NYH + 1) * b / 2,
    y=-(NZH + 1) * c / 2,
    clip_on=False,
    width=0.001,
    head_width=0.5,
    head_length=0.5,
    fc="k",
)
named_arrow(
    ax_3d,
    dx=0,
    dy=c,
    toffset=(-0.5, 0),
    text="$\mathbf{c}$",
    tkwds=dict(
        ha="right",
        va="center",
    ),
    **hole_arrow_kwds,
)
named_arrow(
    ax_3d,
    dx=b,
    dy=0,
    toffset=(0, -0.5),
    text="$\mathbf{b}$",
    tkwds=dict(
        ha="center",
        va="top",
    ),
    **hole_arrow_kwds,
)

# Show element names ----------------------------------------------------------
ax_elements.axis("off")
ax_elements.set_aspect(1)
ax_elements.set_xlim([-1, 1])
ax_elements.set_ylim([-1, 1])
for xy, elem, alpha in zip(
    [(-12, 0), (-4, 0), (4, 0), (12, 0)],
    ["Sn", "Sn (eq)", "Se", "Se (eq)"],
    [1, 0.3, 1, 0.3],
):
    ax_elements.add_patch(
        Circle(
            xy=xy,
            **ATOM_PATCH_PARAMS[elem[0:2]],
            radius=2 * RADII[elem[0:2]],
            clip_on=False,
            alpha=alpha,
        )
    )
    ax_elements.text(x=xy[0], y=1, s=elem, va="bottom", ha="center")

plt.tight_layout()
