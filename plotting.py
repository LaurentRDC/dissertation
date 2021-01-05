from crystals import Crystal
from crystals.affine import change_of_basis

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.patches as mpatches

import numpy as np
from itertools import chain
from functools import partial

# Wrapper around ImageGrid with some parameters fixed
GRID_AXES_PAD = 0.05
GRID_CBAR_PAD = GRID_AXES_PAD
ImageGrid = partial(
    ImageGrid,
    axes_pad=GRID_AXES_PAD,
    cbar_pad=GRID_CBAR_PAD,
    cbar_mode="single",
    cbar_size=0.1,
)

# CONSTANTS -------------------------------------------------------------------

FIGURE_WIDTH = 6 + 3 / 4  # inches

FONTSIZE = 10

# Diffraction patterns are rotated by 8 degrees clockwise from aligned
GRAPHITE_ANGLE = 8  # degrees
GRAPHITE_CAMERA_LENGTH = 0.25  # centi-meters
_peak1 = np.array((754, 905))
_peak2 = np.array((1265, 1318))
GRAPHITE_CENTER = np.array(0.5 * (_peak1 + _peak2), dtype=np.int)

# -----------------------------------------------------------------------------

# Initiate default parameters
plt.rcParams["font.size"] = FONTSIZE


def tag_axis(
    ax,
    text,
    x=0.05,
    y=0.95,
    horizontalalignment="left",
    verticalalignment="top",
):
    """ Tag an axis with some text, e.g. "a)" """
    ax.text(
        x=x,
        y=y,
        s=text,
        color="k",
        bbox={"facecolor": "w", "alpha": 1, "boxstyle": "round"},
        horizontalalignment=horizontalalignment,
        verticalalignment=verticalalignment,
        transform=ax.transAxes,
    )


def draw_hexagon_field(
    ax, radius, crystal, reflections, orientation=np.deg2rad(30), **kwargs
):
    """ Fill the plot with hexagons centered at the Bragg points """
    from_frac = change_of_basis(np.array(crystal.reciprocal_vectors), np.eye(3))
    for refl in reflections:
        xyz = from_frac @ refl
        draw_hexagon(
            ax, center=xyz[0:2], radius=radius, orientation=orientation, **kwargs
        )


def draw_hexagon(ax, center, radius, orientation=np.deg2rad(30), color="w", **kwargs):
    """ Draw a hexagon within an Axes object"""
    hexagon = mpatches.RegularPolygon(
        xy=center,
        numVertices=6,
        radius=radius,
        orientation=orientation,
        facecolor="None",
        fill=False,
        edgecolor=color,
        **kwargs,
    )
    ax.add_patch(hexagon)
