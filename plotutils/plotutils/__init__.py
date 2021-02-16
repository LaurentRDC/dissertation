"""
Package with utilities required to render figures
"""
from functools import partial
from itertools import islice

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from crystals.affine import change_of_basis
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import ImageGrid

from .snse_datasets import DatasetInfo, DatasetInfo200

# CONSTANTS -------------------------------------------------------------------

FIGURE_WIDTH = 6 + 3 / 4  # inches
FONTSIZE = 10

# Diffraction patterns are rotated by 8 degrees clockwise from aligned
GRAPHITE_ANGLE = 8  # degrees
GRAPHITE_CAMERA_LENGTH = 0.25  # meters

SNSE_CAMERA_LENGTH = 0.2939  # meters

# -----------------------------------------------------------------------------

# Wrapper around ImageGrid with some parameters fixed
GRID_AXES_PAD = 0.05
GRID_CBAR_PAD = GRID_AXES_PAD
CBAR_SIZE = 0.1
ImageGrid = partial(
    ImageGrid,
    axes_pad=GRID_AXES_PAD,
    cbar_pad=GRID_CBAR_PAD,
    cbar_mode="single",
    cbar_size=CBAR_SIZE,
)


def named_arrow(ax, x, y, dx, dy, text, toffset=(0, 0), tkwds=dict(), **kwargs):
    """ Draw an arrow annotated with some text """
    ox, oy = toffset
    ax.arrow(x=x, y=y, dx=dx, dy=dy, **kwargs)
    ax.text(x=x + ox + dx / 2, y=y + oy + dy / 2, s=text, **tkwds)


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


def set_height_auto(fig, width):
    """
    Modify the figure height to minimize whitespace, given a width.
    """
    axes = fig.axes

    def get_area(ax):
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        return bbox.width * bbox.height

    largest_ax = sorted(fig.axes, key=get_area)[0]
    bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    factor = width / bbox.width
    fig.set_size_inches(w=width, h=factor * bbox.height)


def draw_hexagon_field(
    ax,
    radius,
    crystal,
    reflections,
    orientation=np.deg2rad(30),
    center=(0, 0),
    **kwargs
):
    """ Fill the plot with hexagons centered at the Bragg points """
    center = np.array(center)
    from_frac = change_of_basis(np.array(crystal.reciprocal_vectors), np.eye(3))
    for refl in reflections:
        xyz = from_frac @ refl
        draw_hexagon(
            ax,
            center=xyz[0:2] + center,
            radius=radius,
            orientation=orientation,
            **kwargs,
        )


def draw_hexagon(
    ax,
    center,
    radius,
    orientation=np.deg2rad(30),
    color="w",
    facecolor="none",
    **kwargs
):
    """ Draw a hexagon within an Axes object"""
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 1
    hexagon = mpatches.RegularPolygon(
        xy=center,
        numVertices=6,
        radius=radius,
        orientation=orientation,
        facecolor=facecolor,
        edgecolor=color,
        **kwargs,
    )
    ax.add_patch(hexagon)


def discrete_colors(num):
    """ Returns a list of discrete colors to plot, for example, various time-traces. """
    return ["goldenrod", "red", "blue", "indigo"][0:num]


def box_errorbars(ax, xdata, ydata, xerr, yerr, colors):
    """
    Create error boxes instead of error bars
    """
    # Create list for all the error patches
    errorboxes = []
    # Loop over data points; create box from errors at each point
    for x, y in zip(xdata, ydata):
        rect = Rectangle((x - xerr, y - yerr), 2 * xerr, 2 * yerr)
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, edgecolor=colors, facecolor=colors, alpha=0.7)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot invisible scatter points to get automatic plot bounds that
    # make sense
    artists = ax.scatter(xdata, ydata, s=0, c=colors)

    return artists
