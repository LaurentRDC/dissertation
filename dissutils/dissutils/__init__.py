"""
Package with utilities required to render figures
"""
from functools import partial
from itertools import islice

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from crystals.affine import change_of_basis
from mpl_toolkits.axes_grid1 import ImageGrid
from skued import spectrum_colors

from .snse import DatasetInfo, DatasetInfo200

# CONSTANTS -------------------------------------------------------------------

LARGE_FIGURE_WIDTH = 6 + 3 / 4  # inches
MEDIUM_FIGURE_WIDTH = 4 + 1 / 4  # inches
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
    edgecolor="k",
):
    """ Tag an axis with some text, e.g. "a)" """
    ax.text(
        x=x,
        y=y,
        s=text,
        color="k",
        bbox={
            "facecolor": "w",
            "alpha": 1,
            "boxstyle": "round",
            "edgecolor": edgecolor,
        },
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
    cmap = plt.get_cmap("inferno")
    mi, ma = 0.11, 0.75

    if num == 1:
        return [cmap(mi)]
    elif num == 2:
        return [cmap(mi), cmap(ma)]

    step = (ma - mi) / (num - 1)
    return [cmap(mi + i * step) for i in range(num)]
