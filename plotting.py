import matplotlib.pyplot as plt
import numpy as np
from itertools import chain

FIGURE_WIDTH = 6 + 3 / 4

FONTSIZE = 8


def tag_axis(
    ax,
    text,
    fontsize=FONTSIZE,
    x=0.05,
    y=0.05,
    horizontalalignment="left",
    verticalalignment="bottom",
):
    """ Tag an axis with some text, e.g. "a)" """
    ax.text(
        x=x,
        y=y,
        s=text,
        color="k",
        bbox={"facecolor": "w", "alpha": 1, "boxstyle": "round"},
        horizontalalignment=horizontalalignment,
        fontsize=fontsize,
        verticalalignment=verticalalignment,
        transform=ax.transAxes,
    )


def normalize_axes_fontsize(ax, fontsize=FONTSIZE):
    """ Set a uniform font size for axes components (e.g. tick labels, axis labels, titles, etc.) """
    items = chain(
        [ax.title, ax.xaxis.label, ax.yaxis.label],
        ax.get_xticklabels(),
        ax.get_yticklabels(),
    )
    for item in items:
        item.set_fontsize(FONTSIZE)
