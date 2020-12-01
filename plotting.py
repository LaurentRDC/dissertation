import matplotlib.pyplot as plt
import numpy as np
from itertools import chain

FIGURE_WIDTH = 6 + 3 / 4

FONTSIZE = 8

# Initiate default parameters
plt.rcParams["font.size"] = FONTSIZE


def tag_axis(
    ax,
    text,
    x=0.02,
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
