"""
Create the title page for the Eisvogel template
"""
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import itertools

random.seed(23)
np.random.seed(23)

# Hexagon radius in mm
RADIUS = 8

# Dimensions of the bounding box of the hexagons
WIDTH = math.sqrt(3) * RADIUS
HEIGHT = 2 * RADIUS

# Page dimensions in mm
# These are for A4 paper
PAGEWIDTH = 210
PAGEHEIGHT = 297
mm_to_in = 0.03937008

# CMAP "cool" also looks, uh, cool
cmap = plt.get_cmap("inferno")

def draw_hexagon(ax, center, color):
    # because the dimensions of the axes have been set in mm,
    # the dimensions of thie hexagon are given in mm as well.
    ax.add_patch(
        mpatches.RegularPolygon(
            xy=center,
            numVertices=6,
            # Due to rounding errors I think, it is better to make the hexagons
            # a little larger than the RADIUS parameter
            radius=RADIUS + 0.2,
            facecolor=color,
            edgecolor="none",
            orientation=0,
            fill=True,
        )
    )

figure, ax = plt.subplots(
    1, 1, figsize=(PAGEWIDTH * mm_to_in, PAGEHEIGHT * mm_to_in), frameon=False
)

# Dimensions of the page in mm
ax.set_xlim([0, PAGEWIDTH])
ax.set_ylim([0, PAGEHEIGHT])

for offset_x, offset_y in [(0, 0), (WIDTH / 2, (3 / 2) * RADIUS)]:
    rows = np.arange(start=offset_x, stop=1.05 * PAGEWIDTH, step=WIDTH)
    columns = np.arange(start=offset_y, stop=1.05 * PAGEHEIGHT, step=3 * RADIUS)
    for x, y in itertools.product(rows, columns):
        color = cmap(
            math.hypot(x, y) / math.hypot(PAGEWIDTH, PAGEHEIGHT) + random.gauss(0, 0.01)
        )
        draw_hexagon(ax, center=(x, y), color=color)

ax.axis("off")
plt.subplots_adjust(top=1, bottom=0, left=0, right=1)
plt.savefig("build/titlepage.pdf", dpi=600)