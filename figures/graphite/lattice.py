""" Visualization of the lattice structure of graphite"""

from itertools import groupby
from math import sqrt

import matplotlib.pyplot as plt
from crystals import Crystal
from matplotlib.patches import Circle

from dissutils import FONTSIZE, discrete_colors, named_arrow

getlast = lambda pos: round(pos[-1], 4)

graphite = Crystal.from_database("C")
positions = [atom.coords_cartesian for atom in graphite.supercell(4, 4, 1)]

layers = list()
for key, layer in groupby(sorted(positions, key=getlast), getlast):
    layers.append(list(layer))  # Store group iterator as a list

fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_aspect("equal")
ax.axis("off")

bot_color, top_color = discrete_colors(2)
for layer, color, radius in zip(layers, [bot_color, top_color], [0.12, 0.08]):

    for (x, y, z) in layer:
        ax.add_patch(
            Circle(xy=(x, y), radius=radius, edgecolor="None", facecolor=color)
        )


a, b, _ = graphite.lattice_vectors
arrow_kwds = dict(
    x=b[0],
    y=b[1],
    length_includes_head=True,
    width=0.001,
    head_width=0.07,
    fc="k",
    clip_on=False,
)
named_arrow(
    ax,
    dx=a[0],
    dy=a[1],
    text=r"$\mathbf{a}_1$",
    tkwds=dict(ha="right", va="top"),
    **arrow_kwds
)
named_arrow(
    ax,
    dx=b[0],
    dy=b[1],
    text=r"$\mathbf{a}_2$",
    tkwds=dict(va="bottom", ha="right"),
    toffset=(-0.05, 0),
    **arrow_kwds
)

txtkwargs = dict(
    fontsize=FONTSIZE + 2, weight="bold", transform=ax.transData, ha="center"
)
ax.text(2 * a[0] / 3, 5.4, s="A", color=top_color, **txtkwargs)
ax.text(4 * a[0] / 3, 5.4, s="B", color=bot_color, **txtkwargs)

ax.set_xlim([-0.3, 5.3])
ax.set_ylim([-0.3, 5.3])
