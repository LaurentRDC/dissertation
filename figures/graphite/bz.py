import matplotlib.pyplot as plt
import numpy as np
from functools import reduce
from math import radians, sqrt
from matplotlib.patches import Circle, RegularPolygon
from plotutils import named_arrow, discrete_colors


def named_arrow(ax, x, y, dx, dy, text, tkwds=dict(), **kwargs):
    """ Draw an arrow annotated with some text """
    ax.arrow(x=x, y=y, dx=dx, dy=dy, **kwargs)
    ax.text(x=x + dx / 2, y=y + dy / 2, s=text, **tkwds)


fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
ax.set_aspect("equal")
ax.axis("off")

# reciprocal lattice vectors
arrow_kwds = dict(
    x=0, y=0, length_includes_head=True, width=0.001, head_width=0.07, fc="k"
)

named_arrow(
    ax,
    dx=-1 / 2,
    dy=-sqrt(3) / 2,
    text=r"$\mathbf{b}_1$",
    tkwds=dict(ha="right", va="bottom"),
    **arrow_kwds
)
named_arrow(
    ax,
    dx=1,
    dy=0,
    text=r"$\mathbf{b}_2$",
    tkwds=dict(va="bottom", ha="center"),
    **arrow_kwds
)

GAMMA = np.array((0, 0))
K = np.array((-1, 0))
M = np.array((0, sqrt(3) / 2))

R60deg = np.array([[1 / 2, -sqrt(3) / 2], [sqrt(3) / 2, 1 / 2]])

hexagon = RegularPolygon(
    xy=(0, 0),
    numVertices=6,
    orientation=radians(30),
    radius=1,
    facecolor="None",
    fill=False,
    edgecolor="k",
)
ax.add_patch(hexagon)
offset = np.array([-0.06, 0.06])
for point, label, sym_color in zip(
    [M, GAMMA, K],
    [r"$\mathbf{M}$", r"$\mathbf{\Gamma}$", r"$\mathbf{K}$"],
    discrete_colors(3),
):
    for n in range(0, 6):
        point_ = reduce(lambda m1, m2: m1 @ m2, n * [R60deg], R60deg) @ point
        ax.add_patch(
            Circle(xy=point_, radius=0.05, edgecolor="None", facecolor=sym_color)
        )

    ax.annotate(
        xy=point,
        ha="right",
        text=label,
        xytext=np.array(point) + offset,
        color=sym_color,
    )
