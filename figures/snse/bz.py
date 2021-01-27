import matplotlib.pyplot as plt
import numpy as np
from functools import reduce
from math import radians, sqrt
from matplotlib.patches import Circle, Rectangle
from plotutils import named_arrow


def named_arrow(ax, x, y, dx, dy, text, tkwds=dict(), **kwargs):
    """ Draw an arrow annotated with some text """
    ax.arrow(x=x, y=y, dx=dx, dy=dy, **kwargs)
    ax.text(x=x + dx / 2, y=y + dy / 2, s=text, **tkwds)


fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
ax.set_aspect("equal")
ax.axis("off")

# reciprocal lattice vectors
arrow_kwds = dict(
    x=0, y=0, length_includes_head=True, width=0.001, head_width=0.07, fc="k"
)

width = 2
height = 1.5

named_arrow(
    ax,
    dx=width / 2,
    dy=0,
    text=r"$\mathbf{b}_3$",
    tkwds=dict(ha="center", va="top"),
    **arrow_kwds
)
named_arrow(
    ax,
    dx=0,
    dy=height / 2,
    text=r"$\mathbf{b}_2$",
    tkwds=dict(va="center", ha="left"),
    **arrow_kwds
)

GAMMA = np.array((0, 0))
Y = np.array((0, 1))
Z = np.array((1, 0))
T = np.array((1, 1))

rect = Rectangle(
    xy=(-width / 2, -height / 2),
    width=width,
    height=height,
    facecolor="None",
    fill=False,
    edgecolor="k",
)
ax.add_patch(rect)

offset = np.array([0.06, 0.06])

for point, label, sym_color in zip(
    [Y, Z],
    [r"$\mathbf{Y}$", r"$\mathbf{Z}$"],
    ["blue", "red"],
):
    for n in range(0, 2):
        point_ = reduce(lambda m1, m2: m1 @ m2, n * [-np.eye(2)], -np.eye(2)) @ point
        p1, p2 = point_
        ax.add_patch(
            Circle(
                xy=(p1 * width / 2, p2 * height / 2),
                radius=0.05,
                edgecolor="None",
                facecolor=sym_color,
            )
        )

    point = point * np.array([width / 2, height / 2])
    ax.annotate(
        xy=point,
        ha="left",
        va="bottom",
        text=label,
        xytext=np.array(point) + offset,
        color=sym_color,
    )

ax.add_patch(Circle(xy=GAMMA, radius=0.05, edgecolor="None", facecolor="k"))

ax.annotate(
    xy=GAMMA,
    ha="right",
    va="top",
    text=r"$\mathbf{\Gamma}$",
    xytext=np.array(GAMMA) - np.abs(offset),
    color="k",
)

color_T = "goldenrod"
R90deg = np.array([[0, -1], [1, 0]])
for n in range(0, 4):
    point_ = reduce(lambda m1, m2: m1 @ m2, n * [R90deg], R90deg) @ T
    p1, p2 = point_
    ax.add_patch(
        Circle(
            xy=(p1 * width / 2, p2 * height / 2),
            radius=0.05,
            edgecolor="None",
            facecolor=color_T,
        )
    )

ax.annotate(
    xy=T,
    ha="left",
    va="bottom",
    text=r"$\mathbf{T}$",
    xytext=np.array(T) * np.array([width / 2, height / 2]) + offset,
    color=color_T,
)
