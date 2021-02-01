import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import itertools as it
from plotutils import FIGURE_WIDTH, named_arrow
from skued import indices_to_text

EWALD_RADIUS = 50
EWALD_RADIUS_XRAY = 5

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 2))

for kx, ky in it.product(range(-5, 6), range(-2, 5)):
    ax.add_patch(
        mpatches.Ellipse(xy=(kx, ky), width=0.07, height=0.3, ec="none", fc="k")
    )

# Ewald spheres
ax.add_patch(
    mpatches.Wedge(
        center=(0, EWALD_RADIUS),
        r=EWALD_RADIUS,
        theta1=260,
        theta2=280,
        fc="none",
        ec="k",
    )
)
ax.add_patch(
    mpatches.Wedge(
        center=(0, EWALD_RADIUS_XRAY),
        r=EWALD_RADIUS_XRAY,
        theta1=180,
        theta2=360,
        fc="none",
        ec="k",
        linestyle="dashed",
    )
)

for (h, k, l) in [(0, 0, 0), (0, 3, 1)]:
    ax.annotate(
        xy=(k, l),
        ha="center",
        va="bottom",
        text=indices_to_text(h, k, l),
        xytext=(k, l + 0.2),
    )

# Lattice vectors
arrow_kwds = dict(
    x=-5, y=-2, length_includes_head=True, width=0.001, head_width=0.1, fc="k"
)
named_arrow(
    ax,
    dx=1,
    dy=0,
    text=r"$\mathbf{b}_2$",
    toffset=(0, -0.1),
    tkwds=dict(va="top", ha="center"),
    **arrow_kwds
)
named_arrow(
    ax,
    dx=0,
    dy=1,
    text=r"$\mathbf{b}_3$",
    toffset=(-0.1, 0),
    tkwds=dict(va="center", ha="right"),
    **arrow_kwds
)

ax.set_xlim([-5.5, 5.5])
ax.set_ylim([-2.5, 4.5])
ax.axis("off")

plt.tight_layout()
