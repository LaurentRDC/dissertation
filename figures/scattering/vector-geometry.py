import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from itertools import product


figure, ax = plt.subplots(1, 1, figsize=(3, 3))

ax.set_xlim([-1.05, 2.05])
ax.set_ylim([-1.05, 2.05])
ax.axis("off")

for x, y in product([-2, -1, 0, 1, 2], repeat=2):
    ax.scatter(x=x, y=y, s=2, c="k")
    ax.add_patch(
        Rectangle(
            xy=(x - 1 / 2, y - 1 / 2),
            width=1,
            height=1,
            facecolor="None",
            fill=False,
            edgecolor="k",
            linestyle=(0, (5, 5)),
            linewidth=0.5,
        )
    )
for x, y in product([0, 1], repeat=2):

    ax.add_patch(
        Rectangle(
            xy=(x - 1 / 2, y - 1 / 2),
            width=1,
            height=1,
            facecolor="None",
            fill=False,
            edgecolor="k",
            linewidth=0.7,
        )
    )

kx, ky = 1.4, 0.7
arrow_kwds = dict(
    arrowstyle="-|>", shrinkA=0.3, shrinkB=0.3, mutation_scale=6, fc="k", ec="k"
)
ax.add_patch(FancyArrowPatch(posA=(0, 0), posB=(1, 1), **arrow_kwds))
ax.add_patch(FancyArrowPatch(posA=(0, 0), posB=(kx, ky), **arrow_kwds))
ax.add_patch(FancyArrowPatch(posA=(1, 1), posB=(kx, ky), **arrow_kwds))

ax.text(x=1 / 2, y=1 / 2 + 0.1, s="$\mathbf{H}_{\mathbf{q}}$", ha="right", va="bottom")
ax.text(x=kx / 2, y=ky / 2 - 0.1, s="$\mathbf{q}$", ha="left", va="top")
ax.text(
    x=1 + (kx - 1) / 2, y=1 + (ky - 1) / 2, s="$\mathbf{k}_0$", ha="left", va="bottom"
)
ax.text(x=0, y=-0.1, s="$(000)$", ha="center", va="top")
ax.text(x=1, y=1.1, s="$(110)$", ha="center", va="bottom")
