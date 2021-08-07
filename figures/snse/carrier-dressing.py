import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Ellipse
import matplotlib.gridspec as gridspec
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis

UNIFORM_COLOR, DISTORT_COLOR = discrete_colors(2)

fig = plt.figure(
    figsize=(MEDIUM_FIGURE_WIDTH, MEDIUM_FIGURE_WIDTH/1.5),
)

gs = gridspec.GridSpec(nrows=2, ncols=2, figure=fig, width_ratios=[0.5, 2])
ax_dos_u = fig.add_subplot(gs[0, 0])  # uniform
ax_dos_d = fig.add_subplot(gs[1, 0], sharex=ax_dos_u)  # distorted
ax_fe = fig.add_subplot(gs[:, 1])  # free energy

for _ax, color, fermi in zip(
    [ax_dos_u, ax_dos_d], [UNIFORM_COLOR, DISTORT_COLOR], [0.75, 0.25]
):
    _ax.set_xlim([0, 0.8])
    _ax.set_ylim([-1, 1])
    _ax.axis("off")

    _ax.axvline(x=0, linewidth=1, clip_on=False, color="k")
    _ax.axhline(y=-1, linewidth=1, clip_on=False, color="k")
    _ax.axhline(y=fermi, linewidth=0.5, linestyle="dashed", color="k")

    # Mimicking yaxis label
    _ax.text(s="$E$", x=-0.05, y=0, rotation=90, ha="right", va="center")

    # Upper band
    _ax.add_patch(Circle(xy=(0, 1.2), radius=0.75, fc="none", ec="k"))

    _x = np.linspace(0, 0.7, num=128)
    lb = 1.2 - np.sqrt(0.75 ** 2 - _x ** 2)
    _ax.fill_between(
        _x, y1=lb, y2=fermi, where=lb < fermi, color=color
    )

    # lower band
    _ax.add_patch(
        Circle(xy=(0, -1.2), radius=0.75, fc=color, ec="k")
    )

ax_dos_u.text(s="$E_F$", x=0.85, y=2*fermi, va="bottom", ha="left", transform=ax_dos_u.transData)

# polaron band
ax_dos_d.add_patch(Ellipse(xy=(0, 0), width=0.6, height=0.25, fc=DISTORT_COLOR, ec="k"))
ax_dos_d.text(s="Polaron", x=0.35, y=0, va="center", ha="left")

# Free energy diagram ---------------------------------------------------------


ax_fe.set_xlim([-1, 2])
ax_fe.set_ylim([0, 1.5])
ax_fe.axis("off")
ax_fe.axvline(x=-1, linewidth=1, color="k", clip_on=False)
ax_fe.axhline(y=0, linewidth=1, color="k", clip_on=False)
ax_fe.yaxis.tick_right()
ax_fe.yaxis.set_label_position("right")

q_uniform = np.linspace(-1, 1, num=128)
q_distorted = np.linspace(0, 2, num=128)
ax_fe.plot(q_uniform, 0.6 * q_uniform ** 2 + 0.75, linestyle="solid", color=UNIFORM_COLOR)
ax_fe.vlines(x=0, ymin=0, ymax=0.75, linestyles="dashed", linewidth=1, colors="k")
ax_fe.plot(
    q_distorted, (q_distorted - 1) ** 2 + 0.25, linestyle="solid", color=DISTORT_COLOR
)
ax_fe.vlines(x=1, ymin=0, ymax=0.25, linestyles="dashed", linewidth=1, colors="k")

# Show phonon dressing
x1 = np.asarray((0, 0.85))
x2 = np.asarray((1, 0.35))
dd = x2 - x1
ax_fe.add_patch(
    FancyArrowPatch(
        posA=x1, posB=x2, arrowstyle="-|>", mutation_scale=10, fc="k", linewidth=2
    )
)
tx, ty = x1 + (x2 - x1) / 2
ax_fe.text(s="phonon\ndressing", x=tx + 0.05, y=ty, ha="left", va="bottom")

# Mimicking axis labels
ax_fe.text(
    s="Lattice distortion coordinate", x=3 / 2 - 1, y=-0.1, va="top", ha="center"
)
ax_fe.text(s="$E$", x=-1.05, y=1, rotation=90, ha="right", va="center")

ax_fe.text(s="No\ndistortion", x=-0.05, y=0, ha='right', va='bottom')
ax_fe.text(s="Polaron\ndistortion", x=1.05, y=0, ha='left', va='bottom')


ax_fe.text(s="Delocalized", x=1, y=1.35, color=UNIFORM_COLOR, ha='center', va='bottom')
ax_fe.text(s="Localized", x=2, y=1.25, color=DISTORT_COLOR, ha='center', va='bottom')

tag_axis(ax_dos_u, 'a)', x=-0.15, horizontalalignment='right')
tag_axis(ax_dos_d, 'b)', x=-0.15, horizontalalignment='right')
tag_axis(ax_fe, 'c)', y=0.975)
