import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from dissutils import (
    LARGE_FIGURE_WIDTH,
    GRID_AXES_PAD,
    CBAR_SIZE,
    discrete_colors,
    tag_axis,
)
from crystals import Crystal
from skued import pelectrostatic, powdersim


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, (3 / 5) * LARGE_FIGURE_WIDTH))
gs = gridspec.GridSpec(nrows=3, ncols=2, figure=figure)
ax_gold = figure.add_subplot(gs[0, 0])
ax_graphite = figure.add_subplot(gs[1, 0], sharex=ax_gold)
ax_m1 = figure.add_subplot(gs[2, 0], sharex=ax_gold)
ax_batio3 = figure.add_subplot(gs[:, 1])

s = np.linspace(2.5, 7, 512)
q = s / (4 * np.pi)

for ax, cryst, color, label in zip(
    [ax_gold, ax_graphite, ax_m1],
    ["Au", "C", "vo2-m1"],
    discrete_colors(3),
    ["a)", "b)", "c)"],
):
    I = powdersim(Crystal.from_database(cryst), s, fwhm_g=0.04, fwhl_l=0.08)
    ax.plot(q, I / I.max(), color=color)
    ax.set_yticks([])
    ax.xaxis.set_visible(False)
    tag_axis(ax, label, x=0.95, y=0.89, horizontalalignment="right")

ax_graphite.set_ylabel("Diffracted intensity [a.u.]")

ax_m1.xaxis.set_visible(True)
ax_m1.set_xlim([q.min(), q.max()])
ax_m1.set_xlabel("Scattering vector $|\mathbf{q}|$ [$1/\AA$]")

# Electrostatic potential
extent = np.linspace(-5, 5, 256)
xx, yy = np.meshgrid(extent, extent)
potential = pelectrostatic(Crystal.from_database("BaTiO3_cubic"), xx, yy)

m = ax_batio3.imshow(
    potential,
    extent=[xx.min(), xx.max(), yy.min(), yy.max()],
    vmax=1500,
    cmap="CMRmap_r",
)
ax_batio3.set_xlabel("Horizontal distance [$\AA$]")
ax_batio3.set_ylabel("Vertical distance [$\AA$]")
ax_batio3.yaxis.tick_right()
ax_batio3.yaxis.set_label_position("right")
tag_axis(ax_batio3, "d)")

cax = make_axes_locatable(ax_batio3).append_axes("top", size=CBAR_SIZE, pad=0.03)
cbar = plt.colorbar(mappable=m, cax=cax, orientation="horizontal")
cbar.set_label("Projected potential [$V \cdot \AA$]")
cax.xaxis.tick_top()
cax.xaxis.set_label_position("top")

plt.subplots_adjust(
    top=0.894, bottom=0.159, left=0.055, right=0.898, hspace=0.0, wspace=0.054
)
