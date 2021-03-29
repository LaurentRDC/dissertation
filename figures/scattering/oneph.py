from pathlib import Path
from crystals import Crystal
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from plotutils import MEDIUM_FIGURE_WIDTH, CBAR_SIZE, draw_hexagon_field

INPUT = Path("data") / "graphite"
DOWNSAMPLING = 4

in_plane_refls = filter(
    lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12)
)
reflections = tuple(in_plane_refls)

figure, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 4))

divider = make_axes_locatable(ax)
cbar_ax = divider.append_axes("top", size=CBAR_SIZE, pad=0.05)

qx = np.load(INPUT / "oneph" / "qx.npy")[::DOWNSAMPLING, ::DOWNSAMPLING]
qy = np.load(INPUT / "oneph" / "qy.npy")[::DOWNSAMPLING, ::DOWNSAMPLING]
bragg_peaks = np.load(INPUT / "oneph" / "bragg_peaks.npy")
cryst = Crystal.from_pwscf(INPUT / "output.out")
image = np.load(INPUT / "oneph" / f"LA_oneph.npy")[::DOWNSAMPLING, ::DOWNSAMPLING]

# Image is scaled so maximum is always 1
m = ax.imshow(
    image / image.max(),
    extent=[qx.min(), qx.max(), qy.min(), qy.max()],
    cmap="CMRmap_r",
    vmin=0,
    vmax=1,
)
ax.scatter(x=bragg_peaks[:, 0], y=bragg_peaks[:, 1], s=1, c="k")

# Bragg peaks might extend beyond the rest of the image
ax.set_xlim([qx.min(), qx.max()])
ax.set_ylim([qy.min(), qy.max()])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

draw_hexagon_field(
    ax=ax,
    radius=1.7,
    crystal=cryst,
    color=(0.7, 0.7, 0.7, 1),  # light gray
    alpha=0.5,
    reflections=reflections,
)

plt.colorbar(mappable=m, cax=cbar_ax, orientation="horizontal")
cbar_ax.xaxis.set_ticks([0, 1])
cbar_ax.xaxis.tick_top()
cbar_ax.xaxis.set_label_position("top")
cbar_ax.set_xlabel(r"$|F_{1\lambda}(\mathbf{q}, \tau<0)|^2$ [a.u.]")

plt.subplots_adjust(bottom=0.01)
