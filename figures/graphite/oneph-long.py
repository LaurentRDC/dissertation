from pathlib import Path
from matplotlib.ticker import FixedFormatter, FixedLocator
from crystals import Crystal
import matplotlib.pyplot as plt
import numpy as np
from plotutils import (
    FIGURE_WIDTH,
    GRAPHITE_ANGLE,
    GRAPHITE_CAMERA_LENGTH,
    GRAPHITE_CENTER,
    ImageGrid,
    draw_hexagon_field,
    tag_axis,
)

INPUT = Path("data") / "graphite"

# Mode ordering of graphite according to the file
# Gra-C_XDM_mode_grid_new2.json
MODE_ORDERING = {
    "LA": 0,
    "TA": 1,
    "ZA": 2,
    "LO1": 3,
    "LO2": 4,
    "LO3": 5,
    "TO1": 6,
    "TO2": 7,
    "TO3": 8,
    "ZO1": 9,
    "ZO2": 10,
    "ZO3": 11,
}
MODES = sorted(MODE_ORDERING.keys())
IN_PLANE_MODES = sorted(set(MODE_ORDERING.keys()) - {"ZA", "ZO1", "ZO2", "ZO3"})

in_plane_refls = filter(
    lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12)
)

reflections = tuple(in_plane_refls)

fig = plt.figure(figsize=(FIGURE_WIDTH, 1.1 * FIGURE_WIDTH))
grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(2, 2),
    cbar_location="top",
)

qx = np.load(INPUT / "oneph" / "qx.npy")
qy = np.load(INPUT / "oneph" / "qy.npy")
bragg_peaks = np.load(INPUT / "oneph" / "bragg_peaks.npy")
cryst = Crystal.from_pwscf(INPUT / "output.out")

# Only Longitudinal modes here
modes = filter(lambda s: s.startswith("L"), IN_PLANE_MODES)
for mode, ax in zip(modes, grid):
    image = np.load(INPUT / "oneph" / f"{mode}.npy")

    # Image is scaled so maximum is always 1
    m = ax.imshow(
        image / image.max(),
        extent=[qx.min(), qx.max(), qy.min(), qy.max()],
        cmap="inferno",
        vmin=0,
        vmax=1,
    )
    ax.scatter(x=bragg_peaks[:, 0], y=bragg_peaks[:, 1], s=1, c="w")

    # Bragg peaks might extend beyond the rest of the image
    ax.set_xlim([qx.min(), qx.max()])
    ax.set_ylim([qy.min(), qy.max()])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    tag_axis(ax, text=f"{mode}")

    draw_hexagon_field(
        ax=ax,
        radius=1.7,
        crystal=cryst,
        color=(0.7, 0.7, 0.7, 1),  # light gray
        linestyle=":",
        reflections=reflections,
    )

cbar = ax.cax.colorbar(
    m,
    ticks=FixedLocator([0, m.get_array().max()]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel(r"$|F_{1j}(\mathbf{q}, \tau=-\infty)|^2$ [a.u.]")

plt.subplots_adjust(bottom=0.01)
