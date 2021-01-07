import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.colorbar import colorbar
from matplotlib.ticker import FixedLocator, FixedFormatter, PercentFormatter
from pathlib import Path
import numpy as np
from crystals import Crystal
from plotutils import (
    ImageGrid,
    FIGURE_WIDTH,
    tag_axis,
    draw_hexagon_field,
    draw_hexagon,
    set_height_auto,
)

INPUT = Path("data") / "graphite" / "debye-waller"

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 1.8))
(ax_rmt, ax_cmp) = ImageGrid(
    fig, 111, nrows_ncols=(1, 2), cbar_mode="each", cbar_location="top"
)

for ax in (ax_rmt, ax_cmp):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

qx = np.load(INPUT / "qx.npy")
qy = np.load(INPUT / "qy.npy")
hot_dw = np.load(INPUT / "hot_dw.npy")
room_temp_dw = np.load(INPUT / "room_temp_dw.npy")

m_rmt = ax_rmt.imshow(
    np.nan_to_num(room_temp_dw / room_temp_dw.max()),
    extent=[qx.min(), qx.max(), qy.max(), qy.min()],
    vmin=0,
    vmax=1,
    cmap="inferno",
)

diff = np.nan_to_num(100 * (hot_dw - room_temp_dw) / room_temp_dw)
m_cmp = ax_cmp.imshow(
    diff,
    extent=[qx.min(), qx.max(), qy.max(), qy.min()],
    cmap="Reds",
)

rmt_cbar = ax_rmt.cax.colorbar(
    mappable=m_rmt,
    ticks=FixedLocator(locs=[0, 1]),
    format=FixedFormatter(["0", r"$1 \leftarrow$      "]),
)
rmt_cbar.ax.set_xlabel("$\sum_s W_s(\mathbf{q})$ [a.u.]")

cmp_cbar = ax_cmp.cax.colorbar(
    mappable=m_cmp, ticks=FixedLocator(locs=[2, 4, 6, 8]), format=PercentFormatter()
)
cmp_cbar.ax.set_xlabel(r"$\Delta \sum_s W_s(\mathbf{q})$ [a.u.]")

tag_axis(ax_rmt, "a)")
tag_axis(ax_cmp, "b)")

plt.subplots_adjust(bottom=0.01)
