from skimage.io import imread
import matplotlib.pyplot as plt
from pathlib import Path

from dissutils import LARGE_FIGURE_WIDTH, ImageGrid, tag_axis

DATADIR = Path("images")
DOWNSAMPLING = 2

fig = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 0.4 * LARGE_FIGURE_WIDTH))
axes = ImageGrid(fig, 111, nrows_ncols=(1, 2), cbar_mode="none")

for (ax, fname, tag) in zip(axes, ["diamond_knife_1.jpg", "diamond_knife_2.jpg"], "ab"):
    im = imread(DATADIR / fname)[::DOWNSAMPLING, ::DOWNSAMPLING, :]

    ax.imshow(im, extent=[0, 1, 1, 0])
    tag_axis(ax, tag + ")")
    ax.axis("off")

# Annotate interesting parts
ax1, ax2 = axes

for text, xy, xytext in zip(
    ["crystal", "knife", "section"],
    [(0.450, 0.270), (0.500, 0.320), (0.5, 0.53)],
    [(0.15, 0.270), (0.8, 0.32), (0.15, 0.47)],
):
    ax1.annotate(
        text,
        xy=xy,
        xytext=xytext,
        arrowprops=dict(facecolor="black", arrowstyle="->"),
        horizontalalignment="center",
        verticalalignment="center",
    )

for text, xy, xytext in zip(
    ["crystal", "knife", "section"],
    [(0.631, 0.332), (0.608, 0.367), (0.365, 0.402)],
    [(0.6, 0.2), (0.4, 0.3), (0.2, 0.2)],
):

    ax2.annotate(
        text,
        xy=xy,
        xytext=xytext,
        arrowprops=dict(facecolor="black", arrowstyle="->"),
        horizontalalignment="center",
        verticalalignment="bottom",
    )

plt.subplots_adjust(bottom=0, top=1)
