from pathlib import Path

import matplotlib.pyplot as plt
from skimage.io import imread

from dissutils import MEDIUM_FIGURE_WIDTH

DATADIR = Path("images")
DOWNSAMPLING = 2

fig, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 0.8 * MEDIUM_FIGURE_WIDTH))

im = imread(DATADIR / "diamond_knife_2.jpg")[::DOWNSAMPLING, ::DOWNSAMPLING, :]

ax.imshow(im, extent=[0, 1, 1, 0])
ax.axis("off")

# Annotate interesting parts
for text, xy, xytext in zip(
    ["crystal", "knife", "section"],
    [(0.631, 0.332), (0.608, 0.367), (0.365, 0.402)],
    [(0.6, 0.2), (0.4, 0.3), (0.2, 0.2)],
):

    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        arrowprops=dict(facecolor="black", arrowstyle="->"),
        horizontalalignment="center",
        verticalalignment="bottom",
    )

plt.subplots_adjust(bottom=0, top=1)
