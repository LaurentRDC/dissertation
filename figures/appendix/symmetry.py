import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from skued import nfold, diffread, autocenter
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedFormatter, FixedLocator
from dissutils import LARGE_FIGURE_WIDTH, ImageGrid, tag_axis

DOWNSAMPLING = 4

im = diffread(Path("data") / "appendix" / "tise2.tif")
mask = np.logical_not(
    diffread(Path("data") / "appendix" / "tise2_mask.tif").astype(bool)
)
r, c = autocenter(im, mask=mask)
sym = nfold(im=im, mod=6, center=(c, r), mask=mask)

im = im[::DOWNSAMPLING, ::DOWNSAMPLING]
mask = mask[::DOWNSAMPLING, ::DOWNSAMPLING]
sym = sym[::DOWNSAMPLING, ::DOWNSAMPLING]

# Zoom on center
l = 689 // DOWNSAMPLING
r, c = int(r) // DOWNSAMPLING, int(c) // DOWNSAMPLING
im = im[r - l : r + l, c - l : c + l]
mask = mask[r - l : r + l, c - l : c + l]
sym = sym[r - l : r + l, c - l : c + l]

fig = plt.figure(figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 2), cbar_location="top")

for ax, image, label in zip(grid, [im, sym], ["a)", "b)"]):
    m = ax.imshow(
        image,
        vmin=0,
        vmax=200,
        cmap="inferno",
    )
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    tag_axis(ax, text=label)

# Show mask as partially-transparent rectangle
grid[0].add_patch(
    Rectangle(
        xy=(598 / DOWNSAMPLING, 0),
        width=(180 / DOWNSAMPLING),
        height=(780 / DOWNSAMPLING),
        ec="k",
        fc="dimgray",
        alpha=0.5,
    )
)

cbar = grid[0].cax.colorbar(
    m,
    ticks=FixedLocator(locs=[0, 200]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")

plt.subplots_adjust(bottom=0.01)
