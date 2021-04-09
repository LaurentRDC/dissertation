import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedFormatter, FixedLocator
import numpy as np
from skued import diffread, align
from skimage.filters import gaussian
from pathlib import Path
from dissutils import MEDIUM_FIGURE_WIDTH, tag_axis, ImageGrid

DOWNSAMPLING = 8

ref = np.clip(diffread(Path("data") / "appendix" / "Cr_1.tif"), 0, 200)
im = np.clip(diffread(Path("data") / "appendix" / "Cr_2.tif"), 0, 200)

mask = np.ones_like(ref, dtype=bool)
mask[0:1200, 950:1250] = False

im = im[::DOWNSAMPLING, ::DOWNSAMPLING]
ref = ref[::DOWNSAMPLING, ::DOWNSAMPLING]
mask = mask[::DOWNSAMPLING, ::DOWNSAMPLING]

shifted = align(image=im, reference=ref, mask=mask)

figure = plt.figure(figsize=(5, MEDIUM_FIGURE_WIDTH))
axes = ImageGrid(
    fig=figure, rect=111, nrows_ncols=(2, 2), cbar_mode="edge", cbar_location="right"
)

diff0 = gaussian(ref - im, sigma=1)
diff0 -= np.mean(diff0)

diff_shifted = gaussian(ref - shifted, sigma=1)
diff_shifted -= np.mean(diff_shifted)

ax1, ax2, ax3, ax4 = axes
map_I = ax1.imshow(ref, vmin=0, vmax=200, cmap="inferno")
ax2.imshow(im, vmin=0, vmax=200, cmap="inferno")
map_diff = ax3.imshow(diff0, vmin=-0.5, vmax=0.5, cmap="RdBu_r")
ax4.imshow(diff_shifted, vmin=-0.5, vmax=0.5, cmap="RdBu_r")


ax1.cax.colorbar(
    mappable=map_I,
    ticks=FixedLocator(locs=[0, 200]),
    format=FixedFormatter(["0", "1"]),
)
ax1.cax.set_ylabel("Intensity [a.u.]")

ax3.cax.colorbar(
    mappable=map_diff,
    ticks=FixedLocator(locs=[-0.4, 0, 0.4]),
    format=FixedFormatter(["-¼%", "0%", "¼%"]),
)
ax3.cax.set_ylabel("Intensity diff. [%]")

# Show mask as partially-transparent rectangle
for ax in [ax1, ax2]:
    ax.add_patch(
        Rectangle(
            xy=(950 / DOWNSAMPLING, 0),
            width=(300 / DOWNSAMPLING),
            height=(1200 / DOWNSAMPLING),
            ec="k",
            fc="dimgray",
            alpha=0.5,
        )
    )

for ax, label in zip([ax1, ax2, ax3, ax4], ["a)", "b)", "c)", "d)"]):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    tag_axis(ax, label)

plt.tight_layout()
