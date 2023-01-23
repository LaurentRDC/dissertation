from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from iris import DiffractionDataset
from matplotlib.patches import Circle
from matplotlib.ticker import FixedFormatter, FixedLocator
from skimage.transform import rotate
from skued import autocenter, diffread

from dissutils import GRAPHITE_ANGLE, LARGE_FIGURE_WIDTH, ImageGrid, tag_axis

DOWNSAMPLING = 8

im_Cr = np.clip(diffread(Path("data") / "appendix" / "Cr_1.tif"), 0, 200)
mask_Cr = np.ones_like(im_Cr, dtype=bool)
mask_Cr[0:1200, 950:1250] = False

with DiffractionDataset(
    Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5", mode="r"
) as dset:
    _c, _r = np.asarray(dset.center).astype(int)
    mask_graphite = rotate(
        dset.valid_mask, angle=GRAPHITE_ANGLE, center=(_c, _r), mode="reflect"
    )
    im_graphite = rotate(
        dset.diff_eq(), angle=GRAPHITE_ANGLE, center=(_c, _r), mode="reflect"
    )

im_Cr = im_Cr[::DOWNSAMPLING, ::DOWNSAMPLING]
mask_Cr = mask_Cr[::DOWNSAMPLING, ::DOWNSAMPLING]

im_graphite = im_graphite[::DOWNSAMPLING, ::DOWNSAMPLING]
mask_graphite = mask_graphite[::DOWNSAMPLING, ::DOWNSAMPLING]

figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2))
axes = ImageGrid(figure, 111, nrows_ncols=(1, 2), cbar_location="top")

for ax, im, mask, label in zip(
    axes, [im_Cr, im_graphite], [mask_Cr, mask_graphite], ["a)", "b)"]
):
    m = ax.imshow(im, vmin=0, vmax=200, cmap="inferno")
    r, c = autocenter(im, mask)

    ax.axhline(y=r, color="w", linestyle="--", linewidth=0.5)
    ax.axvline(x=c, color="w", linestyle="--", linewidth=0.5)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    tag_axis(ax, label)

axes[-1].cax.colorbar(
    mappable=m, ticks=FixedLocator(locs=[0, 200]), format=FixedFormatter(["0", "1"])
)
axes[-1].cax.set_xlabel("Intensity [a.u.]")
axes[-1].cax.xaxis.set_label_position('top')
axes[-1].cax.xaxis.tick_top()

plt.subplots_adjust(bottom=0.01)
