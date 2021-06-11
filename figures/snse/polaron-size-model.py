"""
Visualization of the effect of polaron localization (from large to small)
on the diffuse intensity profile
"""

from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import skued
from crystals import Crystal
from dissutils import (
    LARGE_FIGURE_WIDTH,
    CBAR_SIZE,
    GRID_AXES_PAD,
    tag_axis,
)

DATADIR = Path("data") / "snse"
CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")


def polaron(q, A, rp):
    return A * q * rp ** 3 * np.exp(-(q ** 2 * rp ** 2) / 4)


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 4))
gs = gridspec.GridSpec(2, 2, width_ratios=[2.5, 1])
ax = figure.add_subplot(gs[:, 0])
ax_small = figure.add_subplot(gs[1, 1])
ax_large = figure.add_subplot(gs[0, 1])

cax = make_axes_locatable(ax).append_axes("top", size=CBAR_SIZE, pad=GRID_AXES_PAD)

ks_ = np.linspace(0, 1, num=1024)  # Gamma to T distance is 1 inv Angs
kk, rr = np.meshgrid(ks_, np.linspace(2, 9, num=256))
fwhm = 2 * np.sqrt(2 * np.log(2)) * rr

im = polaron(kk, 1, rr)

# Conserve the spectral weight at each row
weight = np.trapz(y=im, x=kk, axis=1)
im /= weight[:, None]
im /= im.max()

# Visualization in the Brillouin zone
_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters
kx, ky = np.meshgrid(
    np.linspace(-bstar / 2, bstar / 2, num=128),
    np.linspace(-cstar / 2, cstar / 2, num=128),
)
kr = np.hypot(kx, ky)
diff_small = polaron(kr, 1, rr.min())
diff_large = polaron(kr, 1, rr.max())

ax_small.imshow(
    diff_small / diff_small.max(),
    vmin=0,
    vmax=1 / 0.4,
    cmap="inferno",
    extent=[kx.min(), kx.max(), ky.min(), ky.max()],
)
ax_large.imshow(
    diff_large / diff_large.max(),
    vmin=0,
    cmap="inferno",
    extent=[kx.min(), kx.max(), ky.min(), ky.max()],
)

ax_small.set_xlabel("$k_y / b^{*}$")

ax_small.xaxis.set_major_locator(ticker.FixedLocator([kx.min(), 0, kx.max()]))
ax_small.xaxis.set_major_formatter(ticker.FixedFormatter(["-½", "0", "½"]))

for _ax in [ax_small, ax_large]:
    _ax.yaxis.tick_right()
    _ax.yaxis.set_label_position("right")
    _ax.yaxis.set_major_locator(ticker.FixedLocator([ky.min(), 0, ky.max()]))
    _ax.yaxis.set_major_formatter(ticker.FixedFormatter(["-½", "0", "½"]))
    _ax.set_ylabel("$k_z / c^{*}$")

ax_large.xaxis.set_visible(False)

m = ax.imshow(
    im[::-1, :],
    vmin=0,
    vmax=1,
    cmap="inferno",
    aspect="auto",
    extent=[kk.min(), kk.max(), fwhm.min(), fwhm.max()],
)
plt.colorbar(mappable=m, cax=cax, orientation="horizontal")
cax.xaxis.tick_top()
cax.xaxis.set_label_position("top")
cax.set_xlabel("Polaron scattering intensity [a.u.]")

ax.set_ylabel("Polaron FWHM [$\AA$]")
ax.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")

tag_axis(ax, "a)")
tag_axis(ax_large, "b)", x=0.1, y=0.9)
tag_axis(ax_small, "c)", x=0.1, y=0.9)

plt.subplots_adjust(
    top=0.862, bottom=0.161, left=0.098, right=0.897, hspace=0.131, wspace=0.06
)
