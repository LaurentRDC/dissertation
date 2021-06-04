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
import scipy.optimize as opt
from crystals import Crystal
from iris import DiffractionDataset
from dissutils import (
    LARGE_FIGURE_WIDTH,
    CBAR_SIZE,
    GRID_AXES_PAD,
    discrete_colors,
    tag_axis,
)
from dissutils.snse import overnight4
from skimage.filters import gaussian

DATADIR = Path("data") / "snse"
CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")

INDICES_DIFFUSE = [
    (0, -1, 3),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -1, 7),
    (0, -3, 5),
    (0, -5, 3),
]

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, 0, *args, **kwargs)


def polaron(q, A, rp):
    return A * q * rp ** 3 / (1 + (q * rp) ** 2) ** 2


ks, amplitudes, amplitudes_err = np.loadtxt(
    DATADIR / "fast-diffuse.csv", delimiter=",", unpack=True
)

params, pcov = opt.curve_fit(
    polaron,
    ks,
    amplitudes,
    sigma=amplitudes_err,
    bounds=([-np.inf, 0], [np.inf, np.inf]),
)

figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 4))
gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1])
ax = figure.add_subplot(gs[:, 0])
ax_small = figure.add_subplot(gs[1, 1])
ax_large = figure.add_subplot(gs[0, 1])

cax = make_axes_locatable(ax).append_axes("top", size=CBAR_SIZE, pad=GRID_AXES_PAD)

ks_ = np.linspace(0, 1, num=1024)  # Gamma to T distance is 1 inv Angs
kk, rr = np.meshgrid(ks_, np.linspace(0.5, 4, num=256))
fwhm = 2 * np.sqrt(2 * np.log(2)) * rr

im = polaron(kk, params[0], rr)

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
cax.set_xlabel("Diffuse intensity [a.u.]")

ax.set_ylabel("Polaron FWHM [$\AA$]")
ax.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")

tag_axis(ax, "a)")
tag_axis(ax_large, "b)", x=0.1, y=0.9)
tag_axis(ax_small, "c)", x=0.1, y=0.9)

plt.tight_layout()
