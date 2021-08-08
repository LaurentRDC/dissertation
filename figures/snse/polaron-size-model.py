"""
Visualization of the effect of polaron localization (from large to small)
on the diffuse intensity profile
"""

from math import sqrt, log
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter
import skued
from crystals import Crystal
from mpl_toolkits.axes_grid1 import make_axes_locatable

from dissutils import (
    CBAR_SIZE,
    GRID_AXES_PAD,
    LARGE_FIGURE_WIDTH,
    discrete_colors,
    tag_axis,
)

DATADIR = Path("data") / "snse"
CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")


def polaron(q, A, rp):
    return A * q * rp ** 2 * np.exp(-(q ** 2 * rp ** 2) / 2)


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 6))
gs = gridspec.GridSpec(3, 2, width_ratios=[2.5, 1], height_ratios=[1, 1, 1.5])
ax = figure.add_subplot(gs[0:2, 0])
ax_small = figure.add_subplot(gs[1, 1])
ax_large = figure.add_subplot(gs[0, 1])
ax_fits = figure.add_subplot(gs[2, 0])

cax = make_axes_locatable(ax).append_axes("top", size=CBAR_SIZE, pad=GRID_AXES_PAD)

ks_ = np.linspace(0, 1, num=1024)  # Gamma to T distance is 1 inv Angs
kk, rr = np.meshgrid(ks_, np.linspace(1, 7, num=256))
fwhm = 2 * np.sqrt(2 * np.log(2)) * rr

im = polaron(kk, 1, rr)

# Conserve the spectral weight at each row
weight = np.trapz(y=im, x=kk, axis=1)
im /= weight[:, None]
im /= im.max()

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

# Fits ------------------------------------------------------------------------


rf, pf, ef = np.loadtxt(
    DATADIR / "fast-diffuse-profile.csv", delimiter=",", unpack=True
)
pf -= pf.min()
ef /= pf.max()
pf /= pf.max()

rs, ps, es = np.loadtxt(
    DATADIR / "slow-diffuse-profile.csv", delimiter=",", unpack=True
)
ps -= ps.min()
es /= ps.max()
ps /= ps.max()


paramsf, _ = opt.curve_fit(
    polaron,
    rf[rf < 0.32],
    pf[rf < 0.32],
    sigma=ef[rf < 0.32],
    p0=(8e-3, 12),
    bounds=([0, 0], [np.inf, np.inf]),
)

paramss, _ = opt.curve_fit(
    polaron,
    rs,
    ps,
    sigma=es,
    bounds=([-np.inf, 0], [np.inf, np.inf]),
)

# Reversing the order so that the legends are closer to their
# associated data points
for (r, p, e), params, marker, color, label in zip(
    [(rs, ps, es), (rf, pf, ef)],
    [paramss, paramsf],
    ["^", "o"],
    reversed(discrete_colors(2)),
    ["$\\tau=5$ ps", "$\\tau=1$ ps"],
):
    p_smooth = gaussian_filter(p, sigma=3)
    ax_fits.plot(
        r,
        p_smooth,
        marker=marker,
        markerfacecolor="none",
        color=color,
        markersize=6,
        linestyle="None",
        label=label,
    )

    ax_fits.fill_between(
        r,
        y1=p_smooth + e,
        y2=p_smooth - e,
        facecolor=color,
        edgecolor="none",
        linewidth=1,
        alpha=0.2,
    )

    # Show on main plot
    _, rp, *_ = params
    ax.axhline(
        y=2 * sqrt(2 * log(2)) * rp, linestyle="dashed", color="gray", linewidth=1
    )

    r_ = np.linspace(rf.min(), rs.max(), 1024)
    ax_fits.plot(r_, polaron(r_, *params), color=color, linestyle="-")

    ax_fits.set_ylabel("$\Delta I / I_0$ [a.u.]")
    ax_fits.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")

ax_fits.set_xlim([rf.min(), max(ax.get_xlim())])
ax.set_xlim(ax_fits.get_xlim())

# Reciprocal space view -------------------------------------------------------

_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters
kx, ky = np.meshgrid(
    np.linspace(-bstar / 2, bstar / 2, num=128),
    np.linspace(-cstar / 2, cstar / 2, num=128),
)
kr = np.hypot(kx, ky)
diff_small = polaron(kr, 1, paramss[1])
diff_large = polaron(kr, 1, paramsf[1])

ax_small.imshow(
    diff_small / diff_small.max(),
    vmin=0,
    vmax=1,
    cmap="inferno",
    extent=[kx.min(), kx.max(), ky.min(), ky.max()],
)
ax_large.imshow(
    diff_large / diff_large.max(),
    vmin=0,
    vmax=1,
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

ax_fits.legend(
    loc="center", ncol=1, bbox_to_anchor=(0.8, 0.5), edgecolor="none", facecolor="none"
)

tag_axis(ax, "a)", x=0.025, y=0.9625)
tag_axis(ax_fits, "b)", x=0.025)
tag_axis(ax_large, "c)", x=0.1, y=0.9)
tag_axis(ax_small, "d)", x=0.1, y=0.9)


plt.subplots_adjust(
    top=0.91, bottom=0.1, left=0.1, right=0.9, hspace=0.195, wspace=0.105
)
