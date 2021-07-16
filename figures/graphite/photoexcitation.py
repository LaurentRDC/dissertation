"""
Visualization of the effect of photoexcitation on the electronic structure of graphite

Based on Castro-Neto et al, The electronic properties of graphene, Rev Mod Phys (2009)
"""
from math import sqrt

import matplotlib.colors as cl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import FixedFormatter, FixedLocator

from dissutils import MEDIUM_FIGURE_WIDTH, named_arrow


def colormap():
    """Modify the colormap `base` so the last color is completely transparent"""
    purples = plt.get_cmap("Purples")
    yellows = plt.get_cmap("YlOrBr_r")

    colors = [purples(i) for i in np.linspace(0, 1, 128)] + [
        yellows(i) for i in np.linspace(0, 1, 128)
    ]
    colors[-1] = [0, 0, 0, 0]  # Make end color transparent
    return cl.LinearSegmentedColormap.from_list(name="awesome", colors=colors)


# Hopping parameter
t = 2.7  # eV
photon_energy = 1.56  # eV


def E(kx, ky):
    a = 1.42  # Angstroms. carbon - carbon distance
    fk = 2 * np.cos(np.sqrt(3) * ky * a) + 4 * np.cos(np.sqrt(3) * ky * a / 2) * np.cos(
        3 * kx * a / 2
    )

    E_plus = t * np.sqrt(3 + fk)
    E_minus = -t * np.sqrt(3 + fk)
    return E_plus, E_minus


fig, ax1 = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 2.5))

ky = np.linspace(-1.5, 1.5, num=1024)
Eplus, Eminus = E(kx=1.5, ky=ky)

norm = Normalize(vmin=Eminus.min(), vmax=Eplus.max())
scatter_kwargs = dict(cmap=colormap(), norm=norm, s=2)
ax1.scatter(ky, Eplus, c=2 * Eplus, **scatter_kwargs)  # Scale color to get more range
ax1.scatter(ky, Eminus, c=2 * Eminus, **scatter_kwargs)  # Scale color to get more range
ax1.xaxis.set_major_locator(FixedLocator([-1.7 / 2, 0, 1.7 / 2]))
ax1.xaxis.set_major_formatter(
    FixedFormatter([r"$\mathbf{K}$", r"$\mathbf{M}$", r"$\mathbf{K}$"])
)

# draw vertical transition
x = np.argmin(np.abs(Eminus + photon_energy / 2), axis=0)
arrow_kwds = dict(
    length_includes_head=True,
    width=0.001,
    head_width=0.04,
)

named_arrow(
    ax=ax1,
    x=ky[x],
    y=Eminus[x],
    dx=0,
    dy=photon_energy,
    text=r"$\gamma$",
    color="r",
    fc="r",
    tkwds=dict(ha="left", va="top", color="r", fontsize=16),
    toffset=(0.02, 0),
    **arrow_kwds,
)

named_arrow(
    ax=ax1,
    x=ky[x],
    y=np.abs(Eminus[x]),
    dx=2 * ky[len(ky) - x] - 0.02,
    dy=0,
    text=r"$A_1^\prime$",
    fc="k",
    tkwds=dict(ha="center", va="bottom"),
    toffset=(0, 0.02),
    **arrow_kwds,
)

named_arrow(
    ax=ax1,
    x=ky[x],
    y=np.abs(Eminus[x]),
    dx=2 * (ky[np.argmax(Eminus)] - ky[x]) + 0.02,
    dy=0,
    text=r"$E_{2g}$",
    fc="k",
    tkwds=dict(ha="center", va="bottom"),
    toffset=(0, 0.02),
    **arrow_kwds,
)

ax1.set_xlim([-1.2, 1.2])
ax1.set_ylim([-1.2, 1.2])
ax1.axhline(y=0, linestyle="--", linewidth=0.5, color="k")
ax1.set_ylabel(r"$E(\mathbf{k})$ [eV]")

plt.tight_layout()
