"""
Calculation based on:
Hao Zhang and Dmitri V. Talapin, Thermoelectric Tin Selenide: The Beauty of Simplicity
Angew. Chem. Int. Ed. 2014, 53, 2–4
"""
from dissutils import MEDIUM_FIGURE_WIDTH, FONTSIZE, CBAR_SIZE, tag_axis
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

TC = 300


def eta(ZT, TH):
    return ((TH - TC) / TH) * (np.sqrt(1 + ZT) - 1) / (np.sqrt(1 + ZT) + TC / TH)


zt = np.linspace(0, 4, num=128)
th = np.linspace(800, 300, num=128)

zz, tt = np.meshgrid(zt, th)
im = eta(zz, tt)


figure, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))
m = ax.contourf(
    zz,
    tt,
    im,
    cmap="inferno",
    levels=15,
    extent=[zt.min(), zt.max(), tt.min(), tt.max()],
)
ax.contour(m, colors="dimgray", linestyles="dashed", linewidths=0.4)
ax.set_xlabel("$ZT$ [a.u.]")
ax.set_ylabel("$T_H$ [K]")

cax = make_axes_locatable(ax).append_axes("right", size=CBAR_SIZE, pad=0.03)
plt.colorbar(mappable=m, ax=ax, cax=cax, orientation="vertical")
cax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
cax.set_ylabel("Efficiency")

plt.tight_layout()
