"""
Calculation based on:
Hao Zhang and Dmitri V. Talapin, Thermoelectric Tin Selenide: The Beauty of Simplicity
Angew. Chem. Int. Ed. 2014, 53, 2–4
"""
from plotutils import FIGURE_WIDTH, FONTSIZE, CBAR_SIZE, tag_axis, discrete_colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

TC = 300


def eta(ZT, TH):
    return ((TH - TC) / TH) * (np.sqrt(1 + ZT) - 1) / (np.sqrt(1 + ZT) + TC / TH)


zt = np.linspace(0, 3, num=128)
th = np.linspace(800, 300, num=128)

zz, tt = np.meshgrid(zt, th)
im = eta(zz, tt)


figure, ax = plt.subplots(1, 1, figsize=(4, 2))
m = ax.imshow(
    im, cmap="hot", aspect="auto", extent=[zt.min(), zt.max(), tt.min(), tt.max()]
)
ax.set_xlabel("$ZT$ [a.u.]")
ax.set_ylabel("$T_H$ [K]")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size=CBAR_SIZE, pad=0.03)
plt.colorbar(mappable=m, ax=ax, cax=cax, orientation="vertical")
cax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
cax.set_ylabel("Efficiency $\eta$")

plt.tight_layout()
