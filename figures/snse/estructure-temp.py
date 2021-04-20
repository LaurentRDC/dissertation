import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import matplotlib.colors as colors
from dissutils import MEDIUM_FIGURE_WIDTH


def cb1(k, T):
    t = T / 600 - 1
    kc = 0.66 + 0.05 * t
    return 15 * (k - kc) ** 2 + 0.6 + t / 2


def cb2(k, T):
    t = T / 600 - 1
    return 5 * k ** 2 + 0.6 - t / 5


def vb1(k, T):
    t = T / 600 - 1
    kc = 0.66 + 0.05 * t
    return -20 * (k - kc) ** 2 + t / 4 - 0.35


def vb2(k, T):
    t = T / 600
    kc = -0.75

    return -1000 * (k - kc) ** 4 + (15 / t) * (k - kc) ** 2 + t / 10 - 0.4


def gm(k, T):
    t = T / 600
    kc = 0
    return -0.9 * (k - kc) ** 2 - 0.9 + t / 15


figure, ax = plt.subplots(
    1,
    1,
    figsize=(MEDIUM_FIGURE_WIDTH, 2),
)
k = np.linspace(-1, 1, 256)
cmap = plt.get_cmap("hot")
for T in range(800, 300, -25):
    c = cmap((T - 300) / 700)
    ax.plot(k, cb1(k, T), color=c)
    ax.plot(k, vb1(k, T), color=c)
    ax.plot(k, cb2(k, T), color=c)
    ax.plot(k, vb2(k, T), color=c)
    ax.plot(k, gm(k, T), color=c)

ax.axhline(y=0.6, xmin=0.40, xmax=0.92, linestyle="dashed", color="k", linewidth=0.5)
ax.set_ylim([-1, 1.25])
ax.set_xlim([k.min(), k.max()])
ax.xaxis.set_major_locator(ticker.FixedLocator(locs=[-1, 0, 1]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(["$Z$", "$\Gamma$", "$Y$"]))
ax.set_yticks([])
ax.set_ylim([-1, 1.25])

cax = ax.inset_axes([1.02, 0.0, 0.03, 1])
cb = mpl.colorbar.ColorbarBase(
    cax,
    cmap=cmap,
    orientation="vertical",
    ticklocation="right",
    ticks=[300, 800],
    format=ticker.FixedFormatter(["300K", "800K"]),
    norm=mpl.colors.Normalize(vmin=300, vmax=800),
)

ax.axhline(y=0.6, xmin=0.40, xmax=0.92, linestyle="dashed", color="k", linewidth=0.5)

plt.tight_layout()
