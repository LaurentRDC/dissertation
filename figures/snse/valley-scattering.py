import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
from plotutils import tag_axis, discrete_colors

BAND_COLOR, COLOR_TAU1, COLOR_TAU2 = discrete_colors(3)


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


figure, ax = plt.subplots(1, 1, figsize=(3, 2))

k = np.linspace(-1, 1, 256)
ax.plot(k, cb1(k, T=300), color=BAND_COLOR)
ax.plot(k, vb1(k, T=300), color=BAND_COLOR)
ax.plot(k, cb2(k, T=300), color=BAND_COLOR)
ax.plot(k, vb2(k, T=300), color=BAND_COLOR)
ax.plot(k, gm(k, T=300), color=BAND_COLOR)

ax.set_ylim([-1, 1.4])
ax.set_xlim([k.min(), k.max()])
ax.xaxis.set_major_locator(ticker.FixedLocator(locs=[-1, 0, 1]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(["$Z$", "$\Gamma$", "$Y$"]))
ax.yaxis.set_visible(False)

arrow_kwds = dict(
    shrinkA=0, shrinkB=0, mutation_scale=10, clip_on=False, arrowstyle="-|>"
)


ax.add_patch(
    FancyArrowPatch(
        posA=(-0.35, gm(-0.35, 300)),
        posB=(-0.35, cb2(-0.35, 300)),
        ec="r",
        fc="r",
        **arrow_kwds
    )
)
ax.text(
    x=-0.35,
    y=(gm(-0.35, 300) + cb2(-0.35, 300)) / 2,
    s="$\gamma$",
    color="red",
    ha="left",
    va="center",
    transform=ax.transData,
)


# Intravalley relaxation
ax.add_patch(
    FancyArrowPatch(
        posA=(-0.35, cb2(-0.35, 300)),
        posB=(0.2, cb2(0.2, 300)),
        ec=COLOR_TAU1,
        fc=COLOR_TAU1,
        **arrow_kwds
    )
)
ax.text(
    x=0,
    y=(cb2(-0.35, 300) + cb2(0.2, 300)) / 2,
    s="$\\tau_1$",
    color=COLOR_TAU1,
    ha="center",
    va="bottom",
    transform=ax.transData,
)

ax.add_patch(
    FancyArrowPatch(
        posA=(0.2, cb2(0.2, 300)),
        posB=(-0.2, cb2(-0.2, 300)),
        ec=COLOR_TAU1,
        fc=COLOR_TAU1,
        **arrow_kwds
    )
)

ax.add_patch(
    FancyArrowPatch(
        posA=(-0.2, cb2(-0.2, 300)),
        posB=(0, cb2(0, 300)),
        ec=COLOR_TAU1,
        fc=COLOR_TAU1,
        **arrow_kwds
    )
)

# Intervalley relaxation
ax.add_patch(
    FancyArrowPatch(
        posA=(0, cb2(0, 300)),
        posB=(0.66, cb1(0.66, 300)),
        ec=COLOR_TAU2,
        fc=COLOR_TAU2,
        zorder=np.inf,
        **arrow_kwds
    )
)
ax.text(
    x=0.33,
    y=(cb2(0, 300) + cb1(0.66, 300)) / 2,
    s="$\\tau_2$",
    color=COLOR_TAU2,
    ha="center",
    va="bottom",
    transform=ax.transData,
)


ax.set_ylabel("$\epsilon(\mathbf{k})$")

plt.tight_layout()
