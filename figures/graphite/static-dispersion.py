from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.constants import physical_constants

from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, draw_hexagon

HBAR = physical_constants["Planck constant over 2 pi in eV s"][0]
HZ_TO_EV = HBAR * 2 * np.pi

modes = ["LO2", "TA", "TO2", "LA"]

fig, ax = plt.subplots(1, 1, figsize=(LARGE_FIGURE_WIDTH, 3))
ax_thz = ax.twinx()

for mode, color in zip(modes, discrete_colors(len(modes))):
    f = np.load(Path("data") / "graphite" / "static-dispersion" / f"{mode}.npy")[::4]
    # Artifact near K
    k_index = 2 * np.size(f) // 3
    near_k = slice(k_index - 7, k_index + 7)
    f[near_k] = np.clip(
        f[near_k],
        a_min=f[near_k][np.nonzero(f[near_k])].min(),
        a_max=np.inf,
    )
    ax.plot(range(np.size(f)), f * HZ_TO_EV * 1000, color=color)  # meV
    ax_thz.plot(range(np.size(f)), f * 1e-12, color=color)  # THz,

labels = [r"$\mathbf{\Gamma}$", r"$\mathbf{M}$", r"$\mathbf{K}$", r"$\mathbf{\Gamma}$"]
nsteps = np.size(f) / 3
for i in range(1, len(labels)):
    ax.axvline(x=i * nsteps, color="k", linestyle="--", linewidth=0.5)

ax.axhline(y=25, color="k", linestyle="--", linewidth=0.5)  # room temperature

# Draw BZ
hex_radius = 0.4
w, h = fig.get_size_inches()
path_ax = fig.add_axes([0.45, 0.10, 0.15, (w / h) * 0.15])
path_ax.axis("off")
path_ax.set_xlim([-0.5, 0.5])
path_ax.set_ylim([-0.5, 0.5])
draw_hexagon(path_ax, center=(0, 0), radius=hex_radius, color="k", facecolor="w")
for (x, y), label, ha in zip(
    [
        (0, 0),
        (0, hex_radius * sqrt(3) / 2),
        (hex_radius / 2, hex_radius * sqrt(3) / 2),
    ],
    [r"$\mathbf{\Gamma}$", r"$\mathbf{M}$", r"$\mathbf{K}$"],
    ["right", "right", "left"],
):
    path_ax.scatter(x, y, s=5, c="k", zorder=np.inf)
    path_ax.annotate(label, (x, y), va="bottom", ha=ha)
    path_ax.arrow(
        x=0,
        y=0,
        dx=x,
        dy=y,
        linewidth=1,
        color="k",
        length_includes_head=True,
    )

# Annotate the modes
for label, pt in zip(
    ["TA", "LA", "TO", "LO"],
    [(0.14, 0.15), (0.07, 0.35), (0.1, 0.8), (0.25, 0.94)],
):
    ax.annotate(label, pt, color="k", textcoords=ax.transAxes)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(
    ticker.FixedLocator([i * nsteps for i in range(len(labels))])
)

ax.xaxis.set_major_formatter(plt.FixedFormatter(labels))
ax.set_xlim([0, np.size(f)])
ax.set_ylabel("Energy [meV]")
ax.set_ylim([0, 210])

ax_thz.set_ylabel("Frequency [THz]")
ax_thz.yaxis.set_label_position("right")
ax_thz.set_ylim(list(map(lambda i: i * 1e-12 / HZ_TO_EV / 1000, ax.get_ylim())))
