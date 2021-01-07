from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from plotutils import FIGURE_WIDTH, draw_hexagon
from scipy.constants import physical_constants

HBAR = physical_constants["Planck constant over 2 pi in eV s"][0]
HZ_TO_EV = HBAR * 2 * np.pi

modes = ["LA", "LO2", "TA", "TO2"]

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, 3))
ax_thz = ax.twinx()

for mode in modes:
    f = np.load(Path("data") / "graphite" / "static-dispersion" / f"{mode}.npy")
    ax.scatter(x=range(np.size(f)), y=f * HZ_TO_EV * 1000, s=2)  # meV
    ax_thz.scatter(x=range(np.size(f)), y=f * 1e-12, s=0)  # THz,

labels = [r"$\mathbf{\Gamma}$", r"$\mathbf{M}$", r"$\mathbf{K}$", r"$\mathbf{\Gamma}$"]
nsteps = np.size(f) / 3
for i in range(1, len(labels)):
    ax.axvline(x=i * nsteps, color="k", linestyle="--", linewidth=1)

ax.axhline(y=25, color="k", linestyle="--", linewidth=1)  # room temperature

# Draw BZ
# draw_hexagon(ax, center=(0,0), radius=100, color='k', transform=ax.transAxes)
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
