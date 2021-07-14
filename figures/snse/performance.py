import matplotlib.pyplot as plt
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis
from pathlib import Path

DATADIR = Path("data") / "snse"

fig, (ax_ZT, ax_S, ax_klat) = plt.subplots(
    3, 1, figsize=(MEDIUM_FIGURE_WIDTH, 4.5), sharex=True
)

for ax, filepath, ylabel in zip(
    [ax_ZT, ax_S, ax_klat],
    [DATADIR / "zt.csv", DATADIR / "seebeck.csv", DATADIR / "klat.csv"],
    ["$ZT$ [a.u.]", "$S$ [$\mu$V/K]", "$\kappa_l$ [W/m/K]"],
):
    T, c, b, a = np.loadtxt(filepath, delimiter=",", comments="#", unpack=True)
    ax.set_ylabel(ylabel)
    for color, trace, label, marker in zip(
        discrete_colors(3), [a, b, c], ["a-axis", "b-axis", "c-axis"], ["o", "s", "^"]
    ):
        ax.plot(
            T,
            trace,
            label=label,
            color=color,
            linestyle="none",
            marker=marker,
        )

ax_ZT.legend(edgecolor="none", loc="center", bbox_to_anchor=(0.3, 0.6))
ax_klat.set_xlabel("Temperature [K]")

for ax, label, (y, va) in zip(
    [ax_ZT, ax_S, ax_klat],
    ["a)", "b)", "c)"],
    [(0.88, "top"), (0.12, "bottom"), (0.12, "bottom")],
):
    ax.axvline(x=810, color="k", linestyle="--", linewidth=0.5)
    tag_axis(ax, label, y=y, verticalalignment=va)

# Phase transition temperature
_, zt_max = ax_ZT.get_ylim()
ax_ZT.text(x=810, y=1.05 * zt_max, s="$T_c$", ha="center")
offset = 50
ax_ZT.text(x=810 - offset, y=1.05 * zt_max, s="$\leftarrow Pnma$", ha="right")
ax_ZT.text(x=810 + offset, y=1.05 * zt_max, s=r"$Cmcm \rightarrow$", ha="left")

plt.tight_layout()
