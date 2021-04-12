import matplotlib.pyplot as plt
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis
from pathlib import Path

DATADIR = Path("data") / "snse"

fig, (ax_S, ax_klat) = plt.subplots(2, 1, figsize=(MEDIUM_FIGURE_WIDTH, 4), sharex=True)
ax_S.set_ylabel("$s$ [$\mu$V/K]")
ax_S.xaxis.set_visible(False)
ax_klat.set_ylabel("$\kappa_l$ [W/m/K]")

T, S_c, S_b, S_a = np.loadtxt(
    DATADIR / "seebeck.csv", delimiter=",", comments="#", unpack=True
)
for color, trace, label, marker in zip(
    discrete_colors(3), [S_a, S_b, S_c], ["a-axis", "b-axis", "c-axis"], ["o", "s", "^"]
):
    ax_S.plot(
        T,
        trace,
        label=label,
        color=color,
        linestyle="none",
        marker=marker,
    )
ax_S.legend(ncol=3, loc="center", bbox_to_anchor=(0.5, 1.1), edgecolor="none")

T, k_c, k_b, k_a = np.loadtxt(
    DATADIR / "klat.csv", delimiter=",", comments="#", unpack=True
)
for color, trace, label, marker in zip(
    discrete_colors(3), [k_a, k_b, k_c], ["a-axis", "b-axis", "c-axis"], ["o", "s", "^"]
):
    ax_klat.plot(
        T,
        trace,
        label=label,
        color=color,
        linestyle="none",
        marker=marker,
    )

ax_klat.set_xlabel("Temperature [K]")

for ax, label in zip([ax_S, ax_klat], ["a)", "b)"]):
    ax.axvline(x=810, color="k", linestyle="--", linewidth=0.5)
    tag_axis(ax, label, y=0.12, verticalalignment="bottom")

ax_klat.axvline(x=810, color="k", linestyle="--", linewidth=0.5)
ax_klat.text(x=810, y=1.15 * k_c.max(), s="$T_c$", ha="center")
offset = 50
ax_klat.text(x=810 - offset, y=1.15 * k_c.max(), s="$\leftarrow Pnma$", ha="right")
ax_klat.text(x=810 + offset, y=1.15 * k_c.max(), s=r"$Cmcm \rightarrow$", ha="left")

plt.subplots_adjust(
    top=0.896, bottom=0.146, left=0.146, right=0.931, hspace=0.269, wspace=0.2
)
