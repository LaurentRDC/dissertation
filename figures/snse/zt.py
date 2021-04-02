import matplotlib.pyplot as plt
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors
from pathlib import Path

DATADIR = Path("data") / "snse"

fig, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))

T, c, b, a = np.loadtxt(DATADIR / "zt.csv", delimiter=",", comments="#", unpack=True)

for color, axis, label, marker in zip(
    discrete_colors(3), [a, b, c], ["a-axis", "b-axis", "c-axis"], ["o", "s", "^"]
):
    ax.errorbar(
        T,
        axis,
        yerr=0.15 * axis,
        label=label,
        ecolor=color,
        color=color,
        linestyle="none",
        marker=marker,
    )

ax.legend(loc="upper left", edgecolor="none")
ax.set_xlabel("Temperature [K]")
ax.set_ylabel("ZT [a.u.]")
ax.set_ylim([-0.1, 3])
ax.axvline(x=810, color="k", linestyle="--", linewidth=1)
ax.text(x=810, y=3.1, s="$T_c$", ha="center")
offset = 50
ax.text(x=810 - offset, y=3.1, s="$\leftarrow Pnma$", ha="right")
ax.text(x=810 + offset, y=3.1, s=r"$Cmcm \rightarrow$", ha="left")

plt.tight_layout()
