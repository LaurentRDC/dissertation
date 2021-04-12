import matplotlib.pyplot as plt
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors
from pathlib import Path

DATADIR = Path("data") / "snse"

fig, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 2.5))

T, c, b, a = np.loadtxt(DATADIR / "zt.csv", delimiter=",", comments="#", unpack=True)
c /= T
b /= T
a /= T

for color, axis, label, marker in zip(
    discrete_colors(3), [a, b, c], ["a-axis", "b-axis", "c-axis"], ["o", "s", "^"]
):
    ax.errorbar(
        T[T > 300],
        axis[T > 300],
        yerr=0.15 * axis[T > 300],
        label=label,
        ecolor=color,
        color=color,
        linestyle="none",
        marker=marker,
    )

ax.legend(loc="upper left", edgecolor="none")
ax.set_xlabel("Temperature [K]")
ax.set_ylabel("$Z$ [1/K]")
ax.set_ylim([0, 1.2 * b.max()])

ax.axvline(x=810, color="k", linestyle="--", linewidth=0.5)
ax.text(x=810, y=1.25 * b.max(), s="$T_c$", ha="center")
offset = 50
ax.text(x=810 - offset, y=1.25 * b.max(), s="$\leftarrow Pnma$", ha="right")
ax.text(x=810 + offset, y=1.25 * b.max(), s=r"$Cmcm \rightarrow$", ha="left")

plt.tight_layout()
