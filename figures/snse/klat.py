import matplotlib.pyplot as plt
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors
from pathlib import Path

DATADIR = Path("data") / "snse"

fig, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))

T, c, b, a = np.loadtxt(DATADIR / "klat.csv", delimiter=",", comments="#", unpack=True)

for color, axis, label, marker in zip(
    discrete_colors(3), [a, b, c], ["a-axis", "b-axis", "c-axis"], ["o", "s", "^"]
):
    ax.plot(
        T,
        axis,
        label=label,
        color=color,
        linestyle="none",
        marker=marker,
    )

ax.legend(edgecolor="none")
ax.set_xlabel("Temperature [K]")
ax.set_ylabel("$\kappa_l$ [W / m / K]")

ax.axvline(x=810, color="k", linestyle="--", linewidth=1)
ax.text(x=810, y=1.1 * c.max(), s="$T_c$", ha="center")
offset = 50
ax.text(x=810 - offset, y=1.1 * c.max(), s="$\leftarrow Pnma$", ha="right")
ax.text(x=810 + offset, y=1.1 * c.max(), s=r"$Cmcm \rightarrow$", ha="left")

plt.tight_layout()
