from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis

COLOR = discrete_colors(1)[0]

frequencies_fullbw = np.linspace(50e6, 6e9, num=10401) / 1e9  # [GHz]
power_fullbw = np.loadtxt(
    Path("data") / "introduction" / "full_bw_comb_PD_output.csv", skiprows=108
)  # [dBm]

frequencies_3ghz = (
    np.linspace(2995999807.69231, 3000999807.69231, num=10401) / 1e9
)  # [GHz]
power_3ghz = np.loadtxt(
    Path("data") / "introduction" / "5MHz_span_comb_PD_output.csv", skiprows=108
)  # [dBm]

figure, (ax_full, ax_3ghz) = plt.subplots(
    1, 2, figsize=(LARGE_FIGURE_WIDTH, 3), gridspec_kw=dict(width_ratios=[2, 1])
)

for ax in (ax_full, ax_3ghz):
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("Spectral power [dBm]")

ax_full.axvline(x=2.9985, linestyle="dashed", color="k", linewidth=0.5)

ax_full.plot(frequencies_fullbw, power_fullbw, color=COLOR)
ax_full.set_ylim([-75, -15])
ax_full.set_xlim([-0.01, frequencies_fullbw.max()])
ax_full.text(
    x=0.5, y=1.01, ha="center", va="bottom", s="$f_{40}$", transform=ax_full.transAxes
)

maxf = frequencies_3ghz[np.argmax(power_3ghz)]
ax_3ghz.plot(frequencies_3ghz, power_3ghz, color=COLOR)
ax_3ghz.text(
    x=0.5,
    y=1.01,
    ha="center",
    va="bottom",
    s=f"$f_{{40}} = {maxf:.5f}$ GHz",
    transform=ax_3ghz.transAxes,
)

ax_3ghz.set_xlim([maxf - 0.001, maxf + 0.001])
ax_3ghz.yaxis.tick_right()
ax_3ghz.yaxis.set_label_position("right")

tag_axis(ax_full, "a)", x=0.97, y=0.925, horizontalalignment="right")
tag_axis(ax_3ghz, "b)", x=0.075, y=0.925)

plt.tight_layout()
