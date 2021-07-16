from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from skued import baseline_dt, gaussian

from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis

s, intensity = np.load(Path("data") / "appendix" / "rutile_vo2.npy")
q = 4 * np.pi * s

# Double exponential inelastic background and substrate effects
diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
substrate1 = 0.8 * gaussian(s, center=s.mean(), fwhm=s.mean() / 4)
substrate2 = 0.9 * gaussian(s, center=s.mean() / 2.5, fwhm=s.mean() / 4)

signal = intensity + diffuse + substrate1 + substrate2

fig, (ax1, ax2) = plt.subplots(
    nrows=2,
    ncols=1,
    sharex=True,
    figsize=(MEDIUM_FIGURE_WIDTH, MEDIUM_FIGURE_WIDTH),
)

best_baseline = baseline_dt(signal, level=7, max_iter=150, wavelet="qshift3")
for trace, color, label in zip(
    [signal, signal - intensity, best_baseline],
    discrete_colors(3),
    ["signal", "background", "baseline"],
):
    ax1.plot(q, trace, "-", color=color, label=label)

ax1.set_ylim([28, 73])
ax1.legend(edgecolor="none", ncol=3, loc="center", bbox_to_anchor=(0.5, 1.1))

c1, c2 = discrete_colors(2)
ax2.plot(q, intensity, "-", color=c1, label="True intensity")
ax2.plot(
    q,
    signal - baseline_dt(signal, level=7, max_iter=150, wavelet="qshift3"),
    color=c2,
    label="Background-subtracted",
)
ax2.set_xlabel("$|\mathbf{q}|$ ($\AA^{-1}$)")
ax2.set_xlim([0.2 * 4 * np.pi, 0.4 * 4 * np.pi])
ax2.legend(
    edgecolor="none", ncol=2, loc="center", bbox_to_anchor=(0.5, 1.1), framealpha=0
)

for ax, label in zip([ax1, ax2], ["a)", "b)"]):
    ax.set_ylabel("Intensity [a.u.]")
    tag_axis(ax, label, x=0.96, horizontalalignment="right", y=0.9)

plt.tight_layout()
