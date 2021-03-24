import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from plotutils import FIGURE_WIDTH, discrete_colors
from matplotlib.ticker import FixedLocator, FixedFormatter
from skued import gaussian
import numpy as np


fig, ax = plt.subplots(1, 1, sharex=True, figsize=(FIGURE_WIDTH, 2.5))

times = np.linspace(-0.25, 8.84, num=1024)
pump = np.zeros_like(times)
sample = np.zeros_like(times)
probe = np.zeros_like(times)
detector = np.zeros_like(times)

for t in [0, 3, 6, 9]:
    pump += gaussian(times, center=t, fwhm=0.1)
    sample += np.heaviside(times - t, 1 / 2) * np.exp(-(times - t) / 0.8)
    probe += gaussian(times, center=t + 1.5, fwhm=0.2)

for trace, color, label in zip(
    [sample, pump, probe],
    discrete_colors(3),
    ["Sample response", "Pump", "Probe"],
):

    ax.fill_between(times, y1=trace / trace.max(), y2=0, color=color, label=label)

ax.add_patch(
    FancyArrowPatch(
        posA=(0, 1.05),
        posB=(1.5, 1.05),
        arrowstyle="|-|",
        shrinkA=0,
        shrinkB=0,
        mutation_scale=2,
        clip_on=False,
    )
)
ax.text(x=0.75, y=1.1, s="$\\tau$", ha="center")

ax.set_ylim([0, 1.05])
ax.yaxis.set_visible(False)

ax.set_frame_on(False)
ax.set_xlim([times.min(), times.max()])
ax.set_xlabel("Laboratory time [$T$]")
ax.xaxis.set_major_locator(FixedLocator([0, 3, 6, 9]))
ax.xaxis.set_major_formatter(FixedFormatter(["0", "1", "2", "3"]))

ax.legend(ncol=3, loc="center", bbox_to_anchor=(0.5, 1.15), edgecolor="none")

plt.tight_layout()
