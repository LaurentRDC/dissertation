import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from plotutils import FIGURE_WIDTH, discrete_colors
from matplotlib.ticker import FixedLocator, FixedFormatter
from skued import gaussian
import numpy as np

TAU = [0.75, 1.7, 2.5]
PUMPTIMES = [0, 3, 6, 9]

fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, sharex=True, sharey=True, figsize=(FIGURE_WIDTH, 5)
)

times = np.linspace(-0.25, 8.84, num=1024)
pump = np.zeros_like(times)
sample = np.zeros_like(times)
probe1 = np.zeros_like(times)
probe2 = np.zeros_like(times)
probe3 = np.zeros_like(times)
detector = np.zeros_like(times)

for t in PUMPTIMES:
    pump += gaussian(times, center=t, fwhm=0.1)
    sample += np.heaviside(times - t, 1 / 2) * np.exp(-(times - t) / 0.8)
    probe1 += gaussian(times, center=t + TAU[0], fwhm=0.2)
    probe2 += gaussian(times, center=t + TAU[1], fwhm=0.2)
    probe3 += gaussian(times, center=t + TAU[2], fwhm=0.2)


for i, (ax, probe, tau) in enumerate(
    zip([ax1, ax2, ax3], [probe1, probe2, probe3], TAU)
):
    for trace, color, label in zip(
        [sample, pump, probe],
        discrete_colors(3),
        ["Sample response", "Pump", "Probe"],
    ):

        ax.fill_between(times, y1=trace / trace.max(), y2=0, color=color, label=label)

    ax.add_patch(
        FancyArrowPatch(
            posA=(0, 1.05),
            posB=(tau, 1.05),
            arrowstyle="|-|",
            shrinkA=0,
            shrinkB=0,
            mutation_scale=2,
            clip_on=False,
        )
    )
    ax.text(x=tau / 2, y=1.1, s=f"$\\tau_{i+1}$", ha="center")

    ax.set_ylim([0, 1.05])
    ax.yaxis.set_visible(False)

    ax.set_frame_on(False)

ax3.set_xlim([times.min(), times.max()])
ax3.set_xlabel("Laboratory time [$T$]")
ax3.xaxis.set_major_locator(FixedLocator([0, 3, 6, 9]))
ax3.xaxis.set_major_formatter(FixedFormatter(["0", "1", "2", "3"]))

ax1.legend(ncol=3, loc="center", bbox_to_anchor=(0.5, 1.15), edgecolor="none")

plt.tight_layout()
