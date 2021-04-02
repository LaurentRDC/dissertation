import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from dissutils import LARGE_FIGURE_WIDTH, discrete_colors
from matplotlib.ticker import FixedLocator, FixedFormatter
from skued import gaussian, with_irf
import numpy as np

TAU = [0.75, 1.7, 2.5]
PUMPTIMES = [0, 3, 6, 9]

PUMPIRF = 0.1  # "IRF"-equivalent of the pump pulse
PROBEIRF = 0.2


def sample_response(times, tau):
    return np.heaviside(times - tau, 1 / 2) * np.exp(-(times - tau) / 0.8)


fig, ((ax1, ax_dyn1), (ax2, ax_dyn2), (ax3, ax_dyn3)) = plt.subplots(
    3,
    2,
    sharex="col",
    sharey="col",
    figsize=(LARGE_FIGURE_WIDTH, 5),
    gridspec_kw=dict(width_ratios=[3, 1]),
)

times = np.linspace(-0.5, 8.84, num=1024)
dyn_times = np.linspace(-0.5, 3, num=28)

pump = np.zeros_like(times)
sample = np.zeros_like(times)
probe1 = np.zeros_like(times)
probe2 = np.zeros_like(times)
probe3 = np.zeros_like(times)
detector = np.zeros_like(times)

for t in PUMPTIMES:
    pump += gaussian(times, center=t, fwhm=PUMPIRF)
    sample += with_irf(PUMPIRF)(sample_response)(times, t)
    probe1 += gaussian(times, center=t + TAU[0], fwhm=PROBEIRF)
    probe2 += gaussian(times, center=t + TAU[1], fwhm=PROBEIRF)
    probe3 += gaussian(times, center=t + TAU[2], fwhm=PROBEIRF)


for i, (ax, ax_dyn, probe, tau) in enumerate(
    zip([ax1, ax2, ax3], [ax_dyn1, ax_dyn2, ax_dyn3], [probe1, probe2, probe3], TAU)
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
    ax.text(x=tau / 2, y=0.95, s=f"$\\tau_{i+1}$", ha="center", va="top")

    ax_dyn.plot(
        times / 3,
        with_irf(PUMPIRF)(sample_response)(times, 0),
        "-",
        color=discrete_colors(3)[0],
    )
    ax_dyn.plot(
        dyn_times / 3,
        with_irf(PROBEIRF)(sample_response)(dyn_times, 0),
        marker="o",
        markersize=7,
        c="dimgray",
        linestyle="none",
        clip_on=False,
        markerfacecolor="none",
    )
    x = np.asarray([tau / 3])
    y = sample_response(np.asarray([tau]), 0)
    ax_dyn.plot(
        x,
        y,
        marker="o",
        markersize=8,
        c=discrete_colors(3)[-1],
        linestyle="none",
    )
    ax_dyn.text(
        x=x,
        y=y + 0.1,
        s=f"$\\tau_{i+1}$",
        ha="left",
        va="bottom",
        transform=ax_dyn.transData,
    )

    ax.set_ylim([-0.05, 1.05])
    ax_dyn.set_ylim([-0.05, 1.05])

    ax.yaxis.set_visible(False)
    ax_dyn.yaxis.set_visible(False)

    ax.set_frame_on(False)
    ax_dyn.set_frame_on(False)

ax3.set_xlim([times.min(), times.max()])
ax3.set_xlabel("Laboratory time [$T$]")
ax3.xaxis.set_major_locator(FixedLocator([0, 3, 6, 9]))
ax3.xaxis.set_major_formatter(FixedFormatter(["0", "1", "2", "3"]))

ax_dyn3.set_xlim([dyn_times.min() / 3, dyn_times.max() / 3])
ax_dyn3.set_xlabel("Time-delay [T]")

ax1.legend(ncol=3, loc="center", bbox_to_anchor=(0.5, 1.15), edgecolor="none")

plt.subplots_adjust(
    top=0.911, bottom=0.125, left=0.041, right=0.971, hspace=0.142, wspace=0.159
)
