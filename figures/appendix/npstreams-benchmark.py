from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis

NS_COLOR, NP_COLOR = discrete_colors(2)

INPUT = Path("data") / "introduction" / "npstreams-benchmark"

sizes, np_times, ns_times = np.hsplit(np.load(INPUT / "times.npy"), 3)
_, np_memories, ns_memories = np.hsplit(np.load(INPUT / "memory.npy"), 3)
seqsizes, np_seq_times, ns_seq_times = np.hsplit(np.load(INPUT / "seqtimes.npy"), 3)

figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2))
gs = gridspec.GridSpec(nrows=2, ncols=2, figure=figure)
ax1 = figure.add_subplot(gs[0, 0])
ax2 = figure.add_subplot(gs[1, 0], sharex=ax1)
ax3 = figure.add_subplot(gs[:, 1])

for ax in [ax1, ax2]:
    ax.axvline(x=2048, linestyle="dashed", color="k", linewidth=0.5)

for label, trace_t, trace_m, marker, color in zip(
    ["numpy", "npstreams"],
    [np_times, ns_times],
    [np_memories, ns_memories],
    ["x", "^"],
    [NS_COLOR, NP_COLOR],
):

    ax1.plot(
        sizes,
        trace_t * 1e3,  # times in [ms]
        marker=marker,
        linestyle="None",
        color=color,
        label=label,
        markersize=5,
    )
    ax2.plot(
        sizes,
        trace_m,
        color=color,
        marker=marker,
        linestyle="None",
        markersize=5,
    )

ax1.set_ylabel("Time [ms]")
ax1.set_yscale("log")
ax1.legend(
    ncol=2,
    loc="center",
    edgecolor="none",
    bbox_to_anchor=(0.5, 1.15),
    prop=FontProperties(family="monospace"),
)
ax1.xaxis.set_visible(False)


ax2.set_yscale("log")
ax2.set_ylabel("Memory [MB]")

ax2.set_xscale("log", base=2)
ax2.set_xticks([2 ** 9, 2 ** 10, 2 ** 11, 2 ** 12, 2 ** 13])
ax2.set_xlabel("Array size ($n \\times n$)")

ax3.plot(
    seqsizes,
    np_seq_times / ns_seq_times,
    marker="o",
    linestyle="None",
    color=NS_COLOR,
    markersize=3,
)
ax3.set_xlabel("Array sequence length")
ax3.set_ylabel("Speed-up factor [a.u.]")
ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()

tag_axis(ax1, "a)", y=0.9)
tag_axis(ax2, "b)", y=0.9)
tag_axis(ax3, "c)")

plt.subplots_adjust(
    top=0.866, bottom=0.198, left=0.103, right=0.902, hspace=0.06, wspace=0.06
)
