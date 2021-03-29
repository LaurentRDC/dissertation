import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from plotutils import discrete_colors, tag_axis, LARGE_FIGURE_WIDTH

NS_COLOR, NP_COLOR = discrete_colors(2)

INPUT = Path("data") / "introduction" / "npstreams-benchmark"

sizes, np_times, ns_times = np.hsplit(np.load(INPUT / "times.npy"), 3)
_, np_memories, ns_memories = np.hsplit(np.load(INPUT / "memory.npy"), 3)
seqsizes, np_seq_times, ns_seq_times = np.hsplit(np.load(INPUT / "seqtimes.npy"), 3)

fig, (ax1, ax3) = plt.subplots(
    1, 2, figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2)
)
ax2 = ax1.twinx()

ax1.axvline(x=2048, linestyle="dashed", color="k", linewidth=1)
ax1.plot(
    sizes,
    np_times,
    marker="x",
    linestyle="None",
    color=NP_COLOR,
    label="numpy",
    markersize=5,
)
ax1.plot(
    sizes,
    ns_times,
    marker="^",
    linestyle="None",
    color=NS_COLOR,
    label="npstreams",
    markersize=5,
)
ax1.set_ylabel("Wall time [s]")
ax1.set_xlabel("Array size ($n \\times n$)")
ax1.set_yscale("log")
ax1.set_xscale("log", base=2)
ax1.set_xticks([2 ** 9, 2 ** 10, 2 ** 11, 2 ** 12, 2 ** 13])
ax1.legend(ncol=2, loc="center", edgecolor="none", bbox_to_anchor=(0.5, 1.1))

ax2.plot(sizes, np_memories, color=NP_COLOR, linestyle="--")
ax2.plot(sizes, ns_memories, color=NS_COLOR, linestyle="--")
ax2.set_yscale("log")
ax2.set_ylabel("Memory usage [MB]")

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

tag_axis(ax1, "a)")
tag_axis(ax3, "b)")

plt.tight_layout()
