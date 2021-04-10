import numpy as np
import matplotlib.pyplot as plt
from iris import DiffractionDataset
from pathlib import Path
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors

row, col = 1167, 155
DATASET = Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"

with DiffractionDataset(DATASET, mode="r") as dset:
    time_points = dset.time_points
    intensity = dset.time_series(rect=[row - 10, row + 10, col - 10, col + 10])

intensity /= np.mean(intensity[time_points < 0])

intensity = intensity[time_points < 20]
time_points = time_points[time_points < 20]

fig, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 2))

ax.plot(time_points, intensity, ".", color=discrete_colors(1)[0])
ax.axvline(x=0, linestyle="dashed", linewidth=0.5, color="k")
ax.axhline(y=1, linestyle="dashed", linewidth=0.5, color="k")

ax.set_xlim([time_points.min(), time_points.max()])
ax.set_ylim([0.9, 1.01])
ax.set_xlabel("Time-delay [ps]")
ax.set_ylabel("$\Delta I$ [a.u.]")
plt.subplots_adjust(bottom=0.01)
plt.tight_layout()
