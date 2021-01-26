"""
Time-trace of the transient Debye-Waller effect in graphite.
"""
import numpy as np
import matplotlib.pyplot as plt
from skued import autocenter
from iris import DiffractionDataset
from pathlib import Path
from plotutils import FIGURE_WIDTH

row, col = 1167, 155
DATASET = Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5"

with DiffractionDataset(DATASET, mode="r") as dset:
    time_points = dset.time_points
    intensity = dset.time_series(rect=[row - 10, row + 10, col - 10, col + 10])

intensity /= np.mean(intensity[time_points < 0])

intensity = intensity[time_points < 20]
time_points = time_points[time_points < 20]

fig, ax = plt.subplots(1, 1, figsize=(4, 2))

ax.plot(time_points, intensity, ".k")
ax.axvline(x=0, linestyle="dashed", linewidth=1, color="k")
ax.axhline(y=1, linestyle="dashed", linewidth=1, color="k")

ax.set_xlim([time_points.min(), time_points.max()])
ax.set_ylim([0.9, 1.01])
ax.set_xlabel("Time-delay [ps]")
ax.set_ylabel("$\Delta I$ [a.u.]")
plt.subplots_adjust(bottom=0.01)
plt.tight_layout()