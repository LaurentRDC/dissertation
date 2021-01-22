# -*- coding: utf-8 -*-
from pathlib import Path

import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from skued import biexponential, with_irf
from scipy.optimize import curve_fit
import numpy as np
from plotutils import FIGURE_WIDTH

INPUT = Path("data") / "graphite" / "populations"


class PopulationDatabase(h5.File):
    @property
    def modes(self):
        return [
            name.strip("_frequencies")
            for name in self.keys()
            if name.endswith("_frequencies")
        ]

    @property
    def times(self):
        return np.array(self.attrs["times"])

    def frequencies(self, mode):
        return np.array(self[f"{mode}_frequencies"])


fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, 0.75 * FIGURE_WIDTH))
ax_long = inset_axes(ax, width="60%", height="35%", loc="lower right", borderpad=1)


@with_irf(0.3)
def biexp(times, amp1, amp2, tconst1, tconst2, offset):
    return biexponential(
        times,
        tzero=0,
        amp1=amp1,
        amp2=amp2,
        tconst1=tconst1,
        tconst2=tconst2,
        offset=offset,
    )


with PopulationDatabase(INPUT / "population_timeseries.hdf5", mode="r") as dbase:

    mode_energy = dict()
    for mode in dbase.modes:
        population = np.array(dbase[mode])
        total_energy = population * dbase.frequencies(mode)[:, None]
        mode_energy[mode] = np.sum(total_energy, axis=0)

    total_energy = sum(mode_energy[m] for m in mode_energy)
    total_energy /= total_energy.max()
    mode_energy["total"] = total_energy

    for mode, marker, color in zip(
        ["TA", "TO2", "total"], ["x", "o", "s"], ["orange", "blue", "red"]
    ):
        y = mode_energy[mode] / mode_energy[mode].max()
        fit_params, _ = curve_fit(
            biexp,
            xdata=dbase.times,
            ydata=y,
            p0=(1, -1, 1, 1, y.min()),
            bounds=([0, -np.inf, 0, 0, -1], [np.inf, 0, 150, 150, 1]),
        )
        for ax_ in (ax, ax_long):
            ax_.scatter(
                dbase.times,
                y,
                s=5,
                color=color,
                marker=marker,
                label=mode,
            )
            ax_.plot(
                dbase.times, biexp(dbase.times, *fit_params), color=color, linestyle="-"
            )

ax.legend(loc="center", ncol=3, bbox_to_anchor=(0.5, 1.05), bbox_transform=ax.transAxes)
ax.set_ylabel("Energy change [a.u.]")
ax.set_xlim([-5, 30])
ax.set_xlabel("Time-delay [ps]")

ax_long.xaxis.tick_top()
ax_long.xaxis.set_label_position("top")
ax_long.set_xlim([-20, 600])
ax_long.set_ylabel("Energy change [a.u.]")
ax_long.set_xlabel("Time-delay [ps]")
