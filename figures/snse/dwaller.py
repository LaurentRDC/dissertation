from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset
from plotutils.snse_datasets import overnight4

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

INNER_RADIUS = 15

# Determine the peak position and the time-series integration bounding box
INDICES_BRAGG = [(0, 1, 2), (0, -1, -2), (0, 0, 1), (0, 0, -1), (0, 1, 0), (0, -1, 0)]
inferno = plt.get_cmap("inferno")
DWCOLOR, GAMMA_COLOR, AVERAGE_COLOR = inferno(30), inferno(100), inferno(200)

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, *args, **kwargs)


timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    timeseries["debye-waller"] = np.zeros_like(timedelays)
    for h, k, l in INDICES_BRAGG:
        yj, xj = overnight4.miller_to_arrindex(h, k, l)
        q2 = np.linalg.norm(CRYSTAL.scattering_vector((h, k, l))) ** 2
        timeseries["debye-waller"] += (
            dset.time_series_selection(
                skued.DiskSelection(
                    shape=dset.resolution,
                    center=(xj, yj),
                    radius=INNER_RADIUS,
                )
            )
            / q2
        )

timeseries["debye-waller"] /= np.mean(timeseries["debye-waller"][timedelays < 0])

dwparams, dwpcov = opt.curve_fit(
    biexponential,
    timedelays[timedelays < 20],
    timeseries["debye-waller"][timedelays < 20],
    p0=(0, 0.01, 0.02, 0.2, 3.7, timeseries["debye-waller"].min()),
    bounds=([-0.2, -1, -1, 0, 1, 0.9], [0.2, 1, 1, 1, 5, 1.1]),
)

figure, diffuse_ax = plt.subplots(1, 1, figsize=(4.25, 2))
diffuse_ax.axhline(y=1, linestyle="dashed", color="k", linewidth=0.5)
diffuse_ax.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5)

dwfit_curve = biexponential(timedelays, *dwparams)
plot_params = dict(
    color=DWCOLOR,
    markersize=2,
    linestyle="None",
    elinewidth=0.5,
)
diffuse_ax.plot(timedelays, dwfit_curve, linewidth=1, color=DWCOLOR)

diffuse_ax.errorbar(
    x=timedelays,
    y=timeseries["debye-waller"],
    yerr=scipy.stats.sem(timeseries["debye-waller"][timedelays < 0]),
    **plot_params,
)


diffuse_ax.set_xlim([-1.6, 12])
diffuse_ax.set_ylim([0.967, 1.005])

diffuse_ax.set_xlabel("Time-delay [ps]")
diffuse_ax.set_ylabel("$\Delta I/I_0$ [a.u.]")
diffuse_ax.yaxis.set_ticks([0.97, 0.98, 0.99, 1.00])

plt.tight_layout()
