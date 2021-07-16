from math import pi, sin, sqrt
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scipy.constants as constants
import scipy.optimize as opt
import scipy.stats as stats
from crystals import Crystal
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from dissutils import MEDIUM_FIGURE_WIDTH
from dissutils.snse import overnight4, photocarrier_density

DATADIR = Path("data") / "snse"
CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")


def box_errorbars(ax, xdata, ydata, xerr, yerr, colors):
    """
    Create error boxes instead of error bars
    """
    # Create list for all the error patches
    errorboxes = []
    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata, ydata, xerr, yerr):
        rect = Rectangle((x - xe, y - ye / 2), 2 * xe, ye)
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, edgecolor=colors, facecolor=colors)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot invisible scatter points to get automatic plot bounds that
    # make sense
    artists = ax.scatter(xdata, ydata, s=0, c=colors)

    return artists


def transient_u2(timedelays, intensity, q):
    """
    Determine the change in MSD from intensity changes

    Parameters
    ----------
    timedelays : ndarray
        Time-delays [ps]
    intensity : ndarray
        Diffraction intensity, unnormalized [counts]
    q : float
        Amplitude of the scattering vector |q| where `intensity` was generated [Angstroms^-1]

    Returns
    -------
    delta_u2_rms : ndarray
        Transient mean-square-displacement Delta <U^2> [Angstroms^2]
    """
    I = np.array(intensity, copy=True)
    I0 = np.mean(intensity[timedelays < 0])
    deltaI0 = (I - I0) / I0
    u2 = -np.log(1 + deltaI0) / q ** 2

    # Propagation of errors : error in [A ln(1 + x)] is [A x_err/x]
    sigma_I0 = np.max(deltaI0[timedelays < 0]) - np.min(deltaI0[timedelays < 0])
    err = -sigma_I0 / ((q ** 2) * (1 + deltaI0))
    return u2, err


def delta_msd_fluence(fluence):
    filenames = {
        6.6: DATADIR / "debyewaller_6p6mjcm2_022.csv",
        7.9: DATADIR / "debyewaller_7p9mjcm2_022.csv",
        9.5: DATADIR / "debyewaller_9p5mjcm2_022.csv",
        10.7: DATADIR / "debyewaller_10p7mjcm2_022.csv",
        12: DATADIR / "debyewaller_12mjcm2_022.csv",
        13.2: DATADIR / "debyewaller_13p2mjcm2_022.csv",
    }

    filepath = filenames[fluence]
    timedelays, timeseries = np.loadtxt(
        filenames[fluence], skiprows=1, delimiter=",", unpack=True
    )
    # It is assumed that debye-waller curves are taken from the (022) reflection
    q = np.linalg.norm(CRYSTAL.scattering_vector((0, 2, 2)))
    u2, u2err = transient_u2(timedelays, timeseries, q)

    # Factor of 1/sin(45deg) because reflection (022) is not
    # parallel to (001), but 45deg from it.
    return timedelays, u2 / sin(pi / 4), u2err / sin(pi / 4)


figure, ax_displacement = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))
ax_carrdensity = ax_displacement.twiny()

fluences = np.array([6.6, 7.9, 9.5, 10.7, 12, 13.2])
rms = []
errs = []
for fluence in fluences:
    timedelays, delta_msd, delta_msd_err = delta_msd_fluence(fluence)
    rms.append(np.mean(delta_msd[np.logical_and(timedelays > 5, timedelays < 9)]))
    errs.append(np.mean(delta_msd_err[np.logical_and(timedelays > 5, timedelays < 9)]))

rms = np.asarray(rms) * 1e2  # 10^-2 AA^2
errs = np.asarray(errs) * 1e2  # 10^-2 AA^2

# Fit displacement values with linear
sparams, spcov = opt.curve_fit(lambda x, a, b: a * x + b, xdata=fluences, ydata=rms)
slope, intercept = sparams[0], sparams[1]
slope, slope_err = sparams[0], np.sqrt(spcov[0, 0])
intercept, intercept_err = sparams[1], np.sqrt(spcov[1, 1])


cmap = plt.get_cmap("inferno")
colors = [cmap(35 * val + 20) for val in range(len(fluences))]

es = np.linspace(fluences.min() - 0.6, fluences.max() + 0.6, num=1024)
fit_curve = slope * es + intercept
ax_displacement.plot(es, fit_curve, linewidth=1, color="gray", zorder=0)
ax_displacement.fill_between(
    es,
    (slope - slope_err) * es + (intercept - intercept_err),
    (slope + slope_err) * es + (intercept + intercept_err),
    facecolor="gray",
    edgecolor="k",
    linestyle="dashed",
    linewidth=1,
    alpha=0.2,
)

box_errorbars(
    ax_displacement,
    fluences,
    rms,  # 10^-2 AA^2
    xerr=np.full_like(rms, fill_value=0.2),  # mj / cm2
    yerr=errs,  # 10^-2 AA^2
    colors=colors,
)

ax_displacement.set_xlim([es.min(), es.max()])
ax_displacement.set_ylabel("$\\Delta \\langle u_c^2 \\rangle$ [$10^{-2} \AA^2$]")
ax_displacement.set_xlabel("Fluence [mJ cm$^{-2}$]")


carrier_density = photocarrier_density(es) / 1e20
ax_carrdensity.set_xlim([carrier_density.min(), carrier_density.max()])
ax_carrdensity.set_xticks(photocarrier_density(fluences) / 1e20)
ax_carrdensity.xaxis.set_major_formatter(ticker.FuncFormatter(lambda n, _: f"{n:.1f}"))
ax_carrdensity.xaxis.set_label_position("top")
ax_carrdensity.xaxis.tick_top()
ax_carrdensity.set_xlabel("$N_{\gamma}$ [$10^{20}$ cm$^{-3}$]")

plt.tight_layout()
