"""
Determine the change in vibrational amplitude from Debye-Waller
"""

import sys
from math import ceil, cosh, sinh
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np
import scipy.optimize as opt
import scipy.stats
from crystals import Crystal
from scipy.constants import physical_constants
import scipy.constants as constants
from math import floor, log10
from warnings import simplefilter

from plotutils.snse_datasets import overnight4
from plotutils import FIGURE_WIDTH, box_errorbars
import skued
from iris import DiffractionDataset

DATADIR = Path("data") / "snse"
CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")

# SPECTROSCOPIC STUDIES ON LEAD SELENIDE (PbSe) AND TIN SELENIDE(SnSe) THIN FILMS
# Journal of Optoelectronics and Biomedical Materials
ABSORBANCE = 0.075  # at 800nm

SAMPLE_AREA = (50e-6) * (50e-6)  # meters^2
THICKNESS = 30e-9  # meters


def round_sigfig(x, u):
    round_ = lambda n: round(n, -int(floor(log10(abs(u)))))
    return round_(x), round_(u)


def energy_per_cell(fluence, thickness=30e-9):
    """
    Calculate the energy deposited per unit cell.

    Parameters
    ----------
    fluence : float
        Photoexcitation fluence [mJ/cm2]
    thickness : float, optional
        Sample thickness [m]

    Returns
    -------
    energy : float
        Energy deposited per cell [eV]
    """
    fluence_mjm2 = fluence * 10_000  # mJ / m2

    # We assume photoexcitation profile is roughly flat on this sample
    # since the sample is so small
    total_energy = ABSORBANCE * fluence_mjm2 * SAMPLE_AREA  # mJ

    ncells = int(SAMPLE_AREA * thickness / (CRYSTAL.volume * 10 ** (-30)))

    # Conversion to eV deposited
    mJ_per_eV = 1e-3 / constants.electron_volt
    total_energy *= mJ_per_eV
    return total_energy / ncells


def photocarrier_density(energy_per_cell):
    """
    Calculate the photocarrier density deposited in sample, at 800nm.

    Parameters
    ----------
    energy_per_cell : array_like
        Energy per cell [eV / cell]

    Returns
    -------
    carrier_density : float
        Injected photocarrier density [cm^-3]
    """
    energy_per_cell = np.asarray(energy_per_cell)
    Egap = 1.1
    uvol = CRYSTAL.volume * 1e-24  # cm^3
    return (energy_per_cell / Egap) / uvol


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
    normalized = np.array(intensity, copy=True)
    normalized /= np.mean(intensity[timedelays < 0])
    return -3 * np.log(normalized) / (q ** (1 / 2))


def to_frequency(delta_msd):
    """
    Calculate the softening/hardening of mode, assuming no energy changes, i.e. change in amplitude
    is only due to change in frequency, but energy is conserved.

    Parameters
    ----------
    delta_msd : ndarray
        Transient change in RMS motion <u^2> [Ang]

    Returns
    -------
    freq : ndarray
        Fractional change in vibrational frequency [a.u.]
    """
    # We only consider the mean-square displacement for one of the 3N modes
    N = len(CRYSTAL)
    tu = thermal_msd() / (3 * N)
    return 1 / np.sqrt(1 + (delta_msd / tu))


def thermal_msd():
    """ Thermal mean-square-amplitude amplitude for atoms in the "soft" TO mode [Angstroms^2]. """

    def coth(*args):
        return cosh(*args) / sinh(*args)

    amu_to_kg, *_ = physical_constants["atomic mass constant"]
    kb, *_ = physical_constants["Boltzmann constant"]  # J / K
    hbar_SI, *_ = physical_constants["Planck constant over 2 pi"]  # J s

    avg_atom_mass = (
        amu_to_kg * sum(atom.mass for atom in CRYSTAL) / len(CRYSTAL)
    )  # kg/atom

    w0 = 1e12  # Hz
    m2_to_Ang2 = 1e20
    return (
        (hbar_SI / (2 * avg_atom_mass * w0))
        * coth(hbar_SI * w0 / (2 * kb * 300))
        * m2_to_Ang2
    )


def delta_msd_fluence(fluence):
    # if fluence == 2:
    #     return delta_msd_2mjcm2()

    #
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
    return timedelays, transient_u2(timedelays, timeseries, q)


def delta_msd_2mjcm2():
    """
    Compute the delta RMS for a particular fluence

    Returns
    -------
    timedelays : ndarray
    delta_msd : ndarray
    """
    reflections = [(0, 1, 2), (0, -1, 2), (0, 0, 1), (0, 0, -1), (0, 1, 0), (0, -1, 0)]

    with DiffractionDataset(overnight4.path, mode="r") as dset:
        timedelays = dset.time_points

        delta_msd = np.zeros_like(timedelays)
        for (h, k, l) in reflections:
            yj, xj = overnight4.miller_to_arrindex(h, k, l)
            timeseries = dset.time_series_selection(
                skued.DiskSelection(
                    shape=dset.resolution,
                    center=(xj, yj),
                    radius=30,
                )
            )
            q = np.linalg.norm(CRYSTAL.scattering_vector((h, k, l)))
            delta_msd += transient_u2(timedelays, timeseries, q) / len(reflections)

    return timedelays, delta_msd / 15  # thicker sample


figure, ax_softening = plt.subplots(1, 1, figsize=(5, 4))
ax_displacement = ax_softening.inset_axes([0.6, 0.6, 0.38, 0.38])

ax_density = ax_softening.twiny()
ax_density.xaxis.set_label_position("top")
ax_density.xaxis.tick_top()
ax_density.yaxis.set_visible(False)


fluences = np.array([6.6, 7.9, 9.5, 10.7, 12, 13.2])
energies_per_cell = []
rms = []
softenings = []
for fluence in fluences:
    timedelays, delta_msd = delta_msd_fluence(fluence)
    energies_per_cell.append(
        energy_per_cell(fluence=fluence, thickness=30e-9) * 1e3
    )  # eV to meV
    rms_val = np.mean(delta_msd[np.logical_and(timedelays > 5, timedelays < 9)])
    rms.append(rms_val)
    softenings.append(to_frequency(rms_val))

# Fit softening values with 1/x
sparams, spcov = opt.curve_fit(
    lambda x, a, b: a / x + b, xdata=fluences, ydata=softenings
)
slope, intercept = sparams[0], sparams[1]
slope, slope_err = round_sigfig(sparams[0], np.sqrt(spcov[0, 0]))
intercept, intercept_err = round_sigfig(sparams[1], np.sqrt(spcov[1, 1]))


cmap = plt.get_cmap("inferno")
colors = [cmap(35 * val + 20) for val in range(len(fluences))]

es = np.linspace(0.95 * min(fluences), 1.05 * max(fluences), num=1024)
fit_curve = slope / es + intercept
ax_softening.plot(es, fit_curve, linewidth=1, color="gray", zorder=0)
ax_softening.fill_between(
    es,
    (slope - slope_err) / es + (intercept - intercept_err),
    (slope + slope_err) / es + (intercept + intercept_err),
    facecolor="gray",
    edgecolor="k",
    linestyle="dashed",
    linewidth=1,
    alpha=0.2,
)
ax_softening.set_xlim([es.min(), es.max()])

energy_per_cell_err = 5  # mev
fluence_err = 0.2  # mj / cm2
sc = box_errorbars(
    ax_softening,
    fluences,
    softenings,
    xerr=fluence_err,
    yerr=0.007,
    colors=colors,
)

box_errorbars(
    ax_displacement,
    fluences,
    np.array(rms) * 1e2,  # to 10^-2 AA^2
    xerr=fluence_err,
    yerr=0.005 * 1e2,  # to 10^-2 AA^2
    colors=colors,
)

# Plot invisible data to show photocarrier density
densities = photocarrier_density(np.array(energies_per_cell) * 1e-3)  # meV to eV
ax_density.errorbar(
    densities * 1e-20,
    softenings,
    xerr=0.3,  # Error matches 7 meV
    yerr=0.01,
    color="none",
)

# Add colorbar to show fluences
# Note that to match the color of the points,
cax = figure.add_axes([0.19, 0.18, 0.55, 0.02])
cb = mpl.colorbar.ColorbarBase(
    cax,
    cmap=mpl.colors.ListedColormap(colors),
    orientation="horizontal",
    ticklocation="top",
    alpha=0.7,
    ticks=[i + 0.5 for i, _ in enumerate(colors)],
    format=ticker.FixedFormatter([str(ceil(ec)) for ec in energies_per_cell]),
    norm=mpl.colors.Normalize(vmin=0, vmax=len(energies_per_cell)),
)
cb.set_label("$E_c$ [meV]")

ax_displacement.set_ylabel("$\\Delta \\langle u^2 \\rangle$ [$10^{-2} \AA^2$]")
ax_displacement.set_xlabel("Fluence [mJ cm$^{-2}$]")

ax_softening.set_ylabel("$\\omega / \\omega_0$ [a.u.]")
ax_softening.set_xlabel("Fluence [mJ cm$^{-2}$]")
ax_softening.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))

ax_density.set_xlabel("Photocarrier density [$10^{20}$ cm$^{-3}$]")

# Tight layout plot warnings
simplefilter(action="ignore", category=UserWarning)
plt.tight_layout()
