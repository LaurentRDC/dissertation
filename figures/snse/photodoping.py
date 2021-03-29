from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal
import scipy.constants as constants

from plotutils import MEDIUM_FIGURE_WIDTH, discrete_colors

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

SAMPLE_AREA = (50e-6) * (50e-6)  # m^2
SAMPLE_THICKNESS = 45e-9  # m
SAMPLE_VOLUME = SAMPLE_AREA * SAMPLE_THICKNESS  # m^3


def absorbed_energy(fluence, thickness):
    """
    Total amount of absorbed energy

    Parameters
    ----------
    fluence : float
        Photoexcitation fluence [mJ/cm2]

    Returns
    -------
    abs : float
        Absorbed energy [J]

    References
    ----------
    L. Makinistian and E. A. Albanesi: On the band gap location and core spectra
    """
    penetration_depth = 100e-9  # m, from reference
    absorbed_ratio = 1 - np.exp(-thickness / penetration_depth)

    percm2_to_perm2 = 1e4
    mJ_to_J = 1e-3
    radiant_energy = fluence * mJ_to_J * percm2_to_perm2 * SAMPLE_AREA  # J
    return radiant_energy * absorbed_ratio


def photocarrier_density(fluence, thickness):
    """
    Calculate the photocarrier density deposited in sample, at 800nm.

    Parameters
    ----------
    fluence : float
        Photoexcitation fluence [mJ/cm2]

    Returns
    -------
    carrier_density : float
        Injected photocarrier density [10^21 / cm^3]
    """
    energy = absorbed_energy(fluence, thickness)  # J
    total_carriers = energy / (1.55 * constants.eV)
    m3_to_cm3 = 1e6
    return total_carriers / (SAMPLE_VOLUME * m3_to_cm3) / 1e21  # 1/cm^3


figure, carrier_ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 2.5))

fluences = np.linspace(0, 14, 1024)
densities = photocarrier_density(fluences, thickness=SAMPLE_THICKNESS)  # 10^21 / cm^3

carrier_ax.fill_between(
    fluences,
    photocarrier_density(fluences, thickness=SAMPLE_THICKNESS + 5e-9),
    photocarrier_density(fluences, thickness=SAMPLE_THICKNESS - 5e-9),
    facecolor="gray",
    edgecolor="k",
    linestyle="dashed",
    linewidth=1,
    alpha=0.2,
)

carrier_ax.plot(fluences, densities, "-", color=discrete_colors(1)[0])

carrier_ax.set_xlabel("Fluence [mJ/cm$^2$]")
carrier_ax.set_ylabel("$N_{\gamma}$ [$10^{21}$ cm$^{-3}$]")

carrier_ax.grid("on")

carrier_ax.set_xlim([fluences.min(), fluences.max()])
carrier_ax.set_ylim([densities.min(), densities.max()])

plt.tight_layout()
