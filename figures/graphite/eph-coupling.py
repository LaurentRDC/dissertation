from functools import lru_cache, partial
from math import sqrt, pi
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal
from scipy import integrate
from scipy.constants import physical_constants
from scipy.optimize import curve_fit
from scipy.stats import sem
from skued import gaussian, with_irf
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors


INPUT = Path("data") / "graphite"

graphite = Crystal.from_database("C")

h, *_ = physical_constants["molar Planck constant"]
hbar = h / (2 * pi)  # J s / mol
boltzmann, *_ = physical_constants["Boltzmann constant"]  # J / K
avogadro, *_ = physical_constants["Avogadro constant"]  # 1 / mol


@lru_cache(1)
def molar_density(crystal):
    """ Molar density [mol / m**3] """
    mol_in_unitcell = len(crystal) / avogadro  # mol
    vol_of_unitcell = crystal.volume / 1e30  # m**3
    return mol_in_unitcell / vol_of_unitcell


def electron_heat_capacity(temperature):
    """
    Electronic heat capacity for GRAPHITE [J / mol / K].

    References
    ----------
    Nihira PRB 68
    """
    # Nihira 2003, Eq. 43
    heat_capacity = (
        (13.8 * temperature)
        + (1.16e-3 * temperature ** 2)
        + (2.6e-7 * temperature ** 3)
    )  # uJ / mol / K
    heat_capacity *= 1e-6  # J / mol / K
    return heat_capacity


def coth(*args, **kwargs):
    return np.cosh(*args, **kwargs) / np.sinh(*args, **kwargs)


def population(temperature, frequency):
    """
    Determine average phonon population at `temperature`.

    Parameters
    ----------
    temperature : float
        Temperature in Kelvin
    frequency : float
        Phonon frequency in Hertz
    """
    h, *_ = physical_constants["Planck constant over 2 pi in eV s"]
    kb, *_ = physical_constants["Boltzmann constant in eV/K"]
    return 0.5 * coth(h * frequency / (2 * kb * temperature)) - 0.5


def mode_heat_capacity(temperature, frequency):
    """
    Lattice heat capacity for a single phonon mode in GRAPHITE [J / mol / K].

    Parameters
    ----------
    temperature : array-like
        Lattice temperature [K]
    frequency : float
        Mode frequency [Hz]
    """
    # Nihira 2003, eq. 35
    # where frequency distribution function is the delta function at phonon frequency
    # (i.e. no integral)
    omega = 2 * pi * frequency  # rad Hz
    boltzmann_hz_per_K, *_ = physical_constants["Boltzmann constant in Hz/K"]
    operand = hbar * omega / (temperature * boltzmann_hz_per_K)

    heat_capacity = (
        boltzmann * (operand ** 2) * np.exp(operand) / (np.expm1(operand) ** 2)
    )  # J / K
    return heat_capacity * avogadro  # J / mol / K


k_to_heat_capacity = partial(mode_heat_capacity, frequency=38.7e12)


def k_to_sink_heat_capacity(temperature):
    """
    Lattice heat c apacity for all modes associated with the K-TO decay.

    References
    ----------
    Bonini et al., Phonon Anharmonicities in Graphite and Graphene. PRL (2007)
    """
    # The decay from A1' mode is broken down in Fig. 4
    return sum(
        [
            0.09
            * mode_heat_capacity(temperature, 19.35e12)
            * 2,  # Decay into two TA modes
            0.36
            * (
                mode_heat_capacity(temperature, 24.2e12)
                + mode_heat_capacity(temperature, 14.5e12)
            ),
            0.55
            * (
                mode_heat_capacity(temperature, 2.4e12)
                + mode_heat_capacity(temperature, 36.2e12)
                + mode_heat_capacity(temperature, 1.9e12)
                + mode_heat_capacity(temperature, 36.8e12)
            ),
        ]
    )


def pulse_profile(time, fluence, rms):
    """
    Energy deposited by the pulse (assumed gaussian) of the photoexcitation pulse [J / mol / s].

    Parameters
    ----------
    time : array-like
        Time array [ps]
    fluence : float
        Excitation density [mJ/cm^2]
    rms : float
        RMS duration of the pulse.
    """
    # graphite optical depth at 800nm
    optical_penetration_depth = 15e-9  # meters

    # Energy deposited per area
    fluence = 12  #   mJ / cm2
    fluence *= 1e-3  # J / cm2
    fluence *= 1e4  # J / m**2

    # Energy deposited per volme
    excitation = fluence / optical_penetration_depth  # J / m**3

    # Energy deposited per mole
    molar_energy_deposited = excitation / molar_density(graphite)  # J / mol

    # Energy deposited per mole, over time
    return molar_energy_deposited * gaussian(time, center=0, fwhm=rms)


pulse_func = partial(pulse_profile, fluence=12, rms=35e-3)


def system_evolution(T, time, profile, g_ek, g_kl, g_el):
    """
    Evolution of the system over time.

    Parameters
    ----------
    T : np.array, shape (2,)
        Array of the electronic temperature and lattice temperature
    time : float
    profile : callable
        Callable of one argument (time)
    g_ek : float
        Coupling constant between electronic system and K-point phonon
    g_kl : float
        Coupling constant between K-point phonon and rest of lattice
    g_el : float
        Coupling constant between electronic system and rest of lattice.
    """
    Te, Tk, Tl = T[0], T[1], T[2]

    # Energy deposited into the electronic system
    energy_e = -g_ek * (Te - Tk) - g_el * (Te - Tl) + profile(time)

    # Energy deposited into the A1' mode (K-TO)
    energy_k = g_ek * (Te - Tk) - g_kl * (Tk - Tl)

    # Energy deposited by the K-TO mode into the lattice
    energy_l = g_el * (Te - Tl) + g_kl * (Tk - Tl)

    return np.array(
        [
            energy_e / electron_heat_capacity(Te),
            energy_k / k_to_heat_capacity(Tk),
            energy_l / k_to_sink_heat_capacity(Tl),
        ]
    )


with_irf(0.35)


def fit_function(times, amp, g_ek, g_kl, g_el):
    """
    Time-trace representing the dynamics of the K-TO mode population.

    Parameters
    ----------
    times : array_like
        Time points [s]
    amp : float
        Amplitude [a.u.]
    g_ek : float
        Coupling constant between electronic system and K-point phonon [W / mol / K]
    g_kl : float
        Coupling constant between K-point phonon and rest of lattice [W / mol / K]
    g_el : float
        Coupling constant between electronic system and rest of lattice [W / mol / K]
    """
    # NOTE: for some reason, if all initial temperatures are 300K, the model does not evolve.
    #       This is why the initial temperatures are slighly different
    system = partial(
        system_evolution, profile=pulse_func, g_ek=g_ek, g_kl=g_kl, g_el=g_el
    )
    T = integrate.odeint(system, np.array([300, 301, 302]), times)

    return amp * (T[:, 1] - 300)


# -----------------------------------------------------------------------------


def electron_density_of_states(energy):
    """
    Approximate electron density of states of graphite [eV/unitcell].
    From Eq. S5 in ref.

    Parameters
    ----------
    energy : float
        Energy [eV]

    Returns
    -------
    dos : float
        Density of states [1/eV/atom]

    References
    ----------
    [#] A. H. Castro Neto et al., The electronic properties of graphene. Rev. Mod. Phys.
    DOI: 10.1103/RevModPhys.81.109
    """
    # Constants from reference.
    lattice_constant = 1.42  # Angstroms
    in_plane_area = (3 / 2) * np.sqrt(3) * lattice_constant ** 2

    # The estimation of the hopping parameter is provided shown in ref
    # at top of page 114
    # fermi_velocity ~ 10^6 m / s
    hopping_parameter = 2.8  # eV
    fermi_velocity = 3 * hopping_parameter * lattice_constant / 2

    # DOS in [1/eV/atom]
    # Reference has form in [1/eV/unitcell], which is 4 times larger since
    # graphite has 4 unitcell atoms
    return (in_plane_area / (2 * np.pi)) * np.abs(energy) / (fermi_velocity ** 2)


def calculate_ep_coupling():
    """
    calculate e-(K-TO) coupling <g>^2 [eV] from the heat rate G_(e,ph)
    """
    hbar_eV, *_ = physical_constants["Planck constant over 2 pi in eV s"]  # eV s
    Gek = 6.8e17 / molar_density(graphite)  # W / mol / K
    Gek_err = 0.3e16 / molar_density(graphite)  # W / mol / K

    Gel = 8.32e9 / molar_density(graphite)
    Gel_err = 6e15 / molar_density(graphite)

    # Use temperature around peak of K-TO transient intensity (~750fs)
    Ce = electron_heat_capacity(8500)
    Ck = k_to_heat_capacity(4500)
    Cl = k_to_sink_heat_capacity(1000)

    # "harmonic" is the equation that defines 1/tau
    # I call it harmonic because it looks like the equation for
    # current flow in parallel resistors.
    harmonic = sum(
        Gj / Ce - sum(Gi / Cj for Gi in [Gek, Gel])
        for Gj, Cj in zip([Gek, Gel], [Ck, Cl])
    )
    tconst = 1 / harmonic

    # Propagation of error with derivatives
    tconst_err = sqrt(
        abs(Gek_err / (Ce * harmonic ** 2)) ** 2
        + abs(Gel_err / (Ce * harmonic ** 2)) ** 2
    )

    # k-TO mode has energy of 0.2 eV
    dos = electron_density_of_states(1.55 - 0.165)  # eV/unitcell
    g2 = (hbar_eV / tconst) * (1 / (2 * np.pi * dos))
    gerr = (hbar_eV / tconst ** 2) * (1 / (2 * np.pi * dos)) * tconst_err

    print("tau:    {:.2e} ± {:.2e} fs".format(tconst / 1e-15, tconst_err / 1e-15))
    print("<g>^2: ({:.3e} ± {:.0e})".format(g2, gerr), "eV^2")


# ------------------------------------------

# K-TO diffuse intensity dynamics
# Note that the time is given in picoseconds
# Therefore we need to cast to seconds
times, a1prime_dynamics, _ = np.loadtxt(
    INPUT / "a1prime_dynamics_with_errs.csv", delimiter=",", unpack=True, skiprows=1
)

# Error from calculations is ridiculously small
# Instead, we use the standard error in the mean before time zero, where
# a1prime_dynamics should be 0
a1prime_err = np.ones_like(times) * sem(a1prime_dynamics[times < 0])

a1prime_err /= a1prime_dynamics.max()
a1prime_dynamics /= a1prime_dynamics.max()

simulation_times = np.arange(-2, 10, step=1e-3)

# K-TO diffuse intensity dynamics
# Fitting the K-point dynamics
restricted_times = times[np.logical_and(times > -2, times < 10)]
restricted_dynamics = a1prime_dynamics[np.logical_and(times > -2, times < 10)]
restricted_errs = a1prime_err[np.logical_and(times > -2, times < 10)]
best_params, covariance = curve_fit(
    fit_function,
    restricted_times,
    restricted_dynamics,
    sigma=restricted_errs,
    absolute_sigma=True,
    p0=[
        17e-4,
        5e17 / molar_density(graphite) / 1e12,
        1e18 / molar_density(graphite) / 1e12,
        1e16 / molar_density(graphite) / 1e12,
    ],
    bounds=(0, np.inf),
)
fit_amp, fit_ek_coupling, fit_kl_coupling, fit_el_coupling = best_params
_, err_ek_coupling, err_kl_coupling, err_el_coupling = np.sqrt(np.diag(covariance))

# # Additional factor of 1e12 is to convert to W (J/s) from (J/ps)
# print(
#     "Electron-K coupling: ",
#     "({:.8e} ± {:.2e})".format(
#         fit_ek_coupling * molar_density(graphite) * 1e12,
#         err_ek_coupling * molar_density(graphite) * 1e12,
#     ),
#     "W / m3 / K",
# )
# print(
#     "K-lattice coupling: ",
#     "({:.8e} ± {:.2e})".format(
#         fit_kl_coupling * molar_density(graphite) * 1e12,
#         err_kl_coupling * molar_density(graphite) * 1e12,
#     ),
#     "W / m3 / K",
# )
# print(
#     "Electron-lattice coupling: ",
#     "({:.8e} ± {:.2e})".format(
#         fit_el_coupling * molar_density(graphite) * 1e12,
#         err_el_coupling * molar_density(graphite) * 1e12,
#     ),
#     "W / m3 / K",
# )
best_curve = fit_function(simulation_times, *best_params)

# Render simulation using the extracted couplings
example = partial(
    system_evolution,
    profile=pulse_func,
    g_ek=fit_ek_coupling,  # W / mol / K,
    g_kl=fit_kl_coupling,
    g_el=fit_el_coupling,
)
T = integrate.odeint(example, np.array([300, 301, 302]), simulation_times)

electronic_temperature, k_to_temperature, k_to_drain_temperature = (
    T[:, 0],
    T[:, 1],
    T[:, 2],
)

intensity_to_temperature = lambda a: 300 + a / fit_amp
temperature_to_intensity = lambda a: fit_amp * (a - 300)

# Create the figure
# We report times in picoseconds
# Colors are in cycle order (C0, C1, C2, ...)
fig, ax_K = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3))
ax_T = ax_K.twinx()

colors = discrete_colors(2)
ax_K.plot(
    simulation_times,
    best_curve,
    linestyle="solid",
    color=colors[0],
    label="$A_1^{\prime}$ phonon",
)

ax_K.plot(
    simulation_times,
    temperature_to_intensity(k_to_drain_temperature),
    linestyle=(0, (1, 1)),
    color=colors[1],
    label="Lattice drain",
)

# Plot data last so it appears on top of the fit lines
ax_K.errorbar(
    times,
    a1prime_dynamics,
    yerr=a1prime_err,
    color="k",
    marker=".",
    linestyle="None",
    markersize=5,
)

ax_K.set_xlim([simulation_times.min(), simulation_times.max()])

ax_K.set_ylim([-0.25, 1.1])
ax_T.set_ylim(intensity_to_temperature(np.array(ax_K.get_ylim())))

ax_K.axvline(x=0, linestyle="--", linewidth=1, color="k")
ax_K.set_ylabel(r"$\Delta n_{A_1'}(\tau)$ [a.u.]")
ax_T.set_ylabel(r"Temperature [K]")

ax_K.set_xlabel("Delay time [ps]")

ax_K.legend(
    loc="upper right",
    edgecolor="none",
)
plt.tight_layout()
