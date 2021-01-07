import skued
import scipy.constants as constants
from scipy.constants import physical_constants
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from plotutils import FIGURE_WIDTH


def bunch_length(distance, N, spotsize=200e-6, energy=90):
    """
    Determine the bunch length of an electron bunch with `N` electrons
    propagating for `distance` meters. Electrons have `energy` keV. The bunch
    has a transverse radius of `spotsize`.

    Returns
    -------
    length : array_like
        FWHM bunch length in femtosecond.
    """
    # From reference, the FWHM bunch length can be calculated (vs. full bunch length)
    # if the following substitution is made:
    num_electrons = N / 2

    # We cast the 2nd order nonlinear ODE l'' = f(l) into two equations
    #
    #   u' = f(l)
    #   l' = u
    #
    # Let y = [u, l], a vector. Then, the two equations above become:
    #
    #   y' = [f(l), u]
    prefactor = (num_electrons * constants.e ** 2) / (
        constants.m_e * constants.epsilon_0 * np.pi * spotsize ** 2
    )
    f = lambda l: prefactor * (1 - l / np.sqrt(l ** 2 + 4 * spotsize ** 2))

    def system(y, t):
        u, l = y
        return [u, f(l)]  # dy/dt

    # Propagation time is derived from propagation distance
    electron_velocity = skued.electron_velocity(energy) / 1e10  # m/s
    propagation_time = distance / electron_velocity

    # Initial conditions: bunch length is equivalent length to laser pulse (~50 fs),
    # and bunch is not expanding
    # Note that the initial condition for l' is taken from Brad's thesis,
    # where l' = 2 * vz (not that it matters...)
    y0 = [2 * electron_velocity, 500e-15 * electron_velocity]

    sol = odeint(system, y0, propagation_time)
    length = sol[:, 1]  # only care about bunch length
    length /= electron_velocity
    return length


distance = np.linspace(0, 0.05, 1024)  # meters

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, FIGURE_WIDTH))
order_of_magnitudes = [2, 3, 4, 5, 6, 7]
colors = skued.spectrum_colors(len(order_of_magnitudes))

for i, c in zip(order_of_magnitudes, colors):
    num = 10 ** i
    ax.plot(
        distance, bunch_length(distance, N=num), "-", color=c, label=f"$10^{i}$ e/bunch"
    )

# ax.axhline(y=1, linestyle='dashed', linewidth=1, color='k')

ax.set_xlabel("Propagation distance [m]")
ax.set_ylabel("FWHM pulse duration [s]")
plt.legend()
