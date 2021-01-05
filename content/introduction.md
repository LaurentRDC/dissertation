
# Introduction {#sec:introduction}

The basic physical principles of low-energy quantum mechanics are understood reasonably-well, and yet the emergence of new phenomena from the fundamental laws of Physics continues to surprise.

The fields of condensed matter physics and material science tackle a variety of intertwined questions, from fundamental physics to concrete applications.  Given that condensed matter systems are complex, many experimental dimensions are used to distinguish between phenomena. Electron microscopy comes to mind as an example of a technique which can selectively observe physical systems in real-space (via imaging) or in momentum (via diffraction).

In highly-ordered systems (i.e. crystals), where electronic correlation effects are most likely to matter, the time scales associated with fundamental actions are determined by the energy scale around room temperature. Phonons populated at room temperature (\SI{300}{\kelvin} = \SI{25}{\milli\electronvolt}) might have a period shorter than \SI{1}{\pico\second} ($10^{-12}$ \si{\second}).

## Ultrafast laser systems {#sec:laser_sources}

Digital measurement devices are generally limited to nanosecond resolution, and are therefore too slow to observe dynamics in the femtosecond range. Therefore, at the heart of most ultrafast experiment lies an ultrafast laser system. An ultrafast laser system consists of a pulsed source of highly-coherent radiation. 

Two things are required from the laser system: that its pulses be short, and that its pulses have high peak power. The length of laser pulses produced by such systems typically determine the time-resolution of experiments. The high peak power allows for the use of non-linear optical techniques, such as frequency-doubling.

Ultrafast laser systems like those used herein are composed of two optical assemblies. Ultrashort (\SI{<100}{\femto\second}) laser pulses are produced by a mode-locked laser system. These short pulses lack the high peak power, and are amplified via chirped pulse amplification. in a regenerative laser amplifier The resulting pulses are short (\SI{<100}{\femto\second}), while also reaching peak powers of up to \SI{1}{\peta\watt}.

In the following two subsections, a brief description of the two components of the laser system is given. These sections only highlight the basic functionality, while overlooking numerous details for brevity. All bust the most avid reader is encouraged to skip ahead to @sec:spectroscopy.

### Femtosecond mode-locked oscillators

The core working principle of mode-locked oscillators is to fix the phase relationship between allowed amplitude modes within a laser cavity. This phase relationship enforces that constructive interference happen regularly, producing a train of ultrashort pulses. These systems typically use titanium-doped sapphire as the gain medium (Ti-sapphire), because the gain bandwidth of this medium accomodates the intrinsically large bandwidth of short laser pulses (\SI{>45}{\nano\meter}, or \SI{>3}{\tera\hertz}).

@fig:mode-locking shows an example of the first 30 longitudinal modes of a cavity coming together to form a strong pulse. In reality, pulses may combine $10^6$ modes to maximize peak power [@Siegman1986].

```{.python #fig:mode-locking .matplotlib caption="Demonstration of how phase relationship between amplitude modes in a laser cavity can lead to a strong pulse. **a)** The first few longitudinal modes of a laser cavity. **b)** Combination of the first 45 modes of the cavity creates a very strong pulse in the center of the cavity."}
from skued import spectrum_colors

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(FIGURE_WIDTH, 4), sharex=True)

x = np.linspace(-0.5, 0.5, 512)
pulse = np.zeros_like(x)
mode = lambda x, n: np.cos(np.pi*(2*n + 1)*x)

colors = reversed(list(spectrum_colors(4)))
for n, c in enumerate(colors, start=1):
    ax1.plot(x, mode(x,n), color=c)
    pulse += mode(x,n)

for n in range(1, 45):
    pulse += mode(x,n)

ax2.plot(x, pulse/pulse.max(), '-k')

for ax in ax1, ax2:
    ax.xaxis.set_ticks([])

ax1.yaxis.set_ticks([0])
ax2.yaxis.set_ticks([0])

ax2.set_xlabel("Cavity position (a.u.)")
ax1.set_ylabel("E-field amplitude (a.u.)")
ax2.set_ylabel("E-field amplitude (a.u.)")

tag_axis(ax1, "a)", y=0.8)
tag_axis(ax2, "b)", y=0.8)
```

### Regenerative laser amplifiers


## Ultrafast spectroscopic techniques and their limitations {#sec:spectroscopy}

## Electron microscopy {#sec:microscopy}

## Ultrafast electron scattering {#sec:ues}

The extension of electron microscopy into time-domain studies is not new. As early as 1982, Gerard Mourou and Steve Williamson [@Mourou1982] follow the picosecond-scale transformation of an aluminum film following photoexcitation. In this work, the authors prepare \SI{100}{\pico\second} electron bunches using a streak camera. They note that the temporal resolution of the experiment is ultimately limited by the length of the electron bunch:

> *The electron pulse width has been measured by using the camera in the normal streak mode and is found to be ~ 100 ps. This value departs significantly from the 15-ps pulse width expected. The pulse broadening is due to the space-charge effect caused by the relatively high electron flux required to photograph the pattern with our present system.*

The challenges of ultrafast electron scattering on the \SI{100}{\femto\second} time-scale were apparent in the first ultrafast electron scattering experiments. Measurements of atomic structure require high beam brightness (i.e. high signal-to-noise), but dense electron bunches experience space-charge repulsion (i.e. worsening time-resolution). The trade-off between signal-to-noise and time-resolution is represented as the expansion of the bunch length $l$:

$$
\frac{d^2l}{dt^2} = \frac{N e^2}{m_e \epsilon_0 \pi r^2} \left[ 1 - \frac{l}{\sqrt{l^2 + 4 r^2}}\right]
$$ {#eq:bunch_length}

where $N$ is the number of electrons in the bunch, $e$ the electron charge, $m_e$ the electron mass, $\epsilon_0$ the vacuum permittivity, and $r$ is the beam radius[@Siwick2002]. A numerical solution of @eq:bunch_length is presented in @fig:bunch-length.

```{.python #fig:bunch-length .matplotlib caption=""}
import skued
import scipy.constants as constants
from scipy.constants import physical_constants
from scipy.integrate import odeint

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
    num_electrons = N/2

    # We cast the 2nd order nonlinear ODE l'' = f(l) into two equations
    # 
    #   u' = f(l)
    #   l' = u
    #
    # Let y = [u, l], a vector. Then, the two equations above become:
    #
    #   y' = [f(l), u]
    prefactor = (num_electrons * constants.e**2)/(constants.m_e * constants.epsilon_0 * np.pi * spotsize**2)
    f = lambda l: prefactor * (1 - l/np.sqrt(l**2 + 4 * spotsize**2))

    def system(y, t):
        u, l = y
        return [u, f(l)] # dy/dt
    
    # Propagation time is derived from propagation distance
    electron_velocity = skued.electron_velocity(energy) / 1e10 # m/s
    propagation_time = distance / electron_velocity

    # Initial conditions: bunch length is equivalent length to laser pulse (~50 fs),
    # and bunch is not expanding
    # Note that the initial condition for l' is taken from Brad's thesis,
    # where l' = 2 * vz (not that it matters...)
    y0 = [2*electron_velocity, 500e-15 * electron_velocity]
    
    sol = odeint(system, y0, propagation_time)
    length = sol[:,1] # only care about bunch length
    length /= electron_velocity
    return length

distance = np.linspace(0, 0.05, 1024) # meters

fig, ax = plt.subplots(1,1, figsize=(FIGURE_WIDTH,FIGURE_WIDTH))
order_of_magnitudes = [2,3,4,5,6,7]
colors = skued.spectrum_colors(len(order_of_magnitudes))

for i, c in zip(order_of_magnitudes, colors):
    num = 10**i
    ax.plot(distance, bunch_length(distance, N=num), '-', color=c, label=f"$10^{i}$ e/bunch")

#ax.axhline(y=1, linestyle='dashed', linewidth=1, color='k')

ax.set_xlabel('Propagation distance [m]')
ax.set_ylabel('FWHM pulse duration [s]')
plt.legend()
```

### Temporal control of electron bunches {#sec:cavity}

## Experimental apparatus {#sec:experimental_setup}

All ultrafast electron scattering experiments presented in this dissertation were performed with the same experimental apparatus, presented in @fig:apparatus.

![Experimental setup diagram for the ultrafast electron scattering instrument. BS1/2: Beam splitter. HPD: high-bandwidth photodiode. BBO: $\beta$-Barium borate crystal. DM: Dichroic mirror. BS1/2: Beamsplitter. L1/2/3: Focusing lens. EL1/2: Electron lens. APD: alignment photodiode. HPD: High-bandwidth photodiode. DL: Driving loop. PL: Pickup loop.](diagrams/ued_setup.pdf){#fig:apparatus}

The heart of the apparatus consists of the ultrafast laser system. A mode-locked oscillator (Spectra-Physics Tsunami) generates short (\SI{35}{\femto\second}) pulses of \SI{800}{\nano\meter} light, at a pulse repetition rate of \SI{75}{\mega\hertz}. Each oscillator pulse carries \SI{5}{\nano\joule} of energy, which is not enough to generate bright electron pulses required by ultrafast electron scattering experiments. Five percent of the oscillator pulse train is focused down to a \SI{35}{\micro\meter} spot on a photodiode with large bandwidth (\SI{12}{\giga\hertz}), generating a master clock signal with non-negligible spectral power up to \SI{3}{\giga\hertz}. The remaining ninety-five percent of the oscillator pulse train is used to seed a regenerative laser amplifier (Spectra-Physics Spitfire Pro). The oscillator pulses are amplified via chirped-pulse amplification [@Strickland1985], outputting \SI{35}{\femto\second} pulses each carrying \SI{3}{\milli\joule} of energy, at a repetition rate of \SI{1}{\kilo\hertz}. The amplified pulse train is split between two branches: pump and probe. 

On the pump line, the laser pulses propagate to a delay stage (Newport ILS250PP) on which a retroreflector is installed. This delay stage provides a resolution of \SI{2}{\micro\meter}, with total travel of \SI{25}{\centi\meter} (or \SI{50}{\centi\meter} of extra optical path length). These dimensions translate to a time-delay resolution of \SI{6}{\femto\second} and temporal range of \SI{1.6}{\nano\second}. The energy of the pump pulses is then attenuated by using a neutral-density filter, to the desired pump energies. A lens is also used to focus the pump pulses onto the sample. The attenuation and focus of the beam are adjusted to achieve the desired photoexcitation density (or *fluence*) on an area that is much larger than the sample, ensuring uniform illumination conditions. The attenuated and focused pump pulses are routed inside the vacuum chamber through a z-cut quartz (SiO$_2$) window, which is transparent over a large range of wavelengths. Inside the vacuum, two mirrors is positioned so that the pump pulses are transmitted through the sample and routed outside to a photodiode, for alignment purposes.

On the probe line, two barium borate (BBO) crystals are used to generate the third harmonic of \SI{800}{\nano\meter}. A first conversion to \SI{400}{\nano\meter}, and then a second conversion to \SI{266}{\nano\meter}. A calcite crystal is placed between the two BBO crystals to compensate for the wavelength-dependent dispersion in the crystals. This non-linear light conversion is possible because the pulses are extremely energy-dense, resulting in a conversion efficiency of 1%. The ultraviolet (UV) photons are then routed into the electron gun assembly, where they are focused on a copper hemispherical cathode. The work function of copper (\SI{4.5}{\electronvolt}[@Anderson1949]) is slightly lower than the UV photon energy (\SI{4.65}{\electronvolt}), resulting in electrons being freed from the copper cathode with little extra energy. These electrons are accelerated via a static potential of \SI{90}{\kilo\volt} to 53% of the speed of light. After acceleration, a solenoid lens loosely focuses the electron bunches through the \SI{3}{\milli\meter} aperture of a radio-frequency (RF) compression cavity, ensuring that electron bunches do not lose much charge before compression. Electron bunches are compressed using the cavity described in @sec:cavity, amplifying the 40^th^ harmonic of the oscillator reprate (\SI{2.9985}{\giga\hertz}). The RF compression power is tuned so that the temporal focus happens at the sample. After the RF compression cavity, another solenoid lens focuses the transverse profile of the electron bunches onto the detector, so that the diffraction pattern is imaged clearly. Electron bunches are transmitted through the sample at normal incidence. The scattering pattern is collected by an electron camera (Gatan Ultrascan 895), on a cooled charge-coupled detector. The transmitted, unscattered beam is collected by a Faraday cup. A Keithley 6514 electrometer measures the charge on the Faraday cup at a rate of \SI{1}{\kilo\hertz}, giving a rough estimation of the bunch-to-bunch charge.

## Overview of the dissertation {#sec:overview}

## References {.unnumbered}

\printbibliography[heading=none]
