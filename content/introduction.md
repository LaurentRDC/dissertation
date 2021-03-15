
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

```{.python #fig:mode-locking .matplotlib file="figures/introduction/mode-locking.py" caption="Demonstration of how phase relationship between amplitude modes in a laser cavity can lead to a strong pulse. **a)** The first few longitudinal modes of a laser cavity. **b)** Combination of the first 45 modes of the cavity creates a very strong pulse in the center of the cavity."}
```

### Regenerative laser amplifiers


## Ultrafast spectroscopic techniques and their limitations {#sec:spectroscopy}

## Brief history of ultrafast electron scattering {#sec:ues}

The extension of electron microscopy into time-domain studies is not new. This section will be limited to ultrafast electron scattering in the form of diffraction; see the review by King *et al.*[@King2005] for a historical perspective on ultrafast electron microscopy in general. As early as 1982, Gerard Mourou and Steve Williamson [@Mourou1982] followed the picosecond-scale transformation of an aluminum film following photoexcitation. In this work, the authors prepare \SI{100}{\pico\second} electron bunches using a streak camera. They note that the temporal resolution of the experiment is ultimately limited by the length of the electron bunch:

> *The electron pulse width has been measured by using the camera in the normal streak mode and is found to be ~ 100 ps. This value departs significantly from the 15-ps pulse width expected. The pulse broadening is due to the space-charge effect caused by the relatively high electron flux required to photograph the pattern with our present system.*

The challenges of ultrafast electron scattering on the \SI{100}{\femto\second} time-scale were apparent in the first ultrafast electron scattering experiments. Measurements of atomic structure require high beam brightness (i.e. high signal-to-noise), but dense electron bunches experience space-charge repulsion (i.e. worsening time-resolution). The trade-off between signal-to-noise and time-resolution is represented as the expansion of the bunch length $l$:
$$
\frac{d^2l}{dt^2} = \frac{N e^2}{m_e \epsilon_0 \pi r^2} \left[ 1 - \frac{l}{\sqrt{l^2 + 4 r^2}}\right]
$$ {#eq:bunch_length}
where $N$ is the number of electrons in the bunch, $e$ the electron charge, $m_e$ the electron mass, $\epsilon_0$ the vacuum permittivity, and $r$ is the beam radius[@Siwick2002]. The first demonstration of ultrafast electron scattering was done by the same authors[@Mourou1982] to study the ultrafast melting of aluminium[@Williamson1984] with a time-resolution a tens of picoseconds.

The modern form of pump-probe ultrafast electron scattering was brought by Ahmed Zewail. At the time, Zewail was already known as the "Father of Femtochemistry" for the development of ultrafast spectroscopy[@Zewail1988]. Zewail and collaborator proposed to replace the optical probe pulse used in pump-probe femtosecond spectroscopy to an electron pulse[@Williamson1991; @Williamson1997]. What followed was a series of experiments probing the transient molecular structures at the picosecond time-scale following photoexcitation in the gas phase[@Ihee2001; @Ruan2001; @Lobastov2001]. The ability to track the nuclear coordinates directly rather than indirectly represented the perfect complement to Zewail's Nobel prize-winning work on femtochemistry via ultrafast spectroscopy.

It was not until Siwick *et al.*[@Siwick2003; @Siwick2004], however, that the ideas of Mourou and Williamson were applied with sub-picosecond sensitivity (\SI{600}{\femto\second}) to the solid state. This work reached a time-resolution sufficient to address whether or not the ultrafast phase transition from solid to liquid-like was analogous to the thermal phase transition[@Guo2000].

By 2005, the development of the field was booming. At the time, the prevalent sentiment was that the space-charge problem formulated by Mourou and Williamson limited the applicability of ultrafast electron scattering. King *et al.*[@King2005] identify three paths to control space-charge effects:

1. Reduce the electron-emission -- sample distance. By reducing this distance, the space-charge-driven expansion of the electron bunch remains limited[@Siwick2002].

2. Use higher electron energies. For \si{\mega\electronvolt} electrons, the space-charge explosion is reduced significantly[@Gian2002].

3. Reducing space-charge pulse expansion by limiting the energy spread in the electron bunch[@Gallant2000].

### Temporal control of electron bunches {#sec:cavity}

The solution to space-charge grew from the insight of Siwick *et al.* [@Siwick2002], who showed that the electrons in a bunch develop a strong correlation in phase space as the bunch propagates. Space-charge effects effectively slow down the front electrons with respect to the center of charge, while it accelerates the electrons at the front of the bunch. After some propagation time, a linear correlation is established between the relative axial position within the bunch and the relative velocity within the bunch. This is shown in @fig:introduction-bunch-correlation.

![Relative axial velocity vs axial position ($Z$) for all electrons in the pulse at four times ($T$) during its propagation. The pulse is composed of 10 000 electrons at \SI{30}{\kilo\electronvolt}, with an initial bunch length equivalent to \SI{150}{\femto\second}, an initial radius ut \SI{75}{\micro\meter}, and initial beam divergence of \SI{1.5}{\milli\radian}. The width of the velocity distribution becomes larger as the pulse lengthens. The spatial distribution of velocities also evolves. As electrons redistribute themselves inside the packet a linear velocity chirp develops on the electron pulse. For \SI{30}{\kilo\electronvolt} electrons, a propagation time of \SI{2}{\nano\second} is equivalent to a propagation distance of roughly \SI{20}{\centi\meter}. Reproduced with permission from the Americal Institute of Physics.](images/siwick2002_fig6.png){#fig:introduction-bunch-correlation}

It is precisely because the correlation is so strong that a solution can be devised. Siwick, Luiten, and collaborators[@Van2007] showed with advanced charge-particles simulations how an electromagnetic cavity operating in the radio-frequency regime could be used to reverse the correlation between the axial position and velocity of electrons within the bunch. The procedure is shown in @fig:introduction-rf-compression.

![Timeline of the compression of an electron pulse as it propagates to the sample. In the top row, the electron pulse is shown in real-space as it propagates through the compression cavity towards the sample. In the bottom row, the phase-space representation of the bunch is shown like in @fig:introduction-bunch-correlation.](diagrams/rf-compression.pdf){#fig:introduction-rf-compression}

## Experimental apparatus {#sec:experimental_setup}

All ultrafast electron scattering experiments presented in this dissertation were performed with the same experimental apparatus, presented in @fig:apparatus.

![Experimental setup diagram for the ultrafast electron scattering instrument. BS1/2: Beam splitter. HPD: high-bandwidth photodiode. BBO: $\beta$-Barium borate crystal. DM: Dichroic mirror. BS1/2: Beamsplitter. L1/2/3: Focusing lens. EL1/2: Electron lens. APD: alignment photodiode. HPD: High-bandwidth photodiode. DL: Driving loop. PL: Pickup loop.](diagrams/ued_setup.pdf){#fig:apparatus}

The heart of the apparatus consists of the ultrafast laser system. A mode-locked oscillator (Spectra-Physics Tsunami) generates short (\SI{35}{\femto\second}) pulses of \SI{800}{\nano\meter} light, at a pulse repetition rate of \SI{75}{\mega\hertz}. Each oscillator pulse carries \SI{5}{\nano\joule} of energy, which is not enough to generate bright electron pulses required by ultrafast electron scattering experiments. Five percent of the oscillator pulse train is focused down to a \SI{35}{\micro\meter} spot on a photodiode with large bandwidth (\SI{12}{\giga\hertz}), generating a master clock signal with non-negligible spectral power up to \SI{3}{\giga\hertz}. The remaining ninety-five percent of the oscillator pulse train is used to seed a regenerative laser amplifier (Spectra-Physics Spitfire Pro). The oscillator pulses are amplified via chirped-pulse amplification [@Strickland1985], outputting \SI{35}{\femto\second} pulses each carrying \SI{3}{\milli\joule} of energy, at a repetition rate of \SI{1}{\kilo\hertz}. The amplified pulse train is split between two branches: pump and probe. 

On the pump line, the laser pulses propagate to a delay stage (Newport ILS200CCHA) on which a retroreflector is installed. This delay stage provides a resolution of \SI{300}{\nano\meter}, with total travel of \SI{20}{\centi\meter}. These dimensions translate to a time-delay resolution of $\frac{2 \times \SI{300}{\nano\meter}}{c} = \SI{2}{\femto\second}$ and temporal range of $\frac{2 \times \SI{20}{\centi\meter}}{c} =\SI{1.3}{\nano\second}$. The energy of the pump pulses is then attenuated by using a neutral-density filter, to the desired pump energies. A lens is also used to focus the pump pulses onto the sample. The attenuation and focus of the beam are adjusted to achieve the desired photoexcitation density (or *fluence*) on an area that is much larger than the sample, ensuring uniform illumination conditions. The attenuated and focused pump pulses are routed inside the vacuum chamber through a z-cut quartz (SiO$_2$) window, which is transparent over a large range of wavelengths. Inside the vacuum, two mirrors is positioned so that the pump pulses are transmitted through the sample and routed outside to a photodiode, for alignment purposes.

On the probe line, two barium borate (BBO) crystals are used to generate the third harmonic of \SI{800}{\nano\meter}. A first conversion to \SI{400}{\nano\meter}, and then a second conversion to \SI{266}{\nano\meter}. A calcite crystal is placed between the two BBO crystals to compensate for the wavelength-dependent dispersion in the crystals. This non-linear light conversion is possible because the pulses are extremely energy-dense, resulting in a conversion efficiency of 1%. The ultraviolet (UV) photons are then routed into the electron gun assembly, where they are focused on a copper hemispherical cathode. The work function of copper (\SI{4.5}{\electronvolt}[@Anderson1949]) is slightly lower than the UV photon energy (\SI{4.65}{\electronvolt}), resulting in electrons being freed from the copper cathode with little extra energy. These electrons are accelerated via a static potential of \SI{90}{\kilo\volt} to 53% of the speed of light. After acceleration, a solenoid lens loosely focuses the electron bunches through the \SI{3}{\milli\meter} aperture of a radio-frequency (RF) compression cavity, ensuring that electron bunches do not lose much charge before compression. Electron bunches are compressed using the cavity described in @sec:cavity, amplifying the 40^th^ harmonic of the oscillator reprate (\SI{2.9985}{\giga\hertz}). The RF compression power is tuned so that the temporal focus happens at the sample. After the RF compression cavity, another solenoid lens focuses the transverse profile of the electron bunches onto the detector, so that the diffraction pattern is imaged clearly. Electron bunches are transmitted through the sample at normal incidence. The scattering pattern is collected by an electron camera (Gatan Ultrascan 895), on a cooled charge-coupled detector. The transmitted, unscattered beam is collected by a Faraday cup. A Keithley 6514 electrometer measures the charge on the Faraday cup at a rate of \SI{1}{\kilo\hertz}, giving a rough estimation of the bunch-to-bunch charge.

## Overview of the dissertation {#sec:overview}

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
