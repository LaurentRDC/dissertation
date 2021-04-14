
# Introduction {#sec:introduction}

The basic physical principles of low-energy quantum mechanics are understood reasonably-well, and yet the emergence of new phenomena from the fundamental laws of Physics continues to surprise.

The field of condensed matter physics tackles a variety of intertwined questions, from fundamental physics to concrete applications. Given that condensed matter systems are complex, many experimental dimensions are used to distinguish between phenomena. Electron microscopy comes to mind as an example of a technique which can selectively observe physical systems in real-space (via imaging) or in momentum (via diffraction).

TODO: expand here

## Outrunning temperature to resolve atomic dynamics

In highly-ordered systems (i.e. crystals), where electronic correlation effects are most likely to matter, the time scales associated with fundamental actions are determined by the energy scale around room temperature. Phonons populated at room temperature (\SI{300}{\kelvin}) might have a period shorter than \SI{1}{\pico\second} ($10^{-12}$ \si{\second}), as $\omega = k_B T / \hbar = \SI{25}{\milli\electronvolt} / \hbar = \SI{6}{\tera\hertz}$. 

But more important than the absolute energy scale is the energy flow across degrees of freedom. The canonical example of this is conventional superconductivity[@Cooper1956;@Bardeen1957;@Schrieffer1963]. Energy exchange near the Angstrom length scale (\SI{1}{\angstrom} = \SI{1e-10}{\meter}) governs the onset of superconductivity as a coupling between lattice waves and electrons[@Eliashberg1960]. 

In general, the competition between coupling phenomena in a condensed matter system may be revealed as phase transitions. Some of these phases are metastable, hinting at complex networks of energy flow that gives rise to the systems traditionally studied at equilibrium. One such example is an insulator-to-metal phase transition in Pr$_{0.7}$Ca$_{0.3}$MnO$_3$, triggered by directly populating an IR-active phonon mode [@Rini2007]. By photoexciting the sample with the right photons, resistivity in the near-infrared range (\SIrange{0.5}{1.9}{\electronvolt}) drops by 5 orders of magnitude. This state lives for quite a long time in microscopic terms (nanoseconds), and therefore hints at a competition which might also be present at equilibrium.

Another exotic consequence of interactions are charge-density waves phases. These phases are characterized by a spatial distribution of electron density with a characteristic wavelength, which can be commensurate[@Rice1975;@Hellmann2010;@Erasmus2012] or incommensurate[@Mcmillan1976;@Storeck2020] with the atomic lattice dimensions. There is no unifying principle behind the formation of charge-density waves. In the excitonic insulator titanium diselenide, for example, a transition from a commensurate charge-density wave phase below \SI{190}{\kelvin} to a normal phase is accompanied by the freezing of a zone-boundary optical phonon. The appearance of the charge-density wave state is inextricably linked to the order parameter of the phase transition[@Kogar2017;@Kaneko2018;@Otto2020], indicating the importance of interactions between charge-carriers and lattice waves.

On the other hand, ultrafast methodologies are able to *rule out* apparent couplings based on time-scales. Consider vanadium dioxide, a Mott insulator which undergoes a phase transition between semiconducting below \SI{340}{\kelvin} to metallic above[@Morin1959]. This transition in electronic properties is concomitant with a structural phase transition, from a low-temperature monoclinic lattice to a high-temperature rutile lattice[@Goodenough1971]. It has been shown via a combination of structural (ultrafast electron diffraction) and spectroscopic probes (ultrafast infrared reflectivity and time-resolved terahertz spectroscopy) that the structural phase transition is distinct from the insulator-to-metal transition as a metastable metallic state with low-symmetry exists[@Morrison2014;@Otto2019]. In this case, it is the ability to outrun thermalization that revealed how vanadium dioxide's electronic properties and lattice structures are not so coupled after all.

The power of ultrafast methodologies lies in their ability to deconvolve interactions in the time-domain[@Sentef2013]. By definition, measurements at equilibrium are unable to determine energy flow and resulting couplings, as equilibrium is defined by the balanced flow of energy throughout a system. 

### An analogy with Green's functions

Ultrafast methodologies are a way to lift the veil and peer at an aspect of the system's impulse response. In this sense, an analogy with Green's functions is appropriate.

Consider a simple harmonic oscillator, with its position represented by $x$. The behavior of the oscillator is determined by the following differential equation:
$$
    \frac{\partial^2 x}{\partial t ^2} + \omega^2 x = F(t)
$${#eq:introduction-sho}
where $\omega$ is the natural oscillation frequency and is $F(t)$ is an arbitrary driving force. The response of the system $x(t)$ will depend on the driving force $F(t)$; trying to understand this particular system is futile because the intuition might not be transferable to a system with a different driving force $F^\prime(t)$.

Now consider the response of the system to a very specific driving force: the impulse $\delta(t - t_0)$. $\delta( t - t_0)$ is the Dirac delta distribution, the infinite impulse at $t_0$. The response of the system in @eq:introduction-sho to an impulse is special. Assume that the solution to 
$$
    \frac{\partial^2 x}{\partial t ^2} + \omega^2 x = \delta(t - t_0)
$$
is known to be $G(t, t_0)$. Then, it can be shown that the solution to an arbitrary driving force $F(t)$, $x_F(t)$, is given by:
$$
    x_F(t) = \int_{-\infty}^{t}dt_0 ~ G(t, t_0) F(t_0)
$${#eq:introduction-green-sho}

In the general case, consider a physical system whose state is represented by the phase-space vector $\vect{x}$, and its behavior is governed by the linear differential operator $\hat{L}$:
$$
    \hat{L} \vect{x}(t) = F(\vect{x}, t)
$${#eq:introduction-l-system}
where the function $F(\vect{x}, t)$ is a driving force, or *source term*. Given the impulse response of the system, $G(\vect{x}, \vect{x}_0, t, t_0)$, which solves
$$
    \hat{L} ~ G(\vect{x}, \vect{x}_0, t, t_0) = \delta(\vect{x} - \vect{x}_0, t - t_0),
$$
the solution $\vect{x}_F(t)$ can be expressed as a function of $G(\vect{x}, \vect{x}_0, t, t_0)$:
$$
    \vect{x}_F(t) = \int dt_0 d\vect{x}_0 ~ G(\vect{x}, \vect{x}_0, t, t_0) F(\vect{x}_0, t_0)
$$
The impulse response $G$ is called the Green's function. Here's the lesson: the response to the system to any perturbation is *completely determined* by the Green's function.

Why is this relevant? Real physical systems in condensed matter physics are governed by a linear differential operator, the Schrödinger operator:
$$
    \hat{L} = i \hbar \frac{d}{dt} - \hat{H}
$$
where $\hat{H}$ is the Hamiltonian of the (non-driven) system. The response of the total wavefunction $\Psi$ under some external potential $\hat{V}_0$ becomes:
$$
    \hat{L} \Psi = \hat{V}_0 \Psi
$$
which is equivalent to the usual Schrödinger equation $i \hbar \frac{d\Psi}{dt} = (\hat{H} + \hat{V}_0)\Psi$. The Green's function of the Schrödinger operator is sometimes called the propagator[@Sakurai2014Propagators]. In practice, the total wavefunction $\Psi$ is intractable, with more than $10^{23}$ degrees of freedom at the macroscopic scale. Nonetheless, this system is completely determined by the response to an energy impulse, i.e. the Green's function of $\hat{L}$. In fact, the Green's function formalism is at the heart of many-body theory[@Fetter1971GreensFunction;@Sentef2013]. 

Ultrafast measurements allow us to measure an approximation of this impulse response. Distinct excitation conditions and probing techniques explore separate parts of the Hilbert space of possibilities. With the right probes, a shadow of the true impulse response may be assembled to understand important physics which -- while not experimentally-accessible at equilibrium -- still govern macroscopic phenomena. The knowledge of the impulse response of materials can then inform on the behavior of systems under a variety of external conditions $\hat{V}_0$ such as thermal gradients, external electric and magnetic fields, stoichiometric doping, and pressure.

### The reversible pump-probe scheme

The time-resolution with gives experimentalists a glimpse of the impulse response of physical systems is derived from ultrafast laser pulses. These pulses provide the impulse of energy which perturbs the system. However, digital measurement devices are generally limited to nanosecond resolution, and are therefore too slow to observe dynamics in the femtosecond range.

The pump-probe scheme is used to achieve high time-resolution in experiments where detectors are slow. This scheme involves two pulses. One pulse, the *pump pulse*, is an ultrashort laser pulse which impulsively dumps energy in a sample. The sample is pump with a period $T$, which is generally fixed by the laser source. The period $T$ must be great enough to allow the sample to recover. A *probe pulse* is made to interact with the sample after some time-delay $\tau$. This probe pulse might be another ultrashort pulse of light, or electrons in the case of ultrafast electron scattering. The delay can usually be adjusted with great precision by changing the path length of the light used to generate the probe pulse. The probe pulse will encode the state of the sample at time-delay $\tau$ within the width of the probe pulse.

The process will be repeated $n$ times, where $nT$ is a long-enough delay that the detector of the probe pulse can resolve (e.g. \SI{1}{\second}). Then, the time-delay will be adjusted to a new value, and the process is repeated. By scanning over values of $\tau$ in this way, a complete view of the dynamics of the sample following the pump can be assembled. This process is represented in @fig:introduction-pumpprobe.

```{.matplotlib #fig:introduction-pumpprobe file="figures/introduction/pumpprobe.py" caption="Representation of the pump-probe scheme in laboratory time. On the right, dynamics are initiated with a period $T$, and the system is allowed to relax back to its initial state. After a time-delay, the system is probed with a separate probe pulse. The process is repeated multiple times while the detector is recording. On the left, the sampled dynamics are shown. The solid trace is the true response, while the grey circles show the measurements convolved with the instrument time resolution. By scanning the time-delay ($\tau_1$, $\tau_2$, $\tau_3$, ...), a complete picture of sample response can be assembled."}
```

### Exploring the lattice impulse response with ultrafast electron scattering

The impulse response of a macroscopic system is an intractable quantity with more than with $10^{23}$ dimensions. Even in the case of non-interacting systems of electrons and ions, the impulse responses of each subsystem is too large to characterize completely. The problem is greatly reduced by studying highly periodic crystalline systems, where instead the response in the spatial-frequency domain is studied. The spatial-frequency-domain impulse response is reduced in that the degrees-of-freedom which are related by symmetry are equivalent.

The description of the impulse response of ordered systems like crystals is fundamentally-different. The relevant physics of $10^{23}$ particles, $10^{22}$ of which are related by symmetry, is better described using *collective modes* such as lattice waves, also known as phonons. The instantaneous geometry of the crystal lattice determines the allowed phonon modes, and the energy distribution among modes represent a generalization of the idea of temperature at equilibrium. 

While the lattice impulse response in terms of lattice modes represent a large part of the total impulse response of a crystalline system, few techniques are capable of resolving such lattice waves. This is due to two main reasons. First, lattice waves are a mostly low-energy phenomena ($<\SI{200}{\milli\electronvolt}$). Second, lattice waves may not induce a dipole moment, which makes some of the modes optically-dark. At equilibrium, multiple techniques are either sensitive to some modes (e.g. Raman[@Saito2001] and Brillouin[@Bai2018] spectroscopies) or all modes (x-ray[@Holt1999;@Xu2005] and neutron[@Hempelmann2000] thermal diffuse scattering). But at ultrafast time-scales, where modes can be effectively be populated at any wavevector, the choices are much more restrictive. Signal-to-noise in ultrafast electron scattering experiments is now sufficient to measure the effects of lattice waves following impulsive photoexcitation[@Harb2016;@Chase2016;@Waldecker2017;@Stern2018;@RenedeCotret2019;@Otto2020]. 

### Electron-phonon interactions{#sec:introduction-epc}

As stated above, electron-phonon interactions play an important role in many interesting systems. It is therefore no surprise that electron-phonon interactions are central to the understanding of the lattice response to ultrafast photoexcitation, as most of the optical pump energy is absorbed by charge-carriers. The photoexcitation energy must dissipate as lattice waves, as acoustic phonons are some of the lowest energy excitations. Therefore, a complete understanding of systems away from equilibrium inherently involves electron-phonon interactions. In this section, a basic description of electron-phonon interactions is given, and more details are presented when needed in @sec:graphite and @sec:snse.

The Hamiltonian which describes a system of electrons and ions at zero-temperature is given by:
$$
    \hat{H} = \hat{H}_{ph} + \hat{H}_e + \hat{H}_i
$$
where $\hat{H}_{ph}$ is the Hamiltonian of the phonon subsystem, $\hat{H}_e$ is the Hamiltonian of the electronic subsystem, and $\hat{H}_i$ represents the interaction between ions and electrons[@Mahan2000EPInteractions]. The interaction Hamiltonian can be conceptually-expressed as:
$$
    \hat{H}_i = \sum_{j, j^\prime} \hat{V}_i(\vect{r}^e_j - \vect{r}^i_{j^\prime})
$$
where $\vect{r}^e_j$ and $\vect{r}^i_{j^\prime}$ are the positions of electron $j$ and ion $j^\prime$, respectively. The interaction potential $\hat{V}_i$ is kept general for this discussion. Now consider that the atomic positions $\vect{r}^i_{j^\prime}$ are perturbed by phonons, which results in a (small) displacement $\vect{u}_{j^\prime}$ such that $\vect{r}^i_{j^\prime} \to \vect{r}^i_{j^\prime} + \vect{u}_{j^\prime}$. In this scheme, $\vect{r}^i_{j^\prime}$ is now the position of ions at zero-temperature. To first order, the change in interaction strength is given by the gradient of the interaction potential:
$$
    \hat{V}_i(\vect{r}^e_j - \vect{r}^i_{j^\prime} - \vect{u}_{j^\prime}) \approx
        \hat{V}_i(\vect{r}^e_j - \vect{r}^i_{j^\prime}) - \vect{u}_{j^\prime} \cdot \nabla \hat{V}_i(\vect{r}^e_j - \vect{r}^i_{j^\prime})
$$
The first term is the interaction of electrons with ions at their ground-state positions. The second term may be intuitively named the *electron-phonon interaction*:
$$
    \hat{V}_{ep}(\vect{r}) \equiv \sum_{j^\prime} \vect{u}_{j^\prime} \cdot \nabla \hat{V}_i(\vect{r} - \vect{r}^i_{j^\prime})
$$
The associated interaction Hamiltonian is given by:
$$
    \hat{H}_{ep}(\vect{r}) = \int d\vect{r} \hat{\rho}(\vect{r}) \hat{V}_{ep}(\vect{r})
$$
where $\hat{\rho}(\vect{r})$ is the electron density operator. For highly-ordered crystalline systems like the ones studied in this dissertation, it is helpful to express the interaction Hamiltonian in momentum space. In momentum space, $\hat{H}_{ep}$ is expressed most simply in the second quantization framework[@Giustino2017]:
$$
    \hat{H}_{ep} = \sum_{a,b, \lambda} \int \frac{d\vect{k}_e}{(2 \pi)^3} \int \frac{d\vect{k}_p}{(2 \pi)^3} ~ g^{\lambda}_{ab}(\vect{k}_p, \vect{k}_e) ~ \hat{c}^{\dagger}_{a}(\vect{k_e} + \vect{k}_p)\hat{c}_b(\vect{k}_e) \left( \hat{a}_{\lambda}(\vect{k}_p) + \hat{a}^{\dagger}_{\lambda}(-\vect{k}_p)\right)
$$
where $a$ and $b$ label electron bands, $\lambda$ labels phonon modes, $\hat{c}^{\dagger}_{a}$ ($\hat{c}_{a}$) is the creation (annihilation) operator for electrons in band $a$, $\hat{a}^{\dagger}_{\lambda}$ ($\hat{a}_{\lambda}$) is the creation (annihilation) operator for phonon mode $\lambda$, $\vect{k}_e$ and $\vect{k}_p$ are the electronic and phonon wavevectors respectively, and $g$ is the electron-phonon coupling matrix. 

The electron-phonon coupling matrix is effectively impossible to measure at equilibrium. Until recently, the determination of this coupling matrix was the purview of calculations. Ultrafast measurements are able to get at $g$ via the lifetime of excitations[@Brorson1990;@Ferrari2007;@Na2019]. Ultrafast electron diffuse scattering remains the only laboratory-scale technique where electron-phonon coupling matrix elements are accessible with momentum-resolution[@@RenedeCotret2019].

## Brief history of ultrafast electron scattering {#sec:ues}

The extension of electron microscopy into time-domain studies is not new. This section will be limited to ultrafast electron scattering in the form of diffraction; see the review by King *et al.*[@King2005] for a historical perspective on ultrafast electron microscopy in general. As early as 1982, Gerard Mourou and Steve Williamson [@Mourou1982] followed the picosecond-scale transformation of an aluminum film following photoexcitation. In this work, the authors prepare \SI{100}{\pico\second} electron bunches using a streak camera. They note that the temporal resolution of the experiment is ultimately limited by the length of the electron bunch:

> *The electron pulse width has been measured by using the camera in the normal streak mode and is found to be $\sim$\SI{100}{\pico\second}. This value departs significantly from the \SI{15}{\pico\second} pulse width expected. The pulse broadening is due to the space-charge effect caused by the relatively high electron flux required to photograph the pattern with our present system.*

The challenges of ultrafast electron scattering on the \SI{100}{\femto\second} time-scale were apparent in the first ultrafast electron scattering experiments. Measurements of atomic structure require high beam brightness (i.e. high signal-to-noise), but dense electron bunches experience space-charge repulsion (i.e. worsening time-resolution). The trade-off between signal-to-noise and time-resolution is represented as the expansion of the bunch length $l$:
$$
\frac{d^2l}{dt^2} = \frac{N e^2}{m_e \epsilon_0 \pi r^2} \left[ 1 - \frac{l}{\sqrt{l^2 + 4 r^2}}\right]
$$ {#eq:bunch_length}
where $N$ is the number of electrons in the bunch, $e$ the electron charge, $m_e$ the electron mass, $\epsilon_0$ the vacuum permittivity, and $r$ is the beam radius[@Siwick2002]. The first demonstration of ultrafast electron scattering was done by the same authors[@Mourou1982] to study the ultrafast melting of aluminium[@Williamson1984] with a time-resolution a tens of picoseconds.

The modern form of pump-probe ultrafast electron scattering was brought by Ahmed Zewail. At the time, Zewail was already known as the "Father of Femtochemistry" for the development of ultrafast spectroscopy[@Zewail1988]. Zewail and collaborator proposed to replace the optical probe pulse used in pump-probe femtosecond spectroscopy to an electron pulse[@Williamson1991;@Williamson1997]. What followed was a series of experiments probing the transient molecular structures at the picosecond time-scale following photoexcitation in the gas phase[@Ihee2001;@Ruan2001;@Lobastov2001]. The ability to track the nuclear coordinates directly rather than indirectly represented the perfect complement to Zewail's Nobel prize-winning work on femtochemistry via ultrafast spectroscopy.

It was not until Siwick *et al.*[@Siwick2003;@Siwick2004], however, that the ideas of Mourou and Williamson were applied with sub-picosecond sensitivity (\SI{600}{\femto\second}) to the solid state. This work reached a time-resolution sufficient to address whether or not the ultrafast phase transition from solid to liquid-like was analogous to the thermal phase transition[@Guo2000].

By 2005, the prevalent sentiment was that the space-charge problem formulated by Mourou and Williamson limited the applicability of ultrafast electron scattering. King *et al.*[@King2005] identify three paths to control space-charge effects:

1. Reducing the electron-emission -- sample distance. By reducing this distance, the space-charge-driven expansion of the electron bunch remains limited[@Siwick2002].

2. Using higher electron energies. For \si{\mega\electronvolt} electrons, the space-charge explosion is reduced significantly[@Gian2002].

3. Reducing space-charge pulse expansion by limiting the energy spread in the electron bunch[@Gallant2000].

The solution to space-charge came from the insight of Siwick *et al.* [@Siwick2002], who showed that the electrons in a bunch develop a strong correlation in phase space as the bunch propagates. Space-charge effects effectively slow down the front electrons with respect to the center of charge, while it accelerates the electrons at the front of the bunch. After some propagation time, a linear correlation is established between the relative axial position within the bunch and the relative velocity within the bunch. It is precisely because the correlation is so strong that a solution can be devised. Siwick, Luiten, and collaborators[@Van2007] showed with advanced charge-particles simulations how an electromagnetic cavity operating in the radio-frequency (RF) regime could be used to reverse the correlation between the axial position and velocity of electrons within the bunch. The particular implementation in the Siwick research group is presented below in @sec:intro-cavity. The design of a RF compression cavity was validated by Chatelain *et al.*[@Chatelain2012] in 2012.

### The landscape today

At the time of writing, there are many research groups across the world using ultrafast electron scattering techniques with great success. This section is a non-exhaustive enumeration of the field as of the first half of 2021.

First, the most traditional approach of ultrafast electron diffraction in the compact geometry is still alive today. One of the most successful research groups using this tool is led by R. Ernstorfer at the Fritz-Haber Institute in Berlin which has been studying multiple solid-state systems[@Waldecker2016;@Waldecker2017;@Zahn2020]. 

There are also multiple RF-compressed ultrafast electron diffraction instruments in use today. Such an instrument is used in the author's research group led by B. Siwick at McGill University[@Chatelain2014a;@Morrison2014;@Stern2018;@RenedeCotret2019;@Otto2020]. Several instruments are located at the University of Toronto, in the group led by R. J. D. Miller[@Ernstorfer2009;@Gao2013;@Jiang2017;@Jiang2020]. Another instrument based on the design by van Oudheusden *et al.*[@Van2007] has been in use at the University of Nebraska, in the research group led by M. Centurion, with special emphasis on gas-phase electron diffraction[@Yang2015].

Finally, there are a few high-energy (\si{\mega\electronvolt}) accelerator-based ultrafast electron diffraction instruments in operation today. The most well-known is located at the Stanford Linear Accelerator Center (SLAC), part of a group led by X. Wang. Several research groups have been using this instrument as it has turned into a user-facility[@Mo2018;@Sie2019;@Kogar2020]. Another high-energy instrument is located at the University of California, Los Angeles, in a research group led by P. Musumeci[@Maxson2017] which focuses on instrument development.

## Electron bunch compressor{#sec:intro-cavity}

The compression of electron bunches has been the key technological advancement that has allowed the work presented here. It has allowed for high time-resolution while preserving a high bunch charge, which is required to observe the small diffuse signals presented in @sec:graphite and @sec:snse. The electron bunch compression that was used for the experiments presented herein has been the subject of much work and refinement in the past decade. Its function is briefly described in this section, with a particular emphasis on recent developments regarding laser-RF synchronization.

### Bunch correlations

<!-- The figure below is raw latex because an alternate caption without citations needs to be use so that it looks OK in the list of figures-->
\begin{figure}
    \centering
    \includegraphics{images/siwick2002_fig6.png}
    \caption[Electron bunch phase space as it propagates.]{Relative axial velocity vs axial position ($Z$) for all electrons in the pulse at four times ($T$) during its propagation. The pulse is composed of 10 000 electrons at \SI{30}{\kilo\electronvolt}, with an initial bunch length equivalent to \SI{150}{\femto\second}, an initial radius ut \SI{75}{\micro\meter}, and initial beam divergence of \SI{1.5}{\milli\radian}. The width of the velocity distribution becomes larger as the pulse lengthens. The spatial distribution of velocities also evolves. As electrons redistribute themselves inside the packet a linear velocity chirp develops on the electron pulse. For \SI{30}{\kilo\electronvolt} electrons, a propagation time of \SI{2}{\nano\second} is equivalent to a propagation distance of roughly \SI{20}{\centi\meter}. Reused with permission from Siwick \emph{et al.} \autocite{Siwick2002}}
    \label{fig:introduction-bunch-correlation}
\end{figure}

The idea of an electron bunch compressor is to control the electron bunch in phase space. As stated previously, a landmark paper by Siwick *et al.*[@Siwick2002] showed that as the electron bunch propagates, electrons in develop a linear correlation between their position within the bunch and their velocity with respect to the average bunch velocity. This is shown in @fig:introduction-bunch-correlation. The idea behind the RF compressor is to reverse this correlation using a standing radio-frequency wave to slow down the electrons at the front of the bunch, and accelerate the electrons at the back of the bunch, so that there is a single point downstream where the electron bunch is very short. This can be done by using an electromagnetic wave, oriented with the electric field along the propagation axis $z$ (TM$_{010}$ mode), and timed such that the amplitude of the electric field is 0 when the bunch center-of-charge is located in the center of the cavity. This procedure is effectively a rotation in phase-space. The procedure and its effect on the phase-space representation of electron bunches is shown in @fig:introduction-rf-compression. 

![Timeline of the compression of an electron pulse as it propagates to the sample. In the top row, the electron pulse is shown in real-space as it propagates through the compression cavity towards the sample. In the bottom row, the phase-space representation of the bunch is shown like in @fig:introduction-bunch-correlation.](diagrams/rf-compression.pdf){#fig:introduction-rf-compression}

The visualization in @fig:introduction-rf-compression also shows the importance the timing of the compression field. If the field is not timed properly, a net momentum kick will be imparted on the electron bunch, which will change the arrival time of the electrons at the sample. This is equivalent to a *translation* in phase space. This is discussed in @sec:intro-cavity-timing.

### Compressor design

The general idea of phase-space rotation may be realized in multiple ways. Some design considerations that explain the particularities of the RF cavity used in this work are addressed here.

![Machining diagram of the RF compression cavity used in this work.](diagrams/cavity.pdf){#fig:intro-cavity-design}

Consider the choice of the wavelength of the compression electric field. The main requirement here is that the electron bunch experiences a field that varies linearly as it travels through the cavity. Waves in the \SIrange{0.5}{6}{\giga\hertz} range are required given that uncompressed electron bunch with $10^6$ electrons may be $\sim$\SI{100}{\pico\second} long. The choice of \SI{3}{\giga\hertz} was ultimately chosen because of previous work in synchronizing \SI{3}{\giga\hertz} RF cavities to \SI{800}{\nano\meter} laser pulse trains[@Kiewiet2002].

The shape of the cavity is also important. In its simplest form, a "pillbox" cavity (hollow cylinder) of appropriate dimensions would support the right RF compression wave[@Fill2006]. However, pillbox cavities are rather inefficient in terms of field-strength per Watt of input power. The final lobe design reduced input power requirements by 90\% compared to the pillbox geometry[@Van2007]. The design of the cavity as manufactured is shown in @fig:intro-cavity-design, modified from Chatelain[@Chatelain2014].

### Compressor characterization

Initial characterization of the electron compressor in 2012 was performed with a vector network analyzer up to \SI{3}{\giga\hertz} only. The author characterized the cavity again in December 2019 across a much larger bandwidth of \SI{50}{\giga\hertz}, using an Agilent Technologies N5247A vector network analyzer. The transmission spectrum is presented in @fig:introduction-cavity-full-spectrum. The main resonance near \SI{3}{\giga\hertz} is clearly visible.

```{.matplotlib #fig:introduction-cavity-full-spectrum file="figures/introduction/cavity-full-spectrum.py" caption="Characterization of the cavity transmission spectrum at \SI{13}{\celsius}. **a)** Transmitted amplitude near \SI{3}{\giga\hertz}. **b)** Phase of the transmitted signal near \SI{3}{\giga\hertz}. **c)** Amplitude spectrum in a wide band shows multiple resonances up to \SI{12}{\giga\hertz}."}
```

The temperature-dependence of the resonant frequency was also investigated. By changing the temperature of the cavity via liquid cooling, its dimensions are modified slightly due to thermal expansion. The change in the transmission spectrum near \SI{3}{\giga\hertz} is presented in @fig:introduction-cavity-temp-sweep.

```{.matplotlib #fig:introduction-cavity-temp-sweep file="figures/introduction/cavity-temp-sweep.py" caption="Effect of cavity temperature on resonant frequency near \SI{3}{\giga\hertz}. **a)** Transmission function from \SIrange{284}{384}{\kelvin} in steps of \SI{1}{\kelvin}. **b)** Temperature dependence of the resonant frequency shows a linear dependence with a slope of \SI{-52.6 \pm 0.2}{\kilo\hertz\per\kelvin}."}
```

### Driving field generation and timing considerations{#sec:intro-cavity-timing}

The timing of the RF compression field is of paramount importance. As stated previously, jitter or drift in the phase of the field will impart a net momentum kick to electron bunches, which will in turn deteriorate the time-resolution of experiments. The zero-crossing of the sinusoidal compression field needs to be synchronized to the moment the center-of-charge of the bunch passes through the center of the cavity. If the zero-crossing is too early, more of the bunch will be accelerated, resulting in a net momentum kick. If the zero-crossing is too late, more of the bunch will be decelerated, resulting in a net momentum loss. In both cases, the focus (in time) of the compressor will not coincide with the sample position. This is shown in @fig:intro-zero-crossing. 

![Effect of mistiming of electron compression on the electron bunch. The zero-crossing determines the moment where the electric field is 0 in the cavity](diagrams/rf-timing.pdf){#fig:intro-zero-crossing}

Because the electron bunches are photogenerated from a laser pulse train, timing the RF compression field is effectively a problem in synchronizing an electromagnetic wave to laser pulses. The initial synchronization system used until 2017 is described in Kiewiet *et al.*[@Kiewiet2002] In broad strokes, this initial synchronization system generated its own \SI{3}{\giga\hertz} signal from an internal clock. The phase of this signal was adjusted based on the phase of laser pulse train as measured with a photodiode. This approach suffered from relatively large drift over the course of multiple hours.

Since 2017, a direct synchronization approach has been in use. This synchronization system is described in detail in Otto *et al.*[@Otto2017], and a brief summary is given here. There are two parts to this improved synchronization system. First, a \SI{3}{\giga\hertz} signal is generated directly from the \SI{75}{\mega\hertz} laser oscillator pulse train, which is then amplified to drive the electron compressor. Second, the phase of the input to the electron compressor is controlled to eliminate long-term phase drifts. These two aspects of the synchronization system are addressed in their own subsections. 

#### Direct generation of a driving field

High-bandwidth photodiodes ($>$\SI{10}{\giga\hertz}) are now available as commercial products. Their response time is good enough ($\sim$\SI{25}{\pico\second}) that the output spectrum from a laser oscillator pulse train.

```{.matplotlib #fig:intro-sync-spectrum file="figures/introduction/freq-comb.py" caption="Spectrum of a \SI{75}{\mega\hertz} oscillator pulse train on a high-bandwidth (\SI{12.5}{\giga\hertz}) photodiode, with an input power of \SI{20}{\milli\watt}. **a)** Full bandwith measurement shows harmonics up to at least \SI{6}{\giga\hertz}. **b)** Spectrum centered on the 40^th^ harmonic of the oscillator repetition rate."}
```

A master clock signal is generated by splitting 5% of the laser oscillator power (\SI{20}{\milli\watt}) and focusing it on a high-bandwidth photodiode. The harmonics of the \SI{75}{\mega\hertz} pulse train are visible in the power spectrum, as shown in @fig:intro-sync-spectrum; \SI{3}{\giga\hertz} corresponds to the 40^th^ harmonic. A combination of low-phase-noise amplifiers and a band-pass filter is used to isolate a master clock signal at a power of \SI{-2.5}{\dBm} (\SI{0.5}{\milli\watt}). At this stage, a portion of the master clock signal is split off for use in a feedback loop, described in the next section.

The master clock signal is then amplified in a commercial continuous-wave solid-state amplifier, which can provide up to \SI{60}{\deci\bel} of amplification up to \SI{200}{\watt}. Typically, an output power in \SIrange{47.5}{48.5}{\dBm} (\SIrange{56}{71}{\watt}) range was used to drive the electron compressor. To ensure that the cavity is driven at resonance, either the temperature of the cavity (see @fig:introduction-cavity-temp-sweep) or the oscillator pulse train repetition rate can be adjusted.

#### Feedback system to eliminate long-term drift

Ultrafast electron scattering experiments take multiple hours to complete; some of the experiments presented in @sec:snse ran for \SI{72}{\hour}. It is therefore important to ensure that the phase of the compressor driving field (i.e. the timing of the zero-crossing) be stable over this period. Drift in the phase of the driving field can occur due to many factors. One major source of drift is due to temperature changes of the electron compression, which induces thermal expansion/contraction of the RF cavity. The resonance of the cavity is therefore changed, and the driving field becomes detuned from the cavity resonance. Another source of drift is due to amplitude drift in the laser oscillator, which can cause a phase drift in a process called amplitude-phase conversion [@Taylor2011;@Zhang2012].

A solution to the drift problem was implemented based on a feedback loop. The RF cavity was machined with an output port which can be used to measure the RF signal transmitted through the cavity. The phase of this transmitted signal is a measure of how detuned the driving field is compared to the cavity resonance. The phase difference is then measured by comparing the two signals, and a phase shift is applied to minimize the detuning of the driving field.

### Alternative approaches

Electron bunch compression can be achieved in two orthogonal ways. The first method, one form of which was presented above, is to adjust the velocity of electrons with respect to the bunch. The other method is to change the electron path length based on their energy (velocity). The result is the same: the electron bunch is compressed in time at a specific point downstream.

Electron bunch compression via a modification of electron energies must be done with strong electric fields, but electromagnetic cavities are not necessary. Ehberger *et al.*[@Ehberger2019] have demonstrated the compression of electron bunches with \si{\tera\hertz} fields, although with very low bunch charge (\SI{1.6e-7}{\pico\coulomb}, 3 electrons per bunch). This approach is fully-optical, meaning that there is no need to synchronize the laser pulses to electronics. In fact, \si{\tera\hertz} fields can be used to generate, compress, and streak electron bunches[@Ehberger2019;@Matte2021].

Electron bunch compression via path length modulation has been proposed and implemented to varying degrees of success. The advantage of many of these approaches based on electron beam geometry is that laser-electronics synchronization is not required. Geometries include \ang{360} magnet compressors[@Tokita2010], spherical electrostatic capacitor compressors[@Grzelakowski2013], and reflectron-based mirror compressors[@Weber1995;@Kassier2009;@Wang2011]. Most recently, a scheme that combines a magnetic beam separator and a reflectron -- with better experimental alignment and spatial aberration correction -- has been proposed by Mankos *et al.*[@Mankos2017]

## Experimental apparatus {#sec:experimental_setup}

All ultrafast electron scattering experiments presented in this dissertation were performed with the same experimental apparatus, presented in @fig:apparatus.

![Experimental setup diagram for the ultrafast electron scattering instrument. BS1/2: Beam splitter. HPD: high-bandwidth photodiode. BBO: $\beta$-Barium borate crystal. DM: Dichroic mirror. BS1/2: Beamsplitter. L1/2/3: Focusing lens. EL1/2: Electron lens. APD: alignment photodiode. DL: Driving loop. PL: Pickup loop.](diagrams/ued_setup.pdf){#fig:apparatus}

The heart of the apparatus consists of the ultrafast laser system. A mode-locked oscillator (Spectra-Physics Tsunami) generates short (\SI{35}{\femto\second}) pulses of \SI{800}{\nano\meter} light, at a pulse repetition rate of \SI{75}{\mega\hertz}. Each oscillator pulse carries \SI{5}{\nano\joule} of energy, which is not enough to generate bright electron pulses required by ultrafast electron scattering experiments. Five percent of the oscillator pulse train is focused down to a \SI{35}{\micro\meter} spot on a photodiode with large bandwidth (\SI{12}{\giga\hertz}), generating a master clock signal with non-negligible spectral power up to \SI{3}{\giga\hertz}. The remaining ninety-five percent of the oscillator pulse train is used to seed a regenerative laser amplifier (Spectra-Physics Spitfire Pro). The oscillator pulses are amplified via chirped-pulse amplification [@Strickland1985], outputting \SI{35}{\femto\second} pulses each carrying \SI{3}{\milli\joule} of energy, at a repetition rate of \SI{1}{\kilo\hertz}. The amplified pulse train is split between two branches: pump and probe. 

### Pump line

On the pump line, the laser pulses propagate to a delay stage (Newport ILS200CCHA) on which a retroreflector is installed. This delay stage provides a resolution of \SI{300}{\nano\meter}, with total travel of \SI{20}{\centi\meter}. These dimensions translate to a time-delay resolution of \SI{2}{\femto\second} and temporal range of \SI{1.3}{\nano\second}. The energy of the pump pulses is then attenuated by using a neutral-density filter, to the desired pump energies. A lens is also used to focus the pump pulses onto the sample. 

The attenuation and focus of the beam are adjusted to achieve the desired photoexcitation density (or *fluence*) on an area that is much larger than the sample, ensuring uniform illumination conditions. The attenuated and focused pump pulses are routed inside the vacuum chamber through a z-cut quartz (SiO$_2$) window, which is transparent over a large range of wavelengths. Inside the vacuum, two mirrors is positioned so that the pump pulses are transmitted through the sample and routed outside to a photodiode, for alignment purposes.

### Probe line

On the probe line, two barium borate (BBO) crystals are used to generate the third harmonic of \SI{800}{\nano\meter}: a first conversion to \SI{400}{\nano\meter}, and then a second conversion to \SI{266}{\nano\meter}. A calcite crystal is placed between the two BBO crystals to compensate for the wavelength-dependent dispersion in the BBO crystals, and keep the optical pulse short. This non-linear light conversion is possible because the pulses are extremely energy-dense, resulting in a conversion efficiency of 1%. The ultraviolet (UV) photons are then routed into the electron gun assembly, where they are focused on a copper hemispherical cathode. The work function of copper (\SI{4.5}{\electronvolt}[@Anderson1949]) is slightly lower than the UV photon energy (\SI{4.65}{\electronvolt}), resulting in electrons being freed from the copper cathode with little extra energy. 

These electrons are accelerated via a static potential of \SI{90}{\kilo\volt} to 53% of the speed of light. After acceleration, a solenoid lens loosely focuses the electron bunches through the \SI{3}{\milli\meter} aperture of the RF compression cavity described in @sec:intro-cavity, ensuring that electron bunches do not lose much charge before compression. After the RF compression cavity, another solenoid lens focuses the transverse profile of the electron bunches onto the detector, so that the diffraction pattern is imaged clearly. Electron bunches are transmitted through the sample at normal incidence. The scattering pattern is collected by an electron camera (Gatan Ultrascan 895), on a cooled charge-coupled detector. The transmitted, unscattered beam is collected by a Faraday cup. A Keithley 6514 electrometer measures the charge on the Faraday cup at a rate of \SI{1}{\kilo\hertz}, giving a rough estimation of the bunch-to-bunch charge fluctuations.

### Interactive data exploration software {#sec:introduction-data-exploration}

The experimental procedure employed in the Siwick research group acquires data in a repeating fashion, in order to increase the signal-to-noise ratio. The resulting experimental raw data can be hundreds of gigabytes in size. The author has created a software suite that performs data reduction and allows for the interactive exploration of the reduced data, in real-time. While this software suite has had a profound impact on the research presented in this work, the details are not necessary to appreciate the results presented in @sec:graphite and @sec:snse. Interested readers are encouraged to read more in @sec:appendix-software.

## Overview of the dissertation {#sec:overview}

This dissertation has two goals: to push the understanding of nonequilibrium lattice dynamics as measured by ultrafast electron scattering, and to apply this understanding to elucidate the couplings in strongly-correlated materials such as graphite and tin selenide.

To this end, the dissertation starts proper with a description of the theoretical framework of ultrafast electron diffraction and diffuse scattering in @sec:scattering. This chapter holds a special place in the author's heart. At the time of writing, it is the only quantum-mechanical derivation of ultrafast diffuse scattering amplitude for either x-ray or electrons. The thorny subject of multiple electron scattering is also discussed.

Ultrafast electron diffuse scattering measurements in graphite are considered next in @sec:graphite. Graphite is a wonderful benchmark system because of the stiffness of the lattice, which gives diffuse scattering measurements incredibly high contrast. In a way, diffuse scattering in graphite is the most definitive diffuse scattering experiment.

The study of the origin of tin selenide's high thermoelectric performance is the subject of @sec:snse. In this material, ultrafast electron diffraction and diffuse scattering are used in conjunction to understand the physics of the phases of tin selenide, and suggest a pathway to enhance the thermoelectric performance of tin selenide near room temperature.

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
