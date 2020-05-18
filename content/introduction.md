# Introduction {#sec:introduction}

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

```{#fig:mode-locking .matplotlib caption="Demonstration of how phase relationship between amplitude modes in a laser cavity can lead to a strong pulse. **a)** The first six longitudinal modes of a laser cavity. **b)** Combination of the first 30 modes of the cavity creates a very strong pulse in the center of the cavity."}
from skued import spectrum_colors

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(FIGURE_WIDTH, 4), sharex=True)

x = np.linspace(-0.5, 1.5, 256)
pulse = np.zeros_like(x)
mode = lambda x, n: np.sin(np.pi*n*(x-1/2) + np.pi)

colors = reversed(list(spectrum_colors(5)))
for n, c in enumerate(colors, start=1):
    ax1.plot(x, mode(x,n), color=c)
    pulse += mode(x,n)

for n in range(6, 30):
    pulse += mode(x,n)

ax2.plot(x, pulse/pulse.max(), '-k')

for ax in ax1, ax2:
    ax.xaxis.set_ticks([])

ax1.yaxis.set_ticks([-1, 0, 1])
ax2.yaxis.set_ticks([0, 1])

ax2.set_xlabel("Cavity position (a.u.)")
ax1.set_ylabel("E-field amplitude (a.u.)")
ax2.set_ylabel("E-field amplitude (a.u.)")

tag_axis(ax1, "a)")
tag_axis(ax2, "b)")
```

### Regenerative laser amplifiers


## Ultrafast spectroscopic techniques and their limitations {#sec:spectroscopy}

## Electron microscopy {#sec:microscopy}

## Ultrafast electron scattering {#sec:scattering}

## Overview of the dissertation {#sec:overview}

\addcontentsline{toc}{section}{References}
\bibliographystyle{unsrtnatclean}
\bibliography{references}