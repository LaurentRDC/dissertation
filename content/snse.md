
# On the origins of thermoelectricity in tin selenide {#sec:snse}

## High-performance thermoelectric materials

## Tin selenide

```{.matplotlib file="figures/snse/bz.py" caption="In-plane section of the Brillouin zone of SnSe."}
```

### Lattice instability

## Experimental methods

### Sample preparation

Sample preparation involves two steps. Bulk crystals of SnSe were first synthesized, followed by processed by which ultrathin samples were produced.

Bulk SnSe crystals (\SI{20}{\gram}) were synthesized by mixing appropriate ratios of high purity starting materials (Sn chunk, 99.999\%, American Elements, USA and Se shot, 99.999\%, 5N Plus, Canada) in \SI{13}{\milli\meter} diameter quartz tube. The tube was flame-sealed at a residual pressure of $\SI{1e-4}{\mmHg}$, then slowly heated to \SI{1223}{\kelvin} over \SI{10}{\hour}, soaked at this temperature for 6h and subsequently furnace cooled to room temperature. The obtained ingot was crushed into powder and flame-sealed in a quartz tube, which was placed into another, bigger, flame-sealed quartz tube. A crystal with dimensions of $\sim$\SI{13}{\milli\meter} (diameter) $\times$ \SI{20}{\milli\meter} (length) was obtained.

In order to obtain ultrathin samples suitable for ultrafast electron scattering experiments, two methods were tried. Ultramicrotomy yielded suitable samples, but the author wanted to rule out the effect of strain induced by the sample preparation. To this end, a sample was also prepared via mechanical exfoliation.

#### Ultramicrotomy

Six samples were prepared via ultramicrotomy, a sample preparation technique which involves the cutting of samples using a diamond blade. While this technique has mostly been used to prepare organic samples for cryo-electron microscopy, it has been successfully used in the past decade to prepare samples of two-dimensional materials such as 4H-TaSe$_2$[@Erasmus2012], 1T-TaS$_2$[@Eichberger2013]. The author initially tried ultramicrotomy to prepare samples of 1T-TiSe$_2$[@Otto2020].


The stages of sample preparation via ultramicrotomy are shown in @fig:snse-ultramicrotomy. First, a bulk crystal of SnSe is embedded in epoxy, with the cutting plane parallel to the material layers. The surface of the bulk crystal is trimmed with a \ang{45} diamond blade to reveal fresh cutting surface. Sections are then cut with a sharper, \ang{35} diamond blade which then fall in a small water container. Floating sections are then fished out with a carbon-coated TEM grid. 

```{.matplotlib #fig:snse-ultramicrotomy file="figures/snse/sample.py" caption="Stages of sample preparation via ultramicrotome. **a)** Prism of SnSe embedded in epoxy.  **b)** 4x magnification on the cut surface shown in a) shows large crystalline plateaus. **c)** Surface of embedded SnSe prism after trimming with a \ang{45} diamond knife. **d)** \SI{60}{\nano\meter} section of SnSe cut with a \ang{35} diamond knife. "}
```

The images in @fig:snse-ultramicrotomy were provided by H. Gnaegi from Diatome, Ltd, who initially tried to prepare samples for the author. The samples used in this work were prepared by J. Mui from the Facility for Electron Microscopy Research at McGill University. In the end, six samples were produced: three with a thickness of \SI{70}{\nano\meter} and three with a thickness of \SI{90}{\nano\meter}, each with an area of approximately \SI{200 x 200}{\micro\meter}.

#### Mechanical exfoliation

An ultrathin flake of SnSe was prepared via mechanical exfoliation, a procedure that is analogous to the work by Novoselov *et al*[@Novoselov2004]. A small chunk of SnSe was embedded in CrystalBond glue, on a standard \SI{3}{\milli\meter} copper TEM grid with a line spacing of 200 lines per inch. The embedded chunk was then exfoliated repeatedly using ordinary adhesive tape, until the embedded flake was translucent when observed with an optical microscope. The glue was washed away with acetone. This procedure resulted in a sample with an area of approximately \SI{50 x 50}{\micro\meter}.

The sample thickness was determined by taking the ratio of various reflections at thicknesses of \SI{70}{\nano\meter} and \SI{90}{\nano\meter} -- correcting for sample volume and electron flux -- and extrapolating to the intensity of the reflections in the exfoliated sample. Using this procedure, the thickness of the exfoliated sample was estimeted to be $\SI{30 \pm 5}{\nano\meter}$.

### Time-resolved terahertz spectroscopy

![Experimental diagram of the THz time-domain spectroscopy experiments. BS1/2: Beamsplitter. L1/2: Focusing lens. BBO: $\beta$-Barium borate crystal. PM1/2/3/4: Parabolic mirror. APD: Avalanche photodiode.](diagrams/thz_setup.pdf)

### Data acquisition

The experiments presented in this chapter used the same experimental geometry that is presented in @sec:experimental_setup. Ultrashort laser pulses of \SI{1.55}{\electronvolt} light were shone on the sample surface, with an incident angle of \ang{10}, at $t=t_0$, on SnSe samples oriented in the $[100]$ direction. To ensure that the samples had enough time to cool down after every laser shot, the repetition rate of experiments were varied from \SIrange{50}{1000}{\hertz}, but no changes were measured beyond the degradation of signal-to-noise. Therefore, a \SI{1000}{\hertz} repetition rate was used. 

Compressed electron bunches of $10^6$ electrons per bunch were transmitted through the sample at $t=t_0 + \tau$, where the time-delay $\tau$ was scanned from \SIrange{-10}{30}{\pico\second}. The total range of time-delay $\tau$ was limited compared to graphite (\SIrange{-40}{680}{\pico\second}) because the diffuse signals are much weaker in SnSe; limiting the total range of time-delay allowed for more averaging. The flagship measurement presented in this chapter was taken over \SI{72}{\hour}. This was only possible thanks to the advancements in laser-RF synchronization brought by work by the author in Otto *et al*[@Otto2017], which completely eliminated the drift in $t_0$ over the experiment duration.

The samples were photoexcited with a pump spot with a full-width at half-maximum that was at least twice the width of the sample, ensuring nearly-uniform illumination of the sample. The samples were photoexcited with photoexcitation densities ranging from \SIrange{6.6}{13.2}{\milli\joule\per\square\centi\meter}. The absorbed energy will be discussed further in this chapter. The scattering patterns were collected with a Gatan Ultrascan 895 camera, consisting of a 2.54 × 2.54 \si{\square\cm} phosphor screen fiber coupled to a 2048 px × 2048 px charge-coupled detector (CCD) placed \SI{29.39}{\centi\meter} away from the sample. Example static diffraction patterns are shown in 

```{.matplotlib #fig:snse-diff-static file="figures/snse/diff-static.py" caption="Comparison of static diffraction patterns from samples prepared via two techniques. **a)** \SI{90}{\nano\meter}-thick sample prepared via ultramicrotome. **b)** \SI{30}{\nano\meter}-thick samples prepared via mechanical exfoliation."}
```

Contrary to the symmetrization procedure described in @sec:graphite-data-acquisition and shown in @fig:graphite-static, SnSe only has a two-fold symmetry in-plane. Therefore, the data presented in this chapter was *not* symmetrized.

## Ultrafast electron scattering measurements

```{.matplotlib #fig:snse-diffuse file="figures/snse/diffuse.py" caption=""}
```

```{.matplotlib file="figures/snse/highsym.py" caption="Comparison of the diffuse intensity dynamics at various in-plane high-symmetry points."}
```

```{.matplotlib file="figures/snse/widths.py" caption="Dynamics of the width and position of various Bragg peaks following photoexcitation. For every time-delay, Bragg peaks were fit with a Gaussian function. In the right column, the change in full-width at half-maximum $\Delta \sigma$ is shown over time. In the left column, the absolute shift in the center position of the peak $\Delta x_c$ is shown, as a percentage of the average full-width at half-maximum $\bar{\sigma}$. For all plots, the error bars represent the covariance of fit parameter."}
```

### Diffuse scattering measurements

```{.matplotlib file="figures/snse/forbidden.py" caption="Comparison of the relative intensity change $\Delta I / I_0$ at $\Gamma$, for the integration geometry described in @fig:snse-diffuse. **a)** Diffuse intensity dynamics near reflections allowed by the $Pnma$ space group. **b)** Diffuse intensity dynamics near reflections that are forbidden by the $Pnma$ space group. These reflections probably arise from the double scattering of allowed reflections."}
```

### Effect of charge-transfer on diffracted intensity

### Evolution of Bragg peak profile

### Ultrafast phonon softening from Bragg intensity

```{.matplotlib file="figures/snse/softening.py" caption="Measurement of the softening of the TO$_c$ mode, extracted directly from the transient Debye-Waller dynamics, assuming that atomic vibrational amplitude only changes based on TO$_c$ renormalization at early times ($\tau < \SI{1}{\pico\second}$). Color bar shows associated absorbed energy per cell $E_c$. \textbf{Inset} Increase in isotropic mean-square-displacement of all atoms, due to the change in vibrational frequency of the TO$_c$ mode exclusively. Boxes are used to represent error bars along both axes."}
```

## Spectroscopic measurements

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
