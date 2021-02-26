
# On the origins of thermoelectricity in tin selenide {#sec:snse}

## High-performance thermoelectric materials

## Tin selenide

Single-crystal tin selenide at ambient pressure is found in two phases, both of which are layered[@Chattopadhyay1986]. The low-temperature $Pnma$ phase is defined by the lattice vectors $\vec{a}_1 = a ~ \vec{e}_1$, $\vec{a}_2 = b ~ \vec{e}_2$, and $\vec{a}_3 = c ~ \vec{e}_3$ where $a=\SI{11.42}{\angstrom}$, $b=\SI{4.19}{\angstrom}$, and $c=\SI{4.46}{\angstrom}$. The vectors $\set{ \vec{e}_i }$ are understood to be the usual Euclidean vectors. The atomic positions are for the $Pnma$ phase are:
$$
\begin{pmatrix}
    \vec{\text{Se}}_1 \\
    \vec{\text{Se}}_2 \\
    \vec{\text{Se}}_3 \\
    \vec{\text{Se}}_4 \\
    \vec{\text{Sn}}_1 \\
    \vec{\text{Sn}}_2 \\
    \vec{\text{Sn}}_3 \\
    \vec{\text{Sn}}_4 
\end{pmatrix}
=
\begin{pmatrix}
    0.14 & \sfrac{3}{4} & 0.52\\
    0.36 & \sfrac{1}{4} & 0.02\\
    0.64 & \sfrac{3}{4} & 0.98\\
    0.86 & \sfrac{1}{4} & 0.48\\
    0.12 & \sfrac{1}{4} & 0.09\\
    0.38 & \sfrac{3}{4} & 0.59\\
    0.62 & \sfrac{1}{4} & 0.41\\
    0.88 & \sfrac{3}{4} & 0.91
\end{pmatrix}
\begin{pmatrix}
    \vec{a}_1 \\
    \vec{a}_2 \\
    \vec{a}_3
\end{pmatrix}
$$
The lattice vectors for the *conventional cell* of high-temperature phase, $Cmcm$, are given by $\vec{a}_1 = a ~ \vec{e}_1$, $\vec{a}_2 = b ~ \vec{e}_2$, and $\vec{a}_3 = b ~ \vec{e}_3$ where $a=\SI{11.71}{\angstrom}$ and $b=\SI{4.31}{\angstrom}$. The atomic positions in fractional coordinates are:
$$
\begin{pmatrix}
    \vec{\text{Se}}_1 \\
    \vec{\text{Se}}_2 \\
    \vec{\text{Se}}_3 \\
    \vec{\text{Se}}_4 \\
    \vec{\text{Sn}}_1 \\
    \vec{\text{Sn}}_2 \\
    \vec{\text{Sn}}_3 \\
    \vec{\text{Sn}}_4 
\end{pmatrix}
=
\begin{pmatrix}
    0.14 & \sfrac{1}{2}  & \sfrac{3}{4}\\
    0.36 & 0             & \sfrac{1}{4}\\
    0.64 & 0             & \sfrac{3}{4}\\
    0.86 & \sfrac{1}{2}  & \sfrac{1}{4}\\
    0.12 & 0             & \sfrac{1}{4}\\
    0.88 & 0             & \sfrac{3}{4}\\
    0.38 & \sfrac{1}{2}  & \sfrac{3}{4}\\
    0.62 & \sfrac{1}{2}  & \sfrac{1}{4}
\end{pmatrix}
\begin{pmatrix}
    \vec{a}_1 \\
    \vec{a}_2 \\
    \vec{a}_3
\end{pmatrix}
$$
The primitive cell for the $Pnma$ phase is twice the size of the primitice cell for the $Cmcm$ phase.

```{.matplotlib #fig:snse-bz file="figures/snse/bz.py" caption="In-plane section of the Brillouin zone of SnSe."}
```

Experiments presented in this chapter explore the dynamics of SnSe in the $\vec{a}_2$ -- $\vec{a}_3$ plane. The in-plane section of the Brillouin zone is shown in @fig:snse-bz.

### Lattice instability

## Experimental methods

### Sample preparation

Sample preparation involves two steps. Bulk crystals of SnSe were first synthesized, followed by processed by which ultrathin samples were produced.

Bulk SnSe crystals (\SI{20}{\gram}) were synthesized by mixing appropriate ratios of high purity starting materials (Sn chunk, 99.999\%, American Elements, USA and Se shot, 99.999\%, 5N Plus, Canada) in \SI{13}{\milli\meter} diameter quartz tube. The tube was flame-sealed at a residual pressure of $\SI{1e-4}{\mmHg}$, then slowly heated to \SI{1223}{\kelvin} over \SI{10}{\hour}, soaked at this temperature for 6h and subsequently furnace cooled to room temperature. The obtained ingot was crushed into powder and flame-sealed in a quartz tube, which was placed into another, bigger, flame-sealed quartz tube. A crystal with dimensions of $\sim$\SI{13}{\milli\meter} (diameter) $\times$ \SI{20}{\milli\meter} (length) was obtained.

In order to obtain ultrathin samples suitable for ultrafast electron scattering experiments, two methods were tried. Ultramicrotomy yielded suitable samples, but the author wanted to rule out the effect of strain induced by the sample preparation. To this end, a sample was also prepared via mechanical exfoliation.

#### Ultramicrotomy

Six samples were prepared via ultramicrotomy, a sample preparation technique which involves the cutting of samples using a diamond blade. While this technique has mostly been used to prepare organic samples for cryo-electron microscopy, it has been successfully used in the past decade to prepare samples of two-dimensional materials such as 4H-TaSe$_2$[@Erasmus2012], 1T-TaS$_2$[@Eichberger2013]. The author initially tried ultramicrotomy to prepare samples of 1T-TiSe$_2$[@Otto2020].

```{.matplotlib #fig:snse-ultramicrotomy-femr file="figures/snse/ultramicrotome.py" caption="Ultramicrotome setup used to prepare SnSe samples, from two different angles. Crystals embedded in epoxy are mounted, and a diamond knife is used to shave a section. Sections slide down into a boat filled with water, where they are lated fished out. The images were provided by S. K. Sears from McGill University's Facility for Electron Microscopy Research."}
```

First, a bulk crystal of SnSe is embedded in epoxy, with the cutting plane parallel to the material layers. The surface of the bulk crystal is trimmed with a \ang{45} diamond blade to reveal fresh cutting surface. Sections are then cut with a sharper, \ang{35} diamond blade which then fall in a small water container. Floating sections are then fished out with a carbon-coated TEM grid. The assembly is shown in @fig:snse-ultramicrotomy-femr. The effect of cutting on the embedded crystals is shown in @fig:snse-ultramicrotomy. 

```{.matplotlib #fig:snse-ultramicrotomy file="figures/snse/sample.py" caption="Stages of sample preparation via ultramicrotome. **a)** Prism of SnSe embedded in epoxy.  **b)** 4x magnification on the cut surface shown in a) shows large crystalline plateaus. **c)** Surface of embedded SnSe prism after trimming with a \ang{45} diamond knife. **d)** \SI{60}{\nano\meter} section of SnSe cut with a \ang{35} diamond knife. The images were provided by H. Gnaegi from Diatome, Ltd."}
```

The samples used in this work were prepared by J. Mui from the Facility for Electron Microscopy Research (FEMR) at McGill University. In the end, six samples were produced: three with a thickness of \SI{70}{\nano\meter} and three with a thickness of \SI{90}{\nano\meter}, each with an area of approximately \SI{200 x 200}{\micro\meter}. 

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

```{.matplotlib #fig:snse-diffuse file="figures/snse/diffuse.py" caption="Ultrafast electron scattering pattern of SnSe and associated time traces for a sample photoexcited with photoexcitation density of \SI{2}{\milli\joule\per\square\centi\meter}. **a)** Equilibrium scattering pattern of SnSe oriented along the [100] direction. **b)** Integration geometry used to extract time-traces in panels d), differentiating between diffracted intensity (inner disk) and diffuse intensity near $\Gamma$ (outer torus). **c)** Line cut across the horizontal line shown in panel b). The Bragg peak is fit with a Voigt profile (solid black line) with a full-width at half-max of \SI{0.158}{\per\angstrom}, much smaller than the diameter of the inner disk (\SI{0.228}{\per\angstrom}). **d)** Transient ultrafast electron scattering intensity at various points in the Brillouin zone."}
```

Diffuse intensity dynamics for a SnSe sample photoexcited with \SI{1.55}{\electronvolt} light with an excitation density of \SI{2}{\milli\joule\per\square\centi\meter} is shown in @fig:snse-diffuse. The time-resolved dynamics can be separated into three main components that are spatially-separated in reciprocal space. First, on every Bragg peak ($|\vec{k}| < \SI{0.114}{\per\angstrom}$, region (1)), the transient Debye-Waller effect is observed as a decrease in diffracted intensity. The observed change in scattering intensity can be modelled with a biexponential function with a fast time-constant of \SI{300 \pm 100}{\femto\second} and a slow time-constant of \SI{3.7 \pm 0.6}{\pico\second}. In a torus defined by $\SI{0.114}{\per\angstrom} < |\vec{k}| < \SI{0.228}{\per\angstrom}$ (region (2)), the scattering intensity shows completely different dynamics. There, an initial increase is observed with associated time-constant of \SI{300 \pm 100}{\femto\second}, followed by an exponential recovery with time-constant \SI{3.7 \pm 0.6}{\pico\second}. Finally, everywhere else across the Brillouin zone ($|\vec{k}| \geq \SI{0.228}{\per\angstrom}$, region (3)), only a slow increase in scattering intensity is observed, with an associated time-constant of \SI{3.7 \pm 0.6}{\pico\second}. The dynamics across all three regions are identical for every interrogated reflection, with a scaling of the dynamics amplitude with $|\vec{q}|^2$ which is expected from the scaling of the one-phonon structure factors.

A compraison of the slow dynamics associated with the high-symmetry points $Y$, $Z$, and $T$ are shown in @fig:snse-highsym, compared to the average across region (3) shown in @fig:snse-diffuse. Apart from a small variation in amplitude, the diffuse dynamics at $Y$, $Z$, $T$ are identical to the average.

```{.matplotlib #fig:snse-highsym file="figures/snse/highsym.py" caption="Comparison of the diffuse intensity dynamics at various in-plane high-symmetry points."}
```

Effects of photoinduced strain or lattice expansion can be ruled out from the dynamics in region (2). All Bragg peaks were fit with a Gaussian profile at each time-delay, and the changes to their positions and widths were tracked. The results are presented in @fig:snse-peak-profiles for a few representative reflections. This analysis demonstrates that the dynamics associated with region (2) are not due to a modification of the Bragg peak profile, but suggest that the signal might related to diffuse intensity dynamics. 

```{.matplotlib #fig:snse-peak-profiles file="figures/snse/widths.py" caption="Dynamics of the width and position of various Bragg peaks following photoexcitation. For every time-delay, Bragg peaks were fit with a Gaussian function. In the right column, the change in full-width at half-maximum $\Delta \sigma$ is shown over time. In the left column, the absolute shift in the center position of the peak $\Delta x_c$ is shown, as a percentage of the average full-width at half-maximum $\bar{\sigma}$. For all plots, the error bars represent the covariance of fit parameter."}
```

Recall the definition for the diffuse intensity (@eq:scattering-diffuse-intensity):
$$
    I_1(\vec{q}) \propto \sum_{\lambda} \frac{n_{\lambda}(\vec{k}) + 1/2}{\omega_{\lambda}(\vec{k})} |F_{1\lambda}(\vec{q})|^2
$$
This equation implies that an increase in one or more of the $(n_{\lambda}+1/2) / \omega_{\lambda}$ terms at $\vec{k}=\vec{0}$. Contrary to previous work on graphite in @sec:graphite, deconvolving the contribution from all phonon modes is not possible in SnSe because all phonon branches are thermally-occupied at room temperature[@Skelton2016]. The rapidity of the diffuse increase in region (2) suggests impulsive softening of a particular mode ($\downarrow \omega_{\lambda}$), rather than energy transfer via electron-phonon coupling ($\uparrow n_{\lambda}$).

```{.matplotlib #fig:snse-forbidden file="figures/snse/forbidden.py" caption="Comparison of the relative intensity change $\Delta I / I_0$ at $\Gamma$ near reflections that are either allowed (left column) or forbidden (right column) by the $Pnma$ space group. The integration geometry shown in @fig:snse-diffuse is used."}
```

The observation of forbidden reflections in the thicker samples allows to confirm that the dynamics in region (2) are due to changes in diffuse intensity. Normally, only in-plane reflections of the form $\set{ (0kl) \mid k+l \in 2\mathbb{Z}}$ are allowed by the $Pnma$ space group, where $2\mathbb{Z}$ is understood to be the set of even integers. However, forbidden reflections can be observed in thicker samples where double-diffraction can occur. The comparison of the scattering dynamics in region (2) in a thick sample, around reflections that are allowed or forbidden by $Pnma$, are shown in @fig:snse-forbidden. The cross-section for double-diffraction is much higher than the cross-section for a diffraction event *plus* diffuse scattering (see @sec:scattering-multiple). Therefore, the fact that the dynamics in region (2) are not seen around forbidden reflections indicates that those dynamics are due to changes in the diffuse intensity near the $\Gamma$ point.

### Effect of charge-transfer on diffracted intensity

```{.matplotlib file="figures/snse/chargetransfer.py" caption="Expected intensity change $\Delta I_{(hkl)}$, in proportion to the experimental error $\sigma_e$, for the transfer of a single $p$ electron from Se to Sn. $\Delta I_{(hkl)} / \sigma_e > 1$ implies an intensity change that can be measured, indicating that no in-plane reflection is expected to display signals related to charge-transfer. In-plane reflections up to $|\vec{q}| < \SI{10}{\per\angstrom}$ were considered. The experimental error $\sigma_e$ represents the standard error of the variation in diffraction intensity before photoexcitation."}
```

### Evolution of Bragg peak profile

### Ultrafast phonon softening from Bragg intensity


The time-resolved suppression of Bragg intensity due to atomic vibrations is given by the following expression:
\begin{align}
    \frac{I_0(\vec{q},t) - I_0(\vec{q},0)}{I_0(\vec{q},0)} 
        &\equiv \Delta I_0(\vec{q},t) \nonumber \\
        &= \frac{N_c I_e |F|^2}{N_c I_e |F|^2}\frac{\left[e^{-2M(\vec{q},t)} - e^{-2M(\vec{q},0)}\right]}{e^{-2M(\vec{q},0)}} \nonumber \\
        &=e^{-2[M(\vec{q},t) - M(\vec{q},0)]} - 1
\end{align}
Rearranging terms:
$$
    - \frac{1}{2} \ln \left[1 + \Delta I_0(\vec{q},t) \right] = 2 \pi^2 |\vec{q}|^2 \Delta \langle u^2\rangle
$$
The change in root-mean-square vibration $\langle u^2(t)\rangle$ can be used to infer a change in vibrational frequency $\omega_{TO}(t)$ by considering the amount of energy stored in the mode $E_{TO}$, in the classical limit. Provided that the energy is conserved\cite{Kittel1976}:

\begin{align}
    E_{TO} &= 8 \pi^2 m \omega_{TO}^2(t) \langle u^2(t)\rangle \nonumber \\
           &= 8 \pi^2 m \omega_{TO}^2(0) \langle u^2(0)\rangle \nonumber
\end{align}

We can rearrange terms to use our experimental time-series for $\Delta \langle u^2\rangle$:

\begin{equation}
    \frac{\omega_{TO}^2(0)}{\omega_{TO}^2(t)} =  \frac{ \Delta \langle u^2(t)\rangle}{ \langle u^2(0)\rangle} - 1
\end{equation}

```{.matplotlib file="figures/snse/softening.py" caption="Measurement of the softening of the TO$_c$ mode, extracted directly from the transient Debye-Waller dynamics, assuming that atomic vibrational amplitude only changes based on TO$_c$ renormalization at early times ($\tau < \SI{1}{\pico\second}$). Color bar shows associated absorbed energy per cell $E_c$. \textbf{Inset} Increase in isotropic mean-square-displacement of all atoms, due to the change in vibrational frequency of the TO$_c$ mode exclusively. Boxes are used to represent error bars along both axes."}
```

## Spectroscopic measurements

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
