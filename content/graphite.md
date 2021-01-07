
# Momentum-resolved excitation couplings in graphite {#sec:graphite}

Ultrafast electron diffuse scattering is an experimental technique pioneered over the past few years [@Chase2016; @Harb2016; @Waldecker2017]. The author also participated in the early development of the technique with work by M. Stern, L. P. René de Cotret *et al* [@Stern2018]. This early work on graphite was an experimental tour-de-force, but the observations remained qualitative rather than quantitative.

TODO: important differentiator for UEDS vs X-ray is the amount of redundancy. The redundancy leads to the ability to invert the energy-integrated nature of UEDS.

This chapter will detail how to extract the quantitative information encoded in UEDS experiments. The redundancy present in UEDS data will be used to overcome the energy-integration of electron scattering detectors, so that phonon dynamics can be observed across the Brillouin zone.

## Single-crystal graphite

Single-crystal graphite is a two-dimensional material made of a repeating layers of carbon atoms arranged in a hexagonal lattice. Its lattice vectors $\left\{ \vec{a}_i \right\}$ are:
$$
\begin{pmatrix}
    \vec{a}_1 \\
    \vec{a}_2 \\
    \vec{a}_3
\end{pmatrix}
=
\begin{pmatrix}
    \sqrt{3}a  & -a & 0 \\
    0          & 2a & 0 \\
    0          &  0 & c 
\end{pmatrix}
\begin{pmatrix}
    \vec{e}_1 \\
    \vec{e}_2 \\
    \vec{e}_3
\end{pmatrix}
$${#eq:graphite-lattice}
where $a=\SI{1.232}{\angstrom}$, $c=\SI{6.711}{\angstrom}$, and $\left\{ \vec{e}_i \right\}$ are understood to be the usual Euclidean vectors. Carbon atoms $\left\{ \vec{c}_i \right\}$ are positioned as follows:
$$
\begin{pmatrix}
    \vec{c}_1 \\
    \vec{c}_2 \\
    \vec{c}_3 \\
    \vec{c}_4 
\end{pmatrix}
=
\begin{pmatrix}
    0            & 0            & 0            \\
    \sfrac{1}{3} & \sfrac{2}{3} & 0            \\
    0            & 0            & \sfrac{1}{2} \\
    \sfrac{2}{3} & \sfrac{1}{3} & \sfrac{1}{2} \\
\end{pmatrix}
\begin{pmatrix}
    \vec{a}_1 \\
    \vec{a}_2 \\
    \vec{a}_3
\end{pmatrix}
$$

Note carbon atoms $\vec{c}_1$ and $\vec{c}_2$ form a sheet of graphene, while $\vec{c}_3$ and $\vec{c}_4$ form a separate separate sheet of graphene rotated by $\tfrac{\pi}{3}$ (\ang{60}). The stacking of both sublattices along the $\vec{a}_3$ axis confers graphite increased discrete rotational symmetry (6-fold) about the stacking axis, compared to the 3-fold rotational symmetry of graphene. More specifically, the point group for graphite is $6/mmm$, while the Hermann-Mauguin symbol for the space group is $P6_3/mmc$. 

The real-space 6-fold rotational symmetry is mirrored in reciprocal space. The geometry of the Brillouin zone is shown in @fig:graphite-bz, including the high-symmetry points  $\vec{\Gamma}$ (zone-center), $\vec{M}$, and $\vec{K}$.

```{#fig:graphite-bz .matplotlib file="figures/graphite/bz.py" caption="In-plane section of the Brillouin zone of graphite, spanned by reciprocal lattice vectors $\vec{b}_1$ and $\vec{b}_2$. The three classes of high-symmetry points are shown. $\vec{M}$ and $\vec{K}$ have symmetry equivalents on every edge and vertex respectively."}
```

TODO: electronic landscape (Dirac cones)

### Phonon landscape

With four atoms per unit cell, the structure of graphite supports twelve distinct phonon modes. The layered nature of graphite results in a clear distinction between in-plane and out-of-plane phonon modes. These modes are longitudinal modes LA, LO1 - LO3, in-plane transverse modes TA, TO1 - TO3, and out-of-plane traverse modes ZA, ZO1 - ZO3. The phonon dispersion relation for in-plane modes was calculated as described further below in @sec:computational-details, and is shown in @fig:graphite-static-dispersion. The out-of-plane transverse modes are not taken into account as the experiments are not sensitive to them. 

@fig:graphite-static-dispersion reveals why graphite is the *perfect* benchmark material for UEDS. All of the lattice modes are so stiff (high-energy) that, at room temperature, thermally-occupied modes are concentrated at $\vec{\Gamma}$. Therefore, the contrast between static diffuse intensity and diffuse intensity after photoexcitation is bound to be maximal.

TODO: talk about Kohn anomalies, not super visible

```{#fig:graphite-static-dispersion .matplotlib file="figures/graphite/static-dispersion.py" caption="Phonon dispersion relation of graphite for in-plane modes LA, TA, and two-fold degenerate modes LO and TO. The path in reciprocal space is shown in the center. A horizontal dashed line at \SI{25}{\milli\electronvolt} indicates the average energy stored in the phonon modes at room temperature (\SI{300}{\kelvin})."}
```

## The first hundred femtoseconds viewed by trARPES

Electron scattering measurements are not able to isolate to the dynamics of the electronic system. Although UED and UEDS are sensitive to the dynamics of electrons whizzing about in a material, understanding the electronic contribution to measured signals requires some prior knowledge. To understand the ultrafast electron diffuse scattering results presented further, it is important to understand the effect of photoexcitation during first \SI{100}{\femto\second} as observed by time- and angle-resolved photoemission spectroscopy (trARPES).

## Experimental and computational methods

### Sample preparation

Single-crystal flakes of natural graphite \SIrange{10}{90}{\nano\meter} thick were prepared using a mechanical exfoliation procedure analogous to the work by Novoselov *et al* [@Novoselov2004], briefly described here. Thick flakes were embedded in Crystalbond glue on a \SI{3}{\milli\meter} copper TEM grid (200 lines per inch). The embedded flakes are then exfoliated using ordinary adhesive tape. The procedure was repeated until the embedded flakes were translucent when observed under an optical microscope. The glue is then delicately washed away with a solvent. The choice of the solvent is dependent on the glue used; in the present case, acetone was used. Sample thickness has been measured directly using atomic force microscopy characterization (see @fig:graphite-afm). Once a suitable sample has been identified, an aperture was made using aluminum foil to isolate a sample region with uniform thickness. The resulting sample used in this work covered 500 × 500 \si{\square\micro\meter}, with a thickness of \SI{70}{\nano\meter}.

![Characterization of the thickness of potential graphite samples. **Left** optical microscope image of exfoliated graphite. Translucent regions A and B are highlighted as potential samples. **Right** Atomic force microscope measurement shows a sample thickness of \SI{88}{\nano\meter}. Modified from Chatelain [@Chatelain2014].](diagrams/graphite_exfoliation.pdf){#fig:graphite-afm}

### Data acquisition 

The UEDS experiments presented in this chapter made use of the experimental setup presented in @sec:experimental_setup. Ultrashort laser pulses of light  were shone at $t=t_0$ on a thin single-crystal specimen of graphite, oriented in the $\left[001\right]$ direction. Compressed electron bunches with 10^7^ electrons per bunch were transmitted through the sample at $t=t_0 + \tau$. The time-delay $\tau$ was scanned from \SI{-40}{\pico\second} to \SI{680}{\pico\second}.

The interrogated film were pumped with a pump spot of 1 × 1 \si{\square\milli\meter} full-width at half-maximum (FWHM), ensuring nearly uniform illumination of the probed volume. The film was pumped at a fluence of \SI{12}{\milli\joule\per\square\centi\meter}, resulting in an absorbed energy density of \SI{8}{\joule\per\cubic\meter}. The scattering patterns are collected with a Gatan Ultrascan 895 camera: a 2.54 × 2.54 \si{\square\cm} phosphor screen fiber coupled to a 2048 px × 2048 px charge-coupled detector (CCD) placed \SI{25}{\centi\meter} away from the sample. The experiment herein consists of time delays in the range of $\SIrange{-40}{680}{\pico\second}$. Per-pixel scattering intensity fluctuations over laboratory time reveals a transient dynamic range of 1 : 10^8^, allowing the acquisition of diffraction patterns and diffuse scattering patterns simultaneously. A static diffraction pattern is shown in @fig:graphite-static a).

Due to the flatness of the Ewald sphere for \SI{90}{\kilo\electronvolt} electrons, many symmetry-related reflections are visible within each pattern. The information contained in a set of symmetry-equivalent reflections is redundant due to the point-group symmetry of the scattering crystal. As long as the point-group symmetry is not broken by photoexcitation its consequences, it is possible to harness the redundancy to enhance the signal-to-noise ratio of a UEDS dataset. In the case of graphite, no observable symmetry-breaking phenomena is brought on by photoexcitation at \SI{1.55}{\electronvolt} when looking at the raw data. Moreover, trARPES experiments [@Stange2015] do not show the opening of a gap in the electronic band structure, albeit at much lower photoexcitation densities, which would be indicative of point-group symmetry breaking. The point-group of graphite is $6/mmm$, which encompasses 6-fold discrete rotational symmetry in the $\vec{a} \times \vec{b}$ plane. Therefore, specifically in the case of graphite oriented in the $\left[ 001 \right]$ direction, we can safely enhance the diffuse signals by a factor of $\sqrt{6}$ by the use of a six-fold discrete azimuthal average:
$$
    I(\vec{q}, \tau) \to \frac{1}{6} \sum_{n=1}^6 I( \vec{R}(\tfrac{\pi n}{3}) \cdot \vec{q}, \tau)
$$
where $\vec{R}(\theta)$ is the in-plane rotation matrix:
$$
    \vec{R}(\theta) = \begin{pmatrix}
                      \cos{\theta} & -\sin{\theta} & 0\\
                      \sin{\theta} & \cos{\theta}  & 0\\
                      0            & 0             & 1
                     \end{pmatrix}
$$
Throughout the rest of this chapter, "scattering intensity" will imply discrete azimuthal average unless otherwise noted. An example of six-fold averaged diffraction pattern is shown in @fig:graphite-static b).

```{#fig:graphite-static .matplotlib file="figures/graphite/diff-static.py" caption="Static diffraction pattern of graphite. Brillouin zones are shown around each reflection to guide the eye. **a)** static, unprocessed diffraction pattern. **b)** Six-fold discrete azimuthal average of the diffraction pattern in a) results in $\sqrt{6}$ increase in signal-to-noise ratio."}
```

### Computational details {#sec:computational-details}

This section contains the details of the calculations used throughout this chapter, including what is shown in @fig:graphite-static-dispersion. The aim of the computations was to extract the phonon mode frequencies $\left\{ \omega_j(\vec{k})\right\}$ and polarization vectors $\left\{ \vec{e}_{j, s}(\vec{k})\right\}$ that appear in TODO: ONEPH EQUATION.

In order to calculate the force between atoms, the structure of graphite was computed via *structure relaxation*. Structure relaxation was performed using the plane-wave self-consistent field program `PWSCF` from the `QUANTUM ESPRESSO` software suite [@Giannozzi2017]. The structure was relaxed using a 18 × 18 × 10 $\vec{k}$ mesh centered at $\vec{\Gamma}$ -- selected using the Monkhorst-Pack method [@Monkhorst1976] -- and force and energy thresholds of \SI{1e-8}{\rydberg\per\bohr} and \SI{1e-15}{\rydberg} respectively, where \si{\bohr} is the Bohr radius.

The dynamical matrices were computed on a 5 × 5 × 3 $\vec{q}$ mesh using a self-consistency threshold of \SI{1e-18}{\rydberg}. The resulting graphite structure is equivalent to @eq:graphite-lattice, with $a=\SI{1.231}{\angstrom}$ and $c=\SI{6.837}{\angstrom}$.

The phonon mode frequencies $\left\{ \omega_j(\vec{k})\right\}$ and polarization vectors $\left\{ \vec{e}_{j, s}(\vec{k}}\right\)$ were computed using the `PHONON` program, again within the `QUANTUM ESPRESSO` software suite. This calculation made use of the B86b exchange-coupled Perdew-Burke-Ernzerhof generalized-gradient approximation [@Becke1986; @Perdew1996] and the projector augmented-wave method [@Blochl1994]. The cutoff energy of the wave function was set to \SI{100}{\rydberg}, while the cutoff energy for the charge density was set to \SI{1.2e3}{\rydberg}. A Fermi-Dirac smearing of \SI{0.06}{\rydberg} was also applied. To include the dispersion of energy along the stacking axis $\vec{a}_3$, the exchange-hole dipole moment method was used [@Becke2007]. 

#### Clustering of phonon properties into physically-relevant branches

The calculation of phonon properties relevant to his work involves diagonalizing the dynamical matrix. In practice, the force constants that make the dynamical matrix need to be evaluated for every irreducible $\vec{k}$ point, which means that a diagonalization procedure is repeated for every irreducible $\vec{k}$ point. As with most (all?) diagonalization procedures, the eigenvalues and eigenvectors (phonon frequencies and polarization vectors, in the present case) are returned in an order that is not physically-relevant. Therefore, in order to determine what are the phonon properties of a phonon mode $j$, the phonon properties calculated by `PHONON` need to be *clustered*.

As of writing this, there is no component within the `QUANTUM ESPRESSO` software suite that can do this. The scheme described in this section was used instead. The general idea behind the procedure is that phonon properties are continuous. Let $P_{j, \vec{k}_i}$ be the abstract vector representing phonon properties of mode $j$ at one of the irreducible $\vec{k}$ points $\left\{\vec{k}_i\right\}$:
$$
    P_{j, \vec{k}_i} = \begin{bmatrix} \omega_{j, \vec{k}_i} & \vec{e}_{s=1,j,\vec{k}_i} & ... & \vec{e}_{s=M,j,\vec{k}_i} \end{bmatrix}^T
$$
where the index $s$ runs over all $M$ atoms of the unit cell. Define the metric between two such abstract vectors $\vec{P}_{i, \vec{k}}$ and $\vec{P}_{j, \vec{k}'}$ as:
$$
    \left\Vert \vec{P}_{i, \vec{k}} - \vec{P}_{j, \vec{k}'} \right\Vert = |\omega_{i, \vec{k}} - \omega_{j, \vec{k}'}|^2 + \sum_s \left\Vert \vec{e}_{s,i,\vec{k}} - \vec{e}_{s,j,\vec{k}'}\right\Vert
$$

A one-dimensional path $\gamma$ connecting all irreducible points $\left\{\vec{k}_i\right\}$ was defined, starting at $\vec{\Gamma}$. At $\vec{\Gamma}$, polarization vectors are associated with a mode based on geometry. For example, a mode with negligible frequency and polarization vectors at that all point in the same direction physicall corresponds to a longitudinal acoustic mode. The manual assignment for high-frequency optical modes is a bit more arbitrary. Then, following the path $\gamma$, the assignment of phonon branches $i$ at $\gamma(\vec{k} + \vec{\Delta})$ minimizes the quantity $\left\Vert \vec{P}_{i, \vec{k}} - \vec{P}_{j, \vec{k} + \vec{\Delta}} \right\Vert$. 

The procedure described above has been adapted for numerical evaluation and is now part of the `scikit-ued` software package [@RenedeCotret2018]. For the rest of this chapter, it will be assumed that the phonon properties are clustered such that the usual physical phonon branch labels (e.g. TA, LO) are meaningful.

## Diffuse intensity dynamics

The change in scattering intensity $\Delta I(\vec{q}, t=\tau) = I(\vec{q}, \tau) - I(\vec{q}, -\infty)$ for representative time-delays is shown in @fig:graphite-ueds. 

```{#fig:graphite-ueds .matplotlib file="figures/graphite/ueds.py" caption="Change in scattering intensity $\Delta I(\vec{q}, t=\tau) = I(\vec{q}, \tau) - I(\vec{q}, -\infty)$ of photoexcited graphite for a few representative time-delays $\tau$. Hexagonal Brillouin zones are shown on half of the reflections to guide the eye. Scattering patterns show diffuse dynamics in the range of $|\vec{q}| < \SI{12}{\per\angstrom}$. Negative going features (blue) are exclusively due to the transient Debye-Waller effect on the Bragg peaks. All positive changes (red) are dynamics of the diffuse scattering intensity."}
```

First and foremost, note that negative-going features are only visible in the vicinity of Bragg peaks. This is due exclusively to the transient Debye-Waller effect [@Ligges2009]. As photodeposited energy transfers from the electrons to the lattice, average real-space disorder due to phonons lowers the overall symmetry of the lattice, which in turn results in smaller Bragg peaks. This is analogous to the (static) thermal Debye-Waller effect [@Debye1913; @Waller1923]. This leaves only positive-going features everywhere away from Bragg peaks, which must be diffuse in nature.

Second, the structure visible across the Brillouin zone is different for reflections that ought to represent the same physical lattice waves. For example, consider the $(010)$ and $(020)$ reflections, \SI{500}{\femto\second} after photoexcitation, shown in @fig:graphite-ueds-zoomed. Around $(010)$, diffuse intensity increases near two of the $\vec{K}$ points, while diffuse intensity is increasing along steaks between $\vec{\Gamma}$ and $\vec{K}$ near $(020)$. This drastic difference between two physically-equivalent reflections can be explained by the difference in the one-phonon structure factors $|F_{1j}(\vec{q}, \tau)|^2$.

```{#fig:graphite-ueds-zoomed .matplotlib file="figures/graphite/ueds-zoomed.py" caption="Comparison of the diffuse intensity change after \SI{500}{\femto\second} for two Brillouin zones, $(010)$ and $(020)$. The colormap scaling is identical to @fig:graphite-ueds. The difference between those two images can be explained by the difference in one-phonon structure factors $|F_{1j}(\vec{q}, \tau)^2|$."}
```

## The one-phonon structure factors

Recall the definition of diffuse intensity from TODO: ONEPH EQUATION:
$$
    I_1(\vec{q}, \tau) = N_c I_e \sum_j \frac{n_j(\vec{k}, \tau) + \sfrac{1}{2}}{\omega_j(\vec{k}, \tau)} \left| F_{1j}(\vec{q}, \tau)\right|^2
$$
and
$$
    F_{1j}(\vec{q}, \tau) = \sum_s e^{-W_s(\vec{q}, \tau)} \frac{f_s(\vec{q})}{\sqrt{\mu_s}} \left( \vec{q} \cdot \vec{e}_{j,s}(\vec{k})\right) 
$$
where $N_c$ is the number of scattering unit cells, $I_e$ is the intensity of a single scattering event, $s$ are indices that label atoms in the unit cell, $W_s(\vec{q}, \tau)$ is the Debye-Waller factor associated with atom $s$, $f_s(\vec{q})$ are the atomic form factors, $j \in \left\{ 1, 2, ..., 3N\right\}$ runs over the phonon modes, $\left\{ n_j(\vec{k}, \tau)\right\}$ are the phonon populations, and $\left\{ \omega_j(\vec{k}, \tau)\right\}$ are the phonon frequencies.

As shown in @fig:graphite-ueds-zoomed, observed diffuse intensity is the combination of quantities that can be separated in two categories. Some quantities, like phonon frequencies are *local* in reciprocal space (i.e. defined on $\vec{k}$). On the other hand, the one-phonon structure factor is nonlocal in reciprocal space (defined on $\vec{q}$). This results in measurable differences when comparing diffuse intensity at various reflections, which will be used to overcome the energy-integrative nature of UEDS experiments.

In order to recover mode-dependent phonon populations, one-phonon structure factors need to be calculated for every relevant phonon mode. Most importantly, the phonon polarization vectors $\left\{ \vec{e}_{j, s}(\vec{k})\right\}$ need to be determined. These vectors describe the direction of movement for each atom due to a particular lattice wave. Fortunately, polarization vectors are a byproduct of the calculation of the phonon dispersion relation, shown in @fig:graphite-static-dispersion.

### Transient Debye-Waller factors

One more quantity needs to be discussed before one-phonon structure factors can be computed: the Debye-Waller factors $W_s(\vec{q}, \tau)$. The Debye-Waller factors describe the reduction of scattering intensity at $\vec{q}$ due to the effective deformation of the scattering potential of atom $s$ by *all* phonon branches. The general expression for the *anisotropic* Debye-Waller factor for atom $s$ is given by[@Xu2005]:
$$
    W_s(\vec{q}, \tau) = \frac{1}{4 \mu_s} \sum_{s, j} |a_j(\vec{k}, \tau)|^2 |\vec{q} \cdot \vec{e}_{j,s}(\vec{k})|^2
$${#eq:graphite-debye-waller}
where $a_j$ is vibrational amplitude for the phonon mode $j$:
$$
    |a_j(\vec{k}, \tau)|^2 = \frac{2 \hbar}{N_c \omega_j(\vec{k}, \tau)} \left( n_j(\vec{k}, \tau) + \frac{1}{2}\right)
$$
Based on @eq:graphite-debye-waller, the Debye-Waller factor is not sensitive measure to the wavevector-dependent nonequilibrium dynamics of lattice waves because it involves the contribution from *all* lattice waves, integrated over $\vec{k}$. Therefore, phonon population dynamics can only affect the magnitude of Debye-Waller factors. 

The potential time-dependence of the Debye-Waller factors was investigated, via the time-dependence of mode populations. Profoundly non-equilibrium distribution of phonon mode populations were simulated, with all modes populated according to room-temperature average energy except for one mode at high temperature. To account for differences in heat capacities, the maximum temperature for optical modes was \SI{5000}{\kelvin}, while the maximum temperature for acoustic modes was capped at \SI{1500}{\kelvin}. These extreme non-equilibrium scenarios increase the value of $\sum_s W_s(\vec{q})$ by at most 1.5\% for optical modes, and 8\% for acoustic modes. The most significant example is shown in @fig:graphite-debye-waller-la, where the LA mode is populated according to a temperature of \SI{1500}{\kelvin}. Since those fractional changes are constant across $\vec{q}$, wavevector-dependent changes in UEDS signals are not significantly impacted by changes in the Debye-Waller factors $W_s(\vec{q}, \tau)$, and so any time-dependence of the one-phonon structure factors can be ignored to a first approximation.

```{#fig:graphite-debye-waller-la .matplotlib file="figures/graphite/debye-waller.py" caption="Effect of profoundly non-equilibrium distribution of phonon mode populations on the Debye-Waller factors. **a)** Sum of Debye-Waller factors according to the expected phonon distribution at room temperature. **b)** Difference between a) and the sum of Debye-Waller factors where all phonon modes are distributed according to room-temperature, except where the LA mode is populated according to a temperature of \SI{1500}{\kelvin}."}
```

### One-phonon structure factor calculations

The calculation for the one-phonon structure factors $|F_{1j}(\vec{q}, \tau=-\infty)|^2$ was carried out from (TODO: ONEPH EQUATION). As discussed in the previous section, the dynamics in one-phonon structure factors due to photoexcitation via the change in Debye-Waller factors can be neglected in this work. Therefore, the calculations of this section assume room-temperature thermal distribution of lattice waves; this is implied by time-delay $\tau=-\infty$ (pre-photoexcitation). The one-phonon structure factors for longitudinal modes LA, LO1 - LO3 are shown @fig:graphite-oneph-longitudinal. The one-phonon structure factors for transverse modes TA, TO1 - TO3 are shown @fig:graphite-oneph-transverse. Out-of-plane modes yield trivial one-phonon structure factors because $\vec{q} \cdot \vec{e}_{j,s}(\vec{k}) = 0$ for all scattering vectors $\vec{q}$ probed by the experiments presented herein.

```{.matplotlib #fig:graphite-oneph-longitudinal file="figures/graphite/oneph-long.py" caption="Calculated one-phonon structure factors $|F_{1j}(\vec{q}, \tau=-\infty)|^2$ of in-plane longitudinal modes at \SI{300}{\kelvin}, for scattering vectors $\vec{q}$ equivalents to the detector area shown in @fig:graphite-ueds. See @fig:graphite-oneph-transverse for the transverse modes equivalent."}
```

```{.matplotlib #fig:graphite-oneph-transverse file="figures/graphite/oneph-trans.py" caption="Calculated one-phonon structure factors $|F_{1j}(\vec{q}, \tau=-\infty)|^2$ of in-plane transverse modes at \SI{300}{\kelvin}, for scattering vectors $\vec{q}$ equivalents to the detector area shown in @fig:graphite-ueds. See @fig:graphite-oneph-longitudinal for the longitudinal modes equivalent."}
```

The one-phonon structure factors are display complex structure in reciprocal space. The structure for a particular mode can be thought of as selection rule: regions in $\vec{q}$ space where the one-phonon structure factor for mode $j$ is large are regions where phonon mode $j$ contributes importantly to diffuse intensity. The structure of $|F_{1j}|^2$ is most determined by the values of terms like $\left\{ \vec{q} \cdot \vec{e}_{j,s}(\vec{k}) \right\}$. For acoustic modes LA and TA, the structure of the one-phonon structure factor is easier to understand near $\vec{\Gamma}$. For example, $|F_{1j}|^2$ is highest for the LA mode in the radial $\vec{q}$ direction because the polarization vectors of those modes are parallel to $\vec{q}$; this fact defines longitudinal waves. For transverse waves, where atomic motion is perpendicular to the lattice wave propagation direction, $|F_{1j}|^2$ is higher in the azimuthal direction.

An alternative visualization for one-phonon structure factors are weighted dispersion curves, as shown in @fig:graphite-weighted-dispersion.

```{.matplotlib #fig:graphite-weighted-dispersion file="figures/graphite/weighted-dispersion.py" caption="Calculated one-phonon structure factors visualized as weighted dispersion curves for selected in-plane modes. The color saturation of dispersion curves is proportional to $|F_{1j}|^2$ of the associated mode. Equivalent paths in the Brillouin zone around two reflections are shown to highlight the high degree of reciprocal space structure: $(010)$ on the left and $(\bar{1}10)$ on the right. The geometry of the paths with respect to $(000)$ are shown in the inset on the lower left."}
```

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
