
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

### Electronic structure

The electronic structure of graphite has been the subject of much study. Graphite is a very close relative of the Dirac semi-metal graphene; its electronic dispersion is very similar from the point of view of experiments presented in this chapter. The electronic dispersion of graphite the $\pi$ and $\pi^{\star}$ bands can be calculated in the tight-binding framework [@Wallace1947]. Taking into account nearest-neighbor interactions only, the dispersion for both bonding and antibonding bands are given by:
$$
    E(\vec{k})_{\pi^{\star}} = t \sqrt{3 + f(\vec{k})}
$$
$$
    E(\vec{k})_{\pi} = - t \sqrt{3 + f(\vec{k})}
$$
with
\begin{align}
    f(\vec{k}) &= 2  \cos\left[\sqrt{3} a ~ (\vec{k} \cdot \vec{b}_2) \right] \\
               &+ 4  \cos\left[\tfrac{\sqrt{3} a}{2} ~ (\vec{k} \cdot \vec{b}_2) \right] 
                    \cos\left[\tfrac{3 a}{2} ~ (\vec{k} \cdot \vec{b}_1) \right] \nonumber
\end{align}
where $t$ is the tight-binding constant and $a$ is the carbon-carbon distance. The in-plane dispersion for a tight-binding energy of $t=\SI{2.7}{\electronvolt}$ is shown on @fig:graphite-electronic-structure. More complete expressions that are not limited to nearest neighbor interactions are given elsewhere [@Neto2009].

```{.matplotlib #fig:graphite-electronic-structure file="figures/graphite/estructure.py" caption="In-plane electronic dispersion $E(\vec{k})$ for graphite. The in-plane section of the Brillouin zone is shown below. The dotted black line represents the line-cut used for @fig:graphite-photoexcitation."}
```

The nonequilibrium behavior of photoexcited graphite is dominated by the structure of the electronic dispersion near the $\vec{K}$ points, at the corner of the Brillouin zone. These points, where the dispersion of each band meet, are called Dirac points. Near the Dirac points, the dispersion adopts a conical shape, called Dirac cones. It is the Dirac cones that tell us the consequences of photoexcitation with \SI{800}{\nano\meter} light. Photons at these energies (\SI{1.55}{\electronvolt}) can only drive vertical (zero-momentum) transitions near the Dirac cones. As the electron cloud thermalizes, two classes of momentum-conserving decay pathways involving phonons emerge. One such pathway allows for an electron to move across the Dirac cone, emitting an $E_{2g}$ phonon with small wavevector $\vec{k} \sim \vec{0}$. Another pathway allows for an electron to hop onto a neighboring Dirac cone, emitting an $A_1^\prime$ phonon with large wavevector $\vec{k} \sim \vec{K}$. These pathways are shown diagramatically on @fig:graphite-photoexcitation.

```{.matplotlib #fig:graphite-photoexcitation file="figures/graphite/photoexcitation.py" caption="Cut of the electronic dispersion $E(\vec{k})$ along the $\vec{K}$--$\vec{K}$ line (see @fig:graphite-electronic-structure) shows the effects of photoexcitation with \SI{800}{\nano\meter} photons ($\gamma$). The momentum-conserving decay path are highlighted and explained in the text."}
```

### Phonon landscape {#sec:graphite-phonon-landscape}

With four atoms per unit cell, the structure of graphite supports twelve distinct phonon modes. The layered nature of graphite results in a clear distinction between in-plane and out-of-plane phonon modes. These modes are longitudinal modes LA, LO1 - LO3, in-plane transverse modes TA, TO1 - TO3, and out-of-plane traverse modes ZA, ZO1 - ZO3. The phonon dispersion relation for in-plane modes was calculated as described further below in @sec:computational-details, and is shown in @fig:graphite-static-dispersion. The out-of-plane transverse modes are not taken into account as the experiments are not sensitive to them. 

@fig:graphite-static-dispersion reveals why graphite is the *perfect* benchmark material for UEDS. All of the lattice modes are so stiff (high-energy) that, at room temperature, thermally-occupied modes are concentrated at $\vec{\Gamma}$. Therefore, the contrast between static diffuse intensity and diffuse intensity after photoexcitation is bound to be maximal.

TODO: talk about Kohn anomalies, not super visible

```{#fig:graphite-static-dispersion .matplotlib file="figures/graphite/static-dispersion.py" caption="Phonon dispersion relation of graphite for in-plane modes LA, TA, and two-fold degenerate modes LO and TO. The path in reciprocal space is shown in the center. A horizontal dashed line at \SI{25}{\milli\electronvolt} indicates the average energy stored in the phonon modes at room temperature (\SI{300}{\kelvin})."}
```

TODO: The temperature-dependence of the phonon spectrum needs to be addressed, as it will play a factor in the discussion of time-resolved measurements.

## The first hundred femtoseconds viewed by trARPES {#sec:graphite-100-fs}

Electron scattering measurements are not able to isolate to the dynamics of the electronic system. Although UED and UEDS are sensitive to the dynamics of electrons whizzing about in a material, understanding the electronic contribution to measured signals requires some prior knowledge. To understand the ultrafast electron diffuse scattering results presented further, it is important to understand the effect of photoexcitation during first \SI{100}{\femto\second} as observed by time- and angle-resolved photoemission spectroscopy (trARPES).



## Experimental and computational methods

### Sample preparation

Single-crystal flakes of natural graphite \SIrange{10}{90}{\nano\meter} thick were prepared using a mechanical exfoliation procedure analogous to the work by Novoselov *et al* [@Novoselov2004], briefly described here. Thick flakes were embedded in Crystalbond glue on a \SI{3}{\milli\meter} copper TEM grid (200 lines per inch). The embedded flakes are then exfoliated using ordinary adhesive tape. The procedure was repeated until the embedded flakes were translucent when observed under an optical microscope. The glue is then delicately washed away with a solvent. The choice of the solvent is dependent on the glue used; in the present case, acetone was used. Sample thickness has been measured directly using atomic force microscopy characterization (see @fig:graphite-afm). Once a suitable sample has been identified, an aperture was made using aluminum foil to isolate a sample region with uniform thickness. The resulting sample used in this work covered 500 × 500 \si{\square\micro\meter}, with a thickness of \SI{70}{\nano\meter}.

![Characterization of the thickness of potential graphite samples. **Left** optical microscope image of exfoliated graphite. Translucent regions A and B are highlighted as potential samples. **Right** Atomic force microscope measurement shows a sample thickness of \SI{88}{\nano\meter}[^1].](diagrams/graphite_exfoliation.pdf){#fig:graphite-afm}

### Data acquisition {#sec:graphite-data-acquisition}

The UEDS experiments presented in this chapter made use of the experimental setup presented in @sec:experimental_setup. Ultrashort laser pulses of light were shone at $t=t_0$ on a thin single-crystal specimen of graphite, oriented in the $\left[001\right]$ direction. Compressed electron bunches with 10^7^ electrons per bunch were transmitted through the sample at $t=t_0 + \tau$. The time-delay $\tau$ was scanned from \SI{-40}{\pico\second} to \SI{680}{\pico\second}.

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

The phonon mode frequencies $\left\{ \omega_j(\vec{k}) \right\}$ and polarization vectors $\left\{ \vec{e}_{j, s}(\vec{k})\right\}$ were computed using the `PHONON` program, again within the `QUANTUM ESPRESSO` software suite. This calculation made use of the B86b exchange-coupled Perdew-Burke-Ernzerhof generalized-gradient approximation [@Becke1986; @Perdew1996] and the projector augmented-wave method [@Blochl1994]. The cutoff energy of the wave function was set to \SI{100}{\rydberg}, while the cutoff energy for the charge density was set to \SI{1.2e3}{\rydberg}. A Fermi-Dirac smearing of \SI{0.06}{\rydberg} was also applied. To include the dispersion of energy along the stacking axis $\vec{a}_3$, the exchange-hole dipole moment method was used [@Becke2007]. 

#### Clustering of phonon properties into physically-relevant branches

The calculation of phonon properties relevant to his work involves diagonalizing the dynamical matrix. In practice, the force constants that make the dynamical matrix need to be evaluated for every irreducible $\vec{k}$ point, which means that a diagonalization procedure is repeated for every irreducible $\vec{k}$ point. As with most (all?) diagonalization procedures, the eigenvalues and eigenvectors (phonon frequencies and polarization vectors, in the present case) are returned in an order that is not physically-relevant. Therefore, in order to determine what are the phonon properties of a phonon mode $j$, the phonon properties calculated by `PHONON` need to be *clustered*.

As of writing this, there is no component within the `QUANTUM ESPRESSO` software suite that can do this. The scheme described in this section was used instead. The general idea behind the procedure is that phonon properties are continuous. Let $P_{j, \vec{k}_i}$ be the abstract vector representing phonon properties of mode $j$ at one of the irreducible $\vec{k}$ points $\left\{\vec{k}_i\right\}$:
$$
    P_{j, \vec{k}_i} \equiv \begin{bmatrix} \omega_{j, \vec{k}_i} & \vec{e}_{s=1,j,\vec{k}_i} & ... & \vec{e}_{s=M,j,\vec{k}_i} \end{bmatrix}^T
$$
where the index $s$ runs over all $M$ atoms of the unit cell. Define the metric between two such abstract vectors $\vec{P}_{i, \vec{k}}$ and $\vec{P}_{j, \vec{k}'}$ as:
$$
    \left\Vert \vec{P}_{i, \vec{k}} - \vec{P}_{j, \vec{k}'} \right\Vert = |\omega_{i, \vec{k}} - \omega_{j, \vec{k}'}|^2 + \sum_s \left\Vert \vec{e}_{s,i,\vec{k}} - \vec{e}_{s,j,\vec{k}'}\right\Vert
$$
A one-dimensional path $\gamma$ connecting all irreducible points $\left\{\vec{k}_i\right\}$ was defined, starting at $\vec{\Gamma}$. At $\vec{\Gamma}$, polarization vectors are associated with a mode based on geometry. For example, a mode with negligible frequency and polarization vectors at that all point in the same direction physicall corresponds to a longitudinal acoustic mode. The manual assignment for high-frequency optical modes is a bit more arbitrary. Then, following the path $\gamma$, the assignment of phonon branches $i$ at $\gamma(\vec{k} + \vec{\Delta})$ minimizes the quantity $\left\Vert \vec{P}_{i, \vec{k}} - \vec{P}_{j, \vec{k} + \vec{\Delta}} \right\Vert$. 

The procedure described above has been adapted for numerical evaluation and is now part of the `scikit-ued` software package [@RenedeCotret2018]. For the rest of this chapter, it will be assumed that the phonon properties are clustered such that the usual physical phonon branch labels (e.g. TA, LO) are meaningful.

## Diffuse intensity dynamics

The change in scattering intensity $\Delta I(\vec{q}, t=\tau) = I(\vec{q}, \tau) - I(\vec{q}, -\infty)$ for representative time-delays is shown in @fig:graphite-ueds. 

```{#fig:graphite-ueds .matplotlib file="figures/graphite/ueds.py" caption="Change in scattering intensity $\Delta I(\vec{q}, t=\tau) = I(\vec{q}, \tau) - I(\vec{q}, \tau < 0)$ of photoexcited graphite for a few representative time-delays $\tau$. Hexagonal Brillouin zones are shown on half of the reflections to guide the eye. Scattering patterns show diffuse dynamics in the range of $|\vec{q}| < \SI{12}{\per\angstrom}$. Negative going features (blue) are exclusively due to the transient Debye-Waller effect on the Bragg peaks. All positive changes (red) are dynamics of the diffuse scattering intensity."}
```

First and foremost, note that negative-going features are only visible in the vicinity of Bragg peaks. This is due exclusively to the transient Debye-Waller effect [@Ligges2009]. As photodeposited energy transfers from the electrons to the lattice, average real-space disorder due to phonons lowers the overall symmetry of the lattice, which in turn results in smaller Bragg peaks. This is analogous to the (static) thermal Debye-Waller effect [@Debye1913; @Waller1923]. This leaves only positive-going features everywhere away from Bragg peaks, which must be diffuse in nature.

Second, the structure visible across the Brillouin zone is different for reflections that ought to represent the same physical lattice waves. For example, consider the $(010)$ and $(020)$ reflections, \SI{500}{\femto\second} after photoexcitation, shown in @fig:graphite-ueds-zoomed. Around $(010)$, diffuse intensity increases near two of the $\vec{K}$ points, while diffuse intensity is increasing along steaks between $\vec{\Gamma}$ and $\vec{K}$ near $(020)$. This drastic difference between two physically-equivalent reflections can be explained by the difference in the one-phonon structure factors $|F_{1j}(\vec{q}, \tau)|^2$.

```{#fig:graphite-ueds-zoomed .matplotlib file="figures/graphite/ueds-zoomed.py" caption="Comparison of the diffuse intensity change after \SI{100}{\pico\second} for two Brillouin zones, $(010)$ and $(020)$. The colormap scaling is identical to @fig:graphite-ueds. The difference between those two images can be explained by the difference in one-phonon structure factors $|F_{1j}(\vec{q}, \tau)^2|$."}
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

### Transient Debye-Waller factors {#sec:graphite-debye-waller}

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

The calculation for the one-phonon structure factors $|F_{1j}(\vec{q}, \tau<0)|^2$ was carried out from (TODO: ONEPH EQUATION). As discussed in the previous section, the dynamics in one-phonon structure factors due to photoexcitation via the change in Debye-Waller factors can be neglected in this work. Therefore, the calculations of this section assume room-temperature thermal distribution of lattice waves; this is implied by time-delay $\tau<0$ (pre-photoexcitation). The one-phonon structure factors for longitudinal modes LA, LO1 - LO3 are shown @fig:graphite-oneph-longitudinal. The one-phonon structure factors for transverse modes TA, TO1 - TO3 are shown @fig:graphite-oneph-transverse. Out-of-plane modes yield trivial one-phonon structure factors because $\vec{q} \cdot \vec{e}_{j,s}(\vec{k}) = 0$ for all scattering vectors $\vec{q}$ probed by the experiments presented herein.

```{.matplotlib #fig:graphite-oneph-longitudinal file="figures/graphite/oneph-long.py" caption="Calculated one-phonon structure factors $|F_{1j}(\vec{q}, \tau<0)|^2$ of in-plane longitudinal modes at \SI{300}{\kelvin}, for scattering vectors $\vec{q}$ equivalents to the detector area shown in @fig:graphite-ueds. See @fig:graphite-oneph-transverse for the transverse modes equivalent."}
```

```{.matplotlib #fig:graphite-oneph-transverse file="figures/graphite/oneph-trans.py" caption="Calculated one-phonon structure factors $|F_{1j}(\vec{q}, \tau<0)|^2$ of in-plane transverse modes at \SI{300}{\kelvin}, for scattering vectors $\vec{q}$ equivalents to the detector area shown in @fig:graphite-ueds. See @fig:graphite-oneph-longitudinal for the longitudinal modes equivalent."}
```

The one-phonon structure factors are display complex structure in reciprocal space. The structure for a particular mode can be thought of as selection rule: regions in $\vec{q}$ space where the one-phonon structure factor for mode $j$ is large are regions where phonon mode $j$ contributes importantly to diffuse intensity. The structure of $|F_{1j}|^2$ is most determined by the values of terms like $\left\{ \vec{q} \cdot \vec{e}_{j,s}(\vec{k}) \right\}$. For acoustic modes LA and TA, the structure of the one-phonon structure factor is easier to understand near $\vec{\Gamma}$. For example, $|F_{1j}|^2$ is highest for the LA mode in the radial $\vec{q}$ direction because the polarization vectors of those modes are parallel to $\vec{q}$; this fact defines longitudinal waves. For transverse waves, where atomic motion is perpendicular to the lattice wave propagation direction, $|F_{1j}|^2$ is higher in the azimuthal direction.

### Weighted phonon dispersion

An alternative visualization for one-phonon structure factors are weighted dispersion curves, as shown in @fig:graphite-weighted-dispersion. This way of looking at $|F_{1j}|^2$ clearly highlights the difference across Brillouin zones. For example, the difference $|F_{1j=\text{TA}}|^2 - |F_{1j=\text{LA}}|^2$ is either negative ($\vec{q} \approx \vec{\Gamma}_{(010)}$) or positive ($\vec{q} \approx \vec{\Gamma}_{(\bar{1}10)}$). The variation of $|F_{1j}|^2$ explains the diffuse intensity difference shown in @fig:graphite-ueds-zoomed.

```{.matplotlib #fig:graphite-weighted-dispersion file="figures/graphite/weighted-dispersion.py" caption="Calculated one-phonon structure factors visualized as weighted dispersion curves for selected in-plane modes. The color saturation of dispersion curves is proportional to $|F_{1j}|^2$ of the associated mode. Equivalent paths in the Brillouin zone around two reflections are shown to highlight the high degree of reciprocal space structure: $(010)$ on the left and $(\bar{1}10)$ on the right. The geometry of the paths with respect to $(000)$ are shown in the inset on the lower left."}
```

By examining @fig:graphite-weighted-dispersion, it appears that at certain locations, the one-phonon structure factor for one phonon mode dominates; for example, $|F_{1j=\text{LA}}|^2$ near $\vec{q} \approx \vec{\Gamma}_{(010)}$). It might therefore be tempting to attribute the diffuse intensity dynamics at that location exclusively to one mode. If enough such locations in reciprocal space existed, it would be possible to bypass the energy-insensitivity of UEDS and extract mode-dependent phonon dynamics at special points in the Brillouin zone; this idea forms the basis of previous work on graphite [@Stern2018]. However, the information presented by @fig:graphite-weighted-dispersion is incomplete. 

### Relative mode contributions

In order to compare the contribution of a particular phonon mode $j$ on UEDS data, the following weight can be defined:
$$
    W_j(\vec{q}, \tau) \equiv \frac{\left| F_{1j}(\vec{q}, \tau)\right|^2}{\omega_j(\vec{k}, \tau)}
$${#eq:graphite-ueds-weight}
With the weight definition of @eq:graphite-ueds-weight, the change in diffuse intensity (TODO: ONEPH EQUATION) can be re-written as:
\begin{align}
    \Delta I_1(\vec{q}, \tau) &\equiv I_1(\vec{q}, \tau) - I_1(\vec{q}, \tau < 0) \nonumber \\
                              &= \left( n_j(\vec{q}, \tau) + \sfrac{1}{2} \right) W_j(\vec{q}, \tau) - \left( n_j(\vec{q}, \tau < 0) + \sfrac{1}{2} \right) W_j(\vec{q}, \tau < 0) \nonumber \\
                              &\approx \left[ n_j(\vec{q}, \tau) - n_j(\vec{q}, \tau < 0)\right] W_j(\vec{q}, \tau < 0)
\end{align}
In the last line, the quasi-static nature of the one-phonon structure factors $|F_{1j}(\vec{q}, \tau)|^2 \approx |F_{1j}(\vec{q}, \tau < 0)|^2$ was used as discussed in @sec:graphite-debye-waller. The approximation $\omega_j(\vec{k}, \tau) \approx \omega_j(\vec{k}, \tau<0)$, valid for $\vec{k}$ away from $\vec{\Gamma}$, was also used; see @sec:graphite-phonon-landscape for a discussion.

@fig:graphite-oneph-majority shows where, across the reciprocal space area equivalent to @fig:graphite-ueds, a single new phonon from a particular branch contributes more than 50\% or more than 75\% of the diffuse intensity increase. In other words, @fig:graphite-oneph-majority shows where *any* branch $j$, with associated weight $W_j(\vec{q}, \tau<0)$, contributes more than 50\% or 75\% of the total weight $\sum_j W_j(\vec{q}, \tau<0)$. The calculation presented in @fig:graphite-oneph-majority are categorical: dynamics in the transient diffuse intensity **cannot** be associated with the dynamics of a single phonon branch. Furthermore, @fig:graphite-oneph-majority shows that at in-plane high symmetry points $\vec{K}$ and $\vec{M}$, no phonon branch contributes to more than 50\% of the associated dynamics in the diffuse intensity.

In order to access the ultrafast phonon dynamics in a material with UEDS, a more comprehensive analysis of UEDS must be devised -- making good use of one-phonon structure factors.

```{.matplotlib #fig:graphite-oneph-majority file="figures/graphite/oneph-majority.py" caption="Locations in reciprocal space where a single new phonon from a particular branch contributes to more than **a)** 50\% and **b)** 75\% of the associated increase in diffuse intensity."}
```
 
## Phonon spectroscopy across the Brillouin zone{#sec:graphite-ph-spectroscopy}

In this section, the calculation of the one-phonon structure factors will be used to endow UEDS with energy resolution. First, recall that the ultrafast change in diffuse intensity can be expressed as follows:
$$
    \frac{\Delta I(\vec{q}, \tau)}{N_c I_e} = \sum_j \frac{\Delta n_j(\vec{k}, \tau)}{\omega_j(\vec{k}, \tau<0)} \left| F_{1j}(\vec{q}, \tau<0)\right|^2
$${#eq:graphite-ueds-change}
@eq:graphite-ueds-change is approximately valid for $\vec{k}$ away from $\vec{\Gamma}$ (see @sec:graphite-phonon-landscape). Effectively, this restriction is trivial as diffraction peaks drown out the UEDS signals near $\vec{\Gamma}$. 

The sum over all phonon branches $j$ encodes the lack of energy resolution in UEDS experiments. However, by expressing with enough data redundancy, the contribution of every mode $j$ to @eq:graphite-ueds-change can be extracted. Let $\left\{ \vec{H}_1, ..., \vec{H}_M \right\}$ be in-plane reflections. Then, the transient phonon population of mode $j$, $\Delta n_j(\vec{k}, \tau)$, solves the following linear system:
$$
    \vec{I}_{\vec{k}}(\tau) = \vec{F}_{\vec{k}} \vec{D}_{\vec{k}}(\tau)
$${#eq:graphite-vectorized-ueds}
where
\begin{align}
 \vec{I}_{\vec{k}}(\tau) &= 
 	\frac{1}{N_c I_e}
	\begin{bmatrix}
    	\Delta I(\vec{k} + \vec{H}_1, t) & \dots & \Delta I(\vec{k} + \vec{H}_M, t)
	\end{bmatrix}^T \\
\vec{F}_{\vec{k}} &= 
	\begin{bmatrix}
		|F_{11}(\vec{k} + \vec{H}_1, \tau<0)|^2 & \dots  & |F_{1N}(\vec{k} + \vec{H}_1,\tau<0)|^2 \\
		\vdots		   			 			 & \ddots & \vdots					            \\
		|F_{11}(\vec{k} + \vec{H}_M, \tau<0)|^2 & \dots  & |F_{1N}(\vec{k} + \vec{H}_M,\tau<0)|^2
	\end{bmatrix} \\
\vec{D}_{\vec{k}}(\tau) &=
	\begin{bmatrix}
    	\Delta n_1(\vec{k}, \tau)/\omega_1(\vec{k}, \tau<0) & \dots & \Delta n_N(\vec{k}, \tau)/\omega_N(\vec{k}, \tau<0)
	\end{bmatrix}^T
\end{align}

This linear system of equations can be solved numerically provided enough experimental data, i.e. diffuse intensity for at least $M \geq N$ distinct Brillouin zones. The choice to solve for $\vec{D}_{\vec{k}}(\tau)$ rather than $\Delta n_j(\vec{k}, \tau)$ directly, comes down to the degree of confidence that should be placed in the approximations that were made to get to @eq:graphite-ueds-change. Phonon polarization vectors are mostly determined by the symmetries of the crystal, which are fixed throughout the experiments in this chapter. On the other hand, phonon vibrational frequencies in general might be influenced by non-equilibrium carrier and phonon distributions as would be the case with strongly anharmonic crystals (e.g. SnSe, see @sec:snse). Solving for the ratio of phonon population to vibrational frequency is more robust against the weaknesses of this modeling because the one-phonon structure factors only take into account the polarization vectors.

### Numerical procedure

The procedure used to numerically solve for $\vec{D}_{\vec{k}}(\tau)$ is described henceforth.

All visible Brillouin zones in the measurement were used to generate $\vec{I}_{\vec{k}}(\tau)$. Since the scattering patterns have been symmetrized (@fig:graphite-static), one would expect that using the intensity data for reflections related by symmetry would be redundant; however, the following procedure worked better hwn using the entire area of the detector. This is no doubt due to minute misalignment of the scattering patterns and uncertainty in detector position, which are averaged out when using all available data. The Brillouin zones associated with all in-plane reflections $\left\{ \vec{H} \middle| |\vec{H}| \leq \SI{12}{\per\angstrom}\right\}$ were used, for a total of 44 Brillouin zones. Therefore, at every $\vec{k}$ point, $\vec{I}_{\vec{k}}(\tau)$ and $\vec{F}_{\vec{k}}$ have shapes of $(44, 1)$ and $(44, 8)$, respectively.

The possible values for elements of $\vec{D}_{\vec{k}}(\tau)$ were constrained to be non-negative, implying that $\Delta n_j(\vec{k}, \tau) \geq 0 ~ \forall \tau$. In other words, it is assumed that the population of a phonon branch cannot drop below its equilibrium level. While not strictly necessary, it was found to lead to a more stable solution. This constraint allowed for the use of the non-negative approximate matrix inversion algorithm [@Lawson1995] to solve @eq:graphite-vectorized-ueds. This procedure was repeated for every reduced wavevector $\vec{k}$ and time-delay $\vec{\tau}$. Stable solutions were found for $|\vec{k}| \geq \SI{0.45}{\per\angstrom}$, away enough from $\vec{\Gamma}$ and elastic scattering signals.

It must be emphasized that the numerical solution to $\vec{D}_{\vec{k}}(\tau)$ is not the result of fitting (iterative least-squares method), but rather a non-iterative approximate matrix inversion based on linear least-squares. The method described in this section admits no free parameter, other than the phonon vibrational frequencies and polarization vectors.

### Population dynamics{#sec:graphite-pop-dynamics}

The numerical solution for @eq:graphite-vectorized-ueds is shown in @fig:graphite-ph-populations. The transient phonon populations $\left\{ \Delta n_j(\vec{k}, \tau)\right\}$ across the Brillouin zone are shown for three important modes: TO2, TA, and LA. 

```{.matplotlib #fig:graphite-ph-populations file="figures/graphite/decomp.py" caption="Measurement of the change in transient phonon population $\Delta n_j(\vec{k}, \tau)$ following photoexcitation for relevant in-plane modes of graphite across the Brillouin zone. The solution domain is bound by white circle at $|\vec{k}| \leq \SI{0.45}{\per\angstrom}$, and by a solid white hexagon at the Brillouin zone edge. The Brillouin zone midpoint is highlighted with a dashed white hexagon. The location of the $A_1^\prime$ mode is shown in the top row."}
```

While previous work by the author discussed the nonthermal phonon dynamics qualitatively[@Stern2018], @fig:graphite-ph-populations allows to quantitatively determine how the energy deposited in the electronic system flows and thermalizes. A discussion of the observed physical processes is discussed below.

As discussed in @sec:graphite-100-fs, two optical phonons modes are strongly-coupled to the electronic system: $A_1^\prime$ located near the $\vec{K}$ point, and $E_{2g}$, near $\vec{\Gamma}$. $E_{2g}$ is obscured by the elastic signals near $\vec{\Gamma}$; however, its dynamics are accessible via ultrafast spectroscopy [@Kampfrath2005]. On the other hand, $A_1^\prime$ is clearly visible in @fig:graphite-ph-populations. It behaves as expected: a fast early increase in population is seen due to the transfer of energy from the electronic subsystem. By \SI{5}{\pico\second}, energy has already been transferred away through phonon-phonon coupling. 

After a few picoseconds, the strongly-coupled optical modes decay into lower-energy lattice waves. Since one of the primary causes of phonon-phonon scattering is anharmonicity, the decay probabilities are usually computed via DFT [@Bonini2007]. UEDS shows an experimental determination of those decay pathways. The transfer of energy from high-energy optical phonons to lower-energy lattice waves must satisfy the conservation of momentum and energy. This restricts the number of possible decay pathways to mid-Brillouin zone phonons, either directly (from $E_{2g}$) or via Umklapp scattering (from $A_1^\prime$). Specifically, from \SIrange{1.5}{100}{\pico\second}, acoustic phonon population increase is observed along the $\vec{\Gamma}-\vec{M}$ line. This measurement is in accordance with predicted anharmonic decay probabilities from the $E_{2g}$ mode [@Bonini2007]. The small increase of population of the LA mode around \SI{500}{\femto\second} at $\tfrac{1}{2}\vec{K}$ is also a confirmation of predictions by Bonini *et al*[@Bonini2007].

Over longer time-scales (\SI{25}{\pico\second}), the nonthermal TA population has significantly pooled at $\tfrac{1}{3}\vec{M}$ and $\vec{M}$. There are no three-phonon anharmonic decay processes that start in a transverse mode like TA. Allowed *in-plane* interband transitions are of the types $\text{L} \to \text{T} + \text{T}$ or $\text{L} \to \text{L} + \text{T}$, where L and T represent longitudinal and transverse modes respectively [@Lax1981; @Khitun2001]. It is therefore expected *a priori* that a nonthermal TA population remains until decay to out-of-plane phonon modes occur. Predicted decay probabilities also support this observation; both LA and TA modes will favor decay to out-of-plane modes[@Paulatto2013]  which are not visible to UEDS measurements. A special case is visible at $\tfrac{1}{2}\vec{M}$ and $\tfrac{1}{2}\vec{K}$, where predicted anharmonic lifetimes of TA and LA drop significantly due to the activation of Umklapp scattering to out-of-plane acoustic phonons[@Paulatto2013]. The measurements of @fig:graphite-ph-populations corroborate the predictions: energy flow away from TA and LA at the mid-Brillouin zone is fast enough that population never builds up. 

It's worth noting that apart from certain special cases (TODO: multiph THz)

### Long-term decay

It is possible to quantify how long does the non-equilibrium distribution of lattice waves last in graphite. The energy transferred to the sample from photoexcitation, stored in plane in a mode $j$, can be expressed as follows:
$$
    \Delta E(\tau) = \int_{\text{BZ}} d\vec{k} \frac{\Delta n_j(\vec{k}, \tau)}{\omega_j(\vec{k}, \tau)}
$$
where $\int_{\text{BZ}} d\vec{k}$ is understood to be the integral over the in-plane section of the Brillouin zone. Using the population measurements $\left\{ n_j(\vec{k}, \tau)\right\}$ presented in @fig:graphite-ph-populations, the mode-dependent relative change in stored energy (in-plane) can be calculated. The results are shown in @fig:graphite-long-term.

```{.matplotlib #fig:graphite-long-term file="figures/graphite/energy.py" caption="Energy stored in-plane after photoexcitation. **inset** long-term trend shows that thermalization has not yet occurred by \SI{600}{\pico\second}."}
```

The energy trends for modes TA and TO2 are commensurate with the description of @sec:graphite-pop-dynamics. The fascinating aspect of this analysis is revealed by looking at the total energy stored in-plane. Even by \SI{600}{\pico\second}, thermalization of the energy dumped into the sample by photoexcitation has not occurred at all. 

## Mode-projected electron-phonon and phonon-phonon coupling

TODO: make the distinction between G (heat rate?) and g (matrix element) clear

Electron-phonon and phonon-phonon coupling constants describe the strength of the coupling between excitations in a material. A coupling constant $G_{i,j}$ describes the rate of energy transfer from excitation type $i$ to excitation type $j$, if the temperature of one unit of volume was warmed up by \SI{1}{\kelvin}. E-ph and ph-ph coupling is impossible to measure directly at equilibrium. In this section, UEDS measurements of phonon population dynamics will be used to experimentally determine mode-dependent e-ph and ph-ph coupling terms.

The nonequilibrium flow of energy between excitations, both electronic and lattice in nature, has historically been crudely modelled using the *two-temperature model* [@Allen1987]. In summary, this model states that while the concept of temperature does not apply to nonequilibrium situations, the energy distribution of the electrons and lattice waves may be treated separately. In other words, the se[arate thermalization of the electronic and lattice subsystems is much faster than the energy transver between them. It is evident that such a description does not adequately model experimental results from @sec:graphite-ph-spectroscopy for two main reasons. First, the thermalization of decay of optical phonons to acoustic ones overlaps significantly in time with the transfer of energy from the photoexcited electrons to strongly-coupled optical modes. Second, the energy distribution of lattice waves is far from thermal, even at \SI{100}{\pico\second}.

UEDS measurements allow to move beyond the two-temperature approximation as the energy distribution of phonons is known across the Brillouin zone. And yet, theory has not caught up with this opportunity; the best way to extract e-ph and ph-ph coupling values is still based on rate-equations. To this end, the formalism of the two-temperature model needs to be extended.

### The non-thermal lattice model {#sec:graphite-nlm}

The non-thermal lattice model is an extension of the two-temperature model, with the added flexibility that lattice waves need not be thermalized[@Waldecker2016]. It assumes that the energy distribution of each phonon mode admits a thermal description; that is, the energy transfer between phonons of different branches is slower than the scattering of phonons from the same branch. This assumption is a reasonable one, in the light of the results presented in @sec:graphite-pop-dynamics. Within this framework, each phonon branch $j$ can be assigned its own molar heat capacity, $C_{ph, j}$, and temperature, $T_j$:
$$
C_e(T_e) \frac{\partial T_e}{\partial \tau} = \sum_i G_{ep, i}\left[ T_e(\tau) - T_{ph,i}(\tau) \right] + f(\tau)
$${#eq:graphite-nlm-electrons}
$$
\Bigg\{ C_{ph,j}(T_{ph,j}) \frac{\partial T_{ph,j}}{\partial \tau} = 
	\sum_{i\neq j} G_{ep, i} \left[ T_e(\tau)      - T_{ph,i}(\tau) \right] 
				 + G_{pp,ij} \left[ T_{ph,j}(\tau) - T_{ph,i}(\tau) \right] 	\Bigg\}_{j=1}^{N}
$${#eq:graphite-nlm-phonons}
where $f(\tau)$ is the laser pulse profile, and $C_e$ and $T_e$ are the electronic heat capacity and electron temperature, respectively. As discussed in @sec:graphite-100-fs, the use of an electronic temperature is acceptable for $\tau > \SI{150}{\femto\second}$. The coupling constants $G_{ep, i}$ describe the coupling between the electrons and phonon mode $i$, while coupling constants $G_{pp, ij}$ encode the coupling between phonon modes $i$ and $j$. The constants $G_{i,j}$ are related to electron-phonon and phonon-phonon coupling constants. Their relationship is described further below in @sec:graphite-coupling-constants.

Observations of transient phonon populations are more general than mode temperatures. However, in order to make use of the non-thermal lattice model, the mode temperature can be related to transient phonon mode populations via the Bose-Einstein distribution:
$$
    n_j(\vec{k}, \tau) \propto \left[ \exp{\left( \frac{\hbar \omega_j(\vec{k}, \tau<0)}{k_B T_{ph, j}(\tau)} \right)} - 1\right]^{-1}
$$
The above expression can be decomposed with a Laurent series [@Wunsch2005] to extract the quasi-linear relationship between mode population and temperature:
$$
    n_j(\vec{k}, \tau) \propto  \frac{k_B T_{ph, j}(\tau)}{\hbar \omega_j(\vec{k}, \tau<0)} - 1/2 + \mathcal{O}\left( T^{-1}_{ph, j}(\tau) \right)
$${#eq:graphite-population-laurent}
The above holds for appropriately-high mode temperatures. For the remainer of this section, it follows that $\Delta n_j \propto \Delta T_{ph,j}$, where the initial temperature is known to be \SI{300}{\kelvin}.

### The non-thermal lattice model applied to UEDS

Based on the formalism presented in the @sec:graphite-nlm, the couplings to the $A_1^\prime$ phonon mode will be extracted from the TO2 population measurements. 

Let the differential population be $\Delta n_{j=\text{TO2}}(\vec{k} = \vec{K}, \tau) \equiv \Delta n_{A_1^\prime}(\tau)$. $\Delta n_{A_1^\prime}(\tau)$ is obtained by integrating the population of the TO2 mode in a circular arc of radius \SI{0.3}{\per\angstrom} centered at $\vec{k} = \vec{K}$, as shown in @fig:graphite-ph-populations. The heat capacities of the electronic system and every relevant phonon mode must be parametrized. 

#### Heat capacities

The electronic heat capacity is extracted from measurements by Nihira an Iwata [@Nihira2003]:
$$
C_e(T_e) = 13.8 ~ T_e(\tau) + 1.16 \times 10^{-3} ~ T_e^2(\tau) + 2.6 \times 10^{-7}  ~ T_e^3(\tau)
$$

Extracting the lattice specific heat is simplified by the observation that thermal expansion has not occurred on the time-scale of the measurements ($\tau < \SI{680}{\pico\second}$). Bragg peak positions remain static throughout the measurement, as was also reported by Chatelain *et al*[@Chatelain2014a]. The specific heat capacity of each mode can therefore be taken as the heat capacity at constant volume[@Ziman1979]:
$$
C_{ph,j}(T_{ph,j}) = 
k_B \int_0^{\omega_D} d\omega ~ D_j(\omega) 
    \left( 
        \frac{\hbar \omega}{k_B T_{ph,j}} 
    \right)^2 
    \frac{e^{\hbar \omega / k_B T_{ph,j}}}{\left( e^{\hbar \omega / k_B T_{ph,j}} - 1\right)^2}
$$
where $\omega_D$ is the Debye frequency, and $D_j(\omega)$ is the density of states for phonon mode $j$. The momentum-resolution of UEDS allows for a simplification: a single frequency contributes to the density of states, so that $D_j(\omega) \equiv \delta(\omega)$.

#### Reducing complexity with momentum-resolution {#sec:graphite-reducing-complexity}

The number of coupled equations in @eq:graphite-nlm-phonons can be reduced by aggregating lattice heat capacities into two categories: the heat capacity of $A_1^\prime$, $C_{A_1^\prime}$, and the total effective heat capacity of all other relevant modes, defined as $C_l$. The calculation of $C_l$ boils down to adding the contribution of mode that can scatter into $A_1^\prime$ or from it, based on conservation of energy and momentum, weighted by the decay probabilities reported by Bonini *et al* [@Bonini2007]. With this information:
\begin{align}
    C_l &= \frac{9}{100} \left[ C_{ph,j=\text{TA}} + C_{ph,j=\text{TA}} \right] \\
        &+ \frac{36}{100} \left[ C_{ph,j=\text{TA}} + C_{ph,j=\text{LA}}\right] \nonumber \\
        &+ \frac{55}{100} \left[ C_{ph,j=\text{TA}} + C_{ph,j=\text{LA}} + C_{ph,j=\text{LO}} + C_{ph,j=\text{LA}}\right] \nonumber
\end{align}
The system of equations from @sec:graphite-nlm can then be expressed as three equations: the flow of energy in and out of the electronic system, the $A_1^\prime$ mode, and all other modes coupled to the $A_1^\prime$ mode:
$$
\left\{
    \begin{array}{rcl}
        C_e(T_e) \frac{\partial T_e}{\partial \tau} & = & f(\tau) \\ 
                                                    & - & G_{e,A_1^\prime}  ~ \left[ T_e(\tau) - T_{A_1^\prime}(\tau) \right] \\
                                                    & - & G_{e, l}    ~ \left[ T_e(\tau) - T_l(\tau) \right] \\
                                                    $ ~ $ \\
        C_{A_1^\prime}(T_{A_1^\prime}) \frac{\partial T_{A_1^\prime}}{\partial \tau} & = & G_{e,A_1^\prime}  ~ \left[ T_e(\tau) - T_{A_1^\prime}(\tau) \right] \\
                                                      & - & G_{A_1^\prime, l} ~ \left[ T_{A_1^\prime}(\tau) - T_{l}(\tau) \right] \\
                                                      $ ~ $ \\
        C_{l}(T_l) \frac{\partial T_l}{\partial \tau} & = & G_{e,l} ~ \left[ T_e(\tau) - T_l(\tau) \right] \\ 
                                                       & + & G_{A_1^\prime, l} ~ \left[ T_{A_1^\prime}(\tau) - T_{l}(\tau)\right] 
    \end{array}
\right\}
$${#eq:graphite-nlm-system}

### Heat rates {#sec:graphite-eph-solution}

The solution to @eq:graphite-nlm-system was computed using an iterative least-squares method[@Branch1999]. The resulting temperature traces $\left\{ T_e(\tau), T_{A_1^\prime}(\tau), T_{l}(\tau)\right\}$ are shown in [@fig:graphite-eph-coupling]. The coupling constants $\left\{ G_{e,l}, G_{e,A_1^\prime}, G_{A_1^\prime, l}\right\}$ are listed in @tbl:graphite-eph-coupling

```{.matplotlib #fig:graphite-eph-coupling file="figures/graphite/eph-coupling.py" caption="Evolution of the $a_1^\prime$ mode population in graphite after ultrafast photoexcitation. Transient population $\Delta n_{A_1^\prime}(\tau)$ is shown in black (circles). Error bars represent the standard error in the population mean before photoexcitation $\tau < 0$. The biexponential fit to the transient population is shown in pink (solid). The effective temperature of the modes that $A_1^\prime$ can decay into is shown in orange (dotted). **Inset** Temperature dynamics at early times ($\tau < \SI{1}{\pico\second}$) show the thermalization between the electronic system (purple, dashed) and the $A_1^\prime$ mode is very fast, indicative of strong electron-phonon coupling. The traces from the main figure are shown in the inset as well."}
```

\begin{table}
	\centering
	\caption{Coupling strength between electronic system, the $A_1^\prime$ phonon, and the lattice system. Uncertainty is derived from fit covariances.}
	\vspace{2mm}
	\begin{tabular}{c | c}
		~ & Coupling strength [\si{\watt \per \meter \cubed \per \kelvin}] \\ \hline\hline
		$G_{e,A_1^\prime}$  & $(6.8 \pm 0.3) \times 10^{17}$ \\ \hline
		$G_{A_1^\prime, l}$ & $(8.0 \pm 0.5) \times 10^{17}$ \\ \hline
		$G_{e,l}$           & $(0.0 \pm 6.0) \times 10^{15}$ \\ \hline
	\end{tabular} 
	\label{tbl:graphite-eph-coupling}
\end{table}

### Mode-projected electron-phonon coupling matrix elements {#sec:graphite-coupling-constants}

In order to compare to theory and other experiments, the value of electron-phonon coupling matrix-element $g^2_{e, A_1^\prime}$ is calculated from the coupling constants determined from solving @eq:graphite-nlm-system.

In general, the electron-phonon matrix element $g^2_{e,j}(\vec{k})$ between the electronic system and phonon mode $j$ at wavevector $\vec{k}$ is most simply related to the relaxation time $\tau_{e,j}(\vec{k})$ between the two subsystems[@Na2019]:
$$
\frac{\hbar}{\tau_{e,j}(\vec{k})} = 2 \pi \langle g^2_{e,j}(\vec{k}) \rangle D_e(\hbar \omega_{\nu} - \hbar \omega_{j}(\vec{k}))
$$
where $D_e(\epsilon)$ is the electronic density-of-states, $\hbar \omega_{\nu}$ is the optical excitation energy (\SI{1.55}{\electronvolt} in the case of \SI{800}{\nano\meter} light), and $\omega_j(\vec{k})$ is the vibrational frequency of phonon mode $j$ at wavevector $\vec{k}$, as defined previously. Given the nature of the experiments presented here, an approximation to the electronic density of states for graphene close to the Dirac point can be used[@Neto2009]:
$$
    D_e(\epsilon) = \frac{2 A}{\pi} \frac{|\epsilon|}{(\hbar v_{F})^2}
$$
where $A$ is the unit cell area and $v_F = \SI{9.06e5}{\meter \per \second}$ is the Fermi velocity[^2]. It follows that the determination of the mode-dependent electron-phonon coupling matrix element $g^2_{e,j}(\vec{k})$ relies on the calculation of the mode-dependent relaxation time $\tau_{e,j}(\vec{k})$ based on UEDS measurements. 

The calculation for the electron-$A_1^\prime$ coupling matrix element $g^2_{e,A_1^\prime} \equiv g^2_{e,j=\text{TO2}}(\vec{k} \approx \vec{K})$ is demonstrated below. Consider the following sum of @eq:graphite-nlm-system:
\begin{align}
    & \frac{\partial T_e}{\partial \tau} - \sum_j \frac{\partial T_{ph,j}}{\partial \tau} = \label{eq:graphite-tau-1} \\
    & \sum_j \Bigg[ \frac{G_{ep,j}}{C_e} ~ (T_e - T_{ph,j}) \nonumber  -\sum_i \bigg( \frac{G_{ep,i}}{C_{ph,j}} ~ (T_e - T_{ph,i}) + \frac{G_{pp,ij}}{C_{ph,j}} ~ (T_{ph,i} - T_{ph,j}) \bigg) \Bigg] \nonumber
\end{align}
At early times ($< \SI{5}{\pico\second}$), the $A_1^\prime$-phonon coupling $G_{A_1^\prime, l}$ is negligible compared to other coupling constants (@tbl:graphite-eph-coupling). Then, @eq:graphite-tau-1 can be simplified to:
$$
	\frac{\partial T_e}{\partial \tau} - \sum_j \frac{\partial T_{ph,j}}{\partial \tau} = \sum_j \frac{G_{ep,j}}{C_e} (T_e - T_{ph,j}) - \sum_{i,j} \frac{G_{ep,i}}{C_{ph,j}}(T_e - T_{ph,i})
$${#eq:graphite-tau-2}
By performing a substitution $\lambda = T_e - \sum_j T_{ph,j}$, the equation above simplifies to a familiar situation:
$$
    \dot{\lambda}(\tau) - a(\tau) \lambda(\tau) = 0
$${#eq:graphite-tau-3}
where
$$
    a(\tau) = \sum_j \left( \frac{G_{ep,j}}{C_e} - \sum_i \frac{G_{ep,i}}{C_{ph,j}}\right).
$$
The time dependence comes from the time-evolution of the individual temperatures. In the case of phonon temperatures, the phonon population dynamics are directly related to temperature dynamics according to @eq:graphite-population-laurent. @eq:graphite-tau-3 is a separable equation with solution:
$$
\lambda(\tau) = \exp{\int d\tau \left[ a(\tau) \right]}.
$$
For a slow-varying integrand $a(\tau) \approx a$, then $a = 1/\bar{\tau}$, where $\bar{\tau}$ is a compound variable representing the relaxation of the system. This leads to the following form:
$$
\frac{1}{\bar{\tau}} \approx \sum_j \left( \frac{G_{ep,j}}{C_e} - \sum_i \frac{G_{ep,i}}{C_{ph,j}}\right).
\label{eq:graphite-tau-4}
$$
As a specific example, the above expression reduces nicely in the case of the two-temperature model, where all phonon modes are considered to be thermalized with each other, with isochoric heat capacity $C_{ph}$:
$$
\frac{1}{\bar{\tau}} = G_{ep} \left( \frac{1}{C_e} - \frac{1}{C_{ph}}\right)
$$
and we see that $\bar{\tau}$ physically represents the relaxation time of the electronic system into the lattice. @eq:graphite-tau-4 can be thought of as a sum of relaxation times between the electronic subsystem and specific modes $\tau_{e,j}$:
$$
\frac{1}{\tau_{e,j}} = \frac{G_{ep,j}}{C_e} - \sum_i \frac{G_{ep,i}}{C_{ph,j}}.
$$

Using the coupling constants and heat capacities from sec:graphite-eph-solution and @sec:graphite-reducing-complexity respectively, $\tau_{e, A_1^\prime} = \SI{106 \pm 11}{\femto\second}$ and $\langle g^2_{e, A_1^\prime} \rangle = \SI{0.035 \pm 0.001}{\square\electronvolt}$. This value is compared to other references in @tbl:graphite-eph-coupling-comparison.

\begin{table}
	\centering
	\caption{Comparison of measured and calculated values for the electron-phonon coupling matrix element $\langle g^2_{e, A_1^\prime}\rangle$}
	\vspace{2mm}
	\begin{tabular}{l | c | l}
		Source & $\langle g^2_{e, A_1^\prime} \rangle$ [\si{\electronvolt\squared}] & Notes \\ \hline\hline
        This work, René de Cotret \emph{et al}\autocite{RenedeCotret2019} & $0.035 \pm 0.001$ & Experiment \\ \hline
        Johannsen \emph{et al}\autocite{Johannsen2013}                    & $0.033 \pm 0.007$ & Experiment (trARPES) on graphene \\ \hline
        \multirow{2}{*}{Na \emph{et al}\autocite{Na2019}}                 & $0.050 \pm 0.011$ & Experiment (trARPES) \\ \cline{2-3} 
                                                                          & $0.040$           & Theory \\ \hline
	\end{tabular} 
	\label{tbl:graphite-eph-coupling-comparison}
\end{table}

## Conclusion and outlook

TODO: We also note that the procedure presented can be easily extended to (equilibrium) thermal diffuse scattering measurements, where the phonon populations are known at constant  temperature, but the phonon vibrational spectrum is unknown. Therefore, using pre-photoexcitation data of a time-resolved experiment, one could infer the phonon vibrational frequencies, which are then used to determine the change in populations using the measurements after photoexcitation. This scheme only relies on the determination of phonon polarization vectors.

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]

[^1]: Modified from Chatelain [@Chatelain2014].
[^2]: Note that the factor of $\hbar$ has been erroneously ignored by Castro Neto *et al*[@Neto2009].