
# Momentum-resolved excitation couplings in graphite {#sec:graphite}

This chapter is an exploration of ultrafast electron diffuse scattering measurements applied to systems with strong couplings. To this end, graphite is the perfect benchmark system for ultrafast electron diffuse scattering measurements for two major reasons. On the one hand, it has a very stiff lattice in the plane due to inherently strong carbon-carbon bonds and a hexagonal structure[@Bernal1924]. The spectrum of lattice waves in graphite gives ultrafast electron diffuse scattering measurements incredible contrast, as the thermally-occupied phonon modes are localized to the Brillouin zone center[@Tuinstra1970]. On the other hand, much ultrafast work has been done on graphite and graphene. This presents an opportunity to compare and contrast results between spectroscopic probes such as terahertz spectroscopy or time- and angle-resolved photoemission spectroscopy, and ultrafast electron diffuse scattering.

The work in this chapter follows earlier work by the author in Stern *et al.*[@Stern2018] which demonstrated that ultrafast electron diffuse scattering measurements were feasible in graphite. A picture of anharmonic decay from strongly-coupled optical phonon modes was assembled, which meshed very well with the existing literature. However, due to the energy-integrative nature of ultrafast electron diffuse scattering measurements, the analysis remained qualitative. This chapter is focused on the quantitative side of things. It will be demonstrated that the geometry of experiments, combined with the crystal symmetry, can be used to endow ultrafast electron diffuse scattering measurement with energy-resolution. The chapter culminates in the presentation of the time-resolved mode-dependent phonon population *across the Brillouin zone*, describing the energy flow within the graphite lattice.

This chapter is organized as follows. The properties of single-crystal graphite are first described, as well as a brief review of the ultrafast experiments on graphene and graphite. The experimental details are then given, as well as the computational parameters that were used to calculate the in-plane phonon dispersion. The ultrafast electron diffuse scattering measurements are then presented. Next comes the calculation of the geometrical weights known as the one-phonon structure factors (@sec:diffuse-scattering). The high point of this chapter follows, with the presentation of time-resolved mode-dependent phonon populations across the Brillouin zone. Finally, the phenomenological electron-phonon coupling tensor element to the $A_1^\prime$ phonon mode is extracted and compared to other measurements and calculations. 

## Single-crystal graphite

Single-crystal graphite is a two-dimensional material made of a repeating layers of carbon atoms arranged in a hexagonal lattice. Its lattice vectors $\set{ \vect{a}_i }$ are:
$$
\begin{pmatrix}
    \vect{a}_1 \\
    \vect{a}_2 \\
    \vect{a}_3
\end{pmatrix}
=
\begin{pmatrix}
    \sqrt{3}a  & -a & 0 \\
    0          & 2a & 0 \\
    0          &  0 & c 
\end{pmatrix}
\begin{pmatrix}
    \vect{e}_1 \\
    \vect{e}_2 \\
    \vect{e}_3
\end{pmatrix}
\label{eq:graphite-lattice}
$$
where $a=\SI{1.232}{\angstrom}$, $c=\SI{6.711}{\angstrom}$, and $\set{ \vect{e}_i }$ are understood to be the usual Euclidean vectors. Carbon atoms $\set{ \vect{c}_i }$ are positioned as follows:
$$
\begin{pmatrix}
    \vect{c}_1 \\
    \vect{c}_2 \\
    \vect{c}_3 \\
    \vect{c}_4 
\end{pmatrix}
=
\begin{pmatrix}
    0            & 0            & 0            \\
    \sfrac{1}{3} & \sfrac{2}{3} & 0            \\
    0            & 0            & \sfrac{1}{2} \\
    \sfrac{2}{3} & \sfrac{1}{3} & \sfrac{1}{2} \\
\end{pmatrix}
\begin{pmatrix}
    \vect{a}_1 \\
    \vect{a}_2 \\
    \vect{a}_3
\end{pmatrix}
$$
Note that carbon atoms $\vect{c}_1$ and $\vect{c}_2$ form a sheet of graphene, while $\vect{c}_3$ and $\vect{c}_4$ form a separate separate sheet of graphene rotated by $\tfrac{\pi}{3}$ (\ang{60}). The stacking of both sublattices along the $\vect{a}_3$ axis confers graphite increased discrete rotational symmetry (six-fold) about the stacking axis, compared to the 3-fold rotational symmetry of graphene. More specifically, the point group for graphite is $6/mmm$, while the Hermann-Mauguin symbol for the space group is $P6_3/mmc$. 

The real-space six-fold rotational symmetry is mirrored in reciprocal space. The in-plane geometry of the Brillouin zone is shown in @fig:graphite-bz, including the high-symmetry points  $\vect{\Gamma}$ (zone-center), $\vect{M}$, and $\vect{K}$.

```{#fig:graphite-bz .matplotlib file="figures/graphite/bz.py" caption="In-plane section of the Brillouin zone of graphite, spanned by reciprocal lattice vectors $\vect{b}_1$ and $\vect{b}_2$. The three classes of high-symmetry points are shown. $\vect{M}$ and $\vect{K}$ have symmetry equivalents on every edge and vertex respectively."}
```

### Electronic structure

The electronic structure of graphite has been the subject of much study. Graphite is a very close relative of the Dirac semi-metal graphene; its electronic dispersion is very similar from the point of view of experiments presented in this chapter. The electronic dispersion of the $\pi$ and $\pi^{\star}$ bands of graphite can be calculated in the tight-binding framework [@Wallace1947;@Slonczewski1958]. Taking into account nearest-neighbor interactions only, the dispersion for both bonding and antibonding bands are given by:
$$
    E(\vect{k})_{\pi^{\star}} = t \sqrt{3 + f(\vect{k})}
    \label{eq:graphite-estructure-pistar}
$$
$$
    E(\vect{k})_{\pi} = - t \sqrt{3 + f(\vect{k})}
    \label{eq:graphite-estructure-pi}
$$
with
\begin{align}
    f(\vect{k}) &= 2  \cos\left[\sqrt{3} a ~ (\vect{k} \cdot \vect{b}_2) \right] \\
               &+ 4  \cos\left[\tfrac{\sqrt{3} a}{2} ~ (\vect{k} \cdot \vect{b}_2) \right] 
                    \cos\left[\tfrac{3 a}{2} ~ (\vect{k} \cdot \vect{b}_1) \right] \nonumber
\end{align}
where $t$ is the tight-binding constant and $a$ is the carbon-carbon distance. The in-plane dispersion for a tight-binding energy of $t=\SI{2.7}{\electronvolt}$ is shown on @fig:graphite-electronic-structure. More complete expressions that are not limited to nearest neighbor interactions are given elsewhere [@Partoens2006;@Neto2009], but the essential details are captured by @eq:graphite-estructure-pistar and @eq:graphite-estructure-pi.

```{.matplotlib #fig:graphite-electronic-structure file="figures/graphite/estructure.py" caption="In-plane electronic dispersion $E(\vect{k})$ for graphite. The in-plane section of the Brillouin zone is shown below. The dotted black line represents the line-cut used for @fig:graphite-photoexcitation."}
```

The nonequilibrium behavior of photoexcited graphite is dominated by the structure of the electronic dispersion near the $\vect{K}$ points, at the corner of the Brillouin zone. These points, where the dispersion of each band meet, are called Dirac points. In analogy to graphene[@Novoselov2005], the dispersion in graphite adopts a conical shape, called Dirac cones, near the Dirac points[@Zhou2006]. 

### Phonon landscape {#sec:graphite-phonon-landscape}

With four atoms per unit cell, the structure of graphite supports twelve distinct phonon modes. The layered nature of graphite results in a clear distinction between in-plane and out-of-plane phonon modes. These modes are longitudinal modes LA, LO1 - LO3, in-plane transverse modes TA, TO1 - TO3, and out-of-plane traverse modes ZA, ZO1 - ZO3. The phonon dispersion relation for in-plane modes was calculated as described further below in @sec:computational-details, and is shown in @fig:graphite-static-dispersion. The out-of-plane transverse modes are not taken into account as the experiments presented in this chapter are not sensitive to them. 

@fig:graphite-static-dispersion reveals why graphite is the *perfect* benchmark material for ultrafast electron diffuse scattering. All of the lattice modes are so stiff (high-energy) that, at room temperature, thermally-occupied modes are concentrated at $\vect{\Gamma}$. The energy equivalent to room-temperature is shown on @fig:graphite-static-dispersion by a horizontal dashed line at \SI{25}{\milli\electronvolt}. Therefore, the contrast between static diffuse intensity and diffuse intensity after photoexcitation is bound to be maximal.

```{#fig:graphite-static-dispersion .matplotlib file="figures/graphite/static-dispersion.py" caption="Phonon dispersion relation of graphite for in-plane modes LA, TA, and two-fold degenerate modes LO and TO. The path in reciprocal space is shown in the center. A horizontal dashed line at \SI{25}{\milli\electronvolt} indicates the average energy stored in the phonon modes at room temperature (\SI{300}{\kelvin})."}
```

### Kohn anomalies

An important feature of the phonon dispersion of graphite that is not replicated by the simple calculation of @fig:graphite-static-dispersion are *Kohn anomalies* [@Kohn1959]. Kohn anomalies are pronounced dips with a $|E - E_0|$ character in the phonon dispersion of metals. In other words, the derivative of the phononic dispersion curve is discontinuous. In extreme cases, the softening of the phonon dispersion to zero creates a lattice distortion, which is one of the causes of charge-density wave phases[@Williams1974;@Johannes2008;@Otto2021].

Kohn anomalies have been shown to appear in graphite for transverse optical modes at $\vect{\Gamma}$ (mode $E_{2g}$) and $\vect{K}$ (mode $A_1^\prime$) by inelastic x-ray scattering measurements performed by Maultzsch *et al.* [@Maultzsch2004] Theoretical work by Piscanec *et al.* [@Piscanec2004] has linked the slope of these kinks to the strength of the coupling between the electronic system and these specific modes. The calculated band structure including the two Kohn anomalies in graphite is presented in @fig:graphite-kohn-anomalies.

<!-- The figure below is raw latex because an alternate caption without citations needs to be use so that it looks OK in the list of figures-->
\begin{figure}
    \centering
    \includegraphics{images/kohn.pdf}
    \caption[Kohn anomalies in the phonon dispersion of graphene and graphite]{Phonon dispersion of graphene (GE) and graphite (GI), calculated at the experimental and equilibrium lattice spacings ($a_{exp}$ and $a_{th}$). The red straight lines at $\vect{\Gamma}$ and $\vect{K}$ show the slope of the Kohn anomalies. The two lower panels correspond to the dotted windows in the upper panel. The points are theoretical frequencies obtained by direct calculation. Reused with permission from Piscanec \emph{et al.} \autocite{Piscanec2004}}
    \label{fig:graphite-kohn-anomalies}
\end{figure}

### Previous studies of nonequilibrium dynamics in graphite {#sec:graphite-prev-studies}

Ultrafast studies of graphite conceptually start with the work by Kampfrath *et al.* [@Kampfrath2005] The strong coupling between charge carriers and optical phonons $A_1^\prime$ and $E_{2g}$ was demonstrated using time-resolved terahertz spectroscopy by tracking the transient in-plane dielectric function after pumping with short pulses of \SI{1.5}{\electronvolt} photons. The authors reveal that by \SI{500}{\femto\second} after photoexcitation, more than 90% of the absorbed pump energy has been transferred to the lattice. The authors subsequently determine the lifetime of the strongly-coupled optical modes.

The understanding of hot carrier dynamics in allotropes of carbon lept forward significantly following the development of time- and angle-resolved photoemission spectroscopy (trARPES). Early work by Johannsen *et al.*[@Johannsen2013] probed graphene photoexcited with \SI{1.5}{\electronvolt} photons. This experiment showed that for all time-delay after photoexcitation, the energy distribution for electrons near the Dirac cones can be described by Fermi-Dirac statistics. This implies that by \SI{60}{\femto\second} after photoexcitation -- the effective time-resolution of the experiment --, the electrons have thermalized. This estimate has been confirmed by optical pump-probe spectroscopy [@Breusing2009]. Johannsen *et al.*[@Johannsen2013] also determined the process by which energy is transferred to acoustic phonons. According to their modelling, the restricted set of momentum-conserving energy-transfer pathways between the hot electron gas and acoustic phonon can be substantially expanded thanks to defect-mediated interactions, via a process called *supercollision*[@Song2012;@Betz2013].

The experiment by Stange *et al.* [@Stange2015] on graphite showcases some important differences between ultrafast work on graphite and graphene. With a time-resolution of \SI{32}{\femto\second}, it was found that the electronic energy distribution is not well-described by Fermi-Dirac statistics until about \SI{100}{\femto\second}, after which an effective electronic temperature can be defined. The authors also weighed in on the subsequent phonon thermalization mechanism. Contrary to the work by Johannsen *et al.*[@Johannsen2013], Stange *et al.*[@Stange2015] showed that the supercollision model *does not explain* ultrafast thermalization in graphite. Instead, anharmonic coupling between optical and acoustic modes, initially proposed by Yan *et al.*[@Yan2009], was found to be the most important thermalization mechanism. The difference in thermalization process is attributed to the much longer electron mean-free path in graphite (\SIrange{0.2}{10}{\micro\meter})[@Garcia2008] vs. substrate-supported graphene ($<\SI{150}{\nano\meter}$)[@Betz2013].

Ultrafast studies of graphite with structural probes have also been performed. Carbone *et al.*[@Carbone2008] studied the ablation of whole graphene layers from graphite samples photoexcited with an optical density of up to \SI{44.5}{\milli\joule\per\square\centi\meter} via ultrafast electron diffraction. Their results support the existence of strongly-coupled optical phonon modes, as demonstrated by Kampfrath *et al.*[@Kampfrath2005] Ultrafast electron diffraction measurements by Chatelain *et al.*[@Chatelain2014a] also confirm the existence of two strongly-coupled optical modes by observing their coherent generation, indicating that the characteristic energy-transfer time from the hot electron gas is shorter that the phonon period. The dynamics of strongly-coupled optical phonons measured with ultrafast diffraction (electron and X-ray) has been the subject of multiple other studies[@Raman2008;@Schafer2011;@VanMourik2014;@Liang2014;@White2014;@Harb2016;@Stern2018].

### Geometrical interpretation of electron-phonon coupling in photoexcited graphite

```{.matplotlib #fig:graphite-photoexcitation file="figures/graphite/photoexcitation.py" caption="Cut of the electronic dispersion $E(\vect{k})$ along the $\vect{K} - \vect{M} - \vect{K}$ line (see @fig:graphite-electronic-structure) shows the effects of photoexcitation with \SI{1.55}{\electronvolt} photons ($\gamma$). The momentum-conserving decay path are highlighted and explained in the text."}
```

Following the findings presented in the previous section, the effects of electron-phonon coupling in graphite after photoexcitation has been explained intuitively by Kampfrath *et al.*[@Kampfrath2005] The geometric relationship between the electronic structure and the strongly-coupled modes is shown in @fig:graphite-photoexcitation. Specifically, it is the Dirac cones that tell us the consequences of photoexcitation with \SI{800}{\nano\meter} light. Photons can only drive vertical (zero-momentum) transitions near the Dirac cones. As the electron cloud thermalizes, two classes of momentum-conserving decay pathways involving phonons emerge. One such pathway allows for an electron to move across the Dirac cone, emitting an $E_{2g}$ phonon with small wavevector $\vect{k} \sim \vect{0}$. Another pathway allows for an electron to hop onto a neighboring Dirac cone, emitting an $A_1^\prime$ phonon with large wavevector $\vect{k} \sim \vect{K}$.

## Experimental and computational methods

### Sample preparation

Single-crystal flakes of natural graphite \SIrange{10}{90}{\nano\meter} thick were prepared using a mechanical exfoliation procedure analogous to the work by Novoselov *et al.* [@Novoselov2004], briefly described here. Thick flakes were embedded in Crystalbond glue on a \SI{3}{\milli\meter} copper TEM grid (200 lines per inch). The embedded flakes are then exfoliated using ordinary adhesive tape. The procedure was repeated until the embedded flakes were translucent when observed under an optical microscope. The glue is then delicately washed away with a solvent. The choice of the solvent is dependent on the glue used; in the present case, acetone was used. Sample thickness has been measured directly using atomic force microscopy characterization. Once a suitable sample has been identified, an aperture was made using aluminum foil to isolate a sample region with uniform thickness. The resulting sample used in this work covered 500 × 500 \si{\square\micro\meter}, with a thickness of \SI{88}{\nano\meter}.

### Data acquisition {#sec:graphite-data-acquisition}

The ultrafast electron diffuse scattering experiments presented in this chapter made use of the experimental setup presented in @sec:experimental_setup. Ultrashort laser pulses of light were shone at $t=t_0$ on an \SI{88}{\nano\meter}-thick single-crystal specimen of graphite, oriented in the $\left[001\right]$ direction, at a repetition rate of \SI{1}{\kilo\hertz}. Compressed electron bunches with 10^7^ electrons per bunch were transmitted through the sample at $t=t_0 + \tau$. The time-delay $\tau$ was scanned from \SI{-40}{\pico\second} to \SI{680}{\pico\second}.

The interrogated film were pumped with a pump spot of 1 × 1 \si{\square\milli\meter} full-width at half-maximum, ensuring nearly uniform illumination of the probed volume. The film was pumped at a fluence of \SI{12}{\milli\joule\per\square\centi\meter}, resulting in an absorbed energy density of \SI{8}{\joule\per\cubic\meter}. The scattering patterns are collected with a Gatan Ultrascan 895 camera: a 2.54 × 2.54 \si{\square\cm} scintillator fiber-coupled to a 2048 px × 2048 px charge-coupled detector (CCD) placed \SI{25}{\centi\meter} away from the sample. Per-pixel scattering intensity fluctuations over laboratory time reveals a transient dynamic range of 1 : 10^8^, allowing the acquisition of diffraction patterns and diffuse scattering patterns simultaneously. A static diffraction pattern is shown in @fig:graphite-static a).

Due to the flatness of the Ewald sphere for \SI{90}{\kilo\electronvolt} electrons (@sec:scattering-ewald-sphere), many symmetry-related reflections are visible within each pattern. The information contained in a set of symmetry-equivalent reflections is redundant due to the point-group symmetry of the scattering crystal. As long as the point-group symmetry is not broken by photoexcitation, it is possible to harness the redundancy to enhance the signal-to-noise ratio of a ultrafast electron diffuse scattering dataset. In the case of graphite, no observable symmetry-breaking phenomena is brought on by photoexcitation at \SI{1.55}{\electronvolt} when looking at the raw data. Moreover, trARPES experiments [@Stange2015] do not show the opening of a gap in the electronic band structure, albeit at much lower photoexcitation densities, which would be indicative of point-group symmetry breaking. The point-group of graphite is $6/mmm$, which encompasses six-fold discrete rotational symmetry in the $\vect{a}$ -- $\vect{b}$ plane. Therefore, specifically in the case of graphite oriented in the $\left[ 001 \right]$ direction, the diffuse signals can be safely enhanced by a factor of $\sqrt{6}$ by the use of a six-fold discrete azimuthal average:
$$
    I(\vect{q}, \tau) \to \frac{1}{6} \sum_{n=1}^6 I( \vect{R}(\tfrac{\pi n}{3}) \cdot \vect{q}, \tau)
$$
where $\vect{R}(\theta)$ is the in-plane rotation matrix:
$$
    \vect{R}(\theta) = \begin{pmatrix}
                       \cos{\theta} & -\sin{\theta} & 0\\
                       \sin{\theta} & \cos{\theta}  & 0\\
                       0            & 0             & 1
                       \end{pmatrix}
$$
Throughout the rest of this chapter, "scattering intensity" will imply discrete azimuthal average unless otherwise noted. An example of six-fold averaged diffraction pattern is shown in @fig:graphite-static b).

```{.matplotlib #fig:graphite-static file="figures/graphite/diff-static.py" caption="Static diffraction pattern of graphite. **a)** static, unprocessed diffraction pattern. **b)** Six-fold discrete azimuthal average of the diffraction pattern in a) results in $\sqrt{6}$ increase in signal-to-noise ratio. Brillouin zones are shown around each reflection to guide the eye. "}
```

### Computational details {#sec:computational-details}

This section contains the details of the calculations used throughout this chapter, including what is shown in @fig:graphite-static-dispersion. The aim of the computations was to extract the phonon mode frequencies $\set{ \omega_{\lambda}(\vect{k})}$ and polarization vectors $\set{ \vect{e}_{\lambda, s}(\vect{k})}$ that appear in @eq:scattering-one-phonon-structure-factor.

#### Structure-determination{#sec:graphite-structure-relaxation}

In order to calculate the force between atoms, the structure of graphite was computed via *structure relaxation*, performed using the plane-wave self-consistent field program `PWSCF` from the `QUANTUM ESPRESSO` software suite [@Giannozzi2017]. The structure was relaxed using a 18 × 18 × 10 $\vect{k}$ mesh centered at $\vect{\Gamma}$ -- selected using the Monkhorst-Pack method [@Monkhorst1976] -- and force and energy thresholds of \SI{1e-8}{\rydberg\per\bohr} and \SI{1e-15}{\rydberg} respectively, where \si{\bohr} is the Bohr radius. Based on the force constants determined from the relaxed structure, the dynamical matrices were computed on a 5 × 5 × 3 $\vect{q}$ mesh using a self-consistency threshold of \SI{1e-18}{\rydberg}. The resulting graphite structure is equivalent to @eq:graphite-lattice, with $a=\SI{1.231}{\angstrom}$ and $c=\SI{6.837}{\angstrom}$.

#### Phonon properties

From a structure with zero-temperature atomic positions $\set{\vect{r}_s}$ and instantaneous displacements $\set{\vect{u}_s}$, the potential energy $U$ due to ions repulsion is given by:
$$
    U = \frac{1}{2}\sum_{s, s^\prime} \vect{u}_s \cdot \left[ \vect{D}(\vect{r}_s - \vect{r}_{s^\prime}) ~ \vect{u}_{s^\prime} \right]
$$
where the *dynamical matrix* $\vect{D}$ encodes the change in potential energy associated with a small change in distance between ions[@Ashcroft1976DynMatrix]. Given that the displacement vectors $\set{\vect{u}_s}$ can be expressed as a sum of lattice waves (@eq:scattering-displacement), the eigenvalues and eigenvectors of the dynamical matrix are related to vibrational frequencies and polarizations respectively[^dynmat].

The phonon frequencies $\set{ \omega_{\lambda}(\vect{k}) }$ and polarization vectors $\set{ \vect{e}_{\lambda, s}(\vect{k})}$ associated with the structure calculated in @sec:graphite-structure-relaxation were computed using the `PHONON` program, again within the `QUANTUM ESPRESSO` software suite. This calculation made use of the B86b exchange-coupled Perdew-Burke-Ernzerhof generalized-gradient approximation [@Becke1986;@Perdew1996] and the projector augmented-wave method [@Blochl1994]. The cutoff energy of the wave function was set to \SI{100}{\rydberg}, while the cutoff energy for the charge density was set to \SI{1.2e3}{\rydberg}. A Fermi-Dirac smearing of \SI{0.06}{\rydberg} was also applied. To include the dispersion of energy along the stacking axis $\vect{a}_3$, the exchange-hole dipole moment method was used [@Becke2007]. 

[^dynmat]: This is an intuitive picture of the dynamical matrix. In practice, diagonalization is usually performed in reciprocal space. For an example, see Al-Jishi and Dresselhaus[@AlJishi1982].

#### Clustering of phonon properties into physically-relevant branches

As with most diagonalization procedures, the eigenvalues and eigenvectors of the dynamical matrix are returned in an order that is not physically-relevant. Therefore, in order to determine what are the phonon properties of a phonon mode $\lambda$, the phonon properties calculated by `PHONON` need to be assigned a mode, or *clustered*.

As of writing this, there is no component within the `QUANTUM ESPRESSO` software suite that can do this. The scheme described in this section was used instead. The general idea behind the procedure is that phonon properties are continuous. Let $P_{\lambda, \vect{k}_i}$ be the abstract vector representing phonon properties of mode $\lambda$ at one of the irreducible $\vect{k}$ points $\set{\vect{k}_i}$:
$$
    P_{\lambda, \vect{k}_i} \equiv \begin{bmatrix} \omega_{\lambda, \vect{k}_i} & \vect{e}_{s=1,\lambda,\vect{k}_i} & ... & \vect{e}_{s=M,\lambda,\vect{k}_i} \end{bmatrix}^T
$$
where the index $s$ runs over all $M$ atoms of the unit cell. Define the metric between two such abstract vectors $\vect{P}_{\lambda^{\prime}, \vect{k}}$ and $\vect{P}_{\lambda, \vect{k}'}$ as:
$$
    \left\Vert \vect{P}_{\lambda^\prime, \vect{k}} - \vect{P}_{\lambda, \vect{k}'} \right\Vert = |\omega_{\lambda^{\prime}, \vect{k}} - \omega_{\lambda, \vect{k}'}|^2 + \sum_s \left\Vert \vect{e}_{\lambda^{\prime},s,\vect{k}} - \vect{e}_{\lambda,s,\vect{k}'}\right\Vert
$$
A one-dimensional path $\gamma$ connecting all irreducible points $\set{\vect{k}_i}$ was defined, starting at $\vect{\Gamma}$. At $\vect{\Gamma}$, polarization vectors are associated with a mode based on geometry. For example, a mode with negligible frequency and polarization vectors at that all point in the same direction physically corresponds to a longitudinal acoustic mode. The manual assignment for high-frequency optical modes is a bit more arbitrary. Then, following the path $\gamma$, the assignment of phonon branches $\lambda^\prime$ at $\gamma(\vect{k} + \vect{\Delta})$ minimizes the quantity $\left\Vert \vect{P}_{\lambda^{\prime}, \vect{k}} - \vect{P}_{\lambda, \vect{k} + \vect{\Delta}} \right\Vert$. 

The procedure described above has been adapted for numerical evaluation and is now part of the `scikit-ued` software package [@RenedeCotret2018], which is described in more details in @sec:appendix-software. For the rest of this chapter, it will be assumed that the phonon properties are clustered such that the usual physical phonon branch labels (e.g. TA, LO) are meaningful.

## Diffuse intensity dynamics

The change in scattering intensity $\Delta I(\vect{q}, t=\tau) \equiv I(\vect{q}, \tau) - I(\vect{q}, -\infty)$ for representative time-delays is shown in @fig:graphite-ueds. 

```{#fig:graphite-ueds .matplotlib file="figures/graphite/ueds.py" caption="Change in scattering intensity $\Delta I(\vect{q}, t=\tau) \equiv I(\vect{q}, \tau) - I(\vect{q}, \tau < 0)$ of photoexcited graphite for a few representative time-delays $\tau$. Hexagonal Brillouin zones are shown on half of the reflections to guide the eye. Scattering patterns show diffuse scattering in the range of $|\vect{q}| < \SI{12}{\per\angstrom}$. Negative going features (blue) are exclusively due to the transient Debye-Waller effect on the Bragg peaks. All positive changes (red) are diffuse scattering intensity."}
```

First and foremost, note that negative-going features are only visible in the vicinity of Bragg peaks. This is due exclusively to the transient Debye-Waller effect [@Ligges2009]. As photodeposited energy transfers from the electrons to the lattice, average real-space disorder due to phonons lowers the overall periodicity of the lattice, which in turn results in smaller Bragg peaks (@sec:scattering-lattice-waves). This is analogous to the (static) thermal Debye-Waller effect [@Debye1913;@Waller1923]. An example curve for the $(200)$ reflection is shown in @fig:graphite-dw-example. The dynamics associated with the Bragg peaks in photoexcited graphite are discussed in detail elsewhere[@Chatelain2014;@Chatelain2014a]. This leaves only positive-going features everywhere away from Bragg peaks, which are diffuse in nature.

```{.matplotlib #fig:graphite-dw-example file="figures/graphite/dw-example.py" caption="Relative intensity change of the $(200)$ reflection after photoexcitation examplifies the transient Debye-Waller effect."}
```

Second, the structure visible across the Brillouin zone is different for reflections that ought to represent the same physical lattice waves. For example, consider the $(010)$ and $(020)$ reflections, \SI{100}{\pico\second} after photoexcitation, shown in @fig:graphite-ueds-zoomed. Around $(010)$, diffuse intensity increases near two of the $\vect{M}$ points, while diffuse intensity is increasing at two other $\vect{M}$ points near $(020)$. 

```{#fig:graphite-ueds-zoomed .matplotlib file="figures/graphite/ueds-zoomed.py" caption="Comparison of the diffuse intensity change after \SI{100}{\pico\second} for two Brillouin zones, $(010)$ and $(020)$. The colormap scaling is identical to @fig:graphite-ueds. The difference between those two images can be explained by the difference in one-phonon structure factors $|F_{1\lambda}(\vect{q}, \tau)^2|$."}
```

The drastic difference between physically-equivalent reflections can be explained by the difference in the one-phonon structure factors $|F_{1\lambda}(\vect{q}, \tau)|^2$. Therefore, in order to understand the scattering patterns presented in @fig:graphite-ueds, the one-phonon structure factors need to be computed. This is done in the next section.

## The one-phonon structure factors

Recall the definition of diffuse intensity from @eq:scattering-one-phonon-structure-factor:
$$
    I_1(\vect{q}, \tau) = N_c I_e \sum_{\lambda} \frac{n_{\lambda}(\vect{k}, \tau) + \frac{1}{2}}{\omega_{\lambda}(\vect{k}, \tau)} \left| F_{1\lambda}(\vect{q}, \tau)\right|^2
$$
and
$$
    F_{1\lambda}(\vect{q}, \tau) = \sum_s e^{-W_s} \frac{f_s(\vect{q})}{\sqrt{\mu_s}} \left( \vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})\right) 
$$
where $N_c$ is the number of scattering unit cells, $I_e$ is the intensity of a single scattering event, $s$ are indices that label atoms in the unit cell, $W_s$ is the Debye-Waller factor associated with atom $s$, $f_s(\vect{q})$ are the atomic form factors, $j \in \set{ 1, 2, ..., 12}$ runs over the phonon modes, $\set{ n_{\lambda}(\vect{k}, \tau)}$ are the phonon populations, $\set{ \vect{e}_{\lambda,s}(\vect{k}) }$ are the phonon polarization vectors, and $\set{ \omega_{\lambda}(\vect{k}, \tau)}$ are the phonon frequencies.

As shown in @fig:graphite-ueds-zoomed, observed diffuse intensity is the combination of quantities that can be separated in two categories. Some quantities, like phonon frequencies, are *local* in reciprocal space (i.e. defined on $\vect{k}$). On the other hand, the one-phonon structure factors are *nonlocal* in reciprocal space (defined on $\vect{q}$). This results in measurable differences when comparing diffuse intensity at various reflections, which will be used to overcome the energy-integrative nature of ultrafast electron diffuse scattering experiments.

In order to recover mode-dependent phonon populations, one-phonon structure factors need to be calculated for every relevant phonon mode. Most importantly, the phonon polarization vectors $\set{ \vect{e}_{\lambda, s}(\vect{k})}$ need to be determined. These vectors describe the direction of movement for each atom due to a particular lattice wave. Fortunately, polarization vectors are a byproduct of the calculation of the phonon dispersion relation, shown in @fig:graphite-static-dispersion.

### Transient Debye-Waller factors {#sec:graphite-debye-waller}

One more quantity needs to be discussed before one-phonon structure factors can be computed: the Debye-Waller factors $W_s(\vect{q}, \tau)$. The Debye-Waller factors describe the reduction of scattering intensity at $\vect{q}$ due to the effective deformation of the scattering potential of atom $s$ by *all* phonon branches. The general expression for the *anisotropic* Debye-Waller factor for atom $s$ is given by[@Xu2005]:
$$
    W_s(\vect{q}, \tau) = \frac{1}{4 \mu_s} \sum_{s, \lambda} |a_{\lambda}(\vect{k}, \tau)|^2 |\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})|^2
    \label{eq:graphite-debye-waller}
$$
where $a_{\lambda}$ is vibrational amplitude for the phonon mode $\lambda$:
$$
    |a_{\lambda}(\vect{k}, \tau)|^2 = \frac{2 \hbar}{N_c} \left( \frac{n_{\lambda}(\vect{k}, \tau) + \frac{1}{2}}{ \omega_{\lambda}(\vect{k}, \tau)}\right)
$$
Based on @eq:graphite-debye-waller, the Debye-Waller factor is not sensitive measure to the wavevector-dependent nonequilibrium dynamics of lattice waves because it involves the contribution from *all* lattice waves, integrated over $\vect{k}$. Therefore, phonon population dynamics can only affect the magnitude of Debye-Waller factors. 

The potential time-dependence of the Debye-Waller factors was investigated, via the time-dependence of mode populations. Profoundly non-equilibrium distribution of phonon mode populations were simulated, with all modes populated according to room-temperature average energy except for one mode at high temperature. To account for differences in heat capacities, the maximum temperature for optical modes was \SI{5000}{\kelvin}, while the maximum temperature for acoustic modes was capped at \SI{1500}{\kelvin}. These extreme non-equilibrium scenarios increase the value of $\sum_s W_s(\vect{q})$ by at most 1.5\% for optical modes, and 8\% for acoustic modes. Since those fractional changes are constant across $\vect{q}$, wavevector-dependent changes in ultrafast electron diffuse scattering signals are not significantly impacted by changes in the Debye-Waller factors. Any time-dependence of the one-phonon structure factors can be ignored as a first approximation. 

<!-- should I show this plot below? I don't like it. -->
<!-- The most significant example is shown in @fig:graphite-debye-waller-la, where the LA mode is populated according to a temperature of \SI{1500}{\kelvin}. Since those fractional changes are constant across $\vect{q}$, wavevector-dependent changes in ultrafast electron diffuse scattering signals are not significantly impacted by changes in the Debye-Waller factors $W_s(\vect{q}, \tau)$, and so any time-dependence of the one-phonon structure factors can be ignored as a first approximation. 

```{#fig:graphite-debye-waller-la .matplotlib file="figures/graphite/debye-waller.py" caption="Effect of profoundly non-equilibrium distribution of phonon mode populations on the Debye-Waller factors. **a)** Sum of Debye-Waller factors according to the expected phonon distribution at room temperature. **b)** Difference between a) and the sum of Debye-Waller factors where all phonon modes are distributed according to room-temperature, except where the LA mode is populated according to a temperature of \SI{1500}{\kelvin}."}
``` -->

### One-phonon structure factor calculations

The calculation for the one-phonon structure factors $|F_{1\lambda}(\vect{q}, \tau)|^2$ was carried out from @eq:scattering-one-phonon-structure-factor. As discussed in the previous section, the dynamics in one-phonon structure factors due to photoexcitation via the change in Debye-Waller factors can be neglected in this work. Therefore, the calculations of this section assume room-temperature thermal distribution of lattice waves; this is implied by time-delay $\tau<0$ (pre-photoexcitation). The one-phonon structure factors for longitudinal modes LA, LO1 - LO3 are shown @fig:graphite-oneph-longitudinal. The one-phonon structure factors for transverse modes TA, TO1 - TO3 are shown @fig:graphite-oneph-transverse. Out-of-plane modes yield trivial one-phonon structure factors for all scattering vectors $\vect{q}$ probed by the experiments presented herein because $\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}) = 0$.

```{.matplotlib #fig:graphite-oneph-longitudinal file="figures/graphite/oneph-long.py" caption="Calculated one-phonon structure factors $|F_{1\lambda}(\vect{q}, \tau<0)|^2$ of in-plane longitudinal modes at \SI{300}{\kelvin}, for scattering vectors $\vect{q}$ equivalents to the detector area shown in @fig:graphite-ueds. See @fig:graphite-oneph-transverse for the transverse modes equivalent."}
```

```{.matplotlib #fig:graphite-oneph-transverse file="figures/graphite/oneph-trans.py" caption="Calculated one-phonon structure factors $|F_{1\lambda}(\vect{q}, \tau<0)|^2$ of in-plane transverse modes at \SI{300}{\kelvin}, for scattering vectors $\vect{q}$ equivalents to the detector area shown in @fig:graphite-ueds. See @fig:graphite-oneph-longitudinal for the longitudinal modes equivalent."}
```

The one-phonon structure factors display complex structure in reciprocal space. The structure for a particular mode can be thought of as continuous selection rule: regions in $\vect{q}$ space where the one-phonon structure factor for mode $\lambda$ is large are regions where phonon mode $\lambda$ contributes importantly atomic vibrations in reciprocal space. The structure of $|F_{1\lambda}|^2$ is most determined by the values of terms like $\set{ \vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}) }$. For acoustic modes LA and TA, the structure of the one-phonon structure factor is easier to understand near $\vect{\Gamma}$. For example, $|F_{1\lambda}|^2$ is highest for the LA mode in the radial $\vect{q}$ direction because the polarization vectors of those modes are parallel to $\vect{q}$; this fact defines longitudinal waves. For transverse waves, where atomic motion is perpendicular to the lattice wave propagation direction, $|F_{1\lambda}|^2$ is higher in the azimuthal direction.

### Weighted phonon dispersion

An alternative visualization of one-phonon structure factors are weighted dispersion curves, as shown in @fig:graphite-weighted-dispersion. This way of looking at $|F_{1\lambda}|^2$ clearly highlights the difference across Brillouin zones. For example, the difference $|F_{1\lambda=\text{TA}}|^2 - |F_{1\lambda=\text{LA}}|^2$ is either negative for $\vect{q} \approx \vect{\Gamma}_{(010)}$ or positive for $\vect{q} \approx \vect{\Gamma}_{(\bar{1}10)}$, even though the two reflections are adjacent. The variation of $|F_{1\lambda}|^2$ explains the diffuse intensity difference shown in @fig:graphite-ueds-zoomed.

```{.matplotlib #fig:graphite-weighted-dispersion file="figures/graphite/weighted-dispersion.py" caption="Calculated one-phonon structure factors visualized as weighted dispersion curves for selected in-plane modes. The color saturation of dispersion curves is proportional to $|F_{1\lambda}|^2$ of the associated mode. Equivalent paths in the Brillouin zone around two reflections are shown to highlight the high degree of reciprocal space structure: $(010)$ on the left and $(\bar{1}10)$ on the right. The geometry of the paths with respect to $(000)$ are shown in the inset on the lower left."}
```

By examining @fig:graphite-weighted-dispersion, it appears that at certain locations, the one-phonon structure factor for one phonon mode dominates; for example, $|F_{1\lambda=\text{LA}}|^2$ near $\vect{q} \approx \vect{\Gamma}_{(010)}$. It might be tempting to attribute the diffuse intensity dynamics at that location exclusively to one mode. If enough such locations in reciprocal space existed, it would be possible to bypass the energy-insensitivity of ultrafast electron diffuse scattering and extract mode-dependent phonon dynamics at special points in the Brillouin zone; this idea forms the basis of previous work by the author in Stern *et al.*[@Stern2018] The next section describes how crude this intuition really is.

### Relative mode contributions

Are there locations in reciprocal space where the contribution of a particular lattice wave dominates the diffuse intensity? In order to compare the contribution of a particular phonon mode $\lambda$ on ultrafast electron diffuse scattering data, the following weight can be defined:
$$
    W_{\lambda}(\vect{q}, \tau) \equiv \frac{\left| F_{1\lambda}(\vect{q}, \tau)\right|^2}{\omega_{\lambda}(\vect{k}, \tau)}
    \label{eq:graphite-ueds-weight}
$$
With the weight definition of @eq:graphite-ueds-weight, the change in diffuse intensity (@eq:scattering-diffuse-intensity) can be re-written as:
\begin{align}
    \Delta I_1(\vect{q}, \tau) &\equiv I_1(\vect{q}, \tau) - I_1(\vect{q}, \tau < 0) \\
                              &= \left( n_{\lambda}(\vect{q}, \tau) + \sfrac{1}{2} \right) W_{\lambda}(\vect{q}, \tau) - \left( n_{\lambda}(\vect{q}, \tau < 0) + \sfrac{1}{2} \right) W_{\lambda}(\vect{q}, \tau < 0) \nonumber \\
                              &\approx \left[ n_{\lambda}(\vect{q}, \tau) - n_{\lambda}(\vect{q}, \tau < 0)\right] W_{\lambda}(\vect{q}, \tau < 0) \nonumber
\end{align}
In the last line, the quasi-static nature of the one-phonon structure factors $|F_{1\lambda}(\vect{q}, \tau)|^2 \approx |F_{1\lambda}(\vect{q}, \tau < 0)|^2$ was used, as discussed in @sec:graphite-debye-waller. The approximation $\omega_{\lambda}(\vect{k}, \tau) \approx \omega_{\lambda}(\vect{k}, \tau<0)$, valid for $\vect{k}$ away from $\vect{\Gamma}$, was also used; see @sec:graphite-phonon-landscape for a discussion.

@fig:graphite-oneph-majority shows where, across the reciprocal space area equivalent to @fig:graphite-ueds, a single new phonon from a particular branch contributes more than 50\% or more than 75\% of the diffuse intensity increase. In other words, @fig:graphite-oneph-majority shows where *any* branch $\lambda$, with associated weight $W_{\lambda}(\vect{q}, \tau<0)$, contributes more than 50\% or 75\% of the total weight $\sum_{\lambda} W_{\lambda}(\vect{q}, \tau<0)$. The calculation presented in @fig:graphite-oneph-majority is categorical: dynamics in the transient diffuse intensity **cannot** be associated with the dynamics of a single phonon branch. Furthermore, @fig:graphite-oneph-majority shows that at in-plane high symmetry points $\vect{K}$ and $\vect{M}$, no phonon branch contributes to more than 50\% of the associated dynamics in the diffuse intensity.

In order to access the quantitative ultrafast phonon dynamics in a material with ultrafast electron diffuse scattering, a more comprehensive analysis of ultrafast electron diffuse scattering must be devised -- making good use of one-phonon structure factors.

```{.matplotlib #fig:graphite-oneph-majority file="figures/graphite/oneph-majority.py" caption="Locations in reciprocal space where a single new phonon from a particular branch contributes to more than **a)** 50\% and **b)** 75\% of the associated increase in diffuse intensity. $\varnothing$ is used to denote that no phonon mode contributes higher than the threshold."}
```
 
## Phonon spectroscopy across the Brillouin zone{#sec:graphite-ph-spectroscopy}

In this section, the calculation of the one-phonon structure factors will be used to endow ultrafast electron diffuse scattering with energy resolution. An important advantage of ultrafast electron diffuse scattering is that a large number of reflections are simultaneously visible due to the relative flatness of the Ewald sphere (@sec:scattering-ewald-sphere). Compared to ultrafast diffuse x-ray scattering[@Wall2018;@Teitelbaum2018], which is limited to the study of at most a few reflections, the shortcoming of ultrafast electron scattering -- lack of energy-resolution -- can be overcome by the redundancy present in the measurements. 

Recall that the ultrafast change in diffuse intensity can be expressed as follows:
$$
    \frac{\Delta I(\vect{q}, \tau)}{N_c I_e} = \sum_{\lambda} \frac{\Delta n_{\lambda}(\vect{k}, \tau)}{\omega_{\lambda}(\vect{k}, \tau<0)} \left| F_{1\lambda}(\vect{q}, \tau<0)\right|^2
    \label{eq:graphite-ueds-change}
$$
@eq:graphite-ueds-change is approximately valid for $\vect{k}$ away from $\vect{\Gamma}$ (see @sec:graphite-phonon-landscape). Effectively, this restriction is trivial as diffraction peaks drown out the ultrafast electron diffuse scattering signals near $\vect{\Gamma}$. 

The sum over all phonon branches $\lambda$ encodes the lack of energy resolution in ultrafast electron diffuse scattering experiments. However, with enough data redundancy, the contribution of every mode $\lambda$ to @eq:graphite-ueds-change can be extracted. Let $\set{ \vect{H}_1, ..., \vect{H}_M }$ be in-plane reflections that are visible in the experiments. Then, the transient phonon population of mode $\lambda$, $\Delta n_{\lambda}(\vect{k}, \tau)$, solves the following linear system:
$$
    \vect{I}_{\vect{k}}(\tau) = \vect{F}_{\vect{k}} \vect{D}_{\vect{k}}(\tau)
    \label{eq:graphite-vectorized-ueds}
$$
with
\begin{align}
 \vect{I}_{\vect{k}}(\tau) &= 
 	\frac{1}{N_c I_e}
	\begin{bmatrix}
    	\Delta I(\vect{k} + \vect{H}_1, t) & \dots & \Delta I(\vect{k} + \vect{H}_M, t)
	\end{bmatrix}^T \\
\vect{F}_{\vect{k}} &= 
	\begin{bmatrix}
		|F_{11}(\vect{k} + \vect{H}_1, \tau<0)|^2 & \dots  & |F_{1N}(\vect{k} + \vect{H}_1,\tau<0)|^2 \\
		\vdots		   			 			 & \ddots & \vdots					            \\
		|F_{11}(\vect{k} + \vect{H}_M, \tau<0)|^2 & \dots  & |F_{1N}(\vect{k} + \vect{H}_M,\tau<0)|^2
	\end{bmatrix} \\
\vect{D}_{\vect{k}}(\tau) &=
	\begin{bmatrix}
    	\Delta n_1(\vect{k}, \tau)/\omega_1(\vect{k}, \tau<0) & \dots & \Delta n_N(\vect{k}, \tau)/\omega_N(\vect{k}, \tau<0)
	\end{bmatrix}^T
\end{align}
where $N$ runs over the eight in-plane modes. This linear system of equations can be solved numerically provided enough experimental data, i.e. diffuse intensity for at least $M \geq N$ distinct Brillouin zones. Again, this constraint highlights the advantage that ultrafast electron diffuse scattering has over other techniques: the diffuse intensity measured with ultrafast electron diffuse scattering spans a large portion of reciprocal space. 

The choice to solve for $\vect{D}_{\vect{k}}(\tau)$ rather than $\Delta n_{\lambda}(\vect{k}, \tau)$ directly, comes down to the degree of confidence that should be placed in the approximations that were made to get to @eq:graphite-ueds-change. Phonon polarization vectors are mostly determined by the symmetries of the crystal, which are fixed throughout the experiments. On the other hand, phonon vibrational frequencies in general might be influenced by non-equilibrium carrier and phonon distributions as would be the case with strongly anharmonic crystals (see @sec:snse). Solving for the ratio of phonon population to vibrational frequency is more robust against the weaknesses of this modeling because the one-phonon structure factors only take into account the polarization vectors.

### Numerical procedure

The procedure used to numerically solve for $\vect{D}_{\vect{k}}(\tau)$ is described henceforth.

All visible Brillouin zones in the measurement were used to generate $\vect{I}_{\vect{k}}(\tau)$. Since the scattering patterns have been symmetrized (@fig:graphite-static), one would expect that using the intensity data for reflections related by symmetry would be redundant; however, the following procedure worked better when using the entire area of the detector. This is no doubt due to minute misalignment of the scattering patterns and uncertainty in detector position, which are averaged out when using all available data. The Brillouin zones associated with all in-plane reflections $\set{ \vect{H} | |\vect{H}| \leq \SI{12}{\per\angstrom} }$ were used, for a total of 44 Brillouin zones. Therefore, at every $\vect{k}$ point, $\vect{I}_{\vect{k}}(\tau)$ and $\vect{F}_{\vect{k}}$ have matrix dimensions of $(44, 1)$ and $(44, 8)$, respectively.

The possible values for elements of $\vect{D}_{\vect{k}}(\tau)$ were constrained to be non-negative, implying that $\Delta n_{\lambda}(\vect{k}, \tau) \geq 0 ~ \forall \tau$. In other words, it is assumed that the population of a phonon branch cannot drop below its equilibrium level. While not strictly necessary, it was found to lead to a more stable solution. This constraint allowed for the use of the non-negative approximate matrix inversion algorithm [@Lawson1995] to solve @eq:graphite-vectorized-ueds. This procedure was repeated for every reduced wavevector $\vect{k}$ and time-delay $\vect{\tau}$. Stable solutions were found for $|\vect{k}| \geq \SI{0.45}{\per\angstrom}$, away from $\vect{\Gamma}$ and elastic scattering signals.

It must be emphasized that the numerical solution to $\vect{D}_{\vect{k}}(\tau)$ is not the result of fitting (iterative least-squares method), but rather a non-iterative approximate matrix inversion based on linear least-squares. The method described in this section admits no free parameter, other than the phonon vibrational frequencies and polarization vectors.

### Population dynamics{#sec:graphite-pop-dynamics}

The numerical solution for @eq:graphite-vectorized-ueds is shown in @fig:graphite-ph-populations. The transient phonon populations $\set{ \Delta n_{\lambda}(\vect{k}, \tau)}$ across the Brillouin zone are shown for three important modes: TO2, TA, and LA. 

```{.matplotlib #fig:graphite-ph-populations file="figures/graphite/decomp.py" caption="Measurement of the change in transient phonon population $\Delta n_{\lambda}(\vect{k}, \tau)$ following photoexcitation for relevant in-plane modes of graphite across the Brillouin zone. The solution domain is bound by white circle at $|\vect{k}| \leq \SI{0.45}{\per\angstrom}$, and by a solid white hexagon at the Brillouin zone edge. The Brillouin zone midpoint is highlighted with a dashed white hexagon. The location of the $A_1^\prime$ mode is shown in the top row."}
```

While previous work by the author discussed the nonthermal phonon dynamics qualitatively[@Stern2018], @fig:graphite-ph-populations allows to quantitatively determine how the energy deposited in the electronic system flows and thermalizes. A discussion of the observed physical processes is presented below.

As discussed in @sec:graphite-prev-studies, two optical phonons modes are strongly-coupled to the electronic system: $A_1^\prime$ located near the $\vect{K}$ point, and $E_{2g}$, near $\vect{\Gamma}$. $E_{2g}$ is obscured by the elastic signals near $\vect{\Gamma}$; however, its dynamics are accessible via ultrafast spectroscopy [@Kampfrath2005]. On the other hand, $A_1^\prime$ is clearly visible in @fig:graphite-ph-populations. It behaves as expected: a fast early increase in population is seen due to the transfer of energy from the electronic subsystem. By \SI{5}{\pico\second}, energy has already been transferred away through phonon-phonon coupling. 

After a few picoseconds, the strongly-coupled optical modes decay into lower-energy lattice waves. Since one of the primary causes of phonon-phonon scattering is anharmonicity, the decay probabilities are usually computed via density functional theory calculations [@Bonini2007]. Ultrafast electron diffuse scattering shows an experimental determination of those decay pathways. The transfer of energy from high-energy optical phonons to lower-energy lattice waves must satisfy the conservation of momentum and energy. This restricts the number of possible decay pathways to mid-Brillouin zone phonons, either directly (from $E_{2g}$) or via Umklapp scattering[@Peierls1929] (from $A_1^\prime$). Specifically, from \SIrange{1.5}{100}{\pico\second}, acoustic phonon population increase is observed along the $\vect{\Gamma}-\vect{M}$ line. This measurement is in accordance with predicted anharmonic decay probabilities from the $E_{2g}$ mode [@Bonini2007]. The small increase of population of the LA mode around \SI{500}{\femto\second} at $\tfrac{1}{2}\vect{K}$ is also a confirmation of predictions by Bonini *et al.*[@Bonini2007]

Over longer time-scales (\SI{25}{\pico\second}), the nonthermal TA population has significantly pooled at $\tfrac{1}{3}\vect{M}$ and $\vect{M}$. There are no three-phonon anharmonic decay processes that start in a transverse mode like TA. Allowed *in-plane* interband transitions are of the types $\text{L} \to \text{T} + \text{T}$ or $\text{L} \to \text{L} + \text{T}$, where L and T represent longitudinal and transverse modes respectively [@Lax1981;@Khitun2001]. It is therefore expected *a priori* that a nonthermal TA population remains until decay to out-of-plane phonon modes occur. Predicted decay probabilities also support this observation; both LA and TA modes will favor decay to out-of-plane modes[@Paulatto2013]  which are not visible to ultrafast electron diffuse scattering measurements. A special case is visible at $\tfrac{1}{2}\vect{M}$ and $\tfrac{1}{2}\vect{K}$, where predicted anharmonic lifetimes of TA and LA drop significantly due to the activation of Umklapp scattering to out-of-plane acoustic phonons[@Paulatto2013]. The measurements of @fig:graphite-ph-populations corroborate the predictions: energy flow away from TA and LA at the mid-Brillouin zone is fast enough that population never builds up. 

The energy flow that follows photoexcitation described in this section is exactly what one would expect from the *hot phonon model*[@Yan2009;@Scheuch2011]. The measured phonon population dynamics are perfectly in line with the simple three-temperature modelling from Stange *et al.* [@Stange2015] There is no indication in @fig:graphite-ph-populations that defects allow the scattering of electrons and acoustic phonons, contrary to experiments on substrate-supported graphene described in @sec:graphite-prev-studies.

### Long-term decay

```{.matplotlib #fig:graphite-long-term file="figures/graphite/energy.py" caption="Energy stored in-plane after photoexcitation. **Inset** long-term trend shows that thermalization has not yet occurred by \SI{600}{\pico\second}."}
```

It is possible to quantify how long does the non-equilibrium distribution of lattice waves last in graphite. The energy transferred to the sample from photoexcitation, stored in plane in a mode $\lambda$, can be expressed as follows:
$$
    \Delta E_{\lambda}(\tau) = \int_{\text{BZ}} d\vect{k} ~ \Delta n_{\lambda}(\vect{k}, \tau) \hbar \omega_{\lambda}(\vect{k}, \tau)
$$
where $\int_{\text{BZ}} d\vect{k}$ is understood to be the integral over the in-plane section of the Brillouin zone. Using the population measurements $\set{ n_{\lambda}(\vect{k}, \tau)}$ presented in @fig:graphite-ph-populations, the mode-dependent relative change in stored energy (in-plane) can be calculated. The results are shown in @fig:graphite-long-term.

The energy trends for modes TA and TO2 are commensurate with the description of @sec:graphite-pop-dynamics. The fascinating aspect of this analysis is revealed by looking at the total energy stored in-plane. Even by \SI{600}{\pico\second}, thermalization of the energy dumped into the sample by photoexcitation has not occurred at all. 

## Mode-projected excitation couplings

In this section, the phonon population dynamics will be used to extract the electron-phonon coupling matrix element to the $A_1^\prime$ mode, and compare the value to other measurements and calculations.

The nonequilibrium flow of energy between excitations, both electronic and lattice in nature, has historically been crudely modelled using the *two-temperature model* [@Allen1987]. In summary, this model states that while the concept of temperature does not apply to nonequilibrium situations, the energy distribution of the electrons and lattice waves may be treated separately. In other words, the separate thermalization of the electronic and lattice subsystems is much faster than the energy transfer between them. It is evident that such a description does not adequately model experimental results from @sec:graphite-ph-spectroscopy for two main reasons. First, the decay of optical phonons to acoustic ones overlaps significantly in time with the transfer of energy from the photoexcited electrons to strongly-coupled optical modes. Second, the energy distribution of lattice waves is far from thermal, even at \SI{100}{\pico\second}. 

Some authors have extended the two-temperature model based on the details of the system being studied[@Stange2015]. Ultrafast electron diffuse scattering measurements allow to move beyond the two-temperature approximation as the energy distribution of phonons is known across the Brillouin zone. And yet, theory has not caught up with this opportunity. The best way to extract electron-phonon and phonon-phonon couplings is still based on rate-equations. In this section, the extraction of electron-phonon and phonon-phonon couplings from ultrafast electron diffuse scattering measurements will be performed at the $\vect{K}$ point, using an extension of the two-temperature model called the *non-thermal lattice model*[@Waldecker2016]. 

The reader should keep in mind that the procedure presented in this section is a retrofit of the results presented in @sec:graphite-ph-spectroscopy into a more conventional framework. This is the best that can be done until theorist seize the opportunity provided by ultrafast electron diffuse scattering measurements.

### The non-thermal lattice model {#sec:graphite-nlm}

The non-thermal lattice model is an extension of the two-temperature model, with the added flexibility that lattice waves need not be thermalized[@Waldecker2016]. It assumes that the energy distribution of each phonon mode admits a thermal description; that is, the energy transfer between phonons of different branches is slower than the scattering of phonons from the same branch. Within this framework, each phonon branch $\lambda$ can be assigned its own molar heat capacity, $C_{ph, \lambda}$, and temperature, $T_{ph,\lambda}$:
$$
\left\{
    \begin{array}{rcl}
        C_e(T_e) \frac{\partial T_e}{\partial \tau} & = & \sum_{\lambda} G_{ep, \lambda}\left[ T_e(\tau) - T_{ph,\lambda}(\tau) \right] + f(\tau) \\
                                                    & ~ & \\
        C_{ph,\lambda}(T_{ph,\lambda}) \frac{\partial T_{ph,\lambda}}{\partial \tau} & = &\sum_{\lambda^\prime\neq \lambda} G_{ep, \lambda^\prime} \left[ T_e(\tau) - T_{ph,\lambda^\prime}(\tau) \right] \\
				                                                   & + & G_{pp,\lambda \lambda^\prime} \left[ T_{ph,\lambda}(\tau) - T_{ph,\lambda^\prime}(\tau) \right]
    \end{array}
\right\}_{\lambda=1}^{N}
\label{eq:graphite-nlm}
$$
where $f(\tau)$ is the laser pulse profile, and $C_e$ and $T_e$ are the electronic heat capacity and electron temperature, respectively. As discussed in @sec:graphite-prev-studies, the use of an electronic temperature is acceptable for $\tau > \SI{100}{\femto\second}$, in the case of photoexcited graphite specifically. The constants $G_{ep, \lambda}$ describe the rate of energy flow between the electrons and phonon mode $\lambda$, while constants $G_{pp, \lambda \lambda^\prime}$ encode the rate of energy flow between phonon modes $\lambda$ and $\lambda^\prime$. 

A small clarification is needed. A coupling constant $G$ describes the rate of energy transfer from one excitation type to another, if the temperature of one unit of volume was increase by \SI{1}{\kelvin}. This is *not* directly equivalent to the electron-phonon coupling tensor $g$ introduced in @sec:introduction-epc, but it is a helpful intermediate experimental quantity. The relationship between $G$ and $g$ is developed further below in @sec:graphite-coupling-constants.

Observations of transient phonon populations are more general than mode temperatures. However, in order to make use of the non-thermal lattice model, the mode temperature can be related to transient phonon mode populations via the Bose-Einstein distribution[@Bose1924]:
$$
    n_{\lambda}(\vect{k}, \tau) \propto \left[ \exp{\left( \frac{\hbar \omega_{\lambda}(\vect{k}, \tau<0)}{k_B T_{ph, \lambda}(\tau)} \right)} - 1\right]^{-1}
    \label{eq:graphite-population-time-equiv}
$$
The expression for $n_{\lambda}(\vect{k}, \tau)$ makes use of the approximation that phonon vibrational frequencies remain constant through the experiments ($\omega_{\lambda}(\vect{k}, \tau) \approx \omega_{\lambda}(\vect{k}, \tau<0) ~ \forall \tau$), as discussed in @sec:graphite-debye-waller, which is approximately true in the case of graphite away from $\vect{\Gamma}$. @eq:graphite-population-time-equiv can be decomposed with a Laurent series [@Wunsch2005] to extract the quasi-linear relationship between mode population and temperature:
$$
    n_{\lambda}(\vect{k}, \tau) \propto  \frac{k_B T_{ph, \lambda}(\tau)}{\hbar \omega_{\lambda}(\vect{k}, \tau<0)} - 1/2 + \mathcal{O}\left( T^{-1}_{ph, \lambda}(\tau) \right)
    \label{eq:graphite-population-laurent}
$$
The above holds for appropriately-high mode temperatures. It follows that $\Delta n_{\lambda} \propto \Delta T_{ph,\lambda}$:
$$
    \Delta n_{\lambda}(\vect{k}, \tau) \propto \frac{k_B \Delta T_{ph, \lambda}(\tau)}{\hbar \omega_{\lambda}(\vect{k}, \tau<0)}
$$

### The non-thermal lattice model at the $\vect{K}$ point

Based on the formalism presented in the @sec:graphite-nlm, the couplings to the $A_1^\prime$ phonon mode will be extracted from the TO2 population measurements as an example. Let the differential population be $\Delta n_{\lambda=\text{TO2}}(\vect{k} = \vect{K}, \tau) \equiv \Delta n_{A_1^\prime}(\tau)$. $\Delta n_{A_1^\prime}(\tau)$ is obtained by integrating the population of the TO2 mode in a circular arc of radius \SI{0.3}{\per\angstrom} centered at $\vect{k} = \vect{K}$, as shown in @fig:graphite-ph-populations. 

The number of coupled equations in @eq:graphite-nlm can be reduced by aggregating lattice heat capacities into two categories: the heat capacity of $A_1^\prime$, $C_{A_1^\prime}$, and the total effective heat capacity of all other relevant modes, defined as $C_l$. The calculation of $C_l$ boils down to adding the contribution of mode that can scatter into $A_1^\prime$ or from it, based on conservation of energy and momentum, weighted by the decay probabilities reported by Bonini *et al.* [@Bonini2007] With this information:
\begin{align}
    C_l &= \frac{9}{100} \left[ C_{ph,\lambda=\text{TA}} + C_{ph,\lambda=\text{TA}} \right] \\
        &+ \frac{36}{100} \left[ C_{ph,\lambda=\text{TA}} + C_{ph,\lambda=\text{LA}}\right] \nonumber \\
        &+ \frac{55}{100} \left[ C_{ph,\lambda=\text{TA}} + C_{ph,\lambda=\text{LA}} + C_{ph,\lambda=\text{LO}} + C_{ph,\lambda=\text{LA}}\right] \nonumber
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
\label{eq:graphite-nlm-system}
$$

#### Heat capacities

Before @eq:graphite-nlm-system can be solved, the heat capacities of the electronic system and every relevant phonon mode must be parametrized. The electronic heat capacity is extracted from measurements by Nihira and Iwata [@Nihira2003]:
$$
C_e(T_e) = 13.8 ~ T_e(\tau) + 1.16 \times 10^{-3} ~ T_e^2(\tau) + 2.6 \times 10^{-7}  ~ T_e^3(\tau)
$$
Extracting the lattice specific heat is simplified by the observation that thermal expansion has not occurred on the time-scale of the measurements ($\tau < \SI{680}{\pico\second}$). Bragg peak positions remain static throughout the measurement, as was also reported by Chatelain *et al.*[@Chatelain2014a] The specific heat capacity of each mode can therefore be taken as the heat capacity at constant volume:
$$
C_{ph,\lambda}(T_{ph,\lambda}) = 
k_B \int_0^{\omega_D} d\omega ~ D_{\lambda}(\omega) 
    \left( 
        \frac{\hbar \omega}{k_B T_{ph,\lambda}} 
    \right)^2 
    \frac{e^{\hbar \omega / k_B T_{ph,\lambda}}}{\left( e^{\hbar \omega / k_B T_{ph,\lambda}} - 1\right)^2}
$$
where $\omega_D$ is the Debye frequency, and $D_{\lambda}(\omega)$ is the density of states for phonon mode $\lambda$[@Ziman1979]. The momentum-resolution of ultrafast electron diffuse scattering allows for a simplification: a single frequency contributes to the density of states, so that $D_{\lambda}(\omega) \equiv \delta(\omega)$.

### Heat rates {#sec:graphite-eph-solution}

The solution to @eq:graphite-nlm-system using the measurements of the $A_1^\prime$ population was computed using an iterative least-squares method[@Branch1999]. The resulting temperature traces $\set{ T_e(\tau), T_{A_1^\prime}(\tau), T_{l}(\tau)}$ are shown in [@fig:graphite-eph-coupling]. The coupling constants $\set{ G_{e,l}, G_{e,A_1^\prime}, G_{A_1^\prime, l}}$ are listed in @tbl:graphite-eph-coupling

```{.matplotlib #fig:graphite-eph-coupling file="figures/graphite/eph-coupling.py" caption="Evolution of the $A_1^\prime$ mode population in graphite after ultrafast photoexcitation. Transient population $\Delta n_{A_1^\prime}(\tau)$ is shown in black (circles). Error bars represent the standard error in the population mean before photoexcitation $\tau < 0$. The biexponential fit to the transient population is shown in blue (solid). The effective temperature of the modes that $A_1^\prime$ can decay into is shown in orange (dotted)."}
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

### Mode-projected electron-phonon coupling tensor elements {#sec:graphite-coupling-constants}

In order to compare to theory and other experiments, the value of electron-phonon coupling $g$ is calculated from the coupling constants $G$ determined from solving @eq:graphite-nlm-system. Refer to @sec:introduction-epc for an introduction to the electron-phonon coupling tensor.

Photoexcitation with \SI{1.55}{\electronvolt} photons will drive vertical electronic transitions, many of which might couple to a particular phonon mode. Let $\langle \cdots \rangle_{\gamma}$ be the average over all photoexcited electron states. The average electron-phonon coupling $\langle g^2_{e,\lambda}(\vect{k}) \rangle_{\gamma}$ is most simply related to the relaxation time $\tau_{e,\lambda}(\vect{k})$ between the two subsystems[@Sentef2013;@Na2019]:
$$
\frac{1}{\tau_{e,\lambda}(\vect{k})} = \frac{2 \pi}{\hbar} \langle g^2_{e,\lambda}(\vect{k}) \rangle_{\gamma} D_e(\hbar \omega_{\gamma} - \hbar \omega_{\lambda}(\vect{k}))
$${#eq:graphite-scattering-rate}
where $D_e(\epsilon)$ is the electronic density-of-states, $\hbar \omega_{\gamma}$ is the optical excitation energy (\SI{1.55}{\electronvolt} in the case of \SI{800}{\nano\meter} light), and $\omega_{\lambda}(\vect{k})$ is the vibrational frequency of phonon mode $\lambda$ at wavevector $\vect{k}$, as defined previously. 

The scattering rate in terms of heat rate constant $G$ is simply given by the conversion factor between energy and particle number, the heat capacity[@Nihira2003]:
$$
    \tau_{e,A_1^\prime} = \frac{G_{e,A_1^\prime}}{C_e}
$$
with associated error $\sigma_{\tau}$ related to the error in the heat rate $\sigma_G$[@Bevington2003error]:
$$
    \sigma_{\tau} = \left| \frac{\partial \tau}{\partial G}\right| \sigma_G = \frac{C_e}{G^2_{e,A_1^\prime}} \sigma_G
$$
Using the coupling constant $G_{e,A_1^\prime}$ and electronic heat capacity from @sec:graphite-eph-solution, $\tau_{e, A_1^\prime} = \SI{99 \pm 1}{\femto\second}$. Since the electronic heat capacity is highly dependent on temperature, the maximum value of the electronic temperature after photoexcitation was used (\SI{8500}{\kelvin} according to the solution to @eq:graphite-nlm-system).

In order to calculate the average electron-phonon coupling matrix element from the scattering rate, the electronic density-of-states needs to be approximated. Given the nature of the experiments presented here, as shown on @fig:graphite-photoexcitation, an approximation to the electronic density of states for graphite close to the Dirac point can be used[@Partoens2006;@Neto2009]:
$$
    D_e(\epsilon) = \frac{2 A}{\pi} \frac{|\epsilon|}{(\hbar v_{F})^2}
$$
where $A$ is the in-plane unit cell area and $v_F = \SI{9.06e5}{\meter \per \second}$ is the Fermi velocity[^neto-hbar]. Combining this form the the density-of-states with the measurement of the scattering rate in @eq:graphite-scattering-rate gives $\langle g^2_{e, A_1^\prime} \rangle_{\gamma} = \SI{0.032 \pm 0.001}{\square\electronvolt}$. This value compares favourably to calculations and experiments from other works, summarized in @tbl:graphite-eph-coupling-comparison.

\begin{table}
	\centering
	\caption{Comparison of measured and calculated values for the electron-phonon coupling matrix element $\langle g^2_{e, A_1^\prime}\rangle_{\gamma}$}
	\vspace{2mm}
	\begin{tabular}{l | c | l}
		Source & $\langle g^2_{e, A_1^\prime} \rangle_{\gamma}$ [\si{\electronvolt\squared}] & Notes \\ \hline\hline
        This work,                                           & $0.032 \pm 0.001$ & Experiment \\ \hline
        Piscanec \emph{et al.}\autocite{Piscanec2004}        & $< 0.0994$        & Theory (graphene, upper bound) \\ \hline
        Johannsen \emph{et al.}\autocite{Johannsen2013}      & $0.033 \pm 0.007$ & Experiment (trARPES, graphene) \\ \hline
        \multirow{2}{*}{Na \emph{et al.}\autocite{Na2019}}   & $0.050 \pm 0.011$ & Experiment (trARPES) \\ \cline{2-3} 
                                                             & $0.040$           & Theory \\ \hline
	\end{tabular} 
	\label{tbl:graphite-eph-coupling-comparison}
\end{table}

## Conclusion

In this chapter, a clear demonstration of the power of ultrafast electron diffuse scattering was presented. Harnessing the inherent redundancy of electron scattering patterns and prior knowledge of the crystal symmetry, ultrafast electron diffuse scattering measurements were endowed with energy resolution, resulting in a time-, energy-, and momentum-resolved view of phonon dynamics. The measurements were used to understand nonthermal phonon population dynamics. Finally, the phonon population dynamics at the $\vect{K}$ point was retrofitted into the non-thermal lattice model in order to extract mode-dependent electron-phonon and phonon-phonon coupling matrix elements, which were in good agreement with other experiments and calculations.

### Outlook

The procedure to extract phonon populations in @sec:graphite-ph-spectroscopy can be easily extended to other situations. One such situation is the case of thermal diffuse scattering measurements. At constant temperature, the phonon populations are directly related to their vibrational frequency; hence, the phonon band structure could be extracted. Phonon dispersion relations has been extracted before from x-ray diffuse scattering experiments[@Holt1999;@Xu2005], but these schemes are based on iterative (and unstable) fitting procedures. Using electron diffuse scattering instead presents some advantages, namely the inherently stronger scattering cross-section of electrons -- which might be the only way to observe diffuse scattering in monoyalers[@Caruso2021] -- and the relative commodity of electron microscopes, and the ability to directly invert the measurements matrix without fitting. Another situation where the method presented here can be extended is for the case of non-thermal phonon renormalization. For very early time-delays (say, $\tau < \SI{300}{\femto\second}$), the phonon populations may be considered somewhat constant. The dynamics effect on ultrafast electron diffuse scattering signals may be attributed to phonon renormalization, i.e. a change in vibrational frequency. In this case, the change $\Delta \omega_{\lambda}(\vect{k}, \tau)$ can be extracted while keeping the phonon populations fixed. This has been used by the author in the case of photoexcited titanium diselenide, where a change in electronic correlations can change the dielectric screening, thereby hardening a particular phonon mode [@Otto2021].

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]

[^neto-hbar]: Note that the factor of $\hbar$ has been erroneously ignored by Castro Neto *et al.*[@Neto2009]