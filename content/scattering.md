
# The theory of ultrafast electron scattering{#sec:scattering}

The scattering of light or particles by a sample is used by a large class of experimental techniques, dating back a hundred years. Scattering can be broadly defined as the modification of an incoming wave by a potential into outgoing wave, a process which imprints the outgoing wave with some characteristic of the potential. The outgoing wave may lose or gain energy, and its momentum might be changed. When multiple incoming waves are simultaneously used -- forming an incoming wavefront --, the outgoing waves may interfere constructively or destructively. This effect is particularly intense for periodic scattering potentials, for example in crystals.

This chapter will consider the special case of *electron scattering*. In crystals, electrons are scattered by the electrostatic potential of ions and the electronic charge-density. Thanks to improvements to instrument stability (@sec:intro-cavity), as well as advances in data acquisition and data analysis attributable to the author (@sec:introduction-data-exploration), the ultrafast electron scattering instrument used herein can reliably measure the effects of *diffuse scattering*.

In writing this chapter, the author has tried original derivations that emphasize concepts that are important for the remainder of this dissertation. The derivation of ultrafast diffuse scattering intensity is of particular interest, because it is the only full quantum-mechanical treatment relevant to ultrafast electron scattering specifically in the literature today.

## Electron scattering and the Lippmann-Schwinger formalism

Consider an electron wavefunction $\Psi(\vect{x}, t)$. The scattering of $\Psi(\vect{x}, t)$ by an arbitrary potential $V(\vect{x},t)$ is described by the Schrödinger equation:
$$
i \hbar \frac{d}{d t} \Psi(\vect{x}, t) = \left[ \frac{- \hbar^2}{2 m_e} \nabla^2 + V(\vect{x}, t) \right]\Psi(\vect{x}, t).
$${#eq:scattering-schroedinger}

The scattering electrons used in this work have an enormous amount of kinetic energy (\SI{90}{\kilo\electronvolt}). The potential energy from the atomic Coulomb interaction, on the other hand, is much smaller: 
$$
    V = \frac{Z e^2}{4 \pi \epsilon_0 |\vect{x}|}
$${#eq:scattering-electrostatic-potential}
where $\epsilon_0$ is the vacuum dielectric constant, $Z$ the atomic number, $e$ is the unit of charge. For the simple example of graphite ($Z=12$, $|\vect{x}|\approx \SI{1}{\angstrom}$), the potential energy associated with the Coulomb interaction is less than \SI{10}{\electronvolt}. Given that scattering electrons have nearly \SI{100}{\kilo\electronvolt}, the potential $V(\vect{x}, t)$ shall be treated as a perturbation.

### Electrons propagating in free space

In the case of free space ($V(\vect{x}, t) = 0$), [@eq:scattering-schroedinger] reduces to the following equation:
$$
i \hbar \frac{d}{d t} \Psi(\vect{x}, t) = \frac{- \hbar^2}{2 m_e} \nabla^2 \Psi(\vect{x}, t).
$${#eq:scattering-free-electron}
It is instructive to consider the energy eigenfunction satisfying [@eq:scattering-free-electron]. Let the solutions be labeled as $\set{(\Psi_a(\vect{x}, t), \omega_a)}$, where the associated energy eigenvalues are $\omega_a \equiv E_a/\hbar$. The energy eigenfunction can be factorized into spatial and temporal parts as:
$$
\Psi_a(\vect{x}, t) = u_a(\vect{x}) e^{-i \omega_a t}
$${#eq:scattering-free-space}
where $u_a(\vect{x})$ solves the time-independent Schrödinger equation[@Schroedinger1926]:
$$
\hbar \omega_a u_a(\vect{x}) = \frac{- \hbar^2}{2 m_e} \nabla^2 u_a(\vect{x}).
$${#eq:scattering-stationary-schroedinger}
[@eq:scattering-stationary-schroedinger] is an instance of the Helmholtz equation, which can be re-written as:
$$ 
\left[ \nabla^2 + k_a^2\right] u_a(\vect{x}) = 0
$${#eq:scattering-helmholtz}
where $k_a^2 \equiv 2 m_e E_a/\hbar^2$. Physical reasoning and the classical result of the Helmholtz equation reveal that $k_a$ is the *wavevector*, related to the wave momentum as $k_a^2 \equiv \vect{k}_a \cdot \vect{k}_a \equiv \vect{p}_a^2/\hbar^2$. [Eq. @eq:scattering-helmholtz] can be solved using separation of variables, where the solution along linearly-independent directions are considered independently:
$$
u_a(\vect{x})\equiv u_{a,x}(\vect{x} \cdot \hat{\vect{x}}) u_{a,y}(\vect{x} \cdot \hat{\vect{y}}) u_{a,z}(\vect{x} \cdot \hat{\vect{z}}).
$$
[Eq. @eq:scattering-helmholtz] then becomes:
$$
\sum_{i\in\{x,y,z\}}\frac{1}{u_{a,i}(\vect{x} \cdot \hat{\vect{i}})} \frac{d^2 u_{a,i}}{di^2}(\vect{x} \cdot \hat{\vect{i}}) + |\vect{k} \cdot \hat{\vect{i}}| = 0
$$
and the solution can be synthesized as a product of one-dimensional plane wave:
\begin{align}
u_a(\vect{x}) &= \prod_{j \in \{x,y,z\}} A_j e^{i (\vect{k}_a \cdot \hat{\vect{j}}) (\vect{x} \cdot \hat{\vect{j}})} \nonumber\\
             &= A e^{i \vect{k}_a \cdot \vect{x}}
\label{eq:scattering-plane-wave}
\end{align}
for some scalars $\{A_x, A_y, A_z, A\}$. Combining [@eq:scattering-free-space] [@eq:scattering-plane-wave] leads to the energy eigenstate of a free electron propagating in vacuum:
$$
\Psi_a(\vect{x}, t) = A e^{i (\vect{k}_a \cdot \vect{x} - \omega_a t)}
$$
with associated energy eigenvalue $E_a = \hbar^2 \vect{k}_a^2 / 2 m_e=\hbar \omega_a$.

### The Lippmann-Schwinger equation

The problem of scattering with a non-trivial potential $\hat{V}$ will now be considered. This is a fundamental problem in quantum mechanics; approximations appropriate for the research presented in this dissertation will be made to simplify and clarify the presentation. In particular, the following derivations assume that the scattering potential $V(\vect{x})$ is *local*, that is:
$$
    \bra{\vect{x}^\prime} \hat{V} \ket{\vect{x}^{\prime \prime}} = V(\vect{x}^\prime) \delta(\vect{x}^\prime - \vect{x}^{\prime \prime})
$$
This is a reasonable assumption as the electrostatic potential has a $\sfrac{1}{|\vect{x}|}$ dependence (@eq:scattering-electrostatic-potential). Further simplifications can be made if allowed momentum states are assumed to be defined in a (large) cube of size-length $L$, such that
$$
    \braket{\vect{x} | \vect{k}} = \frac{1}{L^{3/2}} e^{i \vect{k} \cdot \vect{x}}
$${#eq:scattering-norm}
and $\vect{k}$ takes discrete values. The calculations will become physical after taking the limit $L \to \infty$. In this approximation, the scattered wave $\braket{\vect{x} | \Psi}$ is given by the *Lippmann-Schwinger* equation[@Lippmann1950]:
$$
    \braket{\vect{x} | \Psi} = \braket{\vect{x} | \vect{k}_i} - \frac{2 m_e}{\hbar^2} \int d^3 x^\prime \frac{e^{i|\vect{k}_i||\vect{x}-\vect{x}^\prime|}}{4 \pi |\vect{x}-\vect{x}^\prime|} V(\vect{x}^\prime) \braket{\vect{x}^\prime | \Psi}
$${#eq:scattering-lippmann-schwinger}
where $\ket{\vect{k}_i}$ represents the initial plane-wave of incoming scattering electrons, and $\ket{\Psi}$ is the final scattered state. Provided that the scattered wavefunction $\braket{\vect{x} | \Psi}$ is measured far from where the scattering potential is localized:
$$
    e^{i|\vect{k}||\vect{x}-\vect{x}^\prime|} \approx e^{ikr}e^{-i\vect{k}_i \cdot \vect{x}^\prime}
$$
where $r \equiv |\vect{x}-\vect{x}^\prime|$ and $k \equiv |\vect{k}|$. @eq:scattering-lippmann-schwinger can then be simplified as:
$$
    \braket{\vect{x} | \Psi} = \braket{\vect{x} | \vect{k}_i} - \frac{m_e}{2 \pi \hbar^2} \frac{e^{i k_f r}}{r}\int d^3 x^\prime e^{-i \vect{k}_f \cdot \vect{x}^\prime} V(\vect{x}^\prime) \braket{\vect{x}^\prime | \Psi}
$${#eq:scattering-lippmann-schwinger-rad}
The canonical description of the Lippmann-Schwinger equation is usually reduced to:
$$
    \braket{\vect{x} | \Psi} = \frac{1}{L^{3/2}} \left[ e^{i \vect{k}_i \cdot \vect{x}} + \frac{e^{i k_f r}}{r} f(\vect{k}_f, \vect{k}_i)\right]
$${#eq:scattering-lippmann-schwinger-general}
where
\begin{align}
    f(\vect{k}_f, \vect{k}_i) & \equiv -\frac{m_e L^3}{2 \pi \hbar^2} \int d\vect{x}^\prime \frac{e^{-i \vect{k}_f \cdot \vect{x}^\prime}}{L^{3/2}} V(\vect{x}^\prime) \braket{\vect{x}^\prime | \Psi} \nonumber \\
                            & =  -\frac{m_e L^3}{2 \pi \hbar^2} \bra{\vect{k}_f} \hat{V} \ket{\Psi}
    \label{eq:scattering-amplitude}
\end{align}
$f(\vect{k}_f, \vect{k})$ is called the *scattering amplitude*. In this notation, the vector $\vect{k}_f$ is a formal variable, and not a known state like $\vect{k}_i$. The form of @eq:scattering-lippmann-schwinger-general complies with intuition: the final scattered wavefunction is composed of an unscattered part ($\braket{\vect{x} | \vect{k}_i} \propto e^{i \vect{k}_i \cdot \vect{x}}$) as well as an outgoing spherical wave with amplitude $f(\vect{k}_f, \vect{k}_i)$ in the $\vect{k}_f$ direction.

The calculation of the scattered wavefunction $\braket{\vect{x}|\Psi}$ has been reduced to the calculation of $\bra{\vect{k}_f} \hat{V} \ket{\Psi}$ for arbitrary values of $\vect{k}_f$. The derivation of an expression for this is beyond the scope of this work, and the final result is stated[^tmatrix]:
$$
    \bra{\vect{k}_f} \hat{V} \ket{\Psi} = \bra{\vect{k}_f} \left[ \sum_{n=1}^{\infty} \hat{V} \left( \frac{1}{E_i - \hat{H}_0 + i \epsilon} \hat{V}\right)^{n-1} \right] \ket{\vect{k}_i}
$${#eq:scattering-potential-decomp}
where $\hat{H}_0$ is the free-space Hamiltonian with eigenvalue $E_i=\hbar^2 |\vect{k}_i|^2/2 m_e$, and $\epsilon$ is a vanishingly small real number. In particular, each term with index $n$ in the sum of @eq:scattering-potential-decomp corresponds to the electron scattering $n$ times[@Feynman1965].

### Measuring the scattered wavefunction

Electron cameras measure the intensity of the wavefunction. In the case of bright-field electron microscopy, where the wavefunction is imaged directly, the intensity on the detector is given by:
$$
    I(\vect{x}) \equiv \left| \braket{\vect{x} | \Psi} \right|^2
$$
In order to sample the scattered wavefunction in reciprocal space, an electron lens can be used to focus the scattered electrons onto the detector. Given that electrons are prepared with definite momentum $\vect{k}_i$, it is trivial to ignore the unscattered part of the wavefunction -- the first term in @eq:scattering-lippmann-schwinger-general -- which is found only at $\vect{k}_f = \vect{k}_i$. Therefore, the intensity away from $\vect{k}_i$ is related only to the scattered wavefunction, $\Psi$:
\begin{align}
I(\vect{k}_f - \vect{k}_i)
    & = \left| \braket{\vect{x} | \Psi} \right|^2 \nonumber \\
    & = \left| \frac{e^{ikr}}{r} f(\vect{k}_f, \vect{k}_i) \right|^2 \nonumber \\
    & = \frac{1}{r^2} |f(\vect{k}_f, \vect{k}_i)|^2
\label{eq:scattering-intensity}
\end{align}
Therefore, the diffracted intensity $I(\vect{k}_f - \vect{k}_i)$ is proportional to the square of the scattering amplitude. 

Note that the factor of $1/r^2$ is generally ignored[@Fultz2002r2]. For the instrument configuration presented in @sec:experimental_setup, $1/r^2$ varies from \SIrange{15.92}{16}{\per\square\meter}, from the corner to the center of the detector respectively. While not an insignificant variation, the experiments presented herein generally measure the relative change in intensity, in which case factors are not important.

## Elastic scattering in a crystal

In this section, the consequences of an electron scattering *once* in crystalline solid will be explored. In this approximation, historically called the *first Born approximation* [@Born1926], only the first term in @eq:scattering-potential-decomp ($n=1$) is considered:
\begin{align}
    f^{(1)}(\vect{k}_f, \vect{k}_i) 
        & = -\frac{m_e L^3}{2 \pi \hbar^2} \bra{\vect{k}_f} \hat{V} \ket{\vect{k}_i} \nonumber \\
        & = -\frac{m_e L^3}{2 \pi \hbar^2} \int d\vect{x}^\prime \bra{\vect{k}_f} V(\vect{x}^\prime) \ket{\vect{x}^\prime}\braket{\vect{x}^\prime | \vect{k}_i} \nonumber \\
        & = -\frac{m_e L^3}{2 \pi \hbar^2} \int d\vect{x}^\prime \frac{e^{i(\vect{k}_i - \vect{k}_f)\cdot \vect{x}^\prime}}{L^3} V(\vect{x}^\prime) \nonumber \\
        & = -\frac{m_e}{2 \pi \hbar^2} \int d\vect{x}^\prime e^{i(\vect{k}_i - \vect{k}_f)\cdot \vect{x}^\prime} V(\vect{x}^\prime)
    \label{eq:scattering-first-born-approx}
\end{align}
where the normalization of @eq:scattering-norm was used. The reader may recognize that the scattering amplitude $f^{(1)}(\vect{k}_f, \vect{k}_i)$ is proportional to the Fourier transform of the scattering potential with respect to $\vect{k}_f - \vect{k}_i \equiv \vect{q}$, the *scattering vector*. If the Fourier transform functional operator is defined as:
$$
    \mathcal{F}\left[ f(\vect{x}) \right] \equiv \ft{f}(\vect{q}) = \frac{1}{2 \pi} \int d\vect{x}^\prime e^{-i \vect{q} \cdot \vect{x}^\prime}f(\vect{x}^\prime).
$$
then @eq:scattering-first-born-approx can be simplified to:
$$
    f^{(1)}(\vect{q}=\vect{k}_f - \vect{k}_i) = -\frac{m_e}{\hbar^2} \ft{V}(\vect{q})
$${#eq:scattering-amplitude-q}

### Scattering potential of a single atom{#sec:affe}

The scattering potential of a single atom is given by:
$$
    V_a(\vect{x}) = -\frac{Z e^2}{|\vect{x}|} + \sum_{i=1}^{Z} \frac{e^2}{|\vect{x} - \vect{x}_i|}
$${#eq:scattering-atom-potential}
where Z is the atomic weight, $\vect{x}$ is a position with respect to the ionic core, and $\vect{x}_i$ is the position of the $i^{\text{th}}$ electron. The potential in @eq:scattering-atom-potential can be calculated from first principles, using relativistic Hartree-Fock calculations[@Fischer1977;@Kirkland2010] to get the real-space electron density (i.e. determining the possible set $\set{\vect{x}_i}$). This is beyond the scope of the present work, and the final result is used here. When discussing electron scattering, the scattering amplitude for a single atom is usually called the *atomic form factors for electron scattering*. To a reasonable degree of accuracy, the atomic form factors for electrons for light atoms are spherically symmetric[@Zheng2009]. They can be parametrized as:
$$
    f_e(\vect{q}) = \sum_{i=1}^3 \frac{a_i}{|\vect{q}|^2 + b_i} + c_i e^{-d_i|\vect{q}|^2}
$${#eq:scattering-affe}
where the constants $\set{ a_i, b_i, c_i, d_i }$ are element-specific fitting parameters which are tabulated in Kirkland [@Kirkland2010]. The associated real-space potential can be calculated via @eq:scattering-amplitude-q. Examples of $f_e(\vect{q})$ and associated $V_a(\vect{x})$ are shown in @fig:scattering-potential for a few elements.

```{.matplotlib #fig:scattering-potential file="figures/scattering/scatt-potential.py" caption="Demonstration of the electrostatic potential of atoms, which scatters electrons. **a)** Radial view of the electrostatic potential in real-space **b)** Radial view of the electrostatic potential in reciprocal space, also known as the *atomic form factor*."}
```

The contribution of individual electronic orbitals to the atomic form factor for electrons is discussed in Zheng *et al.*[@Zheng2009] 

### Scattering potential of a crystal

The scattering potential of a crystalline lattice can be expressed as:
$$
    V_c(\vect{x}) = \sum_i V_{a,i}(\vect{x} - \vect{r}_i)
$$
where  the sum index $i$ runs over atoms in the crystal with positions $\vect{r}_i$. The potential for each atom $V_{a,i}$ is taken to be more general than the single-atom potential of @eq:scattering-atom-potential, in order to model ion-ion interactions (e.g. bonding). Note that for any function $h(\vect{x})$: 
\begin{align}
    \mathcal{F}\left[ h(\vect{x} + \vect{y}) \right] 
        & = \frac{1}{2\pi} \int d\vect{x}^\prime e^{-i\vect{q}\cdot(\vect{x}^\prime + \vect{y})} h(\vect{x}^\prime) \nonumber \\
        & = \frac{e^{-i \vect{q} \cdot \vect{y}}}{2 \pi} \int d\vect{x}^\prime e^{i\vect{q} \vect{x}^\prime} h(\vect{x}^\prime) \nonumber \\
        & = e^{-i \vect{q} \cdot \vect{y}} \mathcal{F}\left[ h(\vect{x}) \right]
    \label{eq:scattering-fourier-shift}
\end{align}
where $\vect{y}$ is some arbitrary translation vector. Therefore, the Fourier transform of the scattering potential of the entire crystal is related to the Fourier transform  the potential its constituent atoms (@eq:scattering-atom-potential) as:
\begin{align}
    \ft{V}_c(\vect{q}) 
        & = \mathcal{F}\left[ \sum_i V_{a,i}(\vect{x} + \vect{r}_i) \right] \nonumber \\
        & = \sum_i \mathcal{F}\left[ V_{a,i}(\vect{x} + \vect{r}_i) \right] \nonumber \\
        & = \sum_i \mathcal{F}\left[ V_{a,i}(\vect{x}) \right] e^{-i \vect{q} \cdot \vect{r}_i} \nonumber \\
        & = \sum_i f_{e,i}(\vect{q}) e^{-i \vect{q} \cdot \vect{r}_i}
    \label{eq:scattering-potential-crystal}
\end{align}

```{.matplotlib #fig:scattering-polonium-example file="figures/scattering/polonium.py" caption="Calculated scattering potential and associated scattering amplitude for an abstract crystal. **a)** Electrostatic potential $V(\vect{x})$ in the $z=0$ plane. The two in-plane lattice vectors $\vect{a}_1$ and $\vect{a}_2$ are shown; lattice vector $\vect{a}_3$ points out of the page. **b)** Scattering amplitude $f(\vect{q})$ associated with the electrostatic potential shown in a). The periodic nature of the potential in real-space creates a structure in reciprocal space called the *reciprocal lattice*."}
``` 

@eq:scattering-potential-crystal has historically been called the *static structure factor*, named thus because the atomic positions are assumed to be fixed. This stands in contrast to the *dynamic structure factor* discussed below in @sec:diffuse-scattering. A visual representation of the scattering potential of a crystal is helpful. Consider the example of a crystal of Polonium, consists in a Pu atom at the every vertex of a rectangular prism of dimensions \SI{5 x 3 x 3}{\angstrom}. The calculated electrostatic potential of this arrangement along the unit cell base is shown in @fig:scattering-polonium-example a). The lattice vectors $\vect{a}_1$ and $\vect{a}_2$ are indicated, with $\vect{a}_3$ being aligned out of the page. The periodic nature of this scattering potential is demonstrated by calculating the resulting scattering amplitude $f_e(\vect{q})$ from @eq:scattering-amplitude-q, which is shown in @fig:scattering-polonium-example b). The periodicity in spatial-frequency-space, also called *reciprocal space*, is evident, and forms a *reciprocal lattice*, formally defined in the next section.

From @eq:scattering-intensity and @eq:scattering-amplitude-q, the measured diffracted intensity is:
$$
    I(\vect{q}) = \frac{m_e^2}{4 \pi^2 \hbar^4 r^2} \left| \sum_i f_{e,i}(\vect{q}) e^{-i \vect{q} \cdot \vect{r}_i} \right|^2
$${#eq:scattering-diffracted-intensity-zero-temp}
This is the standard result for the diffracted intensity being proportional to the square of the static structure factor[@Warren1990intensity;@Kittel1996intensity;@Fultz2002intensity;@Kirkland2010intensity]. 

### The reciprocal lattice

The geometry of reciprocal space and the reciprocal lattice are foundational concepts that drive the understanding of ultrafast electron diffraction.

A perfectly periodic structure in real-space that extends to infinity, with associated lattice vectors $\set{ \vect{a}_i }$, possesses a *dual* lattice in reciprocal space. The lattice vectors $\set{ \vect{b}_i }$ for this reciprocal lattice are defined by the relation
$$
    \vect{a}_i \cdot \vect{b}_j = 2 \pi \delta_{ij}
$$
which leads to the following reciprocal lattice vectors:
$$
    \left\{ \vect{b}_i = 2 \pi \frac{\vect{a}_j \times \vect{a}_k}{\vect{a_i} \cdot (\vect{a_j} \times \vect{a}_k)} ~ \middle| ~ (i,j,k) \in C \right\}
$$
where $C$ is the set of cyclic permutations. For the example of $\alpha$-Pu, where
$$
    \vect{a}_1 = 3.63 ~ \hat{\vect{x}}, ~ \vect{a}_2 = 3.63 ~ \hat{\vect{y}}, ~ \vect{a}_3 = 3.63 ~ \hat{\vect{z}}
$${#eq:scattering-polonium-lattice}
the associated reciprocal lattice vectors are
$$
    \vect{b}_1 = \frac{2\pi}{3.63} ~ \hat{\vect{x}}, ~ \vect{b}_2 = \frac{2\pi}{3.63} ~ \hat{\vect{y}}, ~ \vect{b}_3 = \frac{2\pi}{3.63} ~ \hat{\vect{z}}.
$${#eq:scattering-polonium-recip}
The geometry of vectors in @eq:scattering-polonium-lattice and @eq:scattering-polonium-recip are shown in @fig:scattering-polonium-example. The position of *reciprocal points* $\vect{H}$ -- the location of the fundamental frequencies of the Fourier transform of @eq:scattering-amplitude-q -- is a linear combination of vectors reciprocal basis vectors $\set{ \vect{b}_i }$:
$$
    \vect{H} = h ~ \vect{b}_1 + k ~ \vect{b}_2 + l ~ \vect{b}_3 
$$
where $h$, $k$, and $l$ are all integers. Expressed in the reciprocal basis, reciprocal points are traditionally denoted as $\vect{H} = (hkl)$. The indices $h$, $k$, and $l$ are called *Miller indices*, named for W. H. Miller[@Miller1839]. 

#### Diffraction for large crystals

Now that the reciprocal points $\set{\vect{H}}$ have been introduced, an alternative form of @eq:scattering-diffracted-intensity-zero-temp for large crystals can be written down. Consider the counting of atoms to change from label $i$ to labels $(m, s)$, where $m$ labels unit cells and $s$ labels unit cell atoms. In this case, the atomic positions can be written as:
$$
    \vect{r}_i \equiv \vect{r}_{m,s} = \vect{R}_m + \vect{x}_s
$$
where $\vect{R}_m$ is the position of unit cell $m$, and $\vect{x}_s$ is the position of atom $s$ with respect to the unit cell origin. Then, @eq:scattering-diffracted-intensity-zero-temp becomes:
\begin{align}
    I(\vect{q}) & = \frac{m_e^2}{4 \pi^2 \hbar^4 r^2} \left| \sum_i f_{e,i}(\vect{q}) e^{-i \vect{q} \cdot \vect{r}_i} \right|^2 \nonumber \\
               & = \frac{m_e^2}{4 \pi^2 \hbar^4 r^2} \left| \sum_{m,s} f_{e,s}(\vect{q}) e^{-i \vect{q} \cdot (\vect{R}_m + \vect{x}_s)} \right|^2
\end{align}

Finally, note that because the vectors $\set{ \vect{R}_m }$ are integer multiples of lattice vectors:
$$
    \sum_{m=1}^{N_c} e^{-i \vect{q} \cdot \vect{R}_m} \xrightarrow[]{N_c \to \infty} N_c \sum_{\set{\vect{H}}} \delta(\vect{q} - \vect{H})
$${#eq:scattering-discrete-fourier}
which follows from the definition of the Fourier transform[@Robinson2016DiscreteFourier]. Therefore, for large $N_c$ (large crystals):
$$
    I(\vect{q}) = \frac{N_c^2 m_e^2}{4 \pi^2 \hbar^4 r^2} \left| \sum_{\set{\vect{H}}}\sum_{s} f_{e,s}(\vect{q}) e^{-i \vect{q} \cdot \vect{x}_s} \delta(\vect{q} - \vect{H})\right|^2
$$
This form of @eq:scattering-diffracted-intensity-zero-temp makes it more obvious why, for large periodic structures, the scattering pattern is composed of Bragg peaks at reciprocal points (@fig:scattering-polonium-example).

### Bragg's law

It is possible to *deduce* Bragg's law from @eq:scattering-amplitude-q and our definition of the reciprocal lattice.

By definition, the reciprocal points of the crystal scattering potential, located at the spatial frequencies where the crystal potential strong, form a lattice with basis vectors $\set{ \vect{b}_i }$. Consider an electron with initial wavevector $\vect{k}_i$ that scatters elastically to a final wavevector $\vect{k}_f$. The scattering amplitude for this event, $f(\vect{k}_f, \vect{k}_i)$ is most intense where $\hat{V}(\vect{k}_f - \vect{k}_i)$ is strong; that is, the condition for strong scattering is:
$$
    \vect{k}_f - \vect{k}_i = h ~ \vect{b}_1 + k ~ \vect{b}_2 + l ~ \vect{b}_3 \quad \forall ~ h,k,l \in \mathbb{Z}
$${#eq:scattering-bragg-vector}
This is precisely the *vector* form of Bragg's law [@Warren1990Bragg]. To recover the canonical form of Bragg's law, consider that an electron state with wavevector $\vect{k}$ can be associated with a wavelength of $\lambda = \tfrac{2 \pi}{|\vect{k}|}$. Since for elastic scattering, $|\vect{k}_i| = |\vect{k}_f| = \tfrac{2 \pi}{\lambda}$, @eq:scattering-bragg-vector becomes:
$$
    \frac{2 \pi}{\lambda} \left( \hat{\vect{k}}_f - \hat{\vect{k}}_i\right) = \vect{H}
$${#eq:scattering-bragg-1}
where $\hat{\vect{k}}$ denotes a unit-length vector in the direction of $\vect{k}$. Given that the vectors on both sides of the equation have the same magnitude and direction, the direction of $\vect{H}$ must bisect the angle between $\vect{k}_f$ and $\vect{k}_i$, historically defined as $2\theta$. Taking the amplitude of @eq:scattering-bragg-1:
$$
     \frac{2\pi}{\lambda} \left| \hat{\vect{k}}_f - \hat{\vect{k}}_i\right| = \frac{2 \pi \sin{\theta}}{\lambda} 
$$
and
$$
|\vect{H}| = \frac{1}{d_{hkl}},
$$
which can be combined as
$$
    \frac{4 \pi \sin{\theta}}{\lambda} = \frac{1}{d_{hkl}}.
$${#eq:scattering-bragg-hist}
@eq:scattering-bragg-hist is the historical form of Bragg's law as it relates to polycrystalline diffraction patterns[@Bragg1913]. Note that the vector form of @eq:scattering-bragg-vector is richer than the original form of Bragg's law as it places constraint on the full three-dimensional direction of the scattering vector $\vect{q} = \vect{k}_f - \vect{k}_i$.

### The Ewald sphere{#sec:scattering-ewald-sphere}

Elastic electron scattering, or *electron diffraction*, can be discussed more concretely. Consider an electron initially propagating in the $\hat{\vect{z}}$ direction with wavevector $\vect{k}_i$ that interacts with a scattering potential $\ft{V}(\vect{q})$, and scatters to a final wavevector $\vect{k}_f$. The elastic scattering condition $|\vect{k}_i| = |\vect{k}_f|$ constrains the observation of $\ft{V}(\vect{q})$ to scattering vectors $\vect{q}$ that lie on a sphere of radius $|\vect{q}|=\tfrac{1}{\lambda}$. This sphere is called the *Ewald sphere*[@Ewald1921].

```{.matplotlib #fig:scattering-ewald-sphere file="figures/scattering/ewald.py" caption="Demonstration of the Ewald sphere, a visual representation of the conservation of energy in diffraction. The Fourier transform of the scattering potential from an abstract cubic lattice of side length \SI{5}{\angstrom}, $\ft{V}(\vect{q})$, is shown in the background, with the associated reciprocal lattice vectors $\set{\vect{b}_i}$. The Ewald sphere of radius $\vect{q}$ is shown for two scatterers: electrons (solid) and hard x-ray (dashed)."}
```

The Ewald sphere is a great mental model of the information contained in diffraction patterns. Because diffracting electrons can only sample scattering vectors on the Ewald sphere, any particular measurement of a scattering potential $V(\vect{x})$ is effectively a two-dimensional *slice* of the three-dimensional Fourier transform of $V(\vect{x})$, $\ft{V}(\vect{q})$. This is represented in @fig:scattering-ewald-sphere. In this figure, the potential $\ft{V}(\vect{q})$ for an idealized simple cubic crystal with side-length \SI{5}{\angstrom} is shown in the plane spanned by $\vect{b}_2$ and $\vect{b}_3$. The Ewald spheres associated with \SI{100}{\kilo\electronvolt} electrons (large $|\vect{q}|$) and \SI{13}{\kilo\electronvolt} x-rays (smaller $|\vect{q}|$) are also shown. This electron energy is typical of the work presented in this dissertation, while the x-ray energy is an upper bound on the available energies at the Linac Coherent Light Source as of 2021[@Bostedt2013]. The reciprocal points that intersect the Ewald sphere appear in measurements as diffraction peaks, or Bragg peaks. @fig:scattering-ewald-sphere shows the advantage of electron scattering to study two-dimensional materials: given the proper orientation of the electron beam, a large range of wavevectors can be studied in the plane of interest.

## Multiple scattering of electrons {#sec:scattering-multiple}

Electrons interact quite strongly with matter through the Coulomb interaction. For scattering targets which are thick enough, an electron may scatter more than once before exiting the scattering potential volume. In this section, the scattering of an electron *twice* will be considered. In this case, the second term in @eq:scattering-potential-decomp ($n=2$) is considered:
$$
    f^{(2)}(\vect{k}_f, \vect{k}_i) = -\frac{m_e L^3}{2 \pi \hbar^2} \bra{\vect{k}_f} \hat{V} \frac{1}{E_i - \hat{H}_0 + i \epsilon}\hat{V} \ket{\vect{k}_i}
$$
The calculation of $f^{(2)}(\vect{k}_f, \vect{k}_i)$ involves the insertion of two complete sets of basis states:
\begin{multline}
    \bra{\vect{k}_f} \hat{V} \frac{1}{E_i - \hat{H}_0 + i \epsilon} \hat{V} \ket{\vect{k}_i} = \\ 
        \int d\vect{x}^\prime \int d\vect{x}^{\prime\prime} 
            \braket{\vect{k}_f | \vect{x}^\prime} V(\vect{x}^\prime) \bra{\vect{x}^\prime} \frac{1}{E_i - \hat{H}_0 + i \epsilon} \ket{\vect{x}^{\prime\prime}}V(\vect{x}^{\prime\prime}) \braket{\vect{x}^{\prime\prime} | \vect{k}_i}
\end{multline}
The evaluation of the term $\bra{\vect{x}^\prime} \frac{1}{E_i - \hat{H}_0 + i \epsilon} \ket{\vect{x}^{\prime\prime}}$ naturally happens in the derivation of @eq:scattering-lippmann-schwinger, and so the result is simply stated here:
$$
    \bra{\vect{x}^\prime} \frac{1}{E_i - \hat{H}_0 + i \epsilon} \ket{\vect{x}^{\prime\prime}} = -\frac{m_e}{2 \pi \hbar^2} \frac{e^{i |\vect{k}_i| |\vect{x}^\prime - \vect{x}^{\prime\prime}|}}{|\vect{x}^\prime - \vect{x}^{\prime\prime}|}
$$
It follows that the scattering amplitude for the double-scattering of an electron is given by:
$$
    f^{(2)}(\vect{k}_f, \vect{k}_i) =
        \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vect{x}^\prime \int d\vect{x}^{\prime\prime} 
            e^{-i\vect{k}_f \cdot \vect{x}^\prime} V(\vect{x}^\prime) \left( \frac{e^{i |\vect{k}_i| |\vect{x}^\prime - \vect{x}^{\prime\prime}|}}{|\vect{x}^\prime - \vect{x}^{\prime\prime}|} \right) e^{i\vect{k}_i \cdot \vect{x}^{\prime\prime}} V(\vect{x}^{\prime\prime})
$${#eq:scattering-amplitude-mult}
The structure of @eq:scattering-amplitude-mult informs on the following physical interpretation. Double scattering involves a first scattering at $\vect{x}^{\prime\prime}$ ($V(\vect{x}^{\prime\prime})$), which "radiates" as a spherical wave moving from $\vect{x}^{\prime\prime}$ to $\vect{x}^{\prime}$ ($e^{i |\vect{k}_i| |\vect{x}^\prime - \vect{x}^{\prime\prime}|}/|\vect{x}^\prime - \vect{x}^{\prime\prime}|$), followed by a second scattering at $\vect{x}^\prime$ ($V(\vect{x}^{\prime})$).

### Comparing cross-sections

The computation of the differential scattering cross section for two elastic scattering events is given by:
\begin{align}
    \frac{d\sigma_2}{d\Omega} & = |f^{(2)}(\vect{k}_f, \vect{k}_i)|^2 \nonumber \\
             & = \left| \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vect{x}^\prime \int d\vect{x}^{\prime\prime} 
            e^{-i\vect{k}_f \cdot \vect{x}^\prime} V(\vect{x}^\prime) \left( \frac{e^{i |\vect{k}_i| |\vect{x}^\prime - \vect{x}^{\prime\prime}|}}{|\vect{x}^\prime - \vect{x}^{\prime\prime}|} \right) e^{i\vect{k}_i \cdot \vect{x}^{\prime\prime}} V(\vect{x}^{\prime\prime}) \right|^2
    \label{eq:scattering-cross-sec-2}
\end{align}
The evaluation of @eq:scattering-cross-sec-2 for a realistic potential $V(\vect{x})$ is arduous. However, given that this type of scattering is *undesirable*, getting an upper bound on its scattering cross-section is a worthwhile exercise. 

In electron scattering experiment with thin specimens, the elastic scattering cross section is small enough that an electron is unlikely to scatter twice from the *same atom*. This means that the integral is 0 when $|\vect{x}^\prime - \vect{x}^{\prime\prime}|$ is small. The inner integral over $\vect{x}^{\prime\prime}$ in @eq:scattering-cross-sec-2 can be split:
$$
    \int d\vect{x}^{\prime} \int d\vect{x}^{\prime\prime} \to \int d\vect{x}^{\prime} \left[ \int_{|\vect{x}^\prime - \vect{x}^{\prime\prime}| \leq a} d\vect{x}^{\prime\prime} + \int_{|\vect{x}^\prime - \vect{x}^{\prime\prime}| > a} d\vect{x}^{\prime\prime} \right]
$$
where $a$ is the inter-atomic distance of the crystal. In the approximation described above, where double-scattering is only possible for $|\vect{x}^\prime - \vect{x}^{\prime\prime}| > a$ the first integral over $\vect{x}^{\prime\prime}$ vanishes. Moreover, $|e^{i |\vect{k}_i| |\vect{x}^\prime - \vect{x}^{\prime\prime}|}|/|\vect{x}^\prime - \vect{x}^{\prime\prime}| < 1/a$ for $|\vect{x}^\prime - \vect{x}^{\prime\prime}| > a$. Then:
\begin{align}
    \frac{d\sigma_2}{d\Omega} & < \frac{1}{a} \left| \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vect{x}^\prime \int d\vect{x}^{\prime\prime} 
            e^{-i\vect{k}_f \cdot \vect{x}^\prime} V(\vect{x}^\prime) e^{i\vect{k}_i \cdot \vect{x}^{\prime\prime}} V(\vect{x}^{\prime\prime}) \right|^2 \nonumber \\
            & < \frac{1}{a} \left| \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vect{x}^\prime  
            e^{-i\vect{k}_f \cdot \vect{x}^\prime} V(\vect{x}^\prime) \int d\vect{x}^{\prime\prime} e^{i\vect{k}_i \cdot \vect{x}^{\prime\prime}} V(\vect{x}^{\prime\prime}) \right|^2 \nonumber \\
            & < \frac{1}{a} \left| \frac{m_e}{2 \pi \hbar^2} \int d\vect{x}^\prime  
            e^{-i\vect{k}_f \cdot \vect{x}^\prime} V(\vect{x}^\prime) \right|^2 \left| \frac{m_e}{2 \pi \hbar^2} \int d\vect{x}^{\prime\prime} e^{i\vect{k}_i \cdot \vect{x}^{\prime\prime}} V(\vect{x}^{\prime\prime}) \right|^2 \nonumber \\
            & < \frac{1}{a} \left( \frac{d\sigma_1}{d\Omega}\right)^2
\end{align}
where $d\sigma_1/d\Omega$ is the differential cross-section for a single elastic scattering event[^conj].

## The effect of lattice waves on ultrafast electron scattering{#sec:scattering-lattice-waves}

In this section, the effect of lattice waves on ultrafast electron diffraction measurements will be considered, ending in the derivation of what is known as *diffuse scattering*. 

This section will only consider single-scattering events. A discussion of dynamical diffuse scattering is not necessary to understand the experiments presented in this dissertation, as the samples are not thick enough to display a key signature of multiple diffuse scattering events known as *Kikuchi lines*[@Fultz2013KikuchiLines]. A complete discussion of all orders of *thermal* diffuse scattering is presented in Wang [@Wang1995]. Time-resolved multi-phonon diffuse scattering has also been discussed elsewhere[@Zacharias2021a;@Zacharias2021b].

Consider the vector $\vect{r}_{m,s}$ to be the position of atom $s$ in the unit cell $m$. In this scheme, the indices $s$ run over the size of the unit cell, while the indices $m$ run over the number of unit cells: $1 \leq m \leq N_c$. Due to the presence of lattice waves, the atoms are displaced from their equilibrium positions $\set{\vect{r}_{m,s}}$. Let $\set{\vect{u}_{m,s}}$ be the *displacement vectors* due to lattice waves. Then, the atomic positions can be expressed as $\set{\vect{r}_{m,s} \to \vect{r}_{m,s} + \vect{u}_{m,s}}$. The scattering potential of the crystal (@eq:scattering-potential-crystal) becomes:
\begin{align}
    \ft{V}_c(\vect{q}) & = \mathcal{F}\left[ \sum_m \sum_s V_{a,s}(\vect{x} + \vect{r}_{m,s} + \vect{u}_{m,s}) \right] \nonumber \\
                       & = \sum_m \sum_s f_{e,s}(\vect{q}) e^{-i \vect{q} \cdot \vect{r}_{m,s}} e^{-i \vect{q} \cdot \vect{u}_{m,s}}
    \label{eq:scattering-potential-temp}
\end{align}
Recall from @eq:scattering-intensity that the measurable quantity $|f^{(1)}(\vect{q})|^2$ is proportional to $|\ft{V}(\vect{q})|$:
\begin{align}
    |f^{(1)}(\vect{q})|^2 
        & = \left| -\frac{m_e}{\hbar^2} \ft{V}(\vect{q}) \right|^2 \nonumber \\
        & = \frac{m_e^2}{\hbar^4} \ft{V}(\vect{q}) \ft{V}^{\star}(\vect{q}) \nonumber \\
        & = \frac{m_e^2}{\hbar^4}  
            \left(\sum_m \sum_{s} f_{e,s}(\vect{q}) e^{-i \vect{q} \cdot \vect{r}_{m,s}} e^{-i \vect{q} \cdot \vect{u}_{m,s}} \right)  
            \left(\sum_{m^\prime} \sum_{s^{\prime}} f_{e,s^{\prime}}(\vect{q}) e^{i \vect{q} \cdot \vect{r}_{m^\prime, s^{\prime}}} e^{i \vect{q} \cdot \vect{u}_{m^\prime, s^{\prime}}} \right) \nonumber \\
        & = \frac{m_e^2}{\hbar^4} \sum_{m, m^\prime }\sum_{s, s^{\prime}} f_{e,s}(\vect{q}) f_{e,s^{\prime}}(\vect{q})  
            e^{-i \vect{q} \cdot (\vect{r}_{m,s} - \vect{r}_{m^\prime, s^{\prime}})}   
            e^{-i \vect{q} \cdot \vect{u}_{m,s}} e^{ i \vect{q} \cdot \vect{u}_{m^\prime, s^{\prime}}}
    \label{eq:scattering-average}
\end{align}

The evaluation of the sum in @eq:scattering-average requires some thought. The displacement vectors $\vect{u}_{m,s}$ and $\vect{u}_{m^\prime, s^{\prime}}$ are essentially uncorrelated (for $m \neq m^\prime$) for a large enough crystal. Then, the sum over $m$ and $m^{\prime}$ is equivalent to a thermal average over time. To this end, let us define the average $\langle \cdots \rangle \equiv \frac{1}{N_c}\sum_m \cdots$. @eq:scattering-average can then be expressed as:
$$
    |f^{(1)}(\vect{q})|^2 = \frac{m_e^2}{N_c^2\hbar^4} \sum_{m, m^\prime} \sum_{s,s^{\prime}} f_{e,s}(\vect{q}) f_{e,s^{\prime}}(\vect{q}) e^{-i \vect{q} \cdot (\vect{r}_{m,s} - \vect{r}_{m^\prime, s^{\prime}})}   
            \langle e^{-i \vect{q} \cdot \vect{u}_s} e^{ i \vect{q} \cdot \vect{u}_{s^{\prime}}} \rangle
    \label{eq:scattering-amplitude-average}
$$
where the problem has been reduced to the evaluation of $\langle e^{-i \vect{q} \cdot \vect{u}_s} e^{ i \vect{q} \cdot \vect{u}_{s^{\prime}}} \rangle$. The indices $m$ have been dropped from the displacement vectors as the average is over all values of $m$.

### Quantizing lattice waves

Consider now the description of the displacement vectors as a superposition of lattice waves, or phonons. This is most easily done in the second quantization framework[@Altland2010SecondQuantization;@Giustino2017], which states that atomic displacement can be expressed as the quantum-mechanical operator:
$$
    \hat{\vect{u}}_{m,s} = \sum_{\lambda} \sum_{\set{\vect{k}}} \sqrt{\frac{\hbar}{2 \mu_s N \omega_{\lambda}(\vect{k})}} \left( \hat{a}_{\lambda}(\vect{k})e^{-i\phi_{s,m,\lambda}(\vect{k})} + \hat{a}_{\lambda}^{\dagger}(\vect{k}) e^{i\phi_{s,m,\lambda}(\vect{k})} \right) e^{i \vect{k} \cdot \vect{r}_{m,s}}\vect{e}_{s,\lambda}(\vect{k}) 
$${#eq:scattering-displacement}
where $\set{\lambda}$ label phonon branches, $\mu_s$ is the mass of atom $s$, $N$ is the number of atoms in the crystal, $\omega_{\lambda}(\vect{k})$ is the vibrational frequency of mode $\lambda$ at wavevector $\vect{k}$, $\hat{a}_{\lambda}(\vect{k})$ and $\hat{a}_{\lambda}^{\dagger}(\vect{k})$ are the creation and annihilation operators for the phonon mode $\lambda$, $\phi_{s,m,\lambda}(\vect{k})$ is the phase of a lattice wave, and $\vect{e}_{s,\lambda}(\vect{k})$ is the polarization vector of mode $\lambda$[@Sinha2001]. The expression for $\hat{\vect{u}}_{m,s}$ is the combined effect of all possible phonon modes at the $\vect{r}_{m,s}$ lattice site. The sum $\sum_{\set{\vect{k}}}$ assumes the normalization of @eq:scattering-norm.

The Baker-Campbell-Hausdorff lemma can be used to compute the average $\langle e^{-i \vect{q} \cdot \hat{\vect{u}}_s} e^{ i \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}} \rangle$ [@Campbell1897]. It states that for two operators $\hat{A}$ and $\hat{B}$ with commutator $[\hat{A},\hat{B}]$:
$$
e^{\hat{A}} e^{\hat{B}} = e^{\hat{A} + \hat{B} + \frac{1}{2}[\hat{A},\hat{B}]}
$$
This allows to simplify the average as:
$$
\langle e^{-i \vect{q} \cdot \hat{\vect{u}}_s} e^{ i \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}} \rangle = \langle e^{-i \vect{q} \cdot (\hat{\vect{u}}_s - \hat{\vect{u}}_{s^\prime}) + \frac{1}{2}[\vect{q} \cdot \hat{\vect{u}}_s, \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}]} \rangle
$$

Furthermore, note that $[\hat{a}_{\lambda}(\vect{k}), \hat{a}_{\lambda}^{\dagger}(\vect{k})] = 1$ so that the following simplification is valid:
$$
\langle e^{-i \vect{q} \cdot \hat{\vect{u}}_s} e^{ i \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}} \rangle = \langle e^{-i \vect{q} \cdot (\hat{\vect{u}}_s - \hat{\vect{u}}_{s^\prime})}\rangle \langle e^{\frac{1}{2}[\vect{q} \cdot \hat{\vect{u}}_s, \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}]} \rangle
$$
Finally, for operators $A$ which are a linear combination of position and momentum operators of a harmonic system, $\langle e^A \rangle = e^{\frac{1}{2}\langle A^2 \rangle}$[@Born1941], sometimes known as the Bloch identity. This leads to:
$$
\langle e^{-i \vect{q} \cdot \hat{\vect{u}}_s} e^{ i \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}} \rangle = e^{-\frac{1}{2}\langle (\vect{q} \cdot \hat{\vect{u}}_s)^2\rangle} e^{-\frac{1}{2}\langle (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}})^2\rangle} e^{\langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle}
$$
The terms $e^{-\frac{1}{2}\langle (\vect{q} \cdot \hat{\vect{u}}_s)^2\rangle}$ and $e^{-\frac{1}{2}\langle (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}})^2\rangle}$ are known as the Debye-Waller factors[@Waller1923;@Waller1928], historically defined as:
$$
e^{-\frac{1}{2}\langle (\vect{q} \cdot \hat{\vect{u}}_s)^2\rangle} \equiv e^{-W_s}
$${#eq:scattering-debye-waller}
which means that
$$
\langle e^{-i \vect{q} \cdot \hat{\vect{u}}_s} e^{ i \vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}} \rangle = e^{-W_s} e^{-W_{s^\prime}} e^{\langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle}
$$
For small displacement vectors $\hat{\vect{u}}$, $\vect{q} \cdot \hat{\vect{u}} \leq |\vect{q}| |\hat{\vect{u}}|$ is also small, and so:
\begin{align}
    e^{\langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle} 
        & = 1 + \langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle + \mathcal{O}\left(|\hat{\vect{u}}_s|^2 |\hat{\vect{u}}_{s^{\prime}}|^2 \right) \nonumber \\
        & \approx 1 + \langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle
\end{align}
This approximation limits the final expression to the effects of *one-phonon* scattering. This is a good first approximation for simple crystal structure like graphite and MoS$_2$[@Zacharias2021], but there are reports that compounds with intrinsically-low thermal conductivity -- specifically black Phosphorus -- display a measureable degree of multi-phonon diffuse scattering [@Seiler2021]. Using @eq:scattering-displacement:
\begin{multline}
    (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) = \\
        \frac{\hbar}{2 N} \left( \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})}{\sqrt{\mu_s \omega_{\lambda}(\vect{k})}} \left[ \hat{a}_{\lambda}(\vect{k})e^{-i\phi_{s,m,\lambda}(\vect{k})} + \hat{a}_{\lambda}^{\dagger}(\vect{k}) e^{i\phi_{s,m,\lambda}(\vect{k})} \right] e^{i \vect{k} \cdot \vect{r}_{m,s}} \right) \\
        \left( \sum_{\lambda^\prime} \sum_{\set{\vect{k}^\prime}} \frac{\vect{q} \cdot \vect{e}_{\lambda^\prime,s^\prime}(\vect{k}^\prime)}{\sqrt{\mu_{s^\prime} \omega_{\lambda^\prime}(\vect{k}^\prime)}} \left[ \hat{a}_{\lambda^\prime}(\vect{k}^\prime)e^{-i\phi_{s^\prime,m^\prime,\lambda^\prime}(\vect{k}^\prime)} + \hat{a}_{\lambda^\prime}^{\dagger}(\vect{k}^\prime) e^{i\phi_{s^\prime,m^\prime,\lambda}(\vect{k}^\prime)} \right] e^{i \vect{k}^\prime \cdot \vect{r}_{m^\prime, s^\prime}}\right)
\end{multline}
Equivalently:
\begin{multline}
    \langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle = \\
    \frac{\hbar}{2 N} \sum_{\lambda, \lambda^\prime} \sum_{\set{\vect{k}}, \set{\vect{k}^\prime}} \frac{\left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}) \right) \left(\vect{q} \cdot \vect{e}_{\lambda,s^\prime}(\vect{k})\right)}{\sqrt{\mu_s \mu_{s^\prime} \omega_{\lambda}(\vect{k})\omega_{\lambda^\prime}(\vect{k}^\prime)}} e^{i \vect{k} \cdot \vect{r}_{m,s}} e^{i \vect{k}^\prime \cdot \vect{r}_{m^\prime, s^\prime}} \\
    \left\langle \left[ \hat{a}_{\lambda}(\vect{k})e^{-i\phi_{s,m,\lambda}(\vect{k})} + \hat{a}_{\lambda}^{\dagger}(\vect{k}) e^{i\phi_{s,m,\lambda}(\vect{k})} \right] \left[ \hat{a}_{\lambda^\prime}(\vect{k}^\prime)e^{-i\phi_{s^\prime,m^\prime,\lambda^\prime}(\vect{k}^\prime)} + \hat{a}_{\lambda^\prime}^{\dagger}(\vect{k}^\prime) e^{i\phi_{s^\prime,m^\prime,\lambda}(\vect{k}^\prime)} \right] \right\rangle
\end{multline}

Since the phases $\phi_{s,m,\lambda}(\vect{k})$ are not correlated across unit cells, the cross terms vanish:
\begin{align}
    \langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle 
    = & \frac{\hbar}{2 N} \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{\left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}) \right) \left(\vect{q} \cdot \vect{e}_{\lambda,s^\prime}(\vect{k})\right)}{\omega_{\lambda}(\vect{k})\sqrt{\mu_s \mu_{s^\prime}}} e^{i \vect{k} \cdot \vect{r}_{m,s}} e^{i \vect{k}^\prime \cdot \vect{r}_{m^\prime, s^\prime}} \nonumber \\
      & \left[ \hat{a}_{\lambda}(\vect{k}) \hat{a}_{\lambda}(\vect{k}) 
         + \hat{a}_{\lambda}(\vect{k}) \hat{a}^{\dagger}_{\lambda}(\vect{k}) 
         + \hat{a}^{\dagger}_{\lambda}(\vect{k}) \hat{a}_{\lambda}(\vect{k}) 
         + \hat{a}^{\dagger}_{\lambda}(\vect{k}) \hat{a}^{\dagger}_{\lambda}(\vect{k}) \right] \nonumber \\
    = & \frac{\hbar}{2 N} \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{\left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}) \right) \left(\vect{q} \cdot \vect{e}_{\lambda,s^\prime}(\vect{k})\right)}{\omega_{\lambda}(\vect{k})\sqrt{\mu_s \mu_{s^\prime}}} e^{i \vect{k} \cdot \vect{r}_{m,s}} e^{i \vect{k}^\prime \cdot \vect{r}_{m^\prime, s^\prime}} \left[ 2 \hat{n}_{\lambda}(\vect{k}) + 1\right]
\end{align}
where $\hat{n}_{\lambda}(\vect{k}) \equiv \hat{a}_{\lambda}(\vect{k}) \hat{a}^{\dagger}_{\lambda}(\vect{k}) = \hat{a}^{\dagger}_{\lambda}(\vect{k}) \hat{a}_{\lambda}(\vect{k})  - 1$ is the excitation number operator. Simplifying further:
$$
\langle (\vect{q} \cdot \hat{\vect{u}}_s) ~ (\vect{q} \cdot \hat{\vect{u}}_{s^{\prime}}) \rangle = \frac{\hbar}{N} \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{\hat{n}_{\lambda}(\vect{k}) + 1/2}{\omega_{\lambda}(\vect{k})}\frac{\left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}) \right) \left(\vect{q} \cdot \vect{e}_{\lambda,s^\prime}(\vect{k})\right)}{\sqrt{\mu_s \mu_{s^\prime}}} e^{i \vect{k} \cdot \vect{r}_{m,s}} e^{i \vect{k}^\prime \cdot \vect{r}_{m^\prime, s^\prime}}
$$

### Scattering amplitude

Using the calculation of the previous section, the scattering amplitude can be computed. Since this calculation holds for a prepared initial state $\ket{\vect{k}_i}$ and arbitrary final state $\ket{\vect{k}_f}$ (where $\vect{q} = \vect{k}_f - \vect{k}_i$), the quantities $\hat{\vect{u}}$ and $\hat{n}$ are no longer operators, but observables $\vect{u}$ and $n$ respectively. @eq:scattering-amplitude-average can then be expressed as:
\begin{align}
    |f^{(1)}(\vect{q})|^2 
    = & \frac{m_e^2}{N_c^2 \hbar^4} \sum_{m, m^\prime} \sum_{s,s^{\prime}} f_{e,s}(\vect{q}) f_{e,s^{\prime}}(\vect{q}) e^{-i \vect{q} \cdot (\vect{r}_{m,s} - \vect{r}_{m^\prime, s^{\prime}})} e^{-W_s} e^{-W_{s^\prime}}\left[ 1 + \langle (\vect{q} \cdot \vect{u}_s) ~ (\vect{q} \cdot \vect{u}_{s^{\prime}}) \rangle  \right] \nonumber \\
    = & \frac{m_e^2}{N_c^2 \hbar^4} \left| \sum_m \sum_s f_{e,s}(\vect{q}) e^{-W_s} e^{-i \vect{q} \cdot \vect{r}_{m,s}} \right|^2 \nonumber \\
    + & \frac{m_e^2}{N N_c^2 \hbar^3} \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{n_{\lambda}(\vect{k}) + 1/2}{\omega_{\lambda}(\vect{k})} 
        \left| \sum_m \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})\right) e^{-i (\vect{q} - \vect{k}) \cdot \vect{r}_{m,s}} \right|^2
\end{align}

It is now convenient to express the atomic positions $\vect{r}_{m,s}=\vect{R}_m + \vect{x}_s$, where $\vect{R}_m$ is the absolute position of unit cell $m$, and $\vect{x}_s$ is the position of atom $s$ with respect to the unit cell origin. The above equation becomes:
\begin{align}
    |f^{(1)}(\vect{q})|^2 = & \frac{m_e^2}{N_c^2 \hbar^4} \left| \sum_m \sum_s f_{e,s}(\vect{q}) e^{-W_s} e^{-i \vect{q} \cdot \vect{R}_{m}} e^{-i \vect{q} \cdot \vect{x}_{s}} \right|^2 \nonumber \\
                         + & \frac{m_e^2}{N N_c^2 \hbar^3} \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{n_{\lambda}(\vect{k}) + 1/2}{\omega_{\lambda}(\vect{k})} 
                            \left| \sum_m \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})\right) e^{-i (\vect{q} - \vect{k}) \cdot \vect{R}_{m}} e^{-i (\vect{q} - \vect{k}) \cdot \vect{x}_{s}} \right|^2
\end{align}
The Fourier transform identity can be used (@eq:scattering-discrete-fourier):
$$
    \sum_{m=1}^{N_c} e^{-i \vect{q} \cdot \vect{R}_m} \xrightarrow[]{N_c \to \infty} N_c \sum_{\set{\vect{H}}} \delta(\vect{q} - \vect{H})
$$
The equation for the scattering amplitude becomes:
\begin{align}
    |f^{(1)}(\vect{q})|^2 = & \frac{m_e^2}{\hbar^4} \left| \sum_{\set{\vect{H}}} \sum_s f_{e,s}(\vect{q}) e^{-W_s} e^{-i \vect{q} \cdot \vect{x}_{s}} \delta(\vect{q} - \vect{H}) \right|^2 \nonumber \\
                         + & \frac{m_e^2}{N N_c^2 \hbar^3} \sum_{\lambda} \sum_{\set{\vect{k}}} \frac{n_{\lambda}(\vect{k}) + 1/2}{\omega_{\lambda}(\vect{k})} 
                            \left| \sum_m \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})\right) e^{-i (\vect{q} - \vect{k}) \cdot \vect{R}_{m}} e^{-i (\vect{q} - \vect{k}) \cdot \vect{x}_{s}} \right|^2
\end{align}

The sum over $\vect{k}$ in the second term can be simplified in a similar way. For any $\vect{q}$, there is one reciprocal point which is closest, $\vect{H}_{\vect{q}}$. Define $\vect{k}_0$ to be such that $\vect{q} = \vect{H}_{\vect{q}} + \vect{k}_0$, i.e. $\vect{k}_0$ is constrained to lie in the first Brillouin zone. Then:
$$
    \sum_m e^{-i (\vect{q} - \vect{k}) \cdot \vect{r}_{m,s}} = N_c \delta(\vect{k}-\vect{k}_0) 
$$
The sum $\sum_{\set{\vect{H}}}$ is implicitly contained in the constraints on the values of $\vect{k}_0$. It follows that:
\begin{align}
    |f^{(1)}(\vect{q})|^2 = & \frac{m_e^2}{\hbar^4} \left| \sum_{\set{\vect{H}}} \sum_s f_{e,s}(\vect{q}) e^{-W_s} e^{-i \vect{q} \cdot \vect{x}_{s}} \delta(\vect{q} - \vect{H}) \right|^2 \nonumber \\
                         + &  \frac{m_e^2}{N \hbar^3} \sum_{\lambda} \frac{n_{\lambda}(\vect{k}_0) + 1/2}{\omega_{\lambda}(\vect{k}_0)} 
                            \left| \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}_0)\right) e^{-i \vect{H}_{\vect{q}} \cdot \vect{x}_{s}} \right|^2 
\end{align}

The phase factor $e^{-i \vect{H}_{\vect{q}} \cdot \vect{x}_{s}}$ is a matter of convention[@Xu2005]. In the expression for the displacement vectors in @eq:scattering-displacement, it was assumed that the phonon eigenvectors were not periodic in general, that is:
$$
    \vect{e}_{\lambda,s}(\vect{k} + \vect{k}^\prime) = \vect{e}_{\lambda,s}(\vect{k}) e^{-i \vect{k{^\prime}\cdot \vect{x}_s}}
$$
Henceforth, the phonon eigenvectors are assumed to be periodic, that is, $\vect{e}_{\lambda,s}(\vect{k} + \vect{H}) \equiv \vect{e}_{\lambda,s}(\vect{k})$ for a reciprocal point $\vect{H}$. Then $\vect{e}_{\lambda,s}(\vect{k}_0) e^{-i \vect{H}_{\vect{q}} \cdot \vect{x}_{s}} = \vect{e}_{\lambda,s}(\vect{k}_0 - \vect{H}_{\vect{q}}) = \vect{e}_{\lambda,s}(\vect{k}_0)$. It follows that:
\begin{align}
    |f^{(1)}(\vect{q})|^2 = & \frac{m_e^2}{\hbar^4} \left| \sum_{\set{\vect{H}}} \sum_s f_{e,s}(\vect{q}) e^{-W_s} e^{-i \vect{q} \cdot \vect{x}_{s}} \delta(\vect{q} - \vect{H}) \right|^2 \nonumber \\
                         + &  \frac{m_e^2}{N \hbar^3} \sum_{\lambda} \frac{n_{\lambda}(\vect{k}_0) + 1/2}{\omega_{\lambda}(\vect{k}_0)} 
                            \left| \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k}_0)\right) \right|^2 
    \label{eq:scattering-amplitude-reduced}
\end{align}
A visual representation of the relationship between vectors $\vect{H}_{\vect{q}}$, $\vect{q}$, and $\vect{k_0}$ is shown in @fig:scattering-vector-geometry.

```{.matplotlib #fig:scattering-vector-geometry file="figures/scattering/vector-geometry.py" caption="Geometrical relationship between the scattering vector $\vect{q}$, the reciprocal point closest to $\vect{q}$, $\vect{H}_{\vect{q}}$, and wavevector $\vect{k}_0$ for a hypothetical cubic crystal. The in-plane section of the Brillouin, where $\vect{k}_0$ is confined, is shown as well."}
```

Combining @eq:scattering-amplitude-reduced and @eq:scattering-intensity, the measured intensity is therefore:
$$
    I(\vect{q}) = I_0(\vect{q}) + I_1(\vect{q})
$$
where
$$
    I_0(\vect{q}) = \frac{m_e^2}{r^2 \hbar^4} \left|\sum_{\set{\vect{H}}} \sum_s f_{e,s}(\vect{q}) e^{-W_s} e^{-i \vect{q} \cdot \vect{x}_{s}} \delta(\vect{q} - \vect{H}) \right|^2
$${#eq:scattering-diffracted-intensity-finite-temp}
is the diffracted intensity. The diffracted intensity at finite temperature is equivalent to @eq:scattering-diffracted-intensity-zero-temp, with the substitution $f_{e,s}(\vect{q}) \to f_{e,s}(\vect{q}) e^{-W_s}$. The physical meaning of this substitution is that atomic vibrations due to temperature decrease the periodicity of the lattice, which results in a suppression of the atomic form factor in reciprocal space. The other term, $I_1(\vect{q})$, is known as *first order diffuse scattering*:
$$
    I_1(\vect{q}) = \frac{m_e^2}{r^2 N \hbar^3} \sum_{\lambda} \frac{n_{\lambda}(\vect{k}) + 1/2}{\omega_{\lambda}(\vect{k})} 
                            \left| \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})\right)\right|^2
$${#eq:scattering-diffuse-intensity}

### Diffuse scattering {#sec:diffuse-scattering}

The diffuse intensity can be expressed as:
$$
    I_1(\vect{q}) = I_e \sum_{\lambda} \frac{n_{\lambda}(\vect{k}) + 1/2}{\omega_{\lambda}(\vect{k})} |F_{1\lambda}(\vect{q})|^2
$${#eq:scattering-diffuse-intensity}
where $F_{1\lambda}(\vect{q})$ is known as the *one-phonon structure factor*:
$$
|F_{1\lambda}(\vect{q})|^2 = \left| \sum_s \frac{f_{e,s}(\vect{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})\right)\right|^2
$${#eq:scattering-one-phonon-structure-factor}
named thus in contrast to the static structure factor of @eq:scattering-potential-crystal. 

A few clarifications can be made about diffuse intensity. The diffuse intensity at any scattering vector $\vect{q}$ involves the contribution of all phonon modes $\lambda$. The contribution of each mode can be conceptually separated into two parts. 

The first part of each term in the sum of @eq:scattering-diffuse-intensity, proportional to $(n_{\lambda}(\vect{k}) + 1/2)(\omega_{\lambda}(\vect{k}))$, represents the vibrational amplitude of each mode by analogy with the expression of @eq:scattering-displacement. A higher population ($\uparrow n_{\lambda}$) results in a larger vibrational amplitude because the displacement of atoms is linear in the number of phonons that participate. A lower vibrational frequency ($\downarrow \omega_{\lambda}$) implies a smaller restoring force (in the harmonic oscillator sense), which also intuitively results in a wider vibrational amplitude.

The second part of the summation terms are the one-phonon structure factors. These factors are a geometrical weight ($|F_{1\lambda}(\vect{q})|^2$) which determines if the atomic motion associated with a phonon mode can be captured on the detector. The most important terms to consider are terms of the form $\set{\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})}$. For a phonon polarization which is parallel to the propagation of the scattering electrons, the projection of the polarization onto the detector plane is $0$, and hence the associated diffuse intensity will not contribute. An example of one-phonon structure factor for the in-plane longitudinal acoustic mode of graphite is shown in @fig:scattering-oneph-example. Being a longitudinal mode, the polarization of this mode is in the direction of propagation, which means that terms $\vect{q} \cdot \vect{e}_{\lambda,s}(\vect{k})$ should have a radial character in reciprocal space.

```{.matplotlib #fig:scattering-oneph-example file="figures/scattering/oneph.py" caption="One-phonon structure factor $|F_{1\lambda}(\vect{q})|^2$ for the in-plane longitudinal acoustic mode of graphite. The range of scattering vectors $\vect{q}$ corresponds to the geometry of the instrument described in @sec:experimental_setup."}
```

Diffuse scattering and the effect of one-phonon structure factors is further explored in @sec:graphite.

## Conclusion

In this chapter, the theory of ultrafast electron scattering measurements were presented. First, the scattering of electrons by any potential was considered. This was then applied to the case of elastic scattering of electrons by atoms, and then periodic arrangements of atoms at zero temperature. Finally, the full quantum derivation of ultrafast electron diffuse scattering was presented.

The key takeaway from this chapter is that the Debye-Waller effect and diffuse scattering are two faces of the same coin. Both arise from the same physical phenomenon: atomic vibrations expressed as a superposition of lattice waves. In @sec:graphite, ultrafast measurements in a prototypical benchmark system will show the full power of diffuse scattering measurements, while the link between the Debye-Waller effect and diffuse scattering will play a large role in @sec:snse.

[^tmatrix]: Interested readers are encouraged to peruse Chapter 6 of Sakurai and Napolitano [@Sakurai2014].
[^conj]: Note that because $V(\vect{x})$ is real, $\int d\vect{x}^\prime e^{-i\vect{k} \cdot \vect{x}^\prime} V(\vect{x}) = \int d\vect{x}^\prime e^{i\vect{k} \cdot \vect{x}^\prime} V(\vect{x})$.

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
