
# Electron scattering

The scattering of light or particles by a sample is used by a large class of experimental techniques, dating back a hundred years. Scattering can be broadly defined as the modification of an incoming wave by a potential into outgoing wave, a process which imprints the outgoing wave with some characteristic of the potential. The outgoing wave may lose or gain energy, and its momentum might be changed. When multiple incoming waves are simultaneously used -- forming an incoming wavefront --, the outgoing waves may interfere constructively or destructively. This effect is particularly intense for periodic scattering potentials, for example in crystals.

This chapter will consider the special case of *electron scattering*. In crystals, electrons are scattered by the electrostatic potential due to their charge. The description of X-ray scattering, for example, would be identical provided that the total electron charge-density be considered as the scattering potential, instead of the electrostatic potential.

Specifically, the type of scattering relevant to this work involves the probing of periodic structures, crystals. Thanks to advances in data acquisition and data analysis attributable to the author (see [appendix @sec:appendix]), as well as improvements to instrument stability[@Otto2017], the ultrafast electron scattering instrument used herein can reliably measure the effects of *diffuse scattering*.

In writing this chapter, the author has tried original derivations that emphasize concepts that are important for the remainder of this dissertation. The derivation of ultrafast diffuse scattering intensity is of particular interest, because it is the only full quantum-mechanical treatment relevant to ultrafast electron scattering specifically.

## Elastic electron scattering and the Lippmann-Schwinger formalism

TODO: make it clear that the time-variation of the potential happens on a different scale than the scattering time

Consider an electron wavefunction $\Psi(\vec{x}, t)$. The scattering of $\Psi(\vec{x}, t)$ by an arbitrary potential $V(\vec{x},t)$ is described by the Schrödinger equation:
$$
i \hbar \frac{\partial}{\partial t} \Psi(\vec{x}, t) = \left[ \frac{- \hbar^2}{2 m_e} \nabla^2 + V(\vec{x}, t) \right]\Psi(\vec{x}, t).
$${#eq:scattering-schroedinger}

The scattering electrons have an enormous amount of kinetic energy (\SI{90}{\kilo\electronvolt}). The potential energy from the atomic Coulomb interaction, on the other hand, is much smaller: 
$$
    V = \frac{Z e^2}{4 \pi \epsilon_0 |\vec{x}|}
$$
where $\epsilon_0$ is the vacuum dielectric constant, $Z$ the atomic number, $e$ is the unit of charge. For the simple example of graphite ($Z=12$, $|\vec{x}|\approx \SI{1}{\angstrom}$), the potential energy associated with the Coulomb interaction is less than \SI{10}{\electronvolt}. Therefore, the potential $V(\vec{x}, t)$ shall be treated as a perturbation.

### Electrons propagating in free space

In the case of free space ($V(\vec{x}, t) = 0$), [@eq:scattering-schroedinger] reduces to the following equation:
$$
i \hbar \frac{\partial}{\partial t} \Psi(\vec{x}, t) = \frac{- \hbar^2}{2 m_e} \nabla^2 \Psi(\vec{x}, t).
$${#eq:scattering-free-electron}
It is instructive to consider the energy eigenfunction satisfying [@eq:scattering-free-electron]. Suppose that we can label such solutions as $\set{(\Psi_a(\vec{x}, t), \omega_a)}$, where the associated energy eigenvalues are $\omega_a \equiv E_a/\hbar$. We could then factorize the energy eigenfunction as:
$$
\Psi_a(\vec{x}, t) = u_a(\vec{x}) e^{-i \omega_a t}
$${#eq:scattering-free-space}
where $u_a(\vec{x})$ solves the time-independent Schrödinger equation[@Schroedinger1926]:
$$
\hbar \omega_a u_a(\vec{x}) = \frac{- \hbar^2}{2 m_e} \nabla^2 u_a(\vec{x}).
$${#eq:scattering-stationary-schroedinger}
We note that [@eq:scattering-stationary-schroedinger] is an instance of the Helmholtz equation, which can be re-written as:
$$ 
\left[ \nabla^2 + k_a^2\right] u_a(\vec{x}) = 0
$${#eq:scattering-helmholtz}
where $k_a^2 \equiv 2 m_e E_a/\hbar^2$. We may use physical reasoning and the classical result of the Helmholtz equation to name $k_a$ the **wavevector**, related to momentum as $k_a^2 \equiv \vec{k}_a \cdot \vec{k}_a \equiv \vec{p}_a^2/\hbar^2$. [Eq. @eq:scattering-helmholtz] can be solved using separation of variables, where the solution along linearly-independent directions are considered independently:
$$
u_a(\vec{x})\equiv u_{a,x}(\vec{x} \cdot \hat{\vec{x}}) u_{a,y}(\vec{x} \cdot \hat{\vec{y}}) u_{a,z}(\vec{x} \cdot \hat{\vec{z}}).
$$
[Eq. @eq:scattering-helmholtz] then becomes:
$$
\sum_{i\in\{x,y,z\}}\frac{1}{u_{a,i}(\vec{x} \cdot \hat{\vec{i}})} \frac{d^2 u_{a,i}}{di^2}(\vec{x} \cdot \hat{\vec{i}}) + |\vec{k} \cdot \hat{\vec{i}}| = 0
$$
and the solution can be synthesized as a product of one-dimensional plane wave:
\begin{align}
u_a(\vec{x}) &= \prod_{j \in \{x,y,z\}} A_j e^{i (\vec{k}_a \cdot \hat{\vec{j}}) (\vec{x} \cdot \hat{\vec{j}})} \nonumber\\
             &= A e^{i \vec{k}_a \cdot \vec{x}}
\label{eq:scattering-plane-wave}
\end{align}
for some scalars $\{A_x, A_y, A_z, A\}$. Combining [@eq:scattering-free-space] [@eq:scattering-plane-wave] leads to the energy eigenstate of a free electron propagating in vacuum:
$$
\Psi_a(\vec{x}, t) = A e^{i (\vec{k}_a \cdot \vec{x} - \omega_a t)}
$$
with associated energy eigenvalue $E_a = \hbar^2 \vec{k}_a^2 / 2 m_e=\hbar \omega_a$.

### The Lippmann-Schwinger equation

The problem of scattering with a non-trivial potential will now be considered. This is a fundamental problem in quantum mechanics; approximations appropriate for the research presented in this dissertation will be made to simplify and clarify the presentation. In particular, the following derivations assume that the scattering potential $V(\vec{x})$ is *local*, that is:
$$
    \bra{\vec{x}^\prime} V \ket{\vec{x}^{\prime \prime}} = V(\vec{x}^\prime) \delta(\vec{x}^\prime - \vec{x}^{\prime \prime})
$$
Further simplifications can be made allowed momentum states are assumed to be defined in a (large) cube of size-length $L$, such that
$$
    \braket{\vec{x} | \vec{k}} = \frac{1}{L^{3/2}} e^{i \vec{k} \cdot \vec{x}}
$${#eq:scattering-norm}
and $\vec{k}$ takes discrete values. The calculations will become physical after taking the limit $L \to \infty$.
In this approximation, the scattered wave $\bra{\vec{x} | \Psi}$ is given by the *Lippmann-Schwinger* equation[@Lippmann1950]:
$$
    \braket{\vec{x} | \Psi} = \braket{\vec{x} | \vec{k}_i} - \frac{2 m_e}{\hbar^2} \int d^3 x^\prime \frac{e^{ik|\vec{x}-\vec{x}^\prime|}}{4 \pi |\vec{x}-\vec{x}^\prime|} V(\vec{x}^\prime) \braket{\vec{x}^\prime | \Psi}
$${#eq:scattering-lippmann-schwinger}
Provided that the scattered wavefunction $\braket{\vec{x} | \Psi}$ is measured far from where the scattering potential is localized:
$$
    e^{i|\vec{k}_i||\vec{x}-\vec{x}^\prime|} \approx e^{ikr}e^{-i\vec{k} \cdot \vec{x}^\prime}
$$
where $r \equiv |\vec{x}-\vec{x}^\prime|$ and $k \equiv |\vec{k}_i|$. @eq:scattering-lippmann-schwinger can then be simplified as:
$$
    \braket{\vec{x} | \Psi} = \braket{\vec{x} | \vec{k}_i} - \frac{m_e}{2 \pi \hbar^2} \frac{e^{ikr}}{r}\int d^3 x^\prime e^{-i \vec{k}_f \cdot \vec{x}^\prime} V(\vec{x}^\prime) \braket{\vec{x}^\prime | \Psi}
$${#eq:scattering-lippmann-schwinger-rad}
The canonical description of the Lippmann-Schwinger equation is usually reduced to:
$$
    \braket{\vec{x} | \Psi} = \frac{1}{L^{3/2}} \left[ e^{i \vec{k}_i \cdot \vec{x}} + \frac{e^{ikr}}{r} f(\vec{k}_f, \vec{k}_i)\right]
$${#eq:scattering-lippmann-schwinger-general}
where
\begin{align}
    f(\vec{k}_f, \vec{k}_i) & \equiv -\frac{m_e L^3}{2 \pi \hbar^2} \int d\vec{x}^\prime \frac{e^{-i \vec{k}_f \cdot \vec{x}^\prime}}{L^{3/2}} V(\vec{x}^\prime) \braket{\vec{x}^\prime | \Psi} \nonumber \\
                            & =  -\frac{m_e L^3}{2 \pi \hbar^2} \bra{\vec{k}_f}V\ket{\Psi}
    \label{eq:scattering-amplitude}
\end{align}
$f(\vec{k}_f, \vec{k})$ is called the *scattering amplitude*. The form of @eq:scattering-lippmann-schwinger-general complies with intuition: the incoming state $\ket{\vec{k}_i}$ wave is composed of an unscattered part ($\braket{\vec{x} | \vec{k}_i}$) as well as an outgoing spherical wave with amplitude $f(\vec{k}_f, \vec{k}_i)$.

The calculation of the scattered wavefunction $\braket{\vec{x}|\Psi}$ has been reduced to the calculation of $\bra{\vec{k}_f}V\ket{\Psi}$. The derivation of an expression for this is beyond the scope of this work, and the final result is stated[^tmatrix]:
$$
    \bra{\vec{k}_f}V\ket{\Psi} = \bra{\vec{k}_f} \left[ \sum_{j=0}^{\infty} V \left( \frac{1}{E_i - H_0 + i \epsilon} V\right)^j \right] \ket{\vec{k}_i}
$${#eq:scattering-potential-decomp}
where $H_0$ is the free-space Hamiltonian with eigenvalue $E_i=\hbar^2 \vec{k}_i^2/2 m_e$, and $\epsilon$ is a vanishingly small real number. In particular, each term with index $j$ in the sum of @eq:scattering-potential-decomp corresponds to the electron scattering $j$ times[@Feynman1965].

### Measuring the scattered wavefunction

Electron cameras measure the intensity of the wavefunction. In the case of bright-field electron microscopy, where the wavefunction is imaged directly, the intensity on the detector is given by:
$$
    I(\vec{x}) \equiv \left| \braket{\vec{x} | \Psi} \right|^2
$$
In order to sample the scattered wavefunction in reciprocal space, an electron lens can be used to focus the scattered electrons onto the detector. Given that electrons are prepared with definite momentum $\vec{k}_i$, it is trivial to ignore the unscattered part of the wavefunction -- the first term in @eq:scattering-lippmann-schwinger-general -- which is found only at $|\vec{q}| = 0$. Therefore, the intensity away from $|\vec{q}| = 0$ is related only to the scattered wavefunction, $\Psi_s$:
\begin{align}
I(\vec{q}) 
    & = \left| \braket{\vec{x} | \Psi_s} \right|^2 \nonumber \\
    & = \left| \frac{e^{ikr}}{r} f(\vec{k}_f, \vec{k}_i) \right|^2 \nonumber \\
    & = \frac{1}{r^2} |f(\vec{k}_f, \vec{k}_i)|^2
\label{eq:scattering-intensity}
\end{align}
Therefore, the diffracted intensity $I(\vec{q})$ is proportional to the square of the scattering amplitude. 

Note that the factor of $1/r^2$ is generally ignored[@Fultz2002r2]. For the instrument configuration presented in @sec:experimental_setup, $1/r^2$ varies from \SIrange{15.92}{16}{\per\square\meter}, from the corner to the center of the detector respectively. While not an insignificant variation, the experiments presented herein generally measure the relative change in intensity, in which case factors are not important.

## Elastic scattering in a crystal

In this section, the consequences of an electron scattering *once* in crystalline solid will be explored. In this approximation, historically called the *first Born approximation* [@Born1926], only the first term in @eq:scattering-potential-decomp ($j=0$) is considered:
\begin{align}
    f^{(1)}(\vec{k}_f, \vec{k}_i) & = -\frac{m_e L^3}{2 \pi \hbar^2} \bra{\vec{k}_f} V \ket{\vec{k}_i} \nonumber \\
                                  & = -\frac{m_e L^3}{2 \pi \hbar^2} \int d\vec{x}^\prime \bra{\vec{k}_f} V(\vec{x}^\prime) \ket{\vec{x}^\prime}\braket{\vec{x}^\prime | \vec{k}_i} \nonumber \\
                                  & = -\frac{m_e L^3}{2 \pi \hbar^2} \int d\vec{x}^\prime \frac{e^{i(\vec{k}_i - \vec{k}_f)\cdot \vec{x}^\prime}}{L^3} V(\vec{x}^\prime) \nonumber \\
                                  & = -\frac{m_e}{2 \pi \hbar^2} \int d\vec{x}^\prime e^{i(\vec{k}_i - \vec{k}_f)\cdot \vec{x}^\prime} V(\vec{x}^\prime)
    \label{eq:scattering-first-born-approx}
\end{align}
where the normalization of @eq:scattering-norm was used. The reader may recognize that the scattering amplitude $f^{(1)}(\vec{k}_f, \vec{k}_i)$ is proportional to the Fourier transform of the scattering potential with respect to $\vec{k}_i - \vec{k}_f \equiv \vec{q}$, the *scattering vector*. If the Fourier trnasform functional operator is defined as:
$$
    \mathcal{F}\left[ f(\vec{x}) \right] \equiv \hat{f}(\vec{q}) = \frac{1}{2 \pi} \int d\vec{x}^\prime e^{i \vec{q} \cdot \vec{x}^\prime}f(\vec{x}^\prime).
$$
then @eq:scattering-first-born-approx can be simplified to:
$$
    f^{(1)}(\vec{q}=\vec{k}_f - \vec{k}_i) = -\frac{m_e}{\hbar^2} \hat{V}(\vec{q})
$${#eq:scattering-amplitude-q}

### Scattering potential of a single atom{#sec:affe}

The scattering potential of a single atom is given by:
$$
    V_a(\vec{x}) = -\frac{Z e^2}{|\vec{x}|} + \sum_{i=1}^{Z} \frac{e^2}{|\vec{x} - \vec{x}_i|}
$${#eq:scattering-atom-potential}
where Z is the atomic weight, and $\vec{x}_i$ is the position of the $i^{\text{th}}$ electron. The potential in @eq:scattering-atom-potential can be calculated from first principles, using relativistic Hartree-Fock calculations[@Fischer1977; @Kirkland2010] to get the real-space electron density (i.e. determining the possible set $\set{\vec{x}_i}$). This is beyond the scope of the present work, and the final result is used here. When discussing electron scattering, the scattering amplitude for a single atom is usually called the *atomic form factors for electron scattering*. To a reasonable degree of accuracy, the atomic form factor for electrons are spherically symmetric. They can be parametrized of an atom can be parametrized as:
$$
    f_e(\vec{q}) = \sum_{i=1}^3 \frac{a_i}{|\vec{q}|^2 + b_i} + c_i e^{-d_i|\vec{q}|^2}
$${#eq:scattering-affe}
where the constants $\set{ a_i, b_i, c_i, d_i }$ are element-specific fitting parameters which are tabulated in Kirkland [@Kirkland2010]. The associated real-space potential can be calculated via @eq:scattering-amplitude-q. Examples of $f_e(\vec{q})$ and associated $V_a(\vec{x})$ are shown in @fig:scattering-potential for a few elements.

```{.matplotlib #fig:scattering-potential file="figures/scattering/scatt-potential.py" caption="Demonstration of the electrostatic potential of atoms, which scatters electrons. **a)** Radial view of the electrostatic potential in real-space **b)** Radial view of the electrostatic potential in reciprocal space, also known as the *atomic form factor*."}
```

The contribution of individual electronic orbitals to the atomic form factor for electrons[@Zheng2009] will be discussed in @sec:snse. 

### Scattering potential of a crystal

With the knowledge of the scattering potential of a single atom, the scattering potential of a crystalline lattice can be calculated:
$$
    V_c(\vec{x}) = \sum_i V_a(\vec{x} - \vec{r}_i)
$$
where  the sum index $i$ runs over atoms in the crystal with positions $\vec{r}_i$, and $V_a$ is given by @eq:scattering-atom-potential. Note that for any function $h(\vec{x})$: 
\begin{align}
    \mathcal{F}\left[ h(\vec{x} + \vec{y}) \right] & = \frac{1}{2\pi} \int d\vec{x}^\prime e^{i\vec{q}\cdot(\vec{x}^\prime + \vec{y})} h(\vec{x}^\prime) \nonumber \\
                                                   & = \frac{e^{i \vec{q} \cdot \vec{y}}}{2 \pi} \int d\vec{x}^\prime e^{i\vec{q} \vec{x}^\prime} h(\vec{x}^\prime) \nonumber \\
                                                   & = e^{i \vec{q} \cdot \vec{y}} \mathcal{F}\left[ h(\vec{x}) \right]
\end{align}
where $\vec{y}$ is some arbitrary translation vector. Therefore, the Fourier transform of the scattering potential of the entire crystal is related to the Fourier transform  the potential its constituent atoms (@eq:scattering-affe) as:
\begin{align}
    \hat{V}_c(\vec{q}) & = \mathcal{F}\left[ \sum_i V_a(\vec{x} - \vec{r}_i) \right] \nonumber \\
                       & = \sum_i \mathcal{F}\left[ V_a(\vec{x} - \vec{r}_i) \right] \nonumber \\
                       & = \sum_i \mathcal{F}\left[ V_a(\vec{x}) \right] e^{-i \vec{q} \cdot \vec{r}_i} \nonumber \\
                       & = \sum_i f_{e,i}(\vec{q}) e^{-i \vec{q} \cdot \vec{r}_i}
    \label{eq:scattering-potential-crystal}
\end{align}

@eq:scattering-potential-crystal has historically been called the *static structure factor*, named thus in contrast with the *dynamic structure factor* discussed in @sec:diffuse-scattering. A visual representation of the scattering potential of a crystal is helpful. Consider for moment $\alpha$-polonium ($\alpha$-Pu), one of the simplest crystal structures[@Curie1898], which crystals consists in a Pu atom at the every vertex of a cube with side-length of \SI{3.63}{\angstrom}. The calculated electrostatic potential of this arrangement along the base of a cube is shown in @fig:scattering-polonium-example a). The lattice vectors $\vec{a}_1$ and $\vec{a}_2$ are indicated, with $\vec{a}_3$ being aligned out of the page. The periodic nature of this scattering potential is demonstrated by calculating the resulting scattering amplitude $f(\vec{q})$ from @eq:scattering-amplitude-q, which is shown in @fig:scattering-polonium-example b). The periodicity in spatial-frequency-space, also called *reciprocal space*, is evident, and forms a *reciprocal lattice*, formally defined in the next section.

From @eq:scattering-intensity and @eq:scattering-amplitude-q, the measured diffracted intensity is:
\begin{align}
    I(\vec{q}) & = \frac{m_e^2}{4 \pi^2 \hbar^4 r^2} \left| \sum_i f_{e,i}(\vec{q}) e^{-i \vec{q} \cdot \vec{r}_i} \right|^2 \nonumber \\
               & = \frac{N_c^2 m_e^2}{4 \pi^2 \hbar^4 r^2} \left| \sum_s f_{e,s}(\vec{q}) e^{-i \vec{q} \cdot \vec{r}_s} \right|^2
    \label{eq:scattering-diffracted-intensity-zero-temp}
\end{align}
which is the standard result for the diffracted intensity being proportional to the square of the static structure factor[@Warren1990intensity; @Kittel1996intensity; @Fultz2002intensity; @Kirkland2010intensity]. Note that the sum $\sum_i$ has been replaced with the sum over the unit cell, $N_c \sum_s$.

```{.matplotlib #fig:scattering-polonium-example file="figures/scattering/polonium.py" caption="Calculated scattering potential and associated scattering amplitude for $\alpha$-polonium. **a)** Electrostatic potential $V(\vec{x})$ in the $z=0$ plane. The two in-plane lattice vectors $\vec{a}_1$ and $\vec{a}_2$ are shown; lattice vector $\vec{a}_3$ points out of the page. **b)** Scattering amplitude $f(\vec{q})$ associated with the electrostatic potential shown in a). The periodic nature of the potential in real-space creates a structure in reciprocal space called the *reciprocal lattice*."}
``` 

### The reciprocal lattice

The geometry of reciprocal space and the reciprocal lattice are foundational concepts that drive the understanding of ultrafast electron diffraction.

A perfectly periodic structure in real-space that extends to infinity, with associated lattice vectors $\set{ \vec{a}_i }$, possesses a *dual* lattice in reciprocal space. The lattice vectors $\set{ \vec{b}_i }$ for this reciprocal lattice are defined by the relation
$$
    \vec{a}_i \cdot \vec{b}_j = 2 \pi \delta_{ij}
$$
which leads to the following reciprocal lattice vectors:
$$
    \left\{ \vec{b}_i = 2 \pi \frac{\vec{a}_j \times \vec{a}_k}{\vec{a_i} \cdot (\vec{a_j} \times \vec{a}_k)} ~ \middle| ~ (i,j,k) \in C \right\}
$$
where $C$ is the set of cyclic permutations. For the example of $\alpha$-Pu, where
$$
    \vec{a}_1 = 3.63 ~ \hat{\vec{x}}, ~ \vec{a}_2 = 3.63 ~ \hat{\vec{y}}, ~ \vec{a}_3 = 3.63 ~ \hat{\vec{z}}
$${#eq:scattering-polonium-lattice}
the associated reciprocal lattice vectors are
$$
    \vec{b}_1 = \frac{2\pi}{3.63} ~ \hat{\vec{x}}, ~ \vec{b}_2 = \frac{2\pi}{3.63} ~ \hat{\vec{y}}, ~ \vec{b}_3 = \frac{2\pi}{3.63} ~ \hat{\vec{z}}.
$${#eq:scattering-polonium-recip}
The geometry of vectors in @eq:scattering-polonium-lattice and @eq:scattering-polonium-recip are shown in @fig:scattering-polonium-example. 

The position of *reciprocal points* $\vec{H}$ -- the location of the fundamental frequencies of the Fourier transform of @eq:scattering-amplitude-q -- is a linear combination of
vectors reciprocal basis vectors $\set{ \vec{b}_i }$:
$$
    \vec{H} = h ~ \vec{b}_1 + k ~ \vec{b}_2 + l ~ \vec{b}_3
$$
Expressed in the reciprocal basis, reciprocal points are traditionally denoted as $\vec{H} = (hkl)$. The indices $h$, $k$, and $l$ are called *Miller indices*, named for W. H. Miller[@Miller1839]. 

#### Geometrical interpretation of reciprocal space

The reciprocal lattice is the *dual* of real-space lattice. The duality relationship encodes the idea of orthogonality, that is, reciprocal points are supposed to be in some sense perpendicular to the point on the real lattice. To understand what that means, consider the real-space points $\vec{a}_1/h$ and $\vec{a}_2/k$ defined on the real-space lattice, and the reciprocal point $\vec{H} = (hkl)$:
$$
    \left( \frac{\vec{a}_1}{h} - \frac{\vec{a}_2}{k} \right) \cdot \vec{H} = 0
$$
The vector $\vec{a}_1/h - \vec{a}_2/k$ lies in a plane, and the vector $\vec{H}=(hkl)$ is perpendicular to this plane for any $l$. This example illustrates that reciprocal points $(hkl)$ define planes in real-space. TODO: demonstrate that $|\vec{H}| = 1/d_{hkl}$

### Bragg's law

It is possible to *deduce* Bragg's law from @eq:scattering-amplitude-q and our definition of the reciprocal lattice.

By definition, the reciprocal points of the crystal scattering potential, located at the spatial frequencies where the crystal potential strong, form a lattice with basis vectors $\set{ \vec{b}_i }$. Consider an electron with initial wavevector $\vec{k}_i$ that scatters elastically to a final wavevector $\vec{k}_f$. The scattering amplitude for this event, $f(\vec{k}_f, \vec{k}_i)$ is most intense where $\hat{V}(\vec{k}_f - \vec{k}_i)$ is strong; that is, the condition for strong scattering is:
$$
    \vec{k}_f - \vec{k}_i = h ~ \vec{b}_1 + k ~ \vec{b}_2 + l ~ \vec{b}_3 \quad \forall ~ h,k,l \in \mathbb{Z}
$${#eq:scattering-bragg-vector}
This is precisely the *vector* form of Bragg's law [@Warren1990Bragg]. To recover the canonical form of Bragg's law, consider that an electron state with wavevector $\vec{k}$ can be associated with a wavelength of $\lambda = \tfrac{1}{|\vec{k}|}$. Since for elastic scattering, $|\vec{k}_i| = |\vec{k}_f| = \tfrac{1}{\lambda}$, @eq:scattering-bragg-vector becomes:
$$
    \frac{1}{\lambda} \left( \hat{\vec{k}}_f - \hat{\vec{k}}_i\right) = \vec{H}
$${#eq:scattering-bragg-1}
where $\hat{\vec{k}}$ denotes a unit-length vector in the direction of $\vec{k}$. Given that the vectors on both sides of the equation have the same magnitude and direction, the direction of $\vec{H}$ must bisect the angle between $\vec{k}_f$ and $\vec{k}_i$, historically defined as $2\theta$. Taking the amplitude of @eq:scattering-bragg-1:
$$
     \frac{1}{\lambda} \left| \hat{\vec{k}}_f - \hat{\vec{k}}_i\right| = \frac{2 \sin{\theta}}{\lambda} 
$$
and
$$
|\vec{H}| = \frac{1}{d_{hkl}},
$$
which can be combined as
$$
    \frac{2 \sin{\theta}}{\lambda} = \frac{1}{d_{hkl}}.
$${#eq:scattering-bragg-hist}
@eq:scattering-bragg-hist is the historical form of Bragg's law as it relates to polycrystalline diffraction patterns[@Bragg1913]. Note that the vector form of @eq:scattering-bragg-vector is richer than the original form of Bragg's law as it places constraint on the full three-dimensional direction of the scattering vector $\vec{q} = \vec{k}_f - \vec{k}_i$.

### The Ewald sphere

Elastic electron scattering, or *electron diffraction*, can be discussed more concretely. Consider an electron initially propagating in the $\hat{\vec{z}}$ with wavevector $\vec{k}_i$ direction that interacts with a scattering potential $\hat{V}(\vec{q})$, and scatters to a final wavevector $\vec{k}_f$. The elastic scattering condition $|\vec{k}_i| = |\vec{k}_f|$ constrains the observation of $\hat{V}(\vec{q})$ to scattering vectors $\vec{q}$ that lie on a sphere of radius $|\vec{q}|=\tfrac{1}{\lambda}$. This sphere is called the *Ewald sphere*[@Ewald1921].

```{.matplotlib #fig:scattering-ewald-sphere file="figures/scattering/ewald.py" caption="Demonstration of the Ewald sphere, a visual representation of the conservation of energy in diffraction. The Fourier transform of the scattering potential from an abstract cubic lattice of side length \SI{5}{\angstrom}, $\hat{V}(\vec{q})$, is shown in the background, with the associated reciprocal lattice vectors $\set{\vec{b}_i}$. The Ewald sphere of radius $\vec{q}$ is shown for two scatterers: electrons (solid) and hard x-ray (dashed)."}
```

The Ewald sphere is a great mental model of the information contained in diffraction patterns. Because diffracting electrons can only sample scattering vectors on the Ewald sphere, any particular measurement of a scattering potential $V(\vec{x})$ is effectively a two-dimensional *slice* of the three-dimensional Fourier transform of $V(\vec{x})$, $\hat{V}(\vec{q})$. This is represented in @fig:scattering-ewald-sphere. In this figure, the potential $\hat{V}(\vec{q})$ for an idealized simple cubic crystal with side-length \SI{5}{\angstrom} is shown in the plane spanned by $\vec{b}_2$ and $\vec{b}_3$. The Ewald spheres associated with \SI{100}{\kilo\electronvolt} electrons (large $|\vec{q}|$) and \SI{13}{\kilo\electronvolt} x-rays (smaller $|\vec{q}|$) are also shown. This electron energy is typical of the work presented in this thesis, while the x-ray energy an upper bound on the available energies at the Linac Coherent Light Source as of 2021[@Bostedt2013]. The reciprocal points that intersect the Ewald sphere appear in measurements as diffraction peaks, or Bragg peaks. @fig:scattering-ewald-sphere shows the advantage of electron scattering to study two-dimensional materials: given the proper orientation of the electron beam, a large range of wavevectors can be studied in the plane of interest.

## Multiple scattering of electrons {#sec:scattering-multiple}

Electrons interact quite strongly with matter through the Coulomb interaction. TODO: introduce this more. also talk about how this is not the same as two-beam diffraction because UED bunches are not dense; interference from two separate electrons in less likely than double-scattering by the same electron. 

In this section, the scattering of an electron *twice* will be considered. In this case, the second term in @eq:scattering-potential-decomp ($j=1$) is considered:
$$
    f^{(2)}(\vec{k}_f, \vec{k}_i) = -\frac{m_e L^3}{2 \pi \hbar^2} \bra{\vec{k}_f} V \frac{1}{E_i - H_0 + i \epsilon}V \ket{\vec{k}_i}
$$
The calculation of $f^{(2)}(\vec{k}_f, \vec{k}_i)$ involves the insertion of two complete sets of basis states:
\begin{multline}
    \bra{\vec{k}_f} V \frac{1}{E_i - H_0 + i \epsilon}V \ket{\vec{k}_i} = \\ 
        \int d\vec{x}^\prime \int d\vec{x}^{\prime\prime} 
            \braket{\vec{k}_f | \vec{x}^\prime} V(\vec{x}^\prime) \bra{\vec{x}^\prime} \frac{1}{E_i - H_0 + i \epsilon} \ket{\vec{x}^{\prime\prime}}V(\vec{x}^{\prime\prime}) \braket{\vec{x}^{\prime\prime} | \vec{k}_i}
\end{multline}
The evaluation of the term $\bra{\vec{x}^\prime} \frac{1}{E_i - H_0 + i \epsilon} \ket{\vec{x}^{\prime\prime}}$ naturally happens in the derivation of @eq:scattering-lippmann-schwinger, and so the result is simply stated here:
$$
    \bra{\vec{x}^\prime} \frac{1}{E_i - H_0 + i \epsilon} \ket{\vec{x}^{\prime\prime}} = -\frac{m_e}{2 \pi \hbar^2} \frac{e^{i |\vec{k}_i| |\vec{x}^\prime - \vec{x}^{\prime\prime}|}}{|\vec{x}^\prime - \vec{x}^{\prime\prime}|}
$$
It follows that the scattering amplitude for the double-scattering of an electron is given by:
$$
    f^{(2)}(\vec{k}_f, \vec{k}_i) =
        \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vec{x}^\prime \int d\vec{x}^{\prime\prime} 
            e^{-i\vec{k}_f \cdot \vec{x}^\prime} V(\vec{x}^\prime) \left( \frac{e^{i |\vec{k}_i| |\vec{x}^\prime - \vec{x}^{\prime\prime}|}}{|\vec{x}^\prime - \vec{x}^{\prime\prime}|} \right) e^{i\vec{k}_i \cdot \vec{x}^{\prime\prime}} V(\vec{x}^{\prime\prime})
$${#eq:scattering-amplitude-mult}
The structure of @eq:scattering-amplitude-mult informs on the following physical interpretation. Double scattering involves a first scattering at $\vec{x}^{\prime\prime}$ ($V(\vec{x}^{\prime\prime})$), which "radiates" as a spherical wave moving from $\vec{x}^{\prime\prime}$ to $\vec{x}^{\prime}$ ($e^{i |\vec{k}_i| |\vec{x}^\prime - \vec{x}^{\prime\prime}|}/|\vec{x}^\prime - \vec{x}^{\prime\prime}|$), followed by a second scattering at $\vec{x}^\prime$ ($V(\vec{x}^{\prime})$).

TODO: plot a comparison

### Comparing cross-sections

The computation of the differential scattering cross section for two elastic scattering events is given by:
\begin{align}
    \frac{d\sigma_2}{d\Omega} & = |f^{(2)}(\vec{k}_f, \vec{k}_i)|^2 \nonumber \\
             & = \left| \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vec{x}^\prime \int d\vec{x}^{\prime\prime} 
            e^{-i\vec{k}_f \cdot \vec{x}^\prime} V(\vec{x}^\prime) \left( \frac{e^{i |\vec{k}_i| |\vec{x}^\prime - \vec{x}^{\prime\prime}|}}{|\vec{x}^\prime - \vec{x}^{\prime\prime}|} \right) e^{i\vec{k}_i \cdot \vec{x}^{\prime\prime}} V(\vec{x}^{\prime\prime}) \right|^2
    \label{eq:scattering-cross-sec-2}
\end{align}
The evaluation of @eq:scattering-cross-sec-2 for a realistic potential $V(\vec{x})$ is arduous. However, given that this type of scattering is *undesirable*, getting an upper bound on its scattering cross-section is a worthwhile exercise. 

In electron scattering experiments, the elastic scattering cross section is small enough that an electron is unlikely to scatter twice from the same atom. This means that the integral is 0 when $|\vec{x}^\prime - \vec{x}^{\prime\prime}|$ is small. The inner integral over $\vec{x}^{\prime\prime}$ in @eq:scattering-cross-sec-2 can be split:
$$
    \int d\vec{x}^{\prime} \int d\vec{x}^{\prime\prime} \to \int d\vec{x}^{\prime} \left[ \int_{|\vec{x}^\prime - \vec{x}^{\prime\prime}| \leq a} d\vec{x}^{\prime\prime} + \int_{|\vec{x}^\prime - \vec{x}^{\prime\prime}| > a} d\vec{x}^{\prime\prime} \right]
$$
where $a$ is the inter-atomic distance of the crystal. In the approximation described above, where double-scattering is only possible for $|\vec{x}^\prime - \vec{x}^{\prime\prime}| > a$ the first integral over $\vec{x}^{\prime\prime}$ vanishes. Moreover, $|e^{i |\vec{k}_i| |\vec{x}^\prime - \vec{x}^{\prime\prime}|}|/|\vec{x}^\prime - \vec{x}^{\prime\prime}| < 1/a$ for $|\vec{x}^\prime - \vec{x}^{\prime\prime}| > a$. Then:
\begin{align}
    \frac{d\sigma_2}{d\Omega} & < \frac{1}{a} \left| \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vec{x}^\prime \int d\vec{x}^{\prime\prime} 
            e^{-i\vec{k}_f \cdot \vec{x}^\prime} V(\vec{x}^\prime) e^{i\vec{k}_i \cdot \vec{x}^{\prime\prime}} V(\vec{x}^{\prime\prime}) \right|^2 \nonumber \\
            & < \frac{1}{a} \left| \left( \frac{m_e}{2 \pi \hbar^2} \right)^2 \int d\vec{x}^\prime  
            e^{-i\vec{k}_f \cdot \vec{x}^\prime} V(\vec{x}^\prime) \int d\vec{x}^{\prime\prime} e^{i\vec{k}_i \cdot \vec{x}^{\prime\prime}} V(\vec{x}^{\prime\prime}) \right|^2 \nonumber \\
            & < \frac{1}{a} \left| \frac{m_e}{2 \pi \hbar^2} \int d\vec{x}^\prime  
            e^{-i\vec{k}_f \cdot \vec{x}^\prime} V(\vec{x}^\prime) \right|^2 \left| \frac{m_e}{2 \pi \hbar^2} \int d\vec{x}^{\prime\prime} e^{i\vec{k}_i \cdot \vec{x}^{\prime\prime}} V(\vec{x}^{\prime\prime}) \right|^2 \nonumber \\
            & < \frac{1}{a} \left( \frac{d\sigma_1}{d\Omega}\right)^2
\end{align}
where $d\sigma_1/d\Omega$ is the differential cross-section for a single elastic scattering event[^conj].

## The effect of lattice waves on ultrafast electron scattering

TODO: explain that we only consider single-phonon scattering

In this section, the effect of lattice waves on electron diffraction measurements will be considered, ending in the derivation of what is known as *diffuse scattering*. 

It is helpful to distinguish between atoms of separate unit cells. Consider then the position $\vec{r}_{m,s}$ to be the position of atom $s$ in the unit cell $m$. In this scheme, the indices $s$ run over the size of the unit cell (e.g. $1 \leq s \leq 4$ for graphite), while the indices $m$ run over the number of unit cells: $1 \leq m \leq N_c$. 

Due to the presence of lattice waves, the atoms are vibrating about their equilibrium positions $\set{\vec{r}_{m,s}}$. Let $\set{\vec{u}_{m,s}}$ be the *displacement vectors* due to lattice waves. Then, the atomic positions can be expressed as $\set{\vec{r}_{m,s} \to \vec{r}_{m,s} + \vec{u}_{m,s}}$, where $\set{\vec{r}_{m,s}}$ are the atomic positions as defined by the lattice at zero temperature.Then:
\begin{align}
    \hat{V}_c(\vec{q}) & = \mathcal{F}\left[ \sum_m \sum_s V_a(\vec{x} - \vec{r}_{m,s} - \vec{u}_{m,s}) \right] \nonumber \\
                       & = \sum_m \sum_s f_e(\vec{q}) e^{-i \vec{q} \cdot \vec{r}_{m,s}} e^{-i \vec{q} \cdot \vec{u}_{m,s}}
    \label{eq:scattering-potential-temp}
\end{align}
Recall from @eq:scattering-intensity that the measurable quantity $|f^{(1)}(\vec{q})|^2$ is proportional to $|\hat{V}(\vec{q})|$:
\begin{align}
    |f^{(1)}(\vec{q})|^2 
        & = \left| -\frac{m_e}{\hbar^2} \hat{V}(\vec{q}) \right|^2 \nonumber \\
        & = \frac{m_e^2}{\hbar^4} \hat{V}(\vec{q}) \hat{V}^{\star}(\vec{q}) \nonumber \\
        & = \frac{m_e^2}{\hbar^4}  
            \left(\sum_m \sum_{s} f_{e,s}(\vec{q}) e^{-i \vec{q} \cdot \vec{r}_{m,s}} e^{-i \vec{q} \cdot \vec{u}_{m,s}} \right)  
            \left(\sum_{m^\prime} \sum_{s^{\prime}} f_{e,s^{\prime}}(\vec{q}) e^{i \vec{q} \cdot \vec{r}_{m^\prime, s^{\prime}}} e^{i \vec{q} \cdot \vec{u}_{m^\prime, s^{\prime}}} \right) \nonumber \\
        & = \frac{m_e^2}{\hbar^4} \sum_{m, m^\prime }\sum_{s, s^{\prime}} f_{e,s}(\vec{q}) f_{e,s^{\prime}}(\vec{q})  
            e^{-i \vec{q} \cdot (\vec{r}_{m,s} - \vec{r}_{m^\prime, s^{\prime}})}   
            e^{-i \vec{q} \cdot \vec{u}_{m,s}} e^{ i \vec{q} \cdot \vec{u}_{m^\prime, s^{\prime}}}
    \label{eq:scattering-average}
\end{align}
Note that the term $e^{-i \vec{q} \cdot \vec{u}_{m,s}} e^{ i \vec{q} \cdot \vec{u}_{m^\prime, s^{\prime}}}$ cannot always be simplified because $\set{\vec{u}_{m,s}}$ are quantum-mechanical operators, as will be explicitly stated below. The evaluation of the sum in @eq:scattering-average requires some thought. First, note that the base atomic positions $\set{\vec{r}_{m,s}}$ are independent of the particular unit cell $m$ where atom $s$ is located. The same is not true of the displacement vectors $\set{\vec{u}_{m,s}}$, where $\vec{u}_{m,s}$ and $\vec{u}_{m^\prime, s^{\prime}}$ are essentially uncorrelated. For a large enough crystal, the sum over $m$ and $m^{\prime}$ is equivalent to a thermal average over time. To this end, let us define the average $\langle \cdots \rangle \equiv \frac{1}{N_c}\sum_m \cdots$. @eq:scattering-average can then be expressed as:
$$
    |f^{(1)}(\vec{q})|^2 = N_c^2 \frac{m_e^2}{\hbar^4} \sum_{m, m^\prime} \sum_{s,s^{\prime}} f_{e,s}(\vec{q}) f_{e,s^{\prime}}(\vec{q}) e^{-i \vec{q} \cdot (\vec{r}_{m,s} - \vec{r}_{m^\prime, s^{\prime}})}   
            \langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle
    \label{eq:scattering-amplitude-average}
$$
where the problem is reduced to the evaluation of $\langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle$. 

### Quantizing lattice waves

Consider now the description of the displacement vectors as a superposition of lattice waves, or phonons. In this case, every displacement vector $\vec{u}_s$ can be expressed as:
$$
    \vec{u}_{m,s} = \sum_{\lambda} \sum_{\set{\vec{k}}} \sqrt{\frac{\hbar}{2 \mu_s N \omega_{\lambda}(\vec{k})}} \left( a_{\lambda}(\vec{k})e^{-i\phi_{s,m,\lambda}(\vec{k})} + a_{\lambda}^{\dagger}(\vec{k}) e^{i\phi_{s,m,\lambda}(\vec{k})} \right) \vec{e}_{s,\lambda}(\vec{k}) 
$${#eq:scattering-displacement}
where $\set{\lambda}$ label phonon branches, $\mu_s$ is the reduced mass of atom $s$, $N$ is the number of atoms in the crystal, $\omega_{\lambda}(\vec{k})$ is the vibrational frequency of mode $\lambda$ at wavevector $\vec{k}$, $a_{\lambda}(\vec{k})$ and $a_{\lambda}^{\dagger}(\vec{k})$ are the creation and annihilation operators for phonon mode $\lambda$, $\phi_{s,m,\lambda}(\vec{k})$ is the phase of a lattice wave, and $\vec{e}_{s,\lambda}(\vec{k})$ is the polarization vector of mode $\lambda$[@Sinha2001]. The expression for $\vec{u}_{m,s}$ is the combined effect of all possible phonon modes at the $\vec{r}_{m,s}$ lattice site. 

The Baker-Campbell-Hausdorff lemma can be used to compute the average $\langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle$ [@Hausdorff1906]. It states that for two operators $A$ and $B$ with commutator $[A,B]$:
$$
e^A e^B = e^{A + B + \frac{1}{2}[A,B]}
$$
This allows to simplify the average as:
$$
\langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ -i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle = \langle e^{-i \vec{q} \cdot (\vec{u}_s - \vec{u}_{s^\prime}) + \frac{1}{2}[\vec{q} \cdot \vec{u}_s, \vec{q} \cdot \vec{u}_{s^{\prime}}]} \rangle
$$
Furthermore, note that $[a_{\lambda}(\vec{k}), a_{\lambda}^{\dagger}(\vec{k})] = 1$ so that the following simplification is valid:
$$
\langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ -i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle = \langle e^{-i \vec{q} \cdot (\vec{u}_s - \vec{u}_{s^\prime})}\rangle \langle e^{\frac{1}{2}[\vec{q} \cdot \vec{u}_s, \vec{q} \cdot \vec{u}_{s^{\prime}}]} \rangle
$$
Finally, for operators $A$ which are a linear combination of position and momentum operators of a harmonic system, $\langle e^A \rangle = e^{\frac{1}{2}\langle A^2 \rangle}$[@Born1941]. This leads to:
$$
\langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ -i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle = e^{-\frac{1}{2}\langle (\vec{q} \cdot \vec{u}_s)^2\rangle} e^{-\frac{1}{2}\langle (\vec{q} \cdot \vec{u}_{s^{\prime}})^2\rangle} e^{\langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle}
$$
The terms $e^{-\frac{1}{2}\langle (\vec{q} \cdot \vec{u}_s)^2\rangle}$ and $e^{-\frac{1}{2}\langle (\vec{q} \cdot \vec{u}_{s^{\prime}})^2\rangle}$ are known as the Debye-Waller factors[@Waller1923; @Waller1928], historically defined as:
$$
e^{-\frac{1}{2}\langle (\vec{q} \cdot \vec{u}_s)^2\rangle} \equiv e^{-W_s}
$$
which means that
$$
\langle e^{-i \vec{q} \cdot \vec{u}_s} e^{ -i \vec{q} \cdot \vec{u}_{s^{\prime}}} \rangle = e^{-W_s} e^{-W_{s^\prime}} e^{\langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle}
$$
For small displacement vectors $\vec{u}$, $\vec{q} \cdot \vec{u} \leq |\vec{q}| |\vec{u}|$ is also small, and so:
\begin{align}
    e^{\langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle} 
        & = 1 + \langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle + \mathcal{O}\left(|\vec{u}_s|^2 |\vec{u}_{s^{\prime}}|^2 \right) \nonumber \\
        & \approx 1 + \langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle
\end{align}
Using @eq:scattering-displacement:
\begin{multline}
    (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) = \\
        \frac{\hbar}{2 N} \left( \sum_{\lambda} \sum_{\set{\vec{k}}} \frac{\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k})}{\sqrt{\mu_s \omega_{\lambda}(\vec{k})}} \left[ a_{\lambda}(\vec{k})e^{-i\phi_{s,m,\lambda}(\vec{k})} + a_{\lambda}^{\dagger}(\vec{k}) e^{i\phi_{s,m,\lambda}(\vec{k})} \right]\right) \\
        \left( \sum_{\lambda^\prime} \sum_{\set{\vec{k}^\prime}} \frac{\vec{q} \cdot \vec{e}_{\lambda^\prime,s^\prime}(\vec{k}^\prime)}{\sqrt{\mu_{s^\prime} \omega_{\lambda^\prime}(\vec{k}^\prime)}} \left[ a_{\lambda^\prime}(\vec{k}^\prime)e^{-i\phi_{s^\prime,m^\prime,\lambda^\prime}(\vec{k}^\prime)} + a_{\lambda^\prime}^{\dagger}(\vec{k}^\prime) e^{i\phi_{s^\prime,m^\prime,\lambda}(\vec{k}^\prime)} \right]\right)
\end{multline}
Equivalently:
\begin{multline}
    \langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle = \\
    \frac{\hbar}{2 N} \sum_{\lambda, \lambda^\prime} \sum_{\set{\vec{k}}, \set{\vec{k}^\prime}} \frac{\left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k}) \right) \left(\vec{q} \cdot \vec{e}_{\lambda,s^\prime}(\vec{k})\right)}{\sqrt{\mu_s \mu_{s^\prime} \omega_{\lambda}(\vec{k})\omega_{\lambda^\prime}(\vec{k}^\prime)}} \\
    \left\langle \left[ a_{\lambda}(\vec{k})e^{-i\phi_{s,m,\lambda}(\vec{k})} + a_{\lambda}^{\dagger}(\vec{k}) e^{i\phi_{s,m,\lambda}(\vec{k})} \right] \left[ a_{\lambda^\prime}(\vec{k}^\prime)e^{-i\phi_{s^\prime,m^\prime,\lambda^\prime}(\vec{k}^\prime)} + a_{\lambda^\prime}^{\dagger}(\vec{k}^\prime) e^{i\phi_{s^\prime,m^\prime,\lambda}(\vec{k}^\prime)} \right] \right\rangle
\end{multline}
Since the phases $\phi_{s,m,\lambda}(\vec{k})$ are not correlated across unit cells, the cross terms vanish:
\begin{align}
    \langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle 
    = & \frac{\hbar}{2 N} \sum_{\lambda} \sum_{\set{\vec{k}}} \frac{\left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k}) \right) \left(\vec{q} \cdot \vec{e}_{\lambda,s^\prime}(\vec{k})\right)}{\omega_{\lambda}(\vec{k})\sqrt{\mu_s \mu_{s^\prime}}} \nonumber \\
      & \left[ a_{\lambda}(\vec{k}) a_{\lambda}(\vec{k}) 
         + a_{\lambda}(\vec{k}) a^{\dagger}_{\lambda}(\vec{k}) 
         + a^{\dagger}_{\lambda}(\vec{k}) a_{\lambda}(\vec{k}) 
         + a^{\dagger}_{\lambda}(\vec{k}) a^{\dagger}_{\lambda}(\vec{k}) \right] \nonumber \\
    = & \frac{\hbar}{2 N} \sum_{\lambda} \sum_{\set{\vec{k}}} \frac{\left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k}) \right) \left(\vec{q} \cdot \vec{e}_{\lambda,s^\prime}(\vec{k})\right)}{\omega_{\lambda}(\vec{k})\sqrt{\mu_s \mu_{s^\prime}}} \left[ 2 n_{\lambda}(\vec{k}) + 1\right]
\end{align}
where $n_{\lambda}(\vec{k}) \equiv a_{\lambda}(\vec{k}) a^{\dagger}_{\lambda}(\vec{k}) = a^{\dagger}_{\lambda}(\vec{k}) a_{\lambda}(\vec{k})  - 1$ is the excitation number operator. Simplifying further:
$$
\langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle = \frac{\hbar}{N} \sum_{\lambda} \sum_{\set{\vec{k}}} \frac{n_{\lambda}(\vec{k}) + 1/2}{\omega_{\lambda}(\vec{k})}\frac{\left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k}) \right) \left(\vec{q} \cdot \vec{e}_{\lambda,s^\prime}(\vec{k})\right)}{\sqrt{\mu_s \mu_{s^\prime}}}
$$

### Scattering amplitude

```{.matplotlib #fig:scattering-vector-geometry file="figures/scattering/vector-geometry.py" caption="Diffraction pattern of SnSe (discussed in @sec:snse) showing the geometrical relationship between the scattering vector $\vec{q}$, reciprocal point $\vec{H}$, and wavevector $\vec{k}_0$. The in-plane section of the Brillouin, where $\vec{k}_0$ is confined, is shown as well."}
```

@eq:scattering-amplitude-average can then be expressed as:
\begin{align}
    |f^{(1)}(\vec{q})|^2 = & N_c^2 \frac{m_e^2}{\hbar^4} \sum_{m, m^\prime} \sum_{s,s^{\prime}} f_{e,s}(\vec{q}) f_{e,s^{\prime}}(\vec{q}) e^{-i \vec{q} \cdot (\vec{r}_{m,s} - \vec{r}_{m^\prime, s^{\prime}})} e^{-W_s} e^{-W_{s^\prime}}\left[ 1 + \langle (\vec{q} \cdot \vec{u}_s) ~ (\vec{q} \cdot \vec{u}_{s^{\prime}}) \rangle  \right] \nonumber \\
                         = & N_c^2 \frac{m_e^2}{\hbar^4} \left| \sum_m \sum_s f_{e,s}(\vec{q}) e^{-W_s} e^{-i \vec{q} \cdot \vec{r}_{m,s}} \right|^2 \nonumber \\
                         + & N_c^2 \frac{m_e^2}{N \hbar^3} \sum_{\lambda} \sum_{\set{\vec{k}}} \frac{n_{\lambda}(\vec{k}) + 1/2}{\omega_{\lambda}(\vec{k})} 
                            \left| \sum_m \sum_s \frac{f_{e,s}(\vec{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k})\right) e^{-i \vec{q} \cdot \vec{r}_{m,s}} \right|^2
\end{align}.
Finally, the sum over $\set{\vec{k}}$ can be removed as follows. Let $\vec{q} \equiv \vec{H} + \vec{k}_0$ where $\vec{H}$ is the nearest reciprocal point. Then:
$$
    \sum_m e^{-i \vec{q} \cdot \vec{r}_{m,s}} = \frac{L^3}{(2 \pi)^3} \sum_{\set{\vec{H}}} \delta(\vec{q} - \vec{H} - \vec{k}_0)
$$
where $L$ is some very large length used to discretize the allowed values of $\vec{k}$ (see @eq:scattering-norm). It follows that[@Sinha2001]:
\begin{align}
    |f^{(1)}(\vec{q})|^2 = & N_c^2 \frac{m_e^2}{\hbar^4} \left| \sum_m \sum_s f_{e,s}(\vec{q}) e^{-W_s} e^{-i \vec{q} \cdot \vec{r}_{m,s}} \right|^2 \nonumber \\
                         + &  \frac{L^3 N_c^2m_e^2}{N (2 \pi \hbar)^3} \sum_{\lambda} \frac{n_{\lambda}(\vec{k}_0) + 1/2}{\omega_{\lambda}(\vec{k}_0)} 
                            \left| \sum_s \frac{f_{e,s}(\vec{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k}_0)\right)\right|^2
    \label{eq:scattering-amplitude-reduced}
\end{align}
where $\vec{k}_0$ is understood to be the smallest wavevector such that $\vec{q} = \vec{H} + \vec{k}_0$ for some reciprocal point $\vec{H}$, i.e. $\vec{k}_0$ is constrained to lie in the first Brillouin zone. This implicit condition is equivalent to the sum $\sum_{\set{\vec{H}}}$. A visual representation of the vectors is shown in @fig:scattering-vector-geometry.

Combining @eq:scattering-amplitude-reduced and @eq:scattering-intensity, the measured intensity is therefore:
$$
    I(\vec{q}) = I_0(\vec{q}) + I_1(\vec{q})
$$
where
$$
    I_0(\vec{q}) = \frac{N_c^2 m_e^2}{r^2\hbar^4} \left| \sum_m \sum_s f_{e,s}(\vec{q}) e^{-W_s} e^{-i \vec{q} \cdot \vec{r}_{m,s}} \right|^2
$${#eq:scattering-diffracted-intensity-finite-temp}
is the diffracted intensity, and where
$$
    I_1(\vec{q}) = \frac{N_c^2 I_e}{r^2} \sum_{\lambda} \frac{n_{\lambda}(\vec{k}) + 1/2}{\omega_{\lambda}(\vec{k})} 
                            \left| \sum_s \frac{f_{e,s}(\vec{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k})\right)\right|^2
$${#eq:scattering-diffuse-intensity}
is known as *first order diffuse scattering*, described in the next section. The diffracted intensity at finite temperature is equivalent to @eq:scattering-diffracted-intensity-zero-temp, with the substitution $f_{e,s}(\vec{q}) \to f_{e,s}(\vec{q}) e^{-W_s}$. The physical meaning of this substitution is that atomic vibrations due to temperature decreases the periodicity of the lattice, which results in a suppression of the atomic form factor in reciprocal space.

### Diffuse scattering {#sec:diffuse-scattering}

The diffuse intensity can be expressed as:
$$
    I_1(\vec{q}) = \frac{N_c^2 I_e}{r^2} \sum_{\lambda} \frac{n_{\lambda}(\vec{k}) + 1/2}{\omega_{\lambda}(\vec{k})} |F_{1\lambda}(\vec{q})|^2
$${#eq:scattering-diffuse-intensity}
where $F_{1\lambda}(\vec{q})$ is known as the *one-phonon structure factor*:
$$
|F_{1\lambda}(\vec{q})|^2 = \left| \sum_s \frac{f_{e,s}(\vec{q}) e^{-W_s}}{\sqrt{\mu_s}} \left(\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k})\right)\right|^2
$${#eq:scattering-one-phonon-structure-factor}
named thus in contrast to the static structure factor of @eq:scattering-potential-crystal. 

A few clarifications can be made about diffuse intensity. The diffuse intensity at any scattering vector $\vec{q}$ involves the contribution of all phonon modes $\lambda$. The contribution of each mode can be conceptually separated into two parts: 

1. The vibrational amplitude of each mode is proportional to $(n_{\lambda}(\vec{k}) + 1/2)(\omega_{\lambda}(\vec{k}))$, by analogy with the expression of @eq:scattering-displacement. A higher population ($\uparrow n_{\lambda}$) results in a larger vibrational amplitude because the displacement of atoms is linear in the number of phonons that participate. A lower vibrational frequency ($\downarrow \omega_{\lambda}$) implies a smaller restoring force (in the harmonic oscillator sense), which also intuitively results in a wider vibrational amplitude.

2. A geometrical weight ($|F_{1\lambda}(\vec{q})|^2$) which determines if the atomic motion associated with a phonon mode can be captured on the detector. The most important term to consider are terms of the form $\set{\vec{q} \cdot \vec{e}_{\lambda,s}(\vec{k})}$. For a phonon polarization which is parallel to the propagation of the scattering electrons, the projection of the polarization onto the detector plane is $0$, and hence the associated diffuse intensity will not contribute.

Diffuse scattering and the effect of one-phonon structure factors is further explained in #sec:graphite.


[^tmatrix]: Interested readers are encouraged to peruse chapter 6 of Sakurai and Napolitano [@Sakurai2014].
[^conj]: Note that because $V(\vec{x})$ is real, $\int d\vec{x}^\prime e^{-i\vec{k} \cdot \vec{x}^\prime} V(\vec{x}) = \int d\vec{x}^\prime e^{i\vec{k} \cdot \vec{x}^\prime} V(\vec{x})$.

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
