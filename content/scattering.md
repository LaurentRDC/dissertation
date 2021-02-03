
# Electron scattering

The scattering of light or particles by a sample is used by a large class of experimental techniques, dating back a hundred years. Scattering can be broadly defined as the modification of an incoming wave by a potential into outgoing wave, a process which imprints the outgoing wave with some characteristic of the potential. The outgoing wave may lose or gain energy, and its momentum might be changed. When multiple incoming waves are simultaneously used -- forming an incoming wavefront --, the outgoing waves may interfere constructively or destructively. This effect is particularly intense for periodic scattering potentials, for example in crystals.

This chapter will consider the special case of *electron scattering*. In crystals, electrons are scattered by the electrostatic potential due to their charge. The description of X-ray scattering, for example, would be identical provided that the total electron charge-density be considered as the scattering potential, instead of the electrostatic potential.

Specifically, the type of scattering relevant to this work involves the probing of periodic structures, crystals. Thanks to advances in data acquisition and data analysis attributable to the author (see [appendix @sec:appendix]), as well as improvements to instrument stability[@Otto2017], the ultrafast electron scattering instrument used herein can reliably measure the effects of *diffuse scattering*.

## Single-electron scattering in an arbitrary potential

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

In the case of free space ($V(\vec{x}, t) = 0$), [@eq:scattering-schroedinger] reduces to the following equation:
$$
i \hbar \frac{\partial}{\partial t} \Psi(\vec{x}, t) = \frac{- \hbar^2}{2 m_e} \nabla^2 \Psi(\vec{x}, t).
$${#eq:scattering-free-electron}
It is instructive to consider the energy eigenfunction satisfying [@eq:scattering-free-electron]. Suppose that we can label such solutions as $\left\{(\Psi_a(\vec{x}, t), \omega_a) \right\}$, where the associated energy eigenvalues are $\omega_a \equiv E_a/\hbar$. We could then factorize the energy eigenfunction as:
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
and the solution can be synthetized as a product of one-dimensional plane wave:
\begin{align}
u_a(\vec{x}) &= \prod_{j \in \{x,y,z\}} A_j e^{i (\vec{k}_a \cdot \hat{\vec{j}}) (\vec{x} \cdot \hat{\vec{j}})} \nonumber\\
             &= A e^{i \vec{k}_a \cdot \vec{x}}
\label{eq:scattering-plane-wave}
\end{align}
for some scalars $\{A_x, A_y, A_z, A\}$. Combining [@eq:scattering-free-space] [@eq:scattering-plane-wave] leads to the energy eigenstate of a free electron propagating in vacuum:
$$
\Psi_a(\vec{x}, t) = A e^{i (\vec{k}_a \cdot \vec{x} - \omega_a t)}
$$
with associated energy eigenvalue $E_a = \hbar^2 \vec{k}_a^2 / 2 m_e=\hbar \omega_a$

The arbitrary potential $V(\vec{x}, t)$ allows for electrons to scatter from an initial state $\ket{\vec{k}_i}$ to a final state $\ket{\vec{k}_f}$, where
$$
    \braket{\vec{x} | \vec{k}_a} \propto u_a(\vec{x}) \propto e^{i\vec{k}_a \cdot \vec{x}}
$$

TODO: how to bridge to the end? Lippman-Schwinger equation

$$
    \braket{\vec{x} | \vec{k}_f} = \braket{\vec{x} | \vec{k}_i} + \frac{e^{i |\vec{k}_f| |\vec{x}| }}{|\vec{x}|}f(\vec{k}_f, \vec{k}_i)
$${#eq:scattering-amplitude}
where $f(\vec{k}_f, \vec{k}_i)$ is called the *scattering amplitude*. The form of @eq:scattering-amplitude complies with intuition: the incoming state $\ket{\vec{k}_i}$ wave is composed of an unscattered part ($\braket{\vec{x} | \vec{k}_i}$) as well as an outgoing spherical wave with amplitude $f(\vec{k}_f, \vec{k}_i)$.

The scattering amplitude can be computed most simply by making use of the so-called *first Born approximation* [@Born1926]. This approximation is valid if the scattering potential is weak enough that the electron scatters a single time[@Feynman1965]. In this approximation:
$$
    f(\vec{k}_f, \vec{k}_i) = - \frac{m_e}{2 \pi \hbar^2} \int d\vec{x}^\prime e^{-i (\vec{k}_f - \vec{k}_i) \cdot \vec{x}^\prime} V(\vec{x}^\prime)
$$
The astute reader will have recognized that the scattering amplitude $f(\vec{k}_f, \vec{k}_i)$ is proportional to the Fourier transform of the scattering potential with respect to $\vec{k}_i - \vec{k}_f \equiv \vec{q}$, the *scattering vector*. We can re-express the scattering amplitude as follows:
$$
    f(\vec{q}=\vec{k}_f - \vec{k}_i) = -\frac{m_e}{\hbar^2} \hat{V}(\vec{q})
$${#eq:scattering-amplitude-q}
where the Fourier transform functional operator is defined as
$$
    \mathcal{F}\left[ f(\vec{x}) \right] \equiv \hat{f}(\vec{q}) = \frac{1}{2 \pi} \int d\vec{x}^\prime e^{-i \vec{q} \cdot \vec{x}^\prime}f(\vec{x}^\prime).
$$
Inserting @eq:scattering-amplitude-q in @eq:scattering-amplitude, the scattered wavefunction is given by:
$$
    u_f(\vec{x}) = u_i(\vec{x}) + \frac{e^{i |\vec{k}_f| |\vec{x}| }}{|\vec{x}|}f(\vec{k}_f - \vec{k}_i)
$$
and hence, by @eq:scattering-free-space:
$$
    \Psi_f(\vec{x}, t) = \Psi_i(\vec{x}, t) + \frac{e^{i |\vec{k}_f| |\vec{x}|}e^{ -i\omega_f t }}{|\vec{x}|}f(\vec{k}_f - \vec{k}_i)
$$

### Scattering potential of an atom

The scattering of electrons by atoms will now be considered. The scattering potential of a single atom $V(\vec{x})$ can be calculated from first principles using relativistic Hartree-Fock calculations[@Fischer1977; @Kirkland2010]. This is beyond the scope of the present work, and the final result is used here. 

When discussing electron scattering, the scattering amplitude is usually called the *atomic form factors for electron scattering*. To a reasonable degree of accuracy, the atomic form factor for electrons are spherically symmetric. They can be parametrized of an atom can be parametrized as:
$$
    f(\vec{q}) = \sum_{i=1}^3 \frac{a_i}{|\vec{q}|^2 + b_i} + c_i e^{-d_i|\vec{q}|^2}
$$
where the constants $\left\{ a_i, b_i, c_i, d_i \right\}$ are element-specific fitting parameters which are tabulated in Kirkland [@Kirkland2010]. The associated real-space potential can be calculated via @eq:scattering-amplitude-q. Examples of $f(\vec{q})$ and associated $V(\vec{x})$ are shown in @fig:scattering-potential for a few elements.

```{.matplotlib #fig:scattering-potential file="figures/scattering/scatt-potential.py" caption="Demonstration of the electrostatic potential of atoms, which scatters electrons. **a)** Radial view of the electrostatic potential in real-space **b)** Radial view of the electrostatic potential in reciprocal space, also known as the *atomic form factor*."}
```

The contribution of individual electronic orbitals to the atomic form factor for electrons[@Zheng2009] will be discussed in @sec:snse. 

## Single-electron scattering in a crystal

The main result of the previous section is that single-electron scattering is, to a reasonable approximation, a good measure of the Fourier transform of the scattering potential. It is therefore natural to consider the scattering of an electron by the periodic potential of a crystal. 

The discussion starts with visual representation of the real-space scattering potential of a crystal. Consider for moment $\alpha$-polonium ($\alpha$-Pu), one of the simplest crystal structures[@Curie1898], which crystals consists in a Pu atom at the every vertex of a cube with side-length of \SI{3.63}{\angstrom}. The calculated electrostatic potential of this arrangement along the base of a cube is shown in @fig:scattering-polonium-example a). The lattice vectors $\vec{a}_1$ and $\vec{a}_2$ are indicated, with $\vec{a}_3$ being aligned out of the page. The periodic nature of this scattering potential is demonstrated by calculating the resulting scattering amplitude $f(\vec{q})$ from @eq:scattering-amplitude-q, which is shown in @fig:scattering-polonium-example b). The periodicity in spatial-frequency-space, also called *reciprocal space*, is evident, and forms a *reciprocal lattice*, formally defined in the next section.

```{.matplotlib #fig:scattering-polonium-example file="figures/scattering/polonium.py" caption="Calculated scattering potential and associated scattering amplitude for $\alpha$-polonium. **a)** Electrostatic potential $V(\vec{x})$ in the $z=0$ plane. The two in-plane lattice vectors $\vec{a}_1$ and $\vec{a}_2$ are shown; lattice vector $\vec{a}_3$ points out of the page. **b)** Scattering amplitude $f(\vec{q})$ associated with the electrostatic potential shown in a). The periodic nature of the potential in real-space creates a structure in reciprocal space called the *reciprocal lattice*."}
``` 

### The reciprocal lattice

The geometry of reciprocal space and the reciprocal lattice are foundational concepts that drive the understanding of ultrafast electron diffraction.

A perfectly periodic structure in real-space that extends to infinity, with associated lattice vectors $\left\{ \vec{a}_i \right\}$, possesses a *dual* lattice in reciprocal space. The lattice vectors $\left\{ \vec{b}_i \right\}$ for this reciprocal lattice are defined by the relation
$$
    \vec{a}_i \cdot \vec{b}_j = 2 \pi \delta_ij
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
vectors reciprocal basis vectors $\left\{ \vec{b}_i \right\}$:
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

By definition, the reciprocal points of the crystal scattering potential, located at the spatial frequencies where the crystal potential strong, form a lattice with basis vectors $\left\{ \vec{b}_i \right\}$. Consider an electron with initial wavevector $\vec{k}_i$ that scatters elastically to a final wavevector $\vec{k}_f$. The scattering amplitude for this event, $f(\vec{k}_f, \vec{k}_i)$ is most intense where $\hat{V}(\vec{k}_f - \vec{k}_i)$ is strong; that is, the condition for strong scattering is:
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

```{.matplotlib file="figures/scattering/ewald.py" caption=""}
```



\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
