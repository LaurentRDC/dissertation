
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
    \braket{\vec{x} | \vec{k}_f} = \braket{\vec{x} | \vec{k}_i} + \frac{e^{i |\vec{k}_f| |\vec{x}| }}{|\vec{x}|}f(\vec{k}_i, \vec{k}_f)
$${#eq:scattering-amplitude}
where $f(\vec{k}_i, \vec{k}_f)$ is called the *scattering amplitude*. The form of @eq:scattering-amplitude complies with intuition: the incoming state $\ket{\vec{k}_i}$ wave is composed of an unscattered part ($\braket{\vec{x} | \vec{k}_i}$) as well as an outgoing spherical wave with amplitude $f(\vec{k}_i, \vec{k}_f)$.

The scattering amplitude can be computed most simply by making use of the so-called *first Born approximation* [@Born1926]. In this approximation:
$$
    f(\vec{k}_i, \vec{k}_f) = - \frac{m_e}{2 \pi \hbar^2} \int d\vec{x}^\prime e^{i (\vec{k}_i - \vec{k}_f) \cdot \vec{x}^\prime} V(\vec{x}^\prime)
$$
The astute reader will have recognized that the scattering amplitude $f(\vec{k}_i, \vec{k}_f)$ is proportional to the Fourier transform of the scattering potential with respect to $\vec{k}_i - \vec{k}_f \equiv \vec{q}$, the *scattering vector*. Using the following definition for the Fourier transform functional operator:
$$
    \mathcal{F}\left[ f(\vec{x}) \right] = \frac{1}{2 \pi} \int d\vec{x}^\prime e^{i \vec{q} \cdot \vec{x}^\prime}f(\vec{x}^\prime) \equiv \hat{f}(\vec{q})
$$
, we can re-express the scattering amplitude as follows:
$$
    f(\vec{q}=\vec{k}_i - \vec{k}_f) = -\frac{m_e}{\hbar^2} \hat{V}(\vec{q})
$${#eq:scattering-amplitude-q}
Inserting @eq:scattering-amplitude-q in @eq:scattering-amplitude, the scattered wavefunction is given by:
$$
    u_f(\vec{x}) = u_i(\vec{x}) + \frac{e^{i |\vec{k}_f| |\vec{x}| }}{|\vec{x}|}f(\vec{k}_i - \vec{k}_f)
$$
and hence, by @eq:scattering-free-space:
$$
    \Psi_f(\vec{x}, t) = \Psi_i(\vec{x}, t) + \frac{e^{i |\vec{k}_f| |\vec{x}|}e^{ -i\omega_f t }}{|\vec{x}|}f(\vec{k}_i - \vec{k}_f)
$$

## Single-electron scattering in a periodic potential

```{.matplotlib file="figures/scattering/ewald.py" caption=""}
```

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
