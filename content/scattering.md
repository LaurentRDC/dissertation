
# Electron scattering

The scattering of light or particles by a sample is used by a large class of experimental techniques, dating back a hundred years. Scattering can be broadly defined as the modification of an incoming wave by a potential into outgoing wave, a process which imprints the outgoing wave with some characteristic of the potential. The outgoing wave may lose or gain energy, and its momentum might be changed. When multiple incoming waves are simultaneously used -- forming an incoming wavefront --, the outgoing waves may interfere constructively or destructively. This effect is particularly intense for periodic scattering potentials, for example in crystals.

This chapter will consider the special case of *electron scattering*. In crystals, electrons are scattered by the electrostatic potential due to their charge. The description of X-ray scattering, for example, would be identical provided that the total electron charge-density be considered as the scattering potential, instead of the electrostatic potential.

Specifically, the type of scattering relevant to this work involves the probing of periodic structures, crystals. Thanks to advances in data acquisition and data analysis attributable to the author (see [appendix @sec:appendix]), as well as improvements to instrument stability[@Otto2017], the ultrafast electron scattering instrument used herein can reliably measure the effects of *diffuse scattering*.

## Single-electron scattering in an arbitrary potential

Consider an electron wavefunction $\Psi(\vec{x}, t)$. The scattering of $\Psi(\vec{x}, t)$ by an arbitrary potential $V(\vec{x},t)$ is described by the Schrödinger equation:

$$
i \hbar \frac{\partial}{\partial t} \Psi(\vec{x}, t) = \left[ \frac{- \hbar^2}{2 m_e} \nabla^2 + V(\vec{x}, t) \right]\Psi(\vec{x}, t).
$${#eq:schroedinger}

In the case of free space ($V(\vec{x}, t) = 0$), [@eq:schroedinger] reduces to the following equation:

$$
i \hbar \frac{\partial}{\partial t} \Psi(\vec{x}, t) = \frac{- \hbar^2}{2 m_e} \nabla^2 \Psi(\vec{x}, t).
$${#eq:free-electron}

It is instructive to consider the energy eigenfunction satisfying [@eq:free-electron]. Suppose that we can label such solutions as $\left\{(\Psi_a(\vec{x}, t), \omega_a) \right\}$, where the associated energy eigenvalues are $\omega_a \equiv E_a/\hbar$. We could then factorize the energy eigenfunction as:

$$
\Psi_a(\vec{x}, t) = u_a(\vec{x}) e^{-i \omega_a t}
$$

where $u_a(\vec{x})$ solves the time-independent Schrödinger equation[@Schroedinger1926]:

$$
\hbar \omega_a u_a(\vec{x}) = \frac{- \hbar^2}{2 m_e} \nabla^2 u_a(\vec{x}).
$${#eq:stationary-schroedinger}

We note that [@eq:stationary-schroedinger] is an instance of the Helmholtz equation, which can be re-written as:

$$ 
\nabla^2 u_a(\vec{x}) = -k_a^2 u_a(\vec{x})
$$

where $k_a^2 \equiv 2 m_e E_a/\hbar^2$.

## References {.unnumbered}
\printbibliography[heading=none]
