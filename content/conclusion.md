
# Conclusion {#sec:conclusion}

Ultrafast electron diffuse scattering was the star of this dissertation. It was motivated in @sec:introduction as a way to get at a nearly-inaccessible facet of the impulse response of quantum systems. The theory of ultrafast diffuse scattering was developed in @sec:scattering. It showed that the transient Debye-Waller effect -- often measured in ultrafast electron diffraction -- is intricately linked with diffuse scattering. In @sec:graphite, the perfect benchmark system was studied with an eye for the particular strengths and limitations of diffuse scattering measurements. Those measurements were used to sidestep the ideas of the two-temperature model. Finally, in @sec:snse, ultrafast electron diffraction and diffuse scattering were used in conjunction to elucidate the mystery of thermoelectric performance in one of the best intrinsic thermoelectrics, SnSe. Strong electron-phonon coupling at zone-center suggested that strong carrier-lattice interactions are viable tuning parameters in the optimization of intrinsic thermoelectric performance. 

The last chapter on SnSe is a taste of things to come; the future is bright for ultrafast electron scattering. Following fantastic improvements in electron compression stability, \SI{50}{\femto\second} time-resolution is not far. As new laboratories are built with more stable laser oscillators, further improvements in performance are to be expected, to the point where the observation of direct phonon emission might be realistic. 

The author hopes that in the next 5-10 years, ultrafast electron scattering hardware will be commoditized, at which point it will become a standard tool in the study of low-dimensional materials. Given the unique ability to resolve lattice dynamics across the Brillouin zone, ultrafast electron diffuse scattering is poised to democratize the experimental access to strongly-coupled systems and revolutionize our understanding of strong interactions in functional materials.

## Outlook

Beyond the application of ultrafast electron scattering to more functional materials, the advances in experimental capabilities and sample preparation techniques will allow for the reliable observation of lattice dynamics in monolayers. The inversion-symmetry-breaking in many monolayers suggest an important future directions of research regarding the momentum-resolved measurement of lattice dynamics with ultrafast electron scattering: chiral phonons.

The control of electron flow in semiconductor logic gates is approaching a fundamental limit[@Markov2014]. However, inversion-symmetry breaking in many monolayers have given rise *valleytronics*, the control of electron flow by moving electrons between conduction valleys in reciprocal space[@Xiao2007;@Schaibley2016]. In particular, transition metal dichalcogenides such as MoS$_2$ and WSe$_2$ couple spin and valley properties via spin-orbit coupling, which makes them a rich playground for the light-based control of valley physics[@Mak2010;@Xiao2012;@Cao2012;@Jones2013]. Chiral phonons -- that is, phonons with pseudoangular momentum -- are modes which mediate the intravalley scattering of valley-polarized electrons in valleytronic materials [@Chen2015]. The momentum-resolution of ultrafast electron diffuse scattering would be key in the study of chiral phonons, as in hexagonal lattices they are located at the $K$/$K^\prime$ points[@Zhang2015;@Zhu2018]


```{.matplotlib #fig:conclusion-monolayers file="figures/conclusion/monolayer.py" caption="Diffraction patterns of ultrawide (\SI{250 x 250}{\micro\meter}) monolayers on \SI{10}{\nano\meter}-thick silicon nitride windows. **a)** Monolayer WSe$_2$. **b)** Monolayer MoS$_2$."}
```

Although ultrafast electron scattering interacts very strongly with matter, it was not clear whether diffraction and diffuse scattering measurements were possible on such thin samples. In this regard, the future is looking bright. @fig:conclusion-monolayers shows two diffraction patterns of monolayers (WSe$_2$ and MoS$_2$). These ultrawide monolayers -- (\SI{250 x 250}{\micro\meter}) -- were extracted from a bulk crystal via gold-mediated exfoliation[@Desai2016;@Liu2019;@Liu2020] and deposited on a \SI{10}{\nano\meter}-thick silicon nitride windows. The large area of these monolayers results in a number of diffracting cells that is comparable to typical multilayer samples of smaller area. In fact, the exposition conditions of the diffraction patterns in @fig:conclusion-monolayers (15 000 shots with a bunch charge of \SI{0.16}{\pico\coulomb}) are equivalent to those used in @fig:graphite-static and @fig:snse-diff-static.

There is little doubt that future ultrafast electron scattering measurements will have a large impact on the understanding of valleytronics. TODO: finish this

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]