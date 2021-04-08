\appendix
# Interactive ultrafast electron scattering data exploration {#sec:appendix-software label="Appendix"}

The pursuit of science is accompanied by serious software development in most fields, often in the forms of instrument control, data acquisition, data analysis, and data presentation. One project spearheaded by the author, a graphical program for the interactive exploration of ultrafast electron scattering data [@RenedeCotret2018], has had such a profound impact on the research within the Siwick research group that it deserves its own brief section. 

## Data processing and graphical user interface

Ultrafast electron scattering datasets contain a lot of information, even compared to other time-resolved techniques like optical spectroscopies. The research presented in this dissertation, as well as other research projects within the Siwick group[@Stern2018;@Otto2020], have benefitted from the ability to interactively explore the data. A program was designed and implemented by the author for this very purpose. This program, named `iris`, solves three problems:

1. Data processing: the reduction of raw data, which is often very redundant, into a processed dataset;

2. The interactive exploration of the reduced, processed dataset in real-time in a graphical-user interface;

3. The ability to access the processed data in other computing environments.

### Flexible data reduction

Ultrafast electron scattering experiments within the Siwick research group have always had a similar structure. The experiment is divided into subexperiments called *scans*. One scan involves the acquisition of scattering patterns for all desired time-delays. The scan is then repeated until either some desired number of scans, or until some time has passed. For example, the longest experiments on tin selenide presented in @sec:snse is composed of 66 scans, each of which is a one hour subexperiment.

Over the years, the specific data-acquisition scheme has changed. For example, time-delays started being acquired in random order (to minimize the effect of drifts in the laser power); a diagnostic of the instrument started being performed every 10 minutes. These new additions changed the way the raw data was acquired, but not what the processed data looked like. To abstract the specifics of the data acquisition scheme, `iris` includes a plug-in system. The specifics of the raw data are described in a Python plug-in, which `iris` reads and determines how to reduce the raw data appropriately. This allows the data acquisition scheme to change within the Siwick research group. Another consequence of the plug-in architecture is that other research groups can write plug-ins for `iris`, and use it in combination with their data acquisition scheme.

### Open-source data format

`iris` is designed to store the processed data in an open data format called Hierarchical Data Format version 5 (HDF5). This standardized format possesses three clear advantages when it comes to ultrafast electron scattering datasets. First, data in the form of multi-dimensional arrays can be natively stored in HDF5 files, with control over memory alignment. This means that sections of the data -- most importantly time-traces -- can be read from file at a high rate. Second, metadata of any type can also be stored in the same HDF5 file as the binary data. Third, there are interfaces between the HDF5 format and many programming languages and computing environments: official bindings[@Hdf5Bindings] include C, C++, Fortran, and Java, while third-party bindings are available for Python[@Collette2013], MATLAB[@MatlabHdf5], R[@RHdf5], and many others.

### Interactive exploration

The alignment of HDF5 datasets in memory allows to extract time-series from time-resolved scattering datasets in real-time. This allows datasets reduced by `iris` to be explored interactively. Interactive exploration of time-resolved scattering datasets has profoundly impacted the way research is done in the Siwick research group because the observation of diffuse scattering signals is not limited to Bragg peaks. With time-resolved diffraction experiments, the dynamics -- however complex[@Gao2013] -- can reduced to a set of time-traces, one for each reflection. With diffuse scattering, the space of possibilities is so large that measurements cannot be neatly reduced *a-priori*. Interactive exploration of datasets has been central to the understanding of time-resolved studies on graphite[@Stern2018;@RenedeCotret2019], TiSe$_2$[@Otto2020], and SnSe (@sec:snse).

![Graphical user interface of `iris`. **Background** Interactive exploration of a dataset of photoexcited TiSe$_2$. **Foreground** Interactive exploration of a dataset of photoexcited polycrystalline VO$_2$](images/iris_screen.png){#fig:appendix-iris-screen}

Time-series extraction in the graphical user-interface is shown in @fig:appendix-iris-screen for two types of samples: single-crystal TiSe$_2$ and polycrystalline VO$_2$. Further data transformations are also possible, some of which are described in @sec:appendix-scikit-ued. 

## Streaming data reduction{#sec:appendix-npstreams}

Raw ultrafast electron scattering datasets can reach sizes hundreds of gigabytes. The data reduction operations such as averaging of equivalent scattering patterns is a slow process. Before the author joined the Siwick research group, it used to take multiple hours to finally look at scattering patterns following the end of an experiment.

To effectively reduce the data, a streaming data reduction engine called `npstreams` was developed. `npstreams` extends the de-facto standard array processing library `numpy`[@Van2011;@Harris2020] to work on streams of arrays, that is, sequences of arrays that are not all available in memory at once. Such streams can be created, for example, by loading images one by one from disk; while thousands of images may need to be processed, only 10 or 20 of them might be loaded into memory at any moment. The advantage of streaming lies in the memory savings, which in turns results in much better performance. There are also many streaming algorithms are more performant than the algorithms used by `numpy`[@West1979].

~~~{.matplotlib #fig:introduction-npstreams-benchmark file="figures/appendix/npstreams-benchmark.py" caption="Performance comparison between the npstreams computational engine and the *de-facto* standard numpy at averaging sequences of two-dimensional arrays (representing scattering patterns). **a)** Wall time of averaging a sequence of 10 arrays (data points). The dashed lines represent the maximum memory usage. The vertical line marks the scattering pattern size of the electron camera used in this work. **b)** Speed-up factor of using npstreams vs. numpy to average a sequence of arrays of size $512 \times 512$."}
~~~

The performance comparison between `npstreams` and `numpy` at averaging sequences of images is presented in @fig:introduction-npstreams-benchmark. In @fig:introduction-npstreams-benchmark a), fixed-length sequences of arrays of size $n \times n$ are averaged. This benchmark shows that the time-saving of using `npstreams` is correlated to the memory savings, as expected. @fig:introduction-npstreams-benchmark b), on the other hand, shows how much faster is `npstreams` at averaging variable-length streams of images with a fixed size. For reference, sequences of 100 -- 200 arrays of size $2048 \times 2048$ had to be averaged for the experiments shown in @sec:graphite and @sec:snse. As a final note, the low memory usage of `npstreams` also allows to parallelize the averaging of images at different time-delays, which increases data reduction performance by a factor equal to the number of computer cores. Effectively, the data reduction step in the Siwick research group was improved from multiple hours down to less than five minutes.

## Functions and data structure for ultrafast electron scattering{#sec:appendix-scikit-ued}

The scientific Python community has access to multiple research-oriented packages called *scikits*[@Van2014], which are extensions of the general-purpose SciPy package[@Virtanen2020].

`scikit-ued` is a package which specializes in routines and data structures relevant to ultrafast electron diffraction and related techniques. `scikit-ued` powers most of the functionality of the graphical program `iris`; it is therefore easy for researchers to incorporate the same analysis techniques in their own workflow. The functionality of `scikit-ued` covers a large area of subjects, including but not limited to:

* Baseline-determination;
* Time-series analysis;
* Image processing;
* Diffraction simulation;
* Visualization and plotting;
* Interface to specialized file formats;
* Electron beam properties.

Some important examples are presented below. 

### Baseline-determination

Baseline-determination of polycrystalline and single-crystal diffraction patterns is included in `scikit-ued`. For polycrystalline diffraction patterns, an iterative baseline-removal routine[@RenedeCotret2017] based on the dual-tree complex transform[@Selesnick2005]. This routine allows for the removal of a baseline without any *a-priori* knowledge. Most importantly, the baseline-removal is stable; small changes in the time-resolved data do not result in large spontaneous changes in the baseline. An example of baseline-removal from a static polycrystalline diffraction pattern of rutile VO$_2$ is shown in @fig:appendix-baseline-vo2. This approach has also been extended to optical spectroscopy with great success[@Chevalier2019]. For single-crystal diffraction patterns, a two-dimensional iterative baseline-removal based on the discrete wavelet transform is provided[@Galloway2009].

```{.matplotlib #fig:appendix-baseline-vo2 file="figures/appendix/baseline.py" caption="Example of baseline-determination using the dual-tree complex wavelet transform. **a)** Polycrystalline diffraction of rutile VO$_2$ with an added known background, compared to the calculated baseline. **b)** Comparison of the true intensity (signal without background) and the background-subtracted intensity shows excellent agreement, without any prior knowledge about the background."}
```

### Parsing

`scikit-ued` includes a parser for multiple obscure file formats in which diffraction patterns are regularly found. Supported file formats include Merlin Image Binary (`.mib`), a proprietary image format for Quantum Detectors' MerlinEM direct electron imaging systems[@Maclaren2020], and Digital Micrograph 3 and 4 (`.dm3`, `.dm4`), proprietary image formats for Gatan Digital Micrograph software suite. In addition, all formats supported by `scikit-image`[@Van2014] are also transitively supported, including the Tagged Image File Format (`.tiff`).

### Simulation

Simple kinematic simulations of diffraction patterns and associated quantities are available within `scikit-ued`. These simulations are based on the parametrization of atomic electrostatic potential by Kirkland [@Kirkland2010]. The simulation modes include polycrystalline diffraction, single-crystal diffraction, electrostatic potential in real space, and projected electronstatic potential in real space. Some examples are presented in @fig:appendix-simulation.

```{.matplotlib #fig:appendix-simulation file="figures/appendix/simulation.py" caption="Examples of the simulation capabilities. **a)** Polycrystalline diffraction pattern of gold. **b)** Polycrystalline diffraction of graphite. **c)** Polycrystalline diffraction of M$_1$ vanadium dioxide. **d)** Electrostatic potential of a unit cell of cubic barium titanate, projected down the $c$ axis."}
```

### Image processing

The processing of diffraction patterns largely overlaps with photography and image processing. Many routines pertaining to image processing are available in `scikit-image`[@Van2014]. `scikit-ued` includes a few routines specific to ultrafast electron scattering which are presented below.

#### Image alignment

During data acquisition, the pointing of the electron beam may shift slightly, which appears as a shift in diffraction patterns. This misalignment prevents the reduction of images which should be equivalent (e.g. same time-delay). This problem is hard to solve because some parts of images -- the beam-block, or dead pixels -- don't shift with the electron beam, as they are static with respect to the electron detector. To circumvent this problem, `scikit-ued` includes an alignment procedure based on the masked normalized cross-correlation algorithm[@Padfield2011]. This algorithm extends the cross-correlation in the Fourier domain, which is highly optimized, to images with masked segments. An example of this procedure is shown in @fig:appendix-mnxc.

```{.matplotlib #fig:appendix-mnxc file="figures/appendix/mnxc.py" caption="Scattering pattern alignment based on the masked normalized cross-correlation algorithm. **a)** Reference diffraction pattern of polycrystalline chromium. The static area to be ignored is shown as a light gray overlay. **b)** Shifted diffraction pattern. **c)** Difference between the reference and shifted pattern. **d)** Difference between the reference and shifted pattern after translation."}
```

#### Automatic center-finding

Many data processing steps require the knowledge of the center of the diffraction patterns, such as Bragg peak-finding and symmetrization. `scikit-ued` includes a routine to automatically find the center of diffraction patterns. The routine measures the shift between an image and the centro-inverted image, based on initial guesses of the center. This shift is then used to infer the center of the image in a procedure similar to the previous section. This procedure is a generalization of the approach presented in Liu[@Liu2020], extended to work on any centrosymmetric image. An example of automatic center-finding for both polycrystalline and single-crystal diffraction patterns is shown in @fig:appendix-autocenter.

```{.matplotlib #fig:appendix-autocenter file="figures/appendix/autocenter.py" caption="Automatic determination of the center of diffraction patterns for **a)** polycrystalline chromium and **b)** single-crystal graphite."}
```

#### Symmetrization

The symmetrization of diffraction pattern, first used in @sec:graphite-data-acquisition, is an image processing technique which point-group applies symmetry rules to images to enhance the signal-to-noise ratio. Supported symmetry operations are mirror planes and discrete rotational symmetry. An example of symmetrization based on the $D6h$ point group is shown in @fig:appendix-symmetry.

```{.matplotlib #fig:appendix-symmetry file="figures/appendix/symmetry.py" caption="Symmetrization of a TiSe$_2$ diffraction pattern according to the $D6h$ point group (6-fold rotational symmetry). **a)** Diffraction pattern. The area that bounds the beam-block, which is ignored during symmetrization, is shown as a light-gray overlay. **b)** Symmetrized diffraction pattern."}
```

## Functions and data structures for crystallography

Accessing and manipulating crystal structures is an important component of many analysis workflows in electron microscopy, and ultrafast electron scattering is no different. `crystals` is a Python package designed and implemented by the author which facilitates the handling of crystallographic data. Some features of `crystals` are presented below.

### Parsing of crystal structure

`crystals` includes parsers for many standard crystallographic information files. Crystal structures can be parsed from the *de-facto* standard Crystallography Information File (CIF) format (`.cif`)[@Hester2006;@Bjorkman2011]. Crystal structures can also be downloaded from the Crystallography Open Database [@Gravzulis2009;@Gravzulis2012], which are stored in the CIF format. `crystals` also includes an internal database of simple crystal structures, stored in the CIF format. `crystal` also support downloading and parsing crystal structures from the Protein DataBank[@Berman2000;@Hamelryck2003].

`crystals` also supports the parsing of crystal structures determined by calculations. The conversion to and from the Atomic Simulation Environment[@Larsen2017] is supported. Additionally, the crystal structures resulting from plane-wave self-consistent field calculations from Quantum Espresso[@Giannozzi2017] are parseable by `crystals`. Finally, `crystals` can download and parse structures calculated by the Materials Project[@Hautier2010b;@Jain2013;@Ong2015].

### Representation of crystallographic information

The results of the parsing described in the previous section is stored a fundamental data structure called `Crystal`. A `Crystal` is a container of atomic positions and lattice information from which everything else can be derived. Most of the information in a `Crystal` is calculated on-the-fly, when requested by the user, so that the presented information is always up to date. Below is an example of the information that can be derived from the crystal structure of graphite:

<!--The verbatim code below requires usage of XeLaTeX / LuaLaTeX because of utf-8 symbols-->
```
> python -m crystals info graphite.cif
Crystal object with following unit cell:
    Atom C  @ (0.00, 0.00, 0.25)
    Atom C  @ (0.00, 0.00, 0.75)
    Atom C  @ (0.33, 0.67, 0.25)
    Atom C  @ (0.67, 0.33, 0.75)
Lattice parameters:
    a=2.464Å, b=2.464Å, c=6.711Å
    α=90.000°, β=90.000°, γ=120.000°
Chemical composition:
    C: 100.000%
Symmetry information:
    International symbol
                (short) ..... P6_3/mmc
                (full) ...... P 6_3/m 2/m 2/c
    International number .... 194
    Hermann-Mauguin symbol .. P63/mmc
    Pointgroup .............. 6/mmm
    Hall Number ............. 488
    Centering ............... CenteringType.primitive
```

### Symmetry-determination

The determination of crystal symmetries from atomic positions is included in `crystals`, based on the space-group library `spglib`[@Grosse1999;@Togo2018]. This allows to either check if a known structure possesses the expected symmetries, or to determine the point-group and space-group of a new structure.  

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
