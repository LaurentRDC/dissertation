\appendix
# Interactive ultrafast electron scattering data exploration {#sec:appendix-software label="Appendix"}

TODO: complete this chapter

The pursuit of science is accompanied by serious software development in most fields, often in the forms of instrument control, data acquisition, data analysis, and data presentation. One project spearheaded by the author, a graphical program for the interactive exploration of ultrafast electron scattering data [@RenedeCotret2018], has had such a profound impact on the research within the Siwick research group that it deserves its own brief section. 

## Data processing and graphical user interface

Ultrafast electron scattering datasets contain a lot of information, even compared to other time-resolved techniques like optical spectroscopies. The research presented in this dissertation, as well as other research projects within the Siwick group[@Stern2018; @Otto2020], have benefitted from the ability to interactively explore the data. A program was designed and implemented by the author for this very purpose. This program, named `iris`, solves three problems:

1. Data processing: the reduction of raw data, which is often very redundant, into a processed dataset;

2. The interactive exploration of the reduced, processed dataset in real-time in a graphical-user interface;

3. The ability to access the processed data in other computing environments.

### Flexible data reduction

Ultrafast electron scattering experiments within the Siwick research group have always had a similar structure. The experiment is divided into subexperiments called *scans*. One scan involves the acquisition of scattering patterns for all desired time-delays. The scan is then repeated until either some desired number of scans, or until some time has passed. For example, the longest 

### Interactive exploration

![Graphical user interface of `iris`. **Background** Interactive exploration of a dataset of photoexcited TiSe$_2$. **Foreground** Interactive exploration of a dataset of photoexcited polycrystalline VO$_2$](images/iris_screen.png){#fig:appendix-iris-screen}

### Open-source data format

`h5py`[@Collette2013]

## Streaming data reduction{#sec:appendix-stream-data-reduction}

Raw ultrafast electron scattering datasets can reach sizes hundreds of gigabytes. The data reduction operations such as averaging of equivalent scattering patterns is a slow process. Before the author joined the Siwick research group, it used to take multiple hours to finally look at scattering patterns following the end of an experiment.

To effectively reduce the data, a streaming data reduction engine called `npstreams` was developed. `npstreams` extends the de-facto standard array processing library `numpy`[@Van2011; @Harris2020] to work on streams of arrays, that is, sequences of arrays that are not all available in memory at once. Such streams can be created, for example, by loading images one by one from disk; while thousands of images may need to be processed, only 10 or 20 of them might be loaded into memory at any moment. The advantage of streaming lies in the memory savings, which in turns results in much better performance. There are also many streaming algorithms are more performant than the algorithms used by `numpy`[@West1979].

~~~{.matplotlib #fig:introduction-npstreams-benchmark file="figures/appendix/npstreams-benchmark.py" caption="Performance comparison between the npstreams computational engine and the *de-facto* standard numpy at averaging sequences of two-dimensional arrays (representing scattering patterns). **a)** Wall time of averaging a sequence of 10 arrays (data points). The dashed lines represent the maximum memory usage. The vertical line marks the scattering pattern size of the electron camera used in this work. **b)** Speed-up factor of using npstreams vs. numpy to average a sequence of arrays of size $512 \times 512$."}
~~~

The performance comparison between `npstreams` and `numpy` at averaging sequences of images is presented in @fig:introduction-npstreams-benchmark. In @fig:introduction-npstreams-benchmark a), fixed-length sequences of arrays of size $n \times n$ are averaged. This benchmark shows that the time-saving of using `npstreams` is correlated to the memory savings, as expected. @fig:introduction-npstreams-benchmark b), on the other hand, shows how much faster is `npstreams` at averaging variable-length streams of images with a fixed size. For reference, sequences of 100 -- 200 arrays of size $2048 \times 2048$ had to be averaged for the experiments shown in @sec:graphite and @sec:snse. As a final note, the low memory usage of `npstreams` also allows to parallelize the averaging of images at different time-delays, which increases data reduction performance by a factor equal to the number of computer cores. Effectively, the data reduction step in the Siwick research group was improved from multiple hours down to less than five minutes.

## Functions and data structure for ultrafast electron scattering

The scientific Python community has access to multiple research-oriented packages called *scikits*[@Pedregosa2011; @Van2014], which are extensions of the general-purpose SciPy package[@Virtanen2020].

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

```{.matplotlib file="figures/appendix/baseline.py" caption="Example of baseline-determination using the dual-tree complex wavelet transform. **a)** Polycrystalline diffraction of rutile VO$_2$ with an added known background, compared to the calculated baseline. **b)** Comparison of the true intensity (signal without background) and the background-subtracted intensity shows excellent agreement, without any prior knowledge about the background."}
```

Also applied recently to optical spectroscopy with great success[@Chevalier2019].

### Time-series analysis

### Parsing

### Image processing

#### Image alignment

During data acquisition, the pointing of the electron beam may shift slightly, which appears as a shift in diffraction patterns. This misalignment prevents the reduction of images which should be equivalent (e.g. same time-delay). This problem is hard to solve because some parts of images -- the beam-block, or dead pixels -- don't shift with the electron beam, as they are static with respect to the electron detector. To circumvent this problem, `scikit-ued` includes an alignment procedure based on the masked normalized cross-correlation algorithm[@Padfield2011]. This algorithm extends the cross-correlation in the Fourier domain, which is highly optimized, to images with masked segments. An example of this procedure is shown in @fig:appendix-mnxc.

```{.matplotlib @fig:appendix-mnxc file="figures/appendix/mnxc.py" caption=""}
```

#### Automatic center-finding

[@Liu2020]

```{.matplotlib file="figures/appendix/autocenter.py" caption=""}
```

#### Symmetrization

The symmetrization of diffraction pattern, first used in @sec:graphite-data-acquisition, is an image processing technique which point-group applies symmetry rules to images to enhance the signal-to-noise ratio.

```{.matplotlib file="figures/appendix/symmetry.py" caption="Symmetrization of a TiSe$_2$ diffraction pattern according to the $D6h$ point group (6-fold rotational symmetry). **a)** Diffraction pattern. The area that bounds the beam-block, which is ignored during symmetrization, is shown as a light-gray overlay. **b)** Symmetrized diffraction pattern."}
```

## Functions and data structures for crystallography

[@Togo2018]

Compatible with the Atomic Simulation Environment[@Larsen2017]

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
