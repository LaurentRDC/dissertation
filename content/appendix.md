\appendix
# Interactive ultrafast electron scattering data exploration {#sec:appendix-software label="Appendix"}

The pursuit of science is accompanied by serious software development in most fields, often in the forms of instrument control, data acquisition, data analysis, and data presentation. One project spearheaded by the author, a graphical program for the interactive exploration of ultrafast electron scattering data [@RenedeCotret2018], has had such a profound impact on the research within the Siwick research group that it deserves its own brief section. 

## Graphical user-interface

Ultrafast electron scattering datasets contain a lot of information, even compared to other time-resolved techniques like optical spectroscopies. The research presented in this dissertation, as well as other research projects within the Siwick group[@Stern2018; @Otto2020], have benefitted from the ability to interactively explore the data. A program was designed and implemented by the author for this very purpose, named `iris`. `iris` solves two problems:

1. The reduction of redundant raw data into a processed dataset;

2. The interactive exploration of the reduced, processed dataset in real-time.

### Flexible data reduction

Ultrafast electron scattering experiments within the Siwick research group have always had a similar structure. The experiment is divided into subexperiments called *scans*. One scan involves the acquisition of scattering patterns for all desired time-delays. The scan is then repeated until either some desired number of scans, or until some time has passed. For example, the longest 

### Interactive exploration

## Streaming data reduction{#sec:appendix-stream-data-reduction}

Raw ultrafast electron scattering datasets can reach sizes hundreds of gigabytes. The data reduction operations such as averaging of equivalent scattering patterns is a slow process. Before the author joined the Siwick research group, it used to take multiple hours to finally look at scattering patterns following the end of an experiment.

To effectively reduce the data, a streaming data reduction engine called `npstreams` was developed. `npstreams` extends the de-facto standard array processing library `numpy` to work on streams of arrays, that is, sequences of arrays that are not all available in memory at once. Such streams can be created, for example, by loading images one by one from disk; while thousands of images may need to be processed, only 10 or 20 of them might be loaded into memory at any moment. The advantage of streaming lies in the memory savings, which in turns results in much better performance. There are also many streaming algorithms are more performant than the algorithms used by `numpy`[@West1979].

~~~{.matplotlib #fig:introduction-npstreams-benchmark file="figures/appendix/npstreams-benchmark.py" caption="Performance comparison between the npstreams computational engine and the *de-facto* standard numpy at averaging sequences of two-dimensional arrays (representing scattering patterns). **a)** Wall time of averaging a sequence of 10 arrays (data points). The dashed lines represent the maximum memory usage. The vertical line marks the scattering pattern size of the electron camera used in this work. **b)** Speed-up factor of using npstreams vs. numpy to average a sequence of arrays of size $512 \times 512$."}
~~~

The performance comparison between `npstreams` and `numpy` at averaging sequences of images is presented in @fig:introduction-npstreams-benchmark. In @fig:introduction-npstreams-benchmark a), fixed-length sequences of arrays of size $n \times n$ are averaged. This benchmark shows that the time-saving of using `npstreams` is correlated to the memory savings, as expected. @fig:introduction-npstreams-benchmark b), on the other hand, shows how much faster is `npstreams` at averaging variable-length streams of images with a fixed size. For reference, sequences of 100 -- 200 arrays of size $2048 \times 2048$ had to be averaged for the experiments shown in @sec:graphite and @sec:snse. As a final note, the low memory usage of `npstreams` also allowed to parallelize the averaging of images at different time-delays, which increased data reduction performance by a factor equal to the number of computer cores. Effectively, the data reduction step in the Siwick research group was improved from multiple hours down to less than 5 minutes.

## Functions and data structure for ultrafast electron scattering

### Image alignment

```{.matplotlib file="figures/appendix/mnxc.py" caption=""}
```

```{.matplotlib file="figures/appendix/autocenter.py" caption=""}
```

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
