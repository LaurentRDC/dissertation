# Ultrafast energy flow across the Brillouin zone in two-dimensional materials

## Setting up environment

The following tools are required to render this dissertation:

* `python` 3.7+;
* [`pandoc`](https://pandoc.org) 2.11+;
* [`pandoc-crossref`](https://github.com/lierdakil/pandoc-crossref);
* [`pandoc-plot`](https://github.com/LaurentRDC/pandoc-plot);
* [`inkscape`](https://inkscape.org/) 1.0+
* A LaTeX toolchain, including `pdflatex` and [`biber`](https://sourceforge.net/projects/biblatex-biber/). I tested with `MikTex` (Windows) and `texlive-full` (Ubuntu)

To install the Python dependencies required to render figures:

```bash
python -m pip install -r requirements.txt
```

## Building the dissertation

Once the environment has been set-up, use the `dissc.py` script to build.

```bash
> python dissc.py
usage: dissc [-h] {clean,download-templates,build} ...

Dissertation compiler

positional arguments:
  {clean,download-templates,build}
                        sub-command help
    clean               Clean auxiliary files that are transiently      
                        generated during build.
    download-templates  Download optional template Eisvogel
    build               Build dissertation.

optional arguments:
  -h, --help            show this help message and exit
```

The `build` subcommand has other options:

```bash
> python dissc.py build --help
usage: dissc build [-h] [--style {simple,eisvogel}]

optional arguments:
  -h, --help            show this help message and exit
  --style {simple,eisvogel}
                        Style of dissertation to use. Default is        
                        `eisvogel`.
(dissertation) PS C:\Users\Laurent\OneDrive\McGill\dissertation> 
```

Figures will be automatically generated by [`pandoc-plot`](https://github.com/LaurentRDC/pandoc-plot).

## Mark-up

This sections describes the available mark-up tools that the toolchain supports.

### Citations

Citations go inside square brackets and are separated by semicolons. Each citation must have a key, composed of ‘@’ + the citation identifier from the database, and may optionally have a prefix, a locator, and a suffix. The citation key must begin with a letter, digit, or _, and may contain alphanumerics, _, and internal punctuation characters (:.#$%&-+?<>~/). Here are some examples:

```markdown
Blah blah [see @doe99, pp. 33-35; also @smith04, chap. 1].

Blah blah [@doe99, pp. 33-35, 38-39 and *passim*].

Blah blah [@smith04; @doe99].
```

### Sections, figure, and table labels

Labels based on [`pandoc-crossref`](https://github.com/lierdakil/pandoc-crossref) are possible. To create a label for an item:

```markdown
# Introduction {#sec:intro}

As you can read in the @sec:intro, ...

```

### Equations

```markdown
$$ e = m c ^2 $$ {#eq:label}

Equation [@eq:label] ...
```

For more complex equation forms (e.g. aligned equations), `\label`s can be used. For example:

````markdown
\begin{align}
...
\label{eq:myequation}
\end{align}

You can refer to the equation like so @eq:myequation
````

### Appendices

Because the section label for appendices should be a bit different, use a custom label like so:

```markdown
lorem ipsum (see [appendix @sec:appendix]) ...
```

## Acknowledgement

This repository is partly based on [`pandoc-thesis`](https://github.com/cagix/pandoc-thesis).