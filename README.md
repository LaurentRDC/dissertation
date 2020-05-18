# PhD dissertation repository

## Setting up environment

To ensure that this dissertation content can be reused easily, I make use of conda environments:

```bash
conda env create -f conda_env.yml
```

## Building the dissertation

Once the environment has been set-up, use the `dissc.py` script to build.

```bash
> python dissc.py
usage: dissc [-h] [--simple] [--cleanthesis] [--eisvogel] [--clean] [--download-templates]

Dissertation compiler

optional arguments:
  -h, --help            show this help message and exit
  --simple              Build dissertation in simple book style.
  --cleanthesis         Build dissertation in cleanthesis style.
  --eisvogel            Build dissertation in Eisvogel style.
  --clean               Clean auxiliary files that are transiently generated during build.
  --download-templates  Download optional templates Eisvogel and Cleanthesis.
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

## Aknowledgement

This repository is partly based on [`pandoc-thesis`](https://github.com/cagix/pandoc-thesis).
