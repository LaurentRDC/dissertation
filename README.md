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

## Aknowledgement

This repository is partly based on [`pandoc-thesis`](https://github.com/cagix/pandoc-thesis).

