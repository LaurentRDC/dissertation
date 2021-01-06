"""
Dissertation compiler

This script is meant to be used as a command-line program:

```
python -m dissc --help
```
"""
import contextlib
import argparse
from pathlib import Path
import subprocess
import os
from functools import wraps
import multiprocessing as mp
import sys
from itertools import chain
import warnings
from contextlib import suppress
import shutil
import tempfile
import time
import runpy
import logging

logging.basicConfig(encoding="utf-8", level=logging.INFO)

HERE = Path(os.getcwd())
BUILDDIR = HERE / "build"
BUILDDIR.mkdir(exist_ok=True)

CONTENTDIR = HERE / "content"
PLOTSDIR = BUILDDIR / "plots"
TEMPLATEDIR = HERE / "templates"

PANDOC = "pandoc"
LATEX_ENGINE = "pdflatex"

META = HERE / "metadata.yaml"

SRC = [
    CONTENTDIR / "listoffigures.md",
    CONTENTDIR / "preface.md",
    CONTENTDIR / "introduction.md",
    CONTENTDIR / "scattering.md",
    CONTENTDIR / "graphite.md",
    CONTENTDIR / "tase2.md",
    CONTENTDIR / "snse.md",
    CONTENTDIR / "conclusion.md",
]

BIBFILE = Path("references.bib")

APPENDIX = CONTENTDIR / "appendix.md"

TARGET = HERE / "dissertation.pdf"


## Auxiliary files
## (Do not change!)
TITLEPAGE = "titlepage.tex"
FRONTMATTER = "frontmatter.tex"

TMP1 = BUILDDIR / "__titlepage.filled.tex"
TMP2 = BUILDDIR / "__frontmatter.filled.tex"
TMP = [TMP1, TMP2]

AUX_OPTS = ["--wrap=preserve"]

OPTIONS = [
    "-f markdown+raw_tex"
]  # Some raw tex for \listoffigures macro and siunitx package
OPTIONS += ["--standalone"]

# The order of filters is important!
OPTIONS += ["--filter pandoc-plot"]
OPTIONS += ["-M plot-configuration=plot-config.yml"]

OPTIONS += ["--filter pandoc-crossref"]
OPTIONS += ['-M figPrefix="Figure"']
OPTIONS += ['-M secPrefix="Section"']

# We purposefully bypass pandoc-citeproc because we want
# to have references at the end of each chapter
# This is much easier to do with biblatex.
OPTIONS += ["--biblatex"]

OPTIONS += ["-M lang=en-CA"]
OPTIONS += [f"--metadata-file={META}"]

OPTIONS += [f"--include-in-header=include.tex"]
OPTIONS += [f"--include-in-header={TMP1}"]
OPTIONS += [f"--include-before-body={TMP2}"]

# OPTIONS += ["--listings"]

OPTIONS += ["-V documentclass=scrbook"]
OPTIONS += ["-V papersize=a4"]
OPTIONS += ["-V fontsize=11pt"]

OPTIONS += ["-V classoption:open=right"]
OPTIONS += ["-V classoption:twoside=true"]
OPTIONS += ["-V classoption:cleardoublepage=empty"]
OPTIONS += ["-V classoption:clearpage=empty"]

OPTIONS += ["-V geometry:top=30mm"]
OPTIONS += ["-V geometry:left=25mm"]
OPTIONS += ["-V geometry:bottom=30mm"]
OPTIONS += ["-V geometry:width=150mm"]
OPTIONS += ["-V geometry:bindingoffset=6mm"]

OPTIONS += ["--toc"]
OPTIONS += ["--toc-depth=3"]
OPTIONS += ["--number-sections"]
OPTIONS += ["--top-level-division=chapter"]

OPTIONS += ["-V colorlinks=true"]


## Template variables
EISVOGEL_TEMPLATE = "eisvogel.tex"
EISVOGEL_REPO = "https://github.com/Wandmalfarbe/pandoc-latex-template"
EISVOGEL_VERSION = "27fb7e455536012aa7e92151ffad28ff70986f41"

CLEANTHESIS_TEMPLATE = "cleanthesis.sty"
CLEANTHESIS_REPO = "https://github.com/derric/cleanthesis"
CLEANTHESIS_VERSION = "d89b30e141d2c62ae26bc5e34fe3db515015f258"

parser = argparse.ArgumentParser(prog="dissc", description="Dissertation compiler")

subparsers = parser.add_subparsers(help="sub-command help", dest="command")
parser_clean = subparsers.add_parser(
    "clean", help="Clean auxiliary files that are transiently generated during build."
)
parser_download_tempaltes = subparsers.add_parser(
    "download-templates", help="Download optional templates Eisvogel and Cleanthesis."
)

parser_build = subparsers.add_parser("build", help="Build dissertation.")
parser_build.add_argument(
    "--style",
    action="store",
    help="Style of dissertation to use. Default is `cleanthesis`.",
    choices=["simple", "eisvogel", "cleanthesis"],
    default="cleanthesis",
)


@wraps(subprocess.run)
def run(*args, **kwargs):
    cmd = args[0]
    logging.info("Running " + cmd)
    return subprocess.run(*args, **kwargs)


def runpandoc(options, target, sourcefiles, appendices=None):
    """
    Run pandoc with options and source files

    Parameters
    ----------
    options : List[str]
        List of options, e.g. ["--test true", "--foo=bar"]
    target : path-like
        Target file
    sourcefiles : List[path-like]
        List of source files
    references : path-like
        Reference file
    appendices : path-like
        Appendices
    """
    stringify = lambda xs: " ".join([str(x) for x in xs])
    if appendices is not None:
        sourcefiles += [appendices]
    return run(
        f"pandoc +RTS -N -RTS {stringify(options)} -o {target} {stringify(sourcefiles)}",
        shell=True,
        cwd=HERE,
    )


def render_diagram(source, target):
    """ Render SVG diagram `source` to `target` """
    try:
        run(
            f"inkscape --export-area-page --export-filename {target} {source}",
            shell=True,
            capture_output=True,
            cwd=HERE,
        ).check_returncode()
    except subprocess.CalledProcessError:
        warnings.warn(
            "Rendering of diagrams with inkscape failed. It might not be installed.",
            category=RuntimeWarning,
            stacklevel=0,
        )


def runlatex(source, target):
    """
    Run a full build of latex (i.e. latex, bibtex, 2xlatex)
    """
    target = Path(target)
    run(
        f"{LATEX_ENGINE} -interaction=batchmode -draftmode -aux-directory={BUILDDIR} -job-name={target.stem} {source} "
    )
    run(f"biber --quiet build/{Path(target).stem}")
    run(
        f"{LATEX_ENGINE} -interaction=batchmode -draftmode -aux-directory={BUILDDIR} -job-name={target.stem} {source}"
    )
    run(
        f"{LATEX_ENGINE} -interaction=batchmode -aux-directory={BUILDDIR} -job-name={target.stem} -output-directory={target.parent} {source}"
    )


def build(options, target, sourcefiles, appendices=None):

    runpandoc(
        options=options,
        target=BUILDDIR / "temp.tex",
        sourcefiles=sourcefiles,
        appendices=appendices,
    )

    runlatex(source=BUILDDIR / "temp.tex", target=target)


def download_template_files():
    """ Download template files to download directory """
    for repo, version, template in zip(
        [EISVOGEL_REPO, CLEANTHESIS_REPO],
        [EISVOGEL_VERSION, CLEANTHESIS_VERSION],
        [EISVOGEL_TEMPLATE, CLEANTHESIS_TEMPLATE],
    ):
        TEMPLATEDIR.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=TEMPLATEDIR) as downloaddir:
            run(
                f"git clone --quiet --single-branch --branch master --depth 100 {repo} {downloaddir}",
                shell=True,
            )
            run(f"git checkout --quiet {version}", cwd=downloaddir, shell=True)
            shutil.copy(Path(downloaddir) / template, TEMPLATEDIR / template)
            print("Successfully downloaded", template)


def build_auxiliary(aux_options=AUX_OPTS):
    """
    Build auxiliary files required by other builds.

    Parameters
    ----------
    aux_options : List[str]
        List of options, e.g. ["--test true", "--foo=bar"]
    """
    diagrams = list((HERE / "diagrams").glob("*.svg"))
    targets = [diagram.with_suffix(".pdf") for diagram in diagrams]
    nprocs = min([len(diagrams), os.cpu_count()])
    with mp.Pool(processes=nprocs) as pool:
        pool.starmap(render_diagram, zip(diagrams, targets))

    for template, result in zip([TITLEPAGE, FRONTMATTER], [TMP1, TMP2]):
        runpandoc(
            options=aux_options + [f"--template={template}", f"--metadata-file={META}"],
            target=result,
            sourcefiles=[template],
        )


def build_simple():
    """ Build the dissertation in the simple book style """
    build_auxiliary()
    build(
        options=OPTIONS,
        target=TARGET,
        sourcefiles=SRC,
        appendices=APPENDIX,
    )


def build_cleanthesis():
    """ Build the dissertation in the Cleanthesis style """
    if not (TEMPLATEDIR / CLEANTHESIS_TEMPLATE).exists():
        download_template_files()

    aux_options = AUX_OPTS
    aux_options += [f"-M cleanthesis=true -M cleanthesisbibfile={BIBFILE.stem}"]
    build_auxiliary(aux_options=aux_options)

    options = OPTIONS
    options += [f"--include-in-header=cleanthesis-include.tex"]
    options += aux_options

    build(
        options=options,
        target=TARGET,
        sourcefiles=SRC,
        appendices=APPENDIX,
    )


def build_eisvogel():
    """ Build the dissertation in the Eisvogel style. """
    if not (TEMPLATEDIR / EISVOGEL_TEMPLATE).exists():
        download_template_files()

    # We use run-py as a kind of script-runner
    # because it's too hard to import other python scripts
    # into this one
    if not (BUILDDIR / "titlepage.pdf").exists():
        runpy.run_path(HERE / "scripts" / "mktitlepage.py")

    aux_options = AUX_OPTS
    aux_options += ["-M eisvogel=true"]
    build_auxiliary(aux_options=aux_options)

    options = OPTIONS
    options += ["-V float-placement-figure=ht"]
    options += ["-V listings-no-page-break=true"]
    options += ["-V book=true"]
    options += ["-V toc-own-page=true"]
    options += ["-V titlepage=true"]
    options += ["-V titlepage-text-color=FFFFFF"]
    options += ["-V titlepage-rule-height=0"]
    options += ["-V titlepage-background=build/titlepage.pdf"]
    options += [f"--template={TEMPLATEDIR / EISVOGEL_TEMPLATE}"]

    build(
        options=options,
        target=TARGET,
        sourcefiles=SRC,
        appendices=APPENDIX,
    )


def clean():
    """ Clean generated files """
    shutil.rmtree(BUILDDIR, ignore_errors=True)
    with suppress(FileNotFoundError):
        os.remove(TARGET)


if __name__ == "__main__":

    arguments = parser.parse_args()
    if arguments.command == "build":
        multiplexor = {
            "cleanthesis": build_cleanthesis,
            "simple": build_simple,
            "eisvogel": build_eisvogel,
        }
        multiplexor[arguments.style]()
    elif arguments.command == "clean":
        clean()
    elif arguments.command == "download-templates":
        download_template_files()
    else:
        parser.print_help()
