"""
Dissertation compiler

This script is meant to be used as a command-line program:

```
python -m dissc --help
```
"""
import argparse
from pathlib import Path
from subprocess import run
import os
import sys
from itertools import chain
from contextlib import suppress
import shutil
import tempfile
import time

HERE = Path(__file__).parent
BUILDDIR = HERE / "build"
BUILDDIR.mkdir(exist_ok=True)

CONTENTDIR = HERE / "content"
PLOTSDIR = BUILDDIR / "plots"
TEMPLATEDIR = HERE / "templates"

PANDOC = "pandoc"

META = HERE / "metadata.yaml"

SRC = [
    CONTENTDIR / "listoffigures.md",
    CONTENTDIR / "preface.md",
    CONTENTDIR / "introduction.md",
    CONTENTDIR / "scattering.md",
    CONTENTDIR / "graphite.md",
    CONTENTDIR / "snse.md",
    CONTENTDIR / "tase2.md",
    CONTENTDIR / "conclusion.md",
]

BIBFILE = Path("references.bib")

APPENDIX = CONTENTDIR / "appendix.md"

TARGET = HERE / "dissertation.pdf"


## Auxiliary files
## (Do not change!)
TITLEPAGE = "titlepage.tex"
FRONTMATTER = "frontmatter.tex"
REFERENCES = "references.md"

TMP1 = BUILDDIR / "__titlepage.filled.tex"
TMP2 = BUILDDIR / "__frontmatter.filled.tex"
TMP = [TMP1, TMP2]

AUX_OPTS = ["--wrap=preserve"]

OPTIONS = [
    "-f markdown+raw_tex"
]  # Some raw tex for \listoffigures macro and siunitx package
OPTIONS += ["--pdf-engine=pdflatex"]
OPTIONS += ["--standalone"]

# The order of filters is important!
OPTIONS += ["--filter pandoc-plot"]

OPTIONS += ["--filter pandoc-crossref"]
OPTIONS += ['-M figPrefix="Figure"']

OPTIONS += ["--filter=pandoc-citeproc"]
OPTIONS += [f"-M bibliography={BIBFILE}"]
OPTIONS += ["-M link-citations=true"]

OPTIONS += ["-M lang=en-CA"]
OPTIONS += [f"--metadata-file={META}"]

OPTIONS += [f"--include-in-header=include.tex"]
OPTIONS += [f"--include-in-header={TMP1}"]
OPTIONS += [f"--include-before-body={TMP2}"]

## download from https://www.zotero.org/styles
## cf. https://pandoc.org/MANUAL.html#citations
# OPTIONS                += ["--csl=chicago-author-date-de.csl
# OPTIONS                += ["--csl=chicago-note-bibliography.csl
# OPTIONS                += ["--csl=ieee.csl
# OPTIONS                += ["--csl=oxford-university-press-note.csl

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

## Eisvogel (do not change!)
## https://github.com/Wandmalfarbe/pandoc-latex-template
OPTIONS += ["-V book=true"]
OPTIONS += ["-V titlepage=true"]
OPTIONS += ["-V toc-own-page=true"]


## Template variables
EISVOGEL_TEMPLATE = "eisvogel.tex"
EISVOGEL_REPO = "https://github.com/Wandmalfarbe/pandoc-latex-template"
EISVOGEL_VERSION = "5f740f8de0fb4c96dfb2772ef86e861fd3971654"

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


def runpandoc(options, target, sourcefiles, references=None, appendices=None):
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
    if references is not None:
        sourcefiles += [references]
    if appendices is not None:
        sourcefiles += [appendices]
    return run(
        f"pandoc {stringify(options)} -o {target} {stringify(sourcefiles)}",
        shell=True,
        cwd=HERE,
    )


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
    for template, result in zip([TITLEPAGE, FRONTMATTER], [TMP1, TMP2]):
        runpandoc(
            options=aux_options + [f"--template={template}", f"--metadata-file={META}"],
            target=result,
            sourcefiles=[template],
        )


def build_simple():
    """ Build the dissertation in the simple book style """
    build_auxiliary()
    runpandoc(
        options=OPTIONS,
        target=TARGET,
        sourcefiles=SRC,
        references=REFERENCES,
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

    runpandoc(
        options=options,
        target=TARGET,
        sourcefiles=SRC,
        references=REFERENCES,
        appendices=APPENDIX,
    )


def build_eisvogel():
    """ Build the dissertation in the Eisvogel style. """
    if not (TEMPLATEDIR / EISVOGEL_TEMPLATE).exists():
        download_template_files()

    aux_options = AUX_OPTS
    aux_options += ["-M eisvogel=true"]
    build_auxiliary(aux_options=aux_options)

    options = OPTIONS
    options += ["-V float-placement-figure=htbp"]
    options += ["-V listings-no-page-break=true"]
    options += [f"--template={TEMPLATEDIR / EISVOGEL_TEMPLATE}"]

    runpandoc(
        options=options,
        target=TARGET,
        sourcefiles=SRC,
        references=REFERENCES,
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
