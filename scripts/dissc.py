"""
Dissertation compiler

This script is meant to be used as a command-line program:

```
python -m dissc --help
```
"""
import argparse
import contextlib
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import tempfile
import time
import warnings
from functools import wraps
from pathlib import Path

import termcolor

logging.basicConfig(encoding="utf-8", level=logging.INFO)

HERE = Path(os.getcwd())
BUILDDIR_PDF = HERE / "build"
BUILDDIR_PDF.mkdir(exist_ok=True)

CONTENTDIR = HERE / "content"
TEMPLATEDIR = HERE / "templates"

PANDOC = "pandoc"
LATEX_ENGINE = "lualatex"

META = HERE / "metadata.yaml"
TYPOGRAPHY = HERE / "typography.yaml"

SRC = [
    CONTENTDIR / "preface.md",
    CONTENTDIR / "introduction.md",
    CONTENTDIR / "scattering.md",
    CONTENTDIR / "graphite.md",
    CONTENTDIR / "snse.md",
    CONTENTDIR / "conclusion.md",
    CONTENTDIR / "appendix.md",
]

BIBFILE = Path("references.bib")

TITLEPAGE = "titlepage.tex"
FRONTMATTER = "frontmatter.tex"

TMP1 = BUILDDIR_PDF / "__titlepage.filled.tex"
TMP2 = BUILDDIR_PDF / "__frontmatter.filled.tex"
TMP = [TMP1, TMP2]

AUX_OPTS = ["--wrap=preserve"]

OPTIONS = [
    "-f markdown+raw_tex+latex_macros"
]  # Some raw tex for \listoffigures macro and siunitx package
OPTIONS += ["--standalone"]

# The order of filters is important!
OPTIONS += ["--filter pandoc-plot", "--filter pandoc-crossref"]
OPTIONS += [
    "-M cref:true",
    "-M autoEqnLabels:true",
    "-M plot-configuration=plot-config.yml",
]

OPTIONS += [f"--metadata-file={META}", f"--metadata-file={TYPOGRAPHY}"]

OPTIONS += ["--biblatex", f"-V bibliography={BIBFILE}"]
OPTIONS += [
    "-V biblatexoptions=backend=biber,citestyle=numeric,bibstyle=numeric,refsection=chapter,sorting=none,autocite=superscript,maxnames=99"
]

OPTIONS += [f"--include-in-header=include.tex"]
OPTIONS += [f"--include-in-header={TMP1}"]
OPTIONS += [f"--include-before-body={TMP2}"]

OPTIONS += ["--toc", "--toc-depth=3"]
OPTIONS += ["--number-sections", "--top-level-division=chapter"]

parser = argparse.ArgumentParser(prog="dissc", description="Dissertation compiler")

subparsers = parser.add_subparsers(help="sub-command help", dest="command")
parser_clean = subparsers.add_parser(
    "clean", help="Clean auxiliary files that are transiently generated during build."
)
parser_clean.add_argument(
    "--all",
    action="store_true",
    help="Clean build directory and figure files",
)

parser_build = subparsers.add_parser("build", help="Build dissertation.")
parser_build.add_argument(
    "--style",
    action="store",
    help="Style of dissertation to use. Default is `eisvogel`.",
    choices=["simple", "eisvogel"],
    default="eisvogel",
)
parser_build.add_argument(
    "--print",
    action="store_const",
    const=True,
    help="Build the dissertation for printing media: no hyperlinks, single linespacing, etc.",
    default=False,
)

parser_prereqs = subparsers.add_parser(
    "compute-prerequisites", help="Compute plotting prerequisites from data."
)


def check_for_todo(path):
    """ Scan a file and return if there are any TODOs are left """
    logging.info("Checking for TODOs...")
    with open(path, mode="r") as f:
        for line in f:
            if line.find("TODO"):
                return True
    return False


def precompute_plots():
    """ Compute the plotting prerequisites. """
    for script in [
        "mkdecomp.py",
        "mkdispersion.py",
        "mknpstreams-bench.py",
        "mkoneph.py",
    ]:
        path = HERE / "scripts" / script
        run(f"python -OO {path}")


@wraps(subprocess.run)
def run(cmd, *args, **kwargs):
    logging.info("Running " + cmd)
    cwd = kwargs.pop("cwd", HERE)
    return subprocess.run(cmd, *args, shell=True, cwd=cwd, **kwargs)


def runpandoc(options, target, sourcefiles):
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
    """
    stringify = lambda xs: " ".join([str(x) for x in xs])
    return run(
        f"pandoc +RTS -N -RTS {stringify(options)} -o {target} {stringify(sourcefiles)}",
    ).check_returncode()


def render_diagram(source, target):
    """ Render SVG diagram `source` to `target` """
    try:
        run(
            f"inkscape --export-area-page --export-filename {target} {source}",
        ).check_returncode()
    except subprocess.CalledProcessError:
        warnings.warn(
            "Rendering of diagrams with inkscape failed. It might not be installed.",
            category=RuntimeWarning,
            stacklevel=0,
        )


def runlatex(source):
    """
    Run a full build of latex (i.e. latex, bibtex, 2xlatex)
    """
    # Important: the options -aux-directory is Miktex-only, so I'm not using it
    # so that this script also supports TexLive. Also, confusingly, Texlive uses -jobname
    # and Miktex uses -job-name. Therefore, not using either.
    latex_options = (
        f" -interaction=batchmode -halt-on-error -output-directory={BUILDDIR_PDF}"
    )
    run(f"{LATEX_ENGINE} {latex_options} --draftmode {source}").check_returncode()
    run(f"biber --quiet build/{Path(source).stem}").check_returncode()
    run(f"{LATEX_ENGINE} {latex_options} --draftmode {source}").check_returncode()
    run(f"{LATEX_ENGINE} {latex_options} {source}").check_returncode()


def buildpdf(options, target, sourcefiles):
    # We purposefully bypass pandoc-citeproc because we want
    # to have references at the end of each chapter
    # This is much easier to do with biblatex.
    runpandoc(
        options=options,
        target=BUILDDIR_PDF / target.with_suffix(".tex"),
        sourcefiles=sourcefiles,
    )

    # TODO: also check for undefined references in the log file
    todo_left = check_for_todo(BUILDDIR_PDF / target.with_suffix(".tex"))

    try:
        runlatex(source=BUILDDIR_PDF / target.with_suffix(".tex"))
    except subprocess.CalledProcessError:
        print("--------------------------------")
        print("Error encountered. See log:")
        with open(BUILDDIR_PDF / target.with_suffix(".log"), "r") as f:
            for line in f:
                print(line)
        print("--------------------------------")
    if todo_left:
        print(termcolor.colored("WARNING: There are still TODOs remaining.", "red"))


def build_auxiliary(aux_options=AUX_OPTS, forprint=False):
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


def build_simple(target, forprint=False):
    """ Build the dissertation in the simple book style """
    build_auxiliary(forprint=forprint)

    options = OPTIONS
    if not forprint:
        options += ["-V linestretch=1.5"]

    buildpdf(
        options=options,
        target=target,
        sourcefiles=SRC,
    )


def build_eisvogel(target, forprint=False):
    """ Build the dissertation in the Eisvogel style. """

    # We use run-py as a kind of script-runner
    # because it's too hard to import other python scripts
    # into this one
    if not (BUILDDIR_PDF / "titlepage.pdf").exists():
        run(f"python scripts/mktitlepage.py {BUILDDIR_PDF / 'titlepage.pdf'}")

    aux_options = AUX_OPTS
    aux_options += ["-M eisvogel=true"]
    build_auxiliary(aux_options=aux_options, forprint=forprint)

    options = OPTIONS
    options += ["-V float-placement-figure=htpb"]
    options += ["-V listings-no-page-break=true"]
    options += ["-V book=true"]
    options += ["-V toc-own-page=true"]
    options += ["-V titlepage=true"]
    options += ["-V titlepage-text-color=FFFFFF"]
    options += ["-V titlepage-rule-height=0"]
    # The color DarkViolet is defined in the template
    if not forprint:
        options += ["-V linestretch=1.5"]
        options += ["-V linkcolor=DarkViolet"]
        options += ["-V citecolor=DarkViolet"]
        options += ["-V urlcolor=DarkViolet"]

    # Note that the path separators absolutely needs to be '\'
    titlepage_path = str(BUILDDIR_PDF / "titlepage.pdf").replace("\\", "/")
    options += [f"-V titlepage-background={titlepage_path}"]
    options += [f"--template={TEMPLATEDIR / 'eisvogel.tex'}"]

    buildpdf(
        options=options,
        target=target,
        sourcefiles=SRC,
    )


def clean(full=False):
    """ Clean the build directory. If `full`, delete also the figures cache. """
    if full:
        shutil.rmtree(BUILDDIR_PDF, ignore_errors=True)
        logging.info(f"Removed {BUILDDIR_PDF}")
        return

    files = [entry.path for entry in os.scandir(path=BUILDDIR_PDF) if entry.is_file()]
    for f in files:
        os.remove(f)
        logging.info(f"Removed {f}")


if __name__ == "__main__":

    arguments = parser.parse_args()
    if arguments.command == "build":
        multiplexor = {
            "simple": build_simple,
            "eisvogel": build_eisvogel,
        }
        multiplexor[arguments.style](
            target=Path("dissertation.pdf"), forprint=arguments.print
        )
    elif arguments.command == "clean":
        clean(full=arguments.all)
    elif arguments.command == "compute-prerequisites":
        precompute_plots()
    else:
        parser.print_help()
