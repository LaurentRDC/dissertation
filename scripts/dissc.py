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

LATEX_ENGINE = "lualatex"

CONTENTDIR = HERE / "content"
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

OPTIONS = ["-f markdown+raw_tex+latex_macros"]
OPTIONS += ["--standalone"]

# The order of filters is important!
OPTIONS += [
    "--filter pandoc-plot",
    "--filter scripts/splice.py",
    "--filter pandoc-crossref",
]
OPTIONS += [
    "-M cref:true",
    "-M autoEqnLabels:true",
    "-M plot-configuration=plot-config.yml",
]

OPTIONS += [
    f"--metadata-file={ HERE / 'metadata.yaml' }",  # Document metadata
    f"--metadata-file={ HERE / 'typography.yaml'}",  # Formatting
    f"--metadata-file={ HERE / 'variables.yaml' }",  # Splice variables, eg. fit values
]

OPTIONS += ["--biblatex", f"-V bibliography={BIBFILE}"]
OPTIONS += [
    f"-V biblatexoptions=backend=biber,citestyle=numeric,style=numeric,bibstyle=bibstyle,refsection=chapter,sorting=none,sortcites=true,autocite=superscript,maxnames=99"
]

OPTIONS += [f"--include-in-header=include.tex"]

OPTIONS += ["--toc", "--toc-depth=4"]
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
    "--print",
    action="store_const",
    const=True,
    help="Build the dissertation for printing media: no hyperlinks, single linespacing, etc.",
    default=False,
)

parser_prereqs = subparsers.add_parser(
    "compute-prerequisites", help="Compute plotting prerequisites from data."
)

parser_format = subparsers.add_parser("format", help="Format Python files.")


def check_for_todo(path):
    """Scan a file and return if there are any TODOs are left"""
    logging.info("Checking for TODOs...")
    with open(path, mode="r") as f:
        for line in f:
            if "TODO" in line:
                return True
    return False


def precompute_plots():
    """Compute the plotting prerequisites."""
    logging.info("Running prerequisites")
    for script in (HERE / "scripts" / "prerequisites").glob("*.py"):
        run(f"python -OO {script.relative_to(HERE)}").check_returncode()


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
    """Render SVG diagram `source` to `target`"""
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
    Run a full build of latex (i.e. latex, biber, 2xlatex)
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


def build(target, forprint=False):
    """Build the dissertation from source."""

    # Build diagrams
    diagrams = list((HERE / "diagrams").glob("*.svg"))
    targets = [diagram.with_suffix(".pdf") for diagram in diagrams]
    nprocs = min([len(diagrams), os.cpu_count()])
    with mp.Pool(processes=nprocs) as pool:
        pool.starmap(render_diagram, zip(diagrams, targets))

    # Build titlepage background
    if not (BUILDDIR_PDF / "titlepage.pdf").exists():
        run(f"python scripts/mktitlepage.py {BUILDDIR_PDF / 'titlepage.pdf'}")

    options = OPTIONS
    # The color DarkViolet is defined in the template
    if not forprint:
        options += ["-V linestretch=1.5"]
        options += ["-V linkcolor=DarkViolet"]
        options += ["-V citecolor=DarkViolet"]
        options += ["-V urlcolor=DarkViolet"]

    # Note that the path separators absolutely needs to be '\'
    titlepage_path = str(BUILDDIR_PDF / "titlepage.pdf").replace("\\", "/")
    options += [f"-V titlepage-background={titlepage_path}"]
    options += [f"--template={HERE / 'template.tex'}"]

    # We purposefully bypass pandoc-citeproc because we want
    # to have references at the end of each chapter
    # This is much easier to do with biblatex.
    runpandoc(
        options=options,
        target=BUILDDIR_PDF / target.with_suffix(".tex"),
        sourcefiles=SRC,
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


def clean(full=False):
    """Clean the build directory. If `full`, delete also the figures cache."""
    if full:
        shutil.rmtree(BUILDDIR_PDF, ignore_errors=True)
        logging.info(f"Removed {BUILDDIR_PDF}")
        return

    files = [entry.path for entry in os.scandir(path=BUILDDIR_PDF) if entry.is_file()]
    for f in files:
        with contextlib.suppress(PermissionError):
            os.remove(f)
        logging.info(f"Removed {f}")


def format_scripts():
    """Format all Python files in this directory"""
    run("isort --atomic .")
    run("black --safe .")


if __name__ == "__main__":
    arguments = parser.parse_args()
    if arguments.command == "build":
        build(target=Path("dissertation.pdf"), forprint=arguments.print)
    elif arguments.command == "clean":
        clean(full=arguments.all)
    elif arguments.command == "compute-prerequisites":
        precompute_plots()
    elif arguments.command == "format":
        format_scripts()
    else:
        parser.print_help()
