"""

"""
import argparse
from pathlib import Path
from subprocess import run


WORKDIR = Path(__file__).parent
CONTENTDIR = WORKDIR / "content"

PANDOC = "pandoc"

META = WORKDIR / "metadata.yaml"

SRC = [
    CONTENTDIR / "introduction.md",
    CONTENTDIR / "scattering.md",
    CONTENTDIR / "conclusion.md",
]

BIBFILE = "references.bib"

APPENDIX = CONTENTDIR / "appendix.md"

TARGET = WORKDIR / "dissertation.pdf"


## Auxiliary files
## (Do not change!)
TITLEPAGE = "titlepage.tex"
FRONTMATTER = "frontmatter.tex"
REFERENCES = "references.md"

TMP1 = f"__titlepage.filled.tex"
TMP2 = f"__frontmatter.filled.tex"
TMP = [TMP1, TMP2]


## Pandoc options
AUX_OPTS = ["--wrap=preserve"]

OPTIONS = ["-f markdown"]
OPTIONS += ["--pdf-engine=pdflatex"]
OPTIONS += ["--standalone"]

OPTIONS += ["-M lang=en-CA"]
OPTIONS += [f"--metadata-file={META}"]

OPTIONS += [f"--include-in-header={TMP1}"]
OPTIONS += [f"--include-before-body={TMP2}"]

OPTIONS += ["--filter=pandoc-citeproc"]
OPTIONS += [f"-M bibliography={BIBFILE}"]
OPTIONS += ["-M link-citations=true"]
## download from https://www.zotero.org/styles
## cf. https://pandoc.org/MANUAL.html#citations
# OPTIONS                += ["--csl=chicago-author-date-de.csl
# OPTIONS                += ["--csl=chicago-note-bibliography.csl
# OPTIONS                += ["--csl=ieee.csl
# OPTIONS                += ["--csl=oxford-university-press-note.csl

OPTIONS += ["--listings"]

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

OPTIONS += ["-V colorlinks=true"]

## Eisvogel (do not change!)
## https://github.com/Wandmalfarbe/pandoc-latex-template
OPTIONS += ["-V book=true"]
OPTIONS += ["-V titlepage=true"]
OPTIONS += ["-V toc-own-page=true"]


## Template variables
TEMPLATE_DL_DIR         = WORKDIR / ".tmp_template_dl"

EISVOGEL_TEMPLATE       = "eisvogel.tex"
EISVOGEL_REPO           = "https://github.com/Wandmalfarbe/pandoc-latex-template"
EISVOGEL_VERSION        = "d5155adebf"

CLEANTHESIS_TEMPLATE    = "cleanthesis.sty"
CLEANTHESIS_REPO        = "https://github.com/derric/cleanthesis"
CLEANTHESIS_VERSION     = "c4609c4c70"

TEMPLATE_FILES          = [f"{EISVOGEL_TEMPLATE}",f"{CLEANTHESIS_TEMPLATE}"]


def build_auxiliary():
    """ Build auxiliary files required by other builds """
    for template, result in zip([TITLEPAGE, FRONTMATTER], [TMP1, TMP2]):
        run(
            f"{PANDOC} {' '.join(AUX_OPTS)} --template={template} --metadata-file={META} -o {result} {template}",
            shell=True,
            cwd=WORKDIR
        )

def clean_auxilary():

def build_simple():
    """
    Build the dissertation in the simple book style
    """
    build_auxiliary()

    run(
        f"{PANDOC} {' '.join(OPTIONS)} -o {TARGET} {' '.join([str(path) for path in SRC])} {REFERENCES} {APPENDIX}",
        cwd=WORKDIR,
        shell=True
    )
    
