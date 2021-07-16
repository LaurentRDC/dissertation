"""
Generate the structure required by Rust's mdBook

"""

import argparse
import json
import re
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

Header = namedtuple(typename="Header", field_names="lineno,level,content")

# Matches ATX-style headings
HEADER_REGEX = re.compile(r"^#{1,6}\s+(.)*$")

SPLIT_AT_LEVEL = 2

argparser = argparse.ArgumentParser(description="Split a markdown file into subfiles.")
argparser.add_argument(
    "input_path", metavar="FILE", type=Path, help="Source Markdown file."
)
argparser.add_argument(
    "output_dir",
    metavar="OUTPUT_DIR",
    type=Path,
    help="Directory where to place the subfiles.",
)

if __name__ == "__main__":
    args = argparser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True)

    with open(args.input_path, mode="r", encoding="utf8") as f:
        doc = f.readlines()

    headers_locations = list()
    for lineno, line in enumerate(doc):
        if HEADER_REGEX.match(line):
            level = len(line) - len(line.lstrip("#"))
            if level <= SPLIT_AT_LEVEL:
                # Possibility of classes in markdown header, e.g. "# Preface {#sec:preface .unnumbered)"
                content = line.lstrip("# ").split("{")[0].rstrip(" ")
                # Because certain headers look like so:
                # :::
                # ## Section: ...
                # :::
                #
                # we might have to walk back with an offset
                offset = 0
                if doc[lineno - 1]:
                    offset = -1
                headers_locations.append(
                    Header(lineno + offset, level=level, content=content)
                )

    summary = """
    # Ultrafast energy flow across the Brillouin zone in two-dimensional materials

    [Abstract](abstract.md)
    [Résumé](resume.md)
    [Acknowledgements](acknowledgements.md)
    [Preface](preface.md)

    """
    chapter = 0
    section = 0
    headers_locations = list(sorted(headers_locations, key=lambda h: h.lineno))
    unnumbered = lambda header: header.content in {
        "Abstract",
        "Résumé",
        "Acknowledgements",
        "Preface",
    }
    for hi, header in enumerate(headers_locations):

        if header.level == 1 and not unnumbered(header):
            chapter += 1
            section = 0
        elif header.level == 2 and not unnumbered(header):
            section += 1

        if unnumbered(header):
            fname = output_dir / f"{header.content.replace('é', 'e').lower()}.md"
        else:
            fname = output_dir / f"ch{chapter}_sec{section}.md"

        try:
            part_end = headers_locations[hi + 1].lineno - 1
        except IndexError:
            part_end = -1
        with open(fname, mode="w", encoding="utf8") as f:
            f.writelines(doc[header.lineno : part_end])
