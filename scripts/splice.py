#!/usr/bin/env python

"""
Pandoc filter to allow interpolation of metadata fields
into a document.  %{fields} will be replaced by the field's
value, assuming it is of the type MetaInlines or MetaString.

Modified from the examples in https://github.com/jgm/pandocfilters.
"""

import re

from pandocfilters import Span, Str, attributes, toJSONFilter

pattern = re.compile("%\{(.*)\}$")


def metavars(key, value, format, meta):
    if key == "Str":
        m = pattern.match(value)
        if m:
            field = m.group(1)
            result = meta.get(field, None)
            if result is None:
                raise RuntimeError(f"pandoc-splice: Missing value for key {field}.")
            if "MetaInlines" in result["t"]:
                return Span(
                    attributes({"class": "interpolated", "field": field}), result["c"]
                )
            elif "MetaString" in result["t"]:
                return Str(result["c"])


if __name__ == "__main__":
    toJSONFilter(metavars)
