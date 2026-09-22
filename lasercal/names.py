"""MATLAB-compatible name handling.

The MATLAB tool relies on two behaviours of MATLAB itself:

1. ``readtable`` turns spreadsheet headers into valid MATLAB identifiers
   (``matlab.lang.makeValidName``).  ``7Li`` becomes ``x7Li``; an empty
   header becomes ``Var1`` (1-based column number).
2. The tool cleans RunOrder sample names (``BHVO-2G`` -> ``BHVO_2G``) and
   appends ``_2``, ``_3`` ... to duplicates.

Both are reproduced here so that row / column labels in every XLSX file
match what MATLAB would have written.
"""
from __future__ import annotations

import re
from typing import Iterable, List


def make_valid_name(name, col_number: int = 1) -> str:
    """Approximate ``matlab.lang.makeValidName`` (ReplacementStyle='underscore').

    Rules applied in the order MATLAB documents them:
    * empty / missing -> ``Var<col_number>``
    * whitespace is removed and the following letter capitalised
    * every character that is not a letter, digit or underscore -> ``_``
    * a leading character that is not a letter -> prefix ``x``
    * truncated to 63 characters (MATLAB namelengthmax)
    """
    if name is None:
        return f"Var{col_number}"
    s = str(name)
    if s.strip() == "" or s.lower() == "nan":
        return f"Var{col_number}"
    s = s.strip()
    # remove whitespace, capitalising the next letter
    s = re.sub(r"\s+([a-z])", lambda m: m.group(1).upper(), s)
    s = re.sub(r"\s+", "", s)
    s = re.sub(r"[^A-Za-z0-9_]", "_", s)
    if not re.match(r"[A-Za-z]", s[0]):
        s = "x" + s
    return s[:63]


def mangle_columns(columns: Iterable) -> List[str]:
    """Apply make_valid_name to a whole header row (1-based column numbers)."""
    return [make_valid_name(c, i + 1) for i, c in enumerate(columns)]


def clean_sample_names(sample_names: Iterable[str]) -> List[str]:
    """Reproduce the 'Clean & uniquify names' block of LaserCalTool.

    * a name whose first character parses as a number gets ``S_`` prefixed
    * ' ', '-', '#', '/' -> '_' ; runs of underscores collapsed
    * duplicates get ``_2``, ``_3`` ... appended in order of appearance
    """
    prefixed = []
    for nm in sample_names:
        nm = "" if nm is None else str(nm)
        if nm and _str2double_first_char_is_number(nm[0]):
            nm = "S_" + nm
        for ch in (" ", "-", "#", "/"):
            nm = nm.replace(ch, "_")
        nm = re.sub(r"__+", "_", nm)
        prefixed.append(nm)

    seen = {}
    display = []
    for k in prefixed:
        if k in seen:
            seen[k] += 1
            display.append(f"{k}_{seen[k]}")
        else:
            seen[k] = 1
            display.append(k)
    return display


def _str2double_first_char_is_number(ch: str) -> bool:
    """MATLAB: ``~isnan(str2double(nm(1)))`` for a single character."""
    return ch.isdigit()
