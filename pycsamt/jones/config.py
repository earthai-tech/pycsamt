# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Configuration and canonical constants for the A.G. Jones
J‑format reader/writer (pycsamt.jones).

This module centralizes tokens, regex patterns, field names,
component codes, and small look‑ups that other submodules
(`utils`, `property`, `components`, `j`, etc.) can import.

Notes
-----
The J‑format is defined by A.G. Jones (1994, v2.0). A J file is
organized as:

1. COMMENT BLOCK — lines beginning with ``#`` (ignored by parser).
2. INFORMATION BLOCK — lines beginning with ``>KEY=VALUE``.
3. One or more DATA BLOCKs per site.

Data blocks start with a station line, a data‑type line, the number
of points to follow, then rows whose shape depends on the data type.
See ``FIELDS_R`` and ``FIELDS_TF`` below.

We intentionally keep this module free of I/O logic; only constants
and patterns live here to preserve a clean import graph and make unit
tests deterministic.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Mapping, Pattern, Sequence
import re


# Basic tokens and syntax markers

COMMENT_PREFIX: Final[str] = "#"
INFO_PREFIX: Final[str] = ">"
ENCODING_DEFAULT: Final[str] = "utf-8"
LINE_ENDING: Final[str] = "\n"

# J‑format version we target
JONES_SPEC_VERSION: Final[str] = "2.0"

# Missing/sentinel values used in J‑format
MISSING_FLOAT: Final[float] = -999.0


# Information‑block keys (normalized, upper‑case as they appear in files)
INFO_KEYS: Final[Sequence[str]] = (
    "AZIMUTH",  # site X‑axis azimuth in degrees (TN)
    "LATITUDE", # site latitude in degrees (decimal)
    "LONGITUDE",# site longitude in degrees (decimal)
    "ELEVATION",# site elevation in metres
)

# Normalization map from key -> canonical attribute name for our objects
INFO_KEY_ALIASES: Final[Mapping[str, str]] = {
    "AZIMUTH": "azimuth",
    "LATITUDE": "latitude",
    "LONGITUDE": "longitude",
    "ELEVATION": "elevation",
}


# Data type codes (2nd record) and units tokens
# First character denotes the *kind* of block
KIND_RESISTIVITY: Final[str] = "R"   # rho/phi (or upward‑biased when 'S')
KIND_RESISTIVITY_UP: Final[str] = "S"
KIND_IMPEDANCE: Final[str] = "Z"     # complex TF (or upward‑biased when 'Q')
KIND_IMPEDANCE_UP: Final[str] = "Q"
KIND_SCHMUCKER: Final[str] = "C"     # Schmucker C function (complex TF)
KIND_GDS_TF: Final[str] = "T"        # GDS transfer functions (complex TF)

KIND_COMPLEX_TF: Final[frozenset[str]] = frozenset({
    KIND_IMPEDANCE, KIND_IMPEDANCE_UP, KIND_SCHMUCKER, KIND_GDS_TF
})
KIND_RHO_PHI: Final[frozenset[str]] = frozenset({
    KIND_RESISTIVITY, KIND_RESISTIVITY_UP
})

# Second/third characters denote the *component* or aggregation
COMPONENT_CODES: Final[frozenset[str]] = frozenset(
    {
        "XX", "XY", "YX", "YY",  # impedance tensor entries
        "TE", "TM",                # rotational modes
        "AV", "DE",                # Berdichevsky average / determinant avg
        "ZX", "ZY",                # Tzx, Tzy
    }
)

# Map TE/TM to canonical tensor entries used elsewhere in pycsamt
TE_TM_TO_TENSOR: Final[Mapping[str, str]] = {
    "TE": "RXY",  # E‑polarization
    "TM": "RYX",  # B‑polarization
}

# Units accepted after Z/Q/C/T kinds. Only the keywords are relevant; extra
# words like "UNITS" or punctuation are for readability in the wild.
UNITS_CANONICAL: Final[Mapping[str, str]] = {
    "SI": "ohms",            # impedance in Ω
    "FIELD": "mV/km/nT",     # field units
}

# Liberal matcher to normalize unit tokens ("S.I.", "si", "FIELD UNITS", ...)
_RE_UNIT_TOKEN = r"(?P<unit>SI|S\.?I\.?|FIELD)"  # captured group name is 'unit'

# Column layouts for data rows

FIELDS_R: Final[Sequence[str]] = (
    "period", "rho", "pha", "rhomax", "rhomin",
    "phamax", "phamin", "wrho", "wpha",
)
FIELDS_TF: Final[Sequence[str]] = (
    "period", "real", "imag", "error", "weight",
)

# Rules implied by the specification (applied by the parser):
RULES_R: Final[Mapping[str, str]] = {
    "period_sign": "<0 means input row stores frequency in Hz",
    "rho_reject": "rho < 0 marks a rejected estimate",
    "wrho_reject": "wrho < 0 marks a rejected estimate",
}
RULES_TF: Final[Mapping[str, str]] = {
    "period_sign": "<0 means input row stores frequency in Hz",
}


# Regular expressions used by the lexer

RE_BANNER = re.compile(
    r"""
    ^\s*#\s*WRITTEN\s+BY\s+
    (?P<software>[^:]+):\s*
    (?P<station>\S+)\s+
    (?P<date>\S+)
    (?P<note>.*)$
    """,
    re.IGNORECASE | re.VERBOSE,
)

# Station name line: up to 6 visible ASCII without spaces (case‑insensitive)
# + azimuth 
RE_STATION = re.compile(
    r"""
    ^\s*
    (?P<station>[A-Za-z0-9][A-Za-z0-9._-]{0,31})
    (?:\s+               # optional azimuth on same line (float)
       [+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?
    )?
    \s*$
    """,
    re.VERBOSE,
)

# Data type + optional units, e.g. "RXY", "ZXX SI", "ZXY S.I. UNITS"
RE_DATATYPE_UNITS = re.compile(
    r"""
    ^\s*
    (?P<kind>[RSZQCT])         # R/S/Z/Q/C/T
    \s*                        # (defensive)
    (?P<comp>[A-Za-z]{2})      # XX/XY/YX/YY/TE/TM/AV/DE/ZX/ZY
    (?:\s*
        (?P<unit>SI|S\.?I\.?|FIELD)
        (?:\s*UNITS?)?         # optional 'UNITS' or 'UNIT'
    )?
    \s*$
    """,
    re.IGNORECASE | re.VERBOSE,
)


# Number of rows to follow
RE_NPOINTS: Final[Pattern[str]] = re.compile(r"^\s*(?P<n>\d+)\s*$")

# Flexible numeric field (Fortran‑style E-notation friendly)
NUM: Final[str] = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"

# Row with 9 free‑format floats (R/S kinds)
RE_ROW_R: Final[Pattern[str]] = re.compile(
    rf"^\s*(?P<p>{NUM})\s+(?P<rho>{NUM})\s+(?P<pha>{NUM})\s+"
    rf"(?P<rhomax>{NUM})\s+(?P<rhomin>{NUM})\s+(?P<phamax>{NUM})\s+"
    rf"(?P<phamin>{NUM})\s+(?P<wrho>{NUM})\s+(?P<wpha>{NUM})\s*$"
)

# Row with 5 free‑format floats (Z/Q/C/T kinds)
RE_ROW_TF: Final[Pattern[str]] = re.compile(
    rf"^\s*(?P<p>{NUM})\s+(?P<real>{NUM})\s+(?P<imag>{NUM})\s+"
    rf"(?P<error>{NUM})\s+(?P<weight>{NUM})\s*$"
)

# Information line: ">KEY = value" with optional spaces
RE_INFO: Final[Pattern[str]] = re.compile(
    r"^\s*>\s*(?P<key>[A-Za-z_][A-Za-z0-9_]*)\s*=\s*(?P<val>[^\n\r#]+?)\s*$"
)

# Comment line: leading '#'
RE_COMMENT: Final[Pattern[str]] = re.compile(r"^\s*#")

# Blank/whitespace line
RE_BLANK: Final[Pattern[str]] = re.compile(r"^\s*$")


# Component to tensor index mapping (for assembling 2×2 matrices)

TENSOR_INDEX: Final[Mapping[str, tuple[int, int]]] = {
    "XX": (0, 0),
    "XY": (0, 1),
    "YX": (1, 0),
    "YY": (1, 1),
}

# For convenience, accepted order when constructing averages (RAV/RDE)
CANONICAL_ORDER_R: Final[Sequence[str]] = (
    "RXY", "RYX",  # minimal pair for TE/TM averages
    "RXX", "RYY",  # include diagonals if present
)


@dataclass(frozen=True)
class DTypeSpec:
    """Describe a J data‑block layout.

    Parameters
    ----------
    fields : sequence of str
        Column names for each row in the block.
    regex : Pattern[str]
        Compiled regex to parse a single row.
    notes : str
        Short human‑readable note explaining special rules.
    """

    fields: Sequence[str]
    regex: Pattern[str]
    notes: str = ""

DTYPE_SPECS: Final[Mapping[str, DTypeSpec]] = {
    "R": DTypeSpec(FIELDS_R, RE_ROW_R, "rho/phi rows (9 floats)"),
    "S": DTypeSpec(FIELDS_R, RE_ROW_R, "upward‑biased rho/phi"),
    "Z": DTypeSpec(FIELDS_TF, RE_ROW_TF, "impedance rows (5 floats)"),
    "Q": DTypeSpec(FIELDS_TF, RE_ROW_TF, "upward‑biased impedance"),
    "C": DTypeSpec(FIELDS_TF, RE_ROW_TF, "Schmucker C function"),
    "T": DTypeSpec(FIELDS_TF, RE_ROW_TF, "GDS transfer functions"),
}

# Human‑friendly column descriptions (tooltips/metadata)
COL_DESCR: Final[Mapping[str, str]] = {
    "period": "Period in seconds; if <0, value is frequency (Hz)",
    "rho": "Apparent resistivity (Ω·m); <0 means rejected",
    "pha": "Phase (degrees)",
    "rhomax": "rho + 1σ",
    "rhomin": "rho − 1σ",
    "phamax": "pha + 1σ",
    "phamin": "pha − 1σ",
    "wrho": "Weight for rho; <0 means rejected",
    "wpha": "Weight for phase",
    "real": "Real part of transfer function",
    "imag": "Imaginary part of transfer function",
    "error": "Standard error",
    "weight": "Row weight (may be unreliable in some sources)",
}

# Friendly hints for writer defaults
FLOAT_FORMAT_R: Final[str] = "{:+.6e}"
FLOAT_FORMAT_TF: Final[str] = "{:+.6e}"

__all__ = [
    # tokens
    "COMMENT_PREFIX", "INFO_PREFIX", "ENCODING_DEFAULT", "LINE_ENDING",
    "JONES_SPEC_VERSION", "MISSING_FLOAT",
    # info
    "INFO_KEYS", "INFO_KEY_ALIASES",
    # kinds & components
    "KIND_RESISTIVITY", "KIND_RESISTIVITY_UP", "KIND_IMPEDANCE",
    "KIND_IMPEDANCE_UP", "KIND_SCHMUCKER", "KIND_GDS_TF",
    "KIND_COMPLEX_TF", "KIND_RHO_PHI", "COMPONENT_CODES",
    "TE_TM_TO_TENSOR", "UNITS_CANONICAL",
    # layouts & rules
    "FIELDS_R", "FIELDS_TF", "RULES_R", "RULES_TF",
    # regexes
    "RE_STATION", "RE_DATATYPE_UNITS", "RE_NPOINTS", "RE_ROW_R",
    "RE_ROW_TF", "RE_INFO", "RE_COMMENT", "RE_BLANK",
    # tensor helpers
    "TENSOR_INDEX", "CANONICAL_ORDER_R",
    # specs & descriptions
    "DTypeSpec", "DTYPE_SPECS", "COL_DESCR",
    # output formats
    "FLOAT_FORMAT_R", "FLOAT_FORMAT_TF",
]
