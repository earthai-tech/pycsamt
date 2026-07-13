# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""SEG-EDI schema (keywords, options, ordering).

This module centralizes *tags* (block keywords), *families*,
allowed *options*, and the top-level BNF skeleton for SEG-EDI
files. It is intentionally light-weight and import-safe: there
is no I/O here — just constants and a couple of helpers. Use
it from readers/writers (e.g. :mod:`pycsamt.seg.edi`).

Notes
-----
- Keywords and option names are **UPPERCASE** in the standard.
  We keep tags exactly as they appear (e.g. ``">ZXXR"``).
- Sets are :class:`frozenset` to make them hashable and safe
  for reuse in dataclasses and default arguments.
- The schema aggregates both MT and EMAP blocks, including
  filtered EMAP forms (``F*``) and auxiliary diagnostics.
- The *option* tables here are **non-exhaustive** for some of
  the highly system-dependent blocks (``>TSERIES``, spectra).
  They list the stable core options we validate against; the
  reader may accept extras as opaque ``KEY=VALUE`` pairs.

References
----------
SEG MT/EMAP EDI v1.0 (1991) — consolidated into
:mod:`pycsamt.seg.schema` as constants for reuse.
"""

from __future__ import annotations

from collections.abc import Mapping

# ---------------------------------------------------------------------------
# Top-level file structure (BNF skeleton)
# ---------------------------------------------------------------------------

BNF_EDI_SKELETON: str = (
    "<edi_file> ::=\n"
    "  <head_block>\n"
    "  <info_block> <info_text>\n"
    "  <def_meas_section>\n"
    "  { <tseries_section> | <spectra_section> |\n"
    "    <mt_section> | <emap_section> | <other_section> }\n"
    "  <end_block>\n"
)

# ---------------------------------------------------------------------------
# Canonical section openers / structural tags
# ---------------------------------------------------------------------------

HEAD_TAGS: frozenset[str] = frozenset(
    {
        ">HEAD",
        ">INFO",
        ">FREQ",
        ">ZROT",
        ">RHOROT",
        ">END",
    }
)

SECTION_OPENERS: frozenset[str] = frozenset(
    {
        ">=DEFINEMEAS",
        ">=TSERIESSECT",
        ">=SPECTRASECT",
        ">=MTSECT",
        ">=EMAPSECT",
        ">=OTHERSECT",
    }
)

# Comments ("pseudo-tag" used by the format). Kept for completeness.
COMMENT_TAG: str = ">!"

MEAS_DEF_BLOCKS: frozenset[str] = frozenset({">EMEAS", ">HMEAS"})
TSERIES_BLOCKS: frozenset[str] = frozenset({">TSERIES"})
SPECTRA_BLOCKS: frozenset[str] = frozenset({">SPECTRA"})

# ---------------------------------------------------------------------------
# Option sets for stable blocks (used by validators)
# ---------------------------------------------------------------------------

HEAD_ALLOWED: frozenset[str] = frozenset(
    {
        "DATAID",
        "ACQBY",
        "FILEBY",
        "ACQDATE",
        "ENDDATE",
        "FILEDATE",
        "COUNTRY",
        "STATE",
        "COUNTY",
        "PROSPECT",
        "LOC",
        "LAT",
        "LONG",
        "ELEV",
        "UNITS",
        "STDVERS",
        "PROGVERS",
        "PROGDATE",
        "MAXSECT",
        "BINDATA",
        "EMPTY",
        # Some writers also include non-standard free fields; we ignore.
    }
)
HEAD_REQUIRED: frozenset[str] = frozenset(
    {
        "DATAID",
        "ACQBY",
        "FILEBY",
        "ACQDATE",
        "FILEDATE",
        "STDVERS",
        "PROGVERS",
        "PROGDATE",
    }
)

INFO_ALLOWED: frozenset[str] = frozenset({"MAXINFO"})

FREQ_ALLOWED: frozenset[str] = frozenset({"NFREQ", "ORDER", "CHKSUM"})
FREQ_REQUIRED: frozenset[str] = frozenset({"NFREQ"})
FREQ_ENUMS: Mapping[str, frozenset[str]] = {
    "ORDER": frozenset({"INC", "DEC"}),
}

ZROT_ALLOWED: frozenset[str] = frozenset({"NFREQ", "CHKSUM"})
ZROT_REQUIRED: frozenset[str] = frozenset({})  # dataset length must==NFREQ

RHOROT_ALLOWED: frozenset[str] = frozenset({"NFREQ", "CHKSUM"})
RHOROT_REQUIRED: frozenset[str] = frozenset({})

DEFINEMEAS_ALLOWED: frozenset[str] = frozenset(
    {
        "MAXCHAN",
        "MAXRUN",
        "MAXMEAS",
        "UNITS",
        "REFTYPE",
        "REFLOC",
        "REFLAT",
        "REFLONG",
        "REFELEV",
    }
)
# Required per standard: for measured data REFLAT/REFLONG/REFELEV;
# we treat them as conditionally required in validators.
DEFINEMEAS_REQUIRED: frozenset[str] = frozenset({})

EMEAS_ALLOWED: frozenset[str] = frozenset(
    {
        "ID",
        "CHTYPE",
        "X",
        "Y",
        "Z",
        "X2",
        "Y2",
        "Z2",
        "ACQCHAN",
        "FILTER",
        "SENSOR",
        "GAIN",
        "MEASDATE",
    }
)
EMEAS_REQUIRED_CORE: frozenset[str] = frozenset(
    {
        "ID",
        "CHTYPE",
        "X",
        "Y",
        "X2",
        "Y2",
    }
)
EMEAS_ENUMS: Mapping[str, frozenset[str]] = {
    "CHTYPE": frozenset({"EX", "EY"}),
}

HMEAS_ALLOWED: frozenset[str] = frozenset(
    {
        "ID",
        "CHTYPE",
        "X",
        "Y",
        "Z",
        "AZM",
        "DIP",
        "ACQCHAN",
        "FILTER",
        "SENSOR",
        "GAIN",
        "MEASDATE",
    }
)
HMEAS_REQUIRED_CORE: frozenset[str] = frozenset(
    {
        "ID",
        "CHTYPE",
        "X",
        "Y",
        "AZM",
    }
)
HMEAS_ENUMS: Mapping[str, frozenset[str]] = {
    "CHTYPE": frozenset({"HX", "HY", "HZ"}),
}

TSERIESSECT_ALLOWED: frozenset[str] = frozenset(
    {
        "SECTID",
        "NCHAN",
        "MAXBLKS",
        "CHKSUM",
    }
)
TSERIES_ALLOWED: frozenset[str] = frozenset(
    {
        # Core, kept minimal — systems add many more (e.g. NC, DT,...)
        "SECTID",
        "NPTS",
        "SR",
        "CHKSUM",
    }
)

SPECTRASECT_ALLOWED: frozenset[str] = frozenset(
    {
        "SECTID",
        "NC",
        "MAXBLKS",
        "CHKSUM",
    }
)
SPECTRA_ALLOWED: frozenset[str] = frozenset(
    {
        # Intentional minimal core; writers often include BW/NAV/DT, etc.
        "FREQ",
        "BWIDTH",
        "AVGF",
        "AVGT",
        "NC",
        "CHKSUM",
    }
)

MTSECT_ALLOWED: frozenset[str] = frozenset(
    {
        "SECTID",
        "NFREQ",
        "MAXBLKS",
        "HX",
        "HY",
        "HZ",
        "EX",
        "EY",
        "RX",
        "RY",
    }
)
EMAPSECT_ALLOWED: frozenset[str] = frozenset(
    {
        "SECTID",
        "NFREQ",
        "MAXBLKS",
        "NDIPOLE",
        "TYPE",
        "HX",
        "HY",
        "RX",
        "RY",
        "CHKSUM",
    }
)
OTHERSECT_ALLOWED: frozenset[str] = frozenset({"SECTID", "MAXBLKS"})

# Convenience registry (options we care about by tag)
BLOCK_OPTIONS: Mapping[str, Mapping[str, frozenset[str]]] = {
    ">HEAD": {"allowed": HEAD_ALLOWED, "required": HEAD_REQUIRED},
    ">INFO": {"allowed": INFO_ALLOWED, "required": frozenset()},
    ">FREQ": {"allowed": FREQ_ALLOWED, "required": FREQ_REQUIRED},
    ">ZROT": {"allowed": ZROT_ALLOWED, "required": ZROT_REQUIRED},
    ">RHOROT": {"allowed": RHOROT_ALLOWED, "required": RHOROT_REQUIRED},
    ">=DEFINEMEAS": {
        "allowed": DEFINEMEAS_ALLOWED,
        "required": DEFINEMEAS_REQUIRED,
    },
    ">EMEAS": {"allowed": EMEAS_ALLOWED, "required": EMEAS_REQUIRED_CORE},
    ">HMEAS": {"allowed": HMEAS_ALLOWED, "required": HMEAS_REQUIRED_CORE},
    ">=TSERIESSECT": {
        "allowed": TSERIESSECT_ALLOWED,
        "required": frozenset(),
    },
    ">TSERIES": {"allowed": TSERIES_ALLOWED, "required": frozenset()},
    ">=SPECTRASECT": {
        "allowed": SPECTRASECT_ALLOWED,
        "required": frozenset(),
    },
    ">SPECTRA": {"allowed": SPECTRA_ALLOWED, "required": frozenset()},
    ">=MTSECT": {"allowed": MTSECT_ALLOWED, "required": frozenset()},
    ">=EMAPSECT": {"allowed": EMAPSECT_ALLOWED, "required": frozenset()},
    ">=OTHERSECT": {"allowed": OTHERSECT_ALLOWED, "required": frozenset()},
}

# ---------------------------------------------------------------------------
# Data block keywords by family
# ---------------------------------------------------------------------------

# Impedance (Z) blocks — real/imag, variance, covariance
Z_REAL_IMAG: frozenset[str] = frozenset(
    {
        ">ZXXR",
        ">ZXXI",
        ">ZXYR",
        ">ZXYI",
        ">ZYXR",
        ">ZYXI",
        ">ZYYR",
        ">ZYYI",
    }
)
Z_VAR_COV: frozenset[str] = frozenset(
    {
        ">ZXXR.VAR",
        ">ZXXI.VAR",
        ">ZXX.VAR",
        ">ZXX.COV",
        ">ZXYR.VAR",
        ">ZXYI.VAR",
        ">ZXY.VAR",
        ">ZXY.COV",
        ">ZYXR.VAR",
        ">ZYXI.VAR",
        ">ZYX.VAR",
        ">ZYX.COV",
        ">ZYYR.VAR",
        ">ZYYI.VAR",
        ">ZYY.VAR",
        ">ZYY.COV",
    }
)
# EMAP-filtered impedances
FZ_REAL_IMAG: frozenset[str] = frozenset(
    {
        ">FZXXR",
        ">FZXXI",
        ">FZXYR",
        ">FZXYI",
    }
)

# Apparent resistivity / phase
RHO_BLOCKS: frozenset[str] = frozenset(
    {
        ">RHOXX",
        ">RHOXY",
        ">RHOYX",
        ">RHOYY",
        ">FRHOXX",
        ">FRHOXY",  # filtered (EMAP)
    }
)
RHO_STATS: frozenset[str] = frozenset(
    {
        ">RHOXX.VAR",
        ">RHOXX.ERR",
        ">RHOXX.FIT",
        ">RHOXY.VAR",
        ">RHOXY.ERR",
        ">RHOXY.FIT",
        ">RHOYX.VAR",
        ">RHOYX.ERR",
        ">RHOYX.FIT",
        ">RHOYY.VAR",
        ">RHOYY.ERR",
        ">RHOYY.FIT",
        ">FRHOXX.VAR",
        ">FRHOXX.ERR",
        ">FRHOXX.FIT",
        ">FRHOXY.VAR",
        ">FRHOXY.ERR",
        ">FRHOXY.FIT",
    }
)
PHS_BLOCKS: frozenset[str] = frozenset(
    {
        ">PHSXX",
        ">PHSXY",
        ">PHSYX",
        ">PHSYY",
        ">FPHSXX",
        ">FPHSXY",
    }
)
PHS_STATS: frozenset[str] = frozenset(
    {
        ">PHSXX.VAR",
        ">PHSXX.ERR",
        ">PHSXX.FIT",
        ">PHSXY.VAR",
        ">PHSXY.ERR",
        ">PHSXY.FIT",
        ">PHSYX.VAR",
        ">PHSYX.ERR",
        ">PHSYX.FIT",
        ">PHSYY.VAR",
        ">PHSYY.ERR",
        ">PHSYY.FIT",
        ">FPHSXX.VAR",
        ">FPHSXX.ERR",
        ">FPHSXX.FIT",
        ">FPHSXY.VAR",
        ">FPHSXY.ERR",
        ">FPHSXY.FIT",
    }
)

# Tipper and tipper-related experiment outputs
TIPPER_BLOCKS: frozenset[str] = frozenset(
    {
        ">TIPMAG",
        ">TIPPHS",
    }
)
TIPPER_AUX: frozenset[str] = frozenset(
    {
        ">TIPMAG.ERR",
        ">TIPMAG.FIT",
        ">TIPPHS.FIT",
        ">TXR.EXP",
        ">TXI.EXP",
        ">TXVAR.EXP",
        ">TYR.EXP",
        ">TYI.EXP",
        ">TYVAR.EXP",
    }
)

# Coherency, signal, and noise diagnostics
COH_BLOCKS: frozenset[str] = frozenset(
    {
        ">COH",
        ">EPREDCOH",
        ">HPREDCOH",
        ">SIGAMP",
        ">SIGNOISE",
    }
)

# Strike / skew / ellipticity
GEOM_BLOCKS: frozenset[str] = frozenset(
    {
        ">ZSTRIKE",
        ">ZSKEW",
        ">ZELLIP",
        ">TSTRIKE",
        ">TSKEW",
        ">TELLIP",
    }
)

# Spatial filters (EMAP)
SPATIAL_FILTER_BLOCKS: frozenset[str] = frozenset(
    {
        ">FILWIDTH",
        ">FILANGLE",
        ">EQUIVLEN",
    }
)

# 1D continuous inversion summaries
INV1D_BLOCKS: frozenset[str] = frozenset(
    {
        ">RES1DXX",
        ">DEP1DXX",
        ">RES1DXY",
        ">DEP1DXY",
        ">RES1DYX",
        ">DEP1DYX",
        ">RES1DYY",
        ">DEP1DYY",
        ">FRES1DXX",
        ">FDEP1DXX",
        ">FRES1DXY",
        ">FDEP1DXY",
    }
)

# Union of all known data keywords
DATA_KEYWORDS: frozenset[str] = frozenset().union(
    Z_REAL_IMAG,
    Z_VAR_COV,
    FZ_REAL_IMAG,
    RHO_BLOCKS,
    RHO_STATS,
    PHS_BLOCKS,
    PHS_STATS,
    TIPPER_BLOCKS,
    TIPPER_AUX,
    COH_BLOCKS,
    GEOM_BLOCKS,
    SPATIAL_FILTER_BLOCKS,
    INV1D_BLOCKS,
)

# Everything we recognize as a tag/keyword in EDI files
ALL_KEYWORDS: frozenset[str] = frozenset().union(
    HEAD_TAGS,
    SECTION_OPENERS,
    MEAS_DEF_BLOCKS,
    TSERIES_BLOCKS,
    SPECTRA_BLOCKS,
    DATA_KEYWORDS,
)

# Families (category names used by parsers/writers)
FAMILY_BY_TAG: dict[str, str] = {}
for t in Z_REAL_IMAG | Z_VAR_COV | FZ_REAL_IMAG:
    FAMILY_BY_TAG[t] = "impedance"
for t in RHO_BLOCKS | RHO_STATS:
    FAMILY_BY_TAG[t] = "apparent_resistivity"
for t in PHS_BLOCKS | PHS_STATS:
    FAMILY_BY_TAG[t] = "phase"
for t in TIPPER_BLOCKS | TIPPER_AUX:
    FAMILY_BY_TAG[t] = "tipper"
for t in COH_BLOCKS:
    FAMILY_BY_TAG[t] = "coherency"
for t in GEOM_BLOCKS:
    FAMILY_BY_TAG[t] = "geometry"
for t in SPATIAL_FILTER_BLOCKS:
    FAMILY_BY_TAG[t] = "spatial_filter"
for t in INV1D_BLOCKS:
    FAMILY_BY_TAG[t] = "inversion1d"
for t in MEAS_DEF_BLOCKS:
    FAMILY_BY_TAG[t] = "measurements"
for t in TSERIES_BLOCKS:
    FAMILY_BY_TAG[t] = "time_series"
for t in SPECTRA_BLOCKS:
    FAMILY_BY_TAG[t] = "spectra"
for t in {">=MTSECT"}:
    FAMILY_BY_TAG[t] = "mt_section"
for t in {">=EMAPSECT"}:
    FAMILY_BY_TAG[t] = "emap_section"
for t in {">=TSERIESSECT"}:
    FAMILY_BY_TAG[t] = "tseries_section"
for t in {">=SPECTRASECT"}:
    FAMILY_BY_TAG[t] = "spectra_section"
for t in {">=DEFINEMEAS"}:
    FAMILY_BY_TAG[t] = "define_meas"
for t in {">HEAD", ">INFO", ">FREQ", ">ZROT", ">RHOROT", ">END"}:
    FAMILY_BY_TAG[t] = "structure"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def normalize_tag(tag: str) -> str:
    """Return canonical upper-case tag with leading '>' preserved.

    Accepts things like ``"zxxr"`` or ``">Zxxr"`` and returns
    ``">ZXXR"``. Does **not** try to normalize section openers
    (``">=..."``) beyond case.
    """
    t = tag.strip()
    if not t:
        return t
    if t[0] in ">=":
        head, rest = t[0], t[1:]
        return head + rest.upper()
    if t[0] != ">":
        return ">" + t.upper()
    return ">" + t[1:].upper()


def tag_family(tag: str) -> str | None:
    """Return family/category name for *tag* if known."""
    return FAMILY_BY_TAG.get(normalize_tag(tag))


__all__ = [
    # BNF
    "BNF_EDI_SKELETON",
    # structure
    "HEAD_TAGS",
    "SECTION_OPENERS",
    "COMMENT_TAG",
    # options
    "HEAD_ALLOWED",
    "HEAD_REQUIRED",
    "INFO_ALLOWED",
    "FREQ_ALLOWED",
    "FREQ_REQUIRED",
    "FREQ_ENUMS",
    "ZROT_ALLOWED",
    "ZROT_REQUIRED",
    "RHOROT_ALLOWED",
    "RHOROT_REQUIRED",
    "DEFINEMEAS_ALLOWED",
    "DEFINEMEAS_REQUIRED",
    "EMEAS_ALLOWED",
    "EMEAS_REQUIRED_CORE",
    "EMEAS_ENUMS",
    "HMEAS_ALLOWED",
    "HMEAS_REQUIRED_CORE",
    "HMEAS_ENUMS",
    "TSERIESSECT_ALLOWED",
    "TSERIES_ALLOWED",
    "SPECTRASECT_ALLOWED",
    "SPECTRA_ALLOWED",
    "MTSECT_ALLOWED",
    "EMAPSECT_ALLOWED",
    "OTHERSECT_ALLOWED",
    "BLOCK_OPTIONS",
    # families & keywords
    "Z_REAL_IMAG",
    "Z_VAR_COV",
    "FZ_REAL_IMAG",
    "RHO_BLOCKS",
    "RHO_STATS",
    "PHS_BLOCKS",
    "PHS_STATS",
    "TIPPER_BLOCKS",
    "TIPPER_AUX",
    "COH_BLOCKS",
    "GEOM_BLOCKS",
    "SPATIAL_FILTER_BLOCKS",
    "INV1D_BLOCKS",
    "DATA_KEYWORDS",
    "ALL_KEYWORDS",
    "FAMILY_BY_TAG",
    # helpers
    "normalize_tag",
    "tag_family",
]
