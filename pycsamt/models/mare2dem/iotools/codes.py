# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM data-type code lookup table.

Port of ``m2d_getDataCodeLookupTable.m``.

Each integer code stored in the DATA block of a ``.emdata`` file maps
to a physical quantity (component name) and a representation (e.g. Real,
Imaginary, Amplitude, Phase, ApRes, log10(ApRes)).  Use
:func:`code_label` for a human-readable string and :func:`is_mt_code`
/ :func:`is_csem_code` to classify a code.
"""

from __future__ import annotations

__all__ = [
    "DATA_CODES",
    "code_label",
    "code_component",
    "code_representation",
    "is_mt_code",
    "is_csem_code",
]

# ---------------------------------------------------------------------------
# Code table: {int_code: (component, representation)}
# Ported directly from getDataCodeLookupTable.m
# ---------------------------------------------------------------------------

DATA_CODES: dict[int, tuple[str, str]] = {
    # ---- CSEM E-field (real / imaginary) ----
    1:  ("Ex",       "Real"),
    2:  ("Ex",       "Imaginary"),
    3:  ("Ey",       "Real"),
    4:  ("Ey",       "Imaginary"),
    5:  ("Ez",       "Real"),
    6:  ("Ez",       "Imaginary"),
    # ---- CSEM B-field (real / imaginary) ----
    11: ("Bx",       "Real"),
    12: ("Bx",       "Imaginary"),
    13: ("By",       "Real"),
    14: ("By",       "Imaginary"),
    15: ("Bz",       "Real"),
    16: ("Bz",       "Imaginary"),
    # ---- CSEM E-field (amplitude / phase / log10 amplitude) ----
    21: ("Ex",       "Amplitude"),
    22: ("Ex",       "Phase"),
    23: ("Ey",       "Amplitude"),
    24: ("Ey",       "Phase"),
    25: ("Ez",       "Amplitude"),
    26: ("Ez",       "Phase"),
    27: ("Ex",       "Log10 Amplitude"),
    28: ("Ey",       "Log10 Amplitude"),
    29: ("Ez",       "Log10 Amplitude"),
    # ---- CSEM B-field (amplitude / phase / log10 amplitude) ----
    31: ("Bx",       "Amplitude"),
    32: ("Bx",       "Phase"),
    33: ("By",       "Amplitude"),
    34: ("By",       "Phase"),
    35: ("Bz",       "Amplitude"),
    36: ("Bz",       "Phase"),
    37: ("Bx",       "Log10 Amplitude"),
    38: ("By",       "Log10 Amplitude"),
    39: ("Bz",       "Log10 Amplitude"),
    # ---- CSEM polarization ellipse ----
    41: ("Ep",       "PEmax"),
    42: ("Ep",       "PEmin"),
    43: ("Bp",       "PEmax"),
    44: ("Bp",       "PEmin"),
    # ---- MT impedance (apparent resistivity / phase) ----
    103: ("Zxy (TE)",  "ApRes"),
    104: ("Zxy (TE)",  "Phase"),
    105: ("Zyx (TM)",  "ApRes"),
    106: ("Zyx (TM)",  "Phase"),
    # ---- MT impedance (real / imaginary) ----
    113: ("Zxy (TE)",  "Real"),
    114: ("Zxy (TE)",  "Imaginary"),
    115: ("Zyx (TM)",  "Real"),
    116: ("Zyx (TM)",  "Imaginary"),
    # ---- MT log10 apparent resistivity ----
    123: ("Zxy (TE)",  "log10(ApRes)"),
    125: ("Zyx (TM)",  "log10(ApRes)"),
    # ---- MT determinant ----
    109: ("Det |Z|",   "ApRes"),
    110: ("Det |Z|",   "Phase"),
    129: ("Det |Z|",   "log10(ApRes)"),
    # ---- MT tipper ----
    133: ("Mzy (TE)",  "Real"),
    134: ("Mzy (TE)",  "Imaginary"),
    135: ("Mzy (TE)",  "Amplitude"),
    136: ("Mzy (TE)",  "Phase"),
    # ---- MT EM fields (real / imaginary) ----
    151: ("Ex",        "Real"),
    152: ("Ex",        "Imaginary"),
    153: ("Ey",        "Real"),
    154: ("Ey",        "Imaginary"),
    155: ("Ez",        "Real"),
    156: ("Ez",        "Imaginary"),
    161: ("Hx",        "Real"),
    162: ("Hx",        "Imaginary"),
    163: ("Hy",        "Real"),
    164: ("Hy",        "Imaginary"),
    165: ("Hz",        "Real"),
    166: ("Hz",        "Imaginary"),
}

# Integer code sets for fast classification
_CSEM_CODES: frozenset[int] = frozenset(range(1, 50)) | frozenset(range(151, 167))
_MT_CODES: frozenset[int] = frozenset(range(100, 140))


def code_label(code: int) -> str:
    """Return ``"Component — Representation"`` for an integer data code.

    Parameters
    ----------
    code : int
        Integer data-type code from the DATA block of a ``.emdata`` file.

    Returns
    -------
    str
        Human-readable label, e.g. ``"Zxy (TE) — Phase"`` or
        ``"Unknown (code=7)"``.

    Examples
    --------
    >>> code_label(104)
    'Zxy (TE) — Phase'
    >>> code_label(27)
    'Ex — Log10 Amplitude'
    """
    if code in DATA_CODES:
        comp, rep = DATA_CODES[code]
        return f"{comp} — {rep}"
    return f"Unknown (code={code})"


def code_component(code: int) -> str:
    """Return the EM component name for *code*, or ``""`` if unknown."""
    return DATA_CODES[code][0] if code in DATA_CODES else ""


def code_representation(code: int) -> str:
    """Return the representation name for *code*, or ``""`` if unknown."""
    return DATA_CODES[code][1] if code in DATA_CODES else ""


def is_mt_code(code: int) -> bool:
    """Return ``True`` when *code* belongs to the MT data type range."""
    return code in _MT_CODES


def is_csem_code(code: int) -> bool:
    """Return ``True`` when *code* belongs to the CSEM data type range."""
    return code in _CSEM_CODES
