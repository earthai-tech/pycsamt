# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""File-format validators for ModEM files.

Each ``is_*`` function accepts a path and returns ``True`` if the file
matches the expected ModEM format signature, ``False`` otherwise
(including when the path does not exist).
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

PathLike = Union[str, Path]

__all__ = [
    "is_data_file",
    "is_model_file",
    "is_model_2d_file",
    "is_model_3d_file",
    "is_covariance_file",
    "is_control_file",
    "is_log_file",
    "ModEmFileType",
    "detect_file_type",
]

# Known component type tokens that appear after '>' in data files
_DATA_COMPONENT_TYPES = frozenset({
    "FULL_IMPEDANCE",
    "OFF_DIAGONAL_IMPEDANCE",
    "DETERMINANT_IMPEDANCE",
    "TE_IMPEDANCE",
    "TM_IMPEDANCE",
    "FULL_VERTICAL_COMPONENTS",
    "PHASE_TENSOR",
})


class ModEmFileType:
    DATA        = "data"
    MODEL_2D    = "model_2d"
    MODEL_3D    = "model_3d"
    COVARIANCE  = "covariance"
    CONTROL     = "control"
    LOG         = "log"
    UNKNOWN     = "unknown"


def _head(path: PathLike, n_lines: int = 20) -> list[str]:
    """Return up to *n_lines* non-empty stripped lines from *path*."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    lines: list[str] = []
    with p.open("r", errors="replace") as fh:
        for raw in fh:
            s = raw.strip()
            if s:
                lines.append(s)
            if len(lines) >= n_lines:
                break
    return lines


def is_data_file(path: PathLike) -> bool:
    """Return True if *path* is a ModEM data file.

    Signature: a ``>`` line whose token matches a known component type,
    or the comment header contains ``Period(s) Code GG_Lat``.
    """
    try:
        lines = _head(path, 20)
    except (FileNotFoundError, OSError):
        return False
    for ln in lines:
        if ln.startswith(">"):
            token = ln[1:].strip().upper()
            if token in _DATA_COMPONENT_TYPES:
                return True
        if "PERIOD(S)" in ln.upper() and "CODE" in ln.upper() and "GG_LAT" in ln.upper():
            return True
    return False


def _parse_model_header(lines: list[str]):
    """Return (n_ints, loge_present) for first non-comment line with ≥2 ints."""
    for ln in lines:
        if ln.startswith("#"):
            continue
        parts = ln.split()
        ints = []
        for p in parts:
            try:
                ints.append(int(p))
            except ValueError:
                break
        loge = "LOGE" in ln.upper() or "LOG10" in ln.upper() or "LINEAR" in ln.upper()
        if len(ints) >= 2:
            return len(ints), loge
    return 0, False


def is_model_2d_file(path: PathLike) -> bool:
    """Return True if *path* is a ModEM 2D model file (N_x N_z LOGE header)."""
    try:
        lines = _head(path, 5)
    except (FileNotFoundError, OSError):
        return False
    n_ints, loge = _parse_model_header(lines)
    return n_ints == 2 and loge


def is_model_3d_file(path: PathLike) -> bool:
    """Return True if *path* is a ModEM 3D model file (N_x N_y N_z LOGE header)."""
    try:
        lines = _head(path, 5)
    except (FileNotFoundError, OSError):
        return False
    n_ints, loge = _parse_model_header(lines)
    return n_ints >= 3 and loge


def is_model_file(path: PathLike) -> bool:
    """Return True if *path* is any ModEM model file (2D or 3D)."""
    return is_model_2d_file(path) or is_model_3d_file(path)


def is_covariance_file(path: PathLike) -> bool:
    """Return True if *path* is a ModEM covariance file."""
    try:
        lines = _head(path, 10)
    except (FileNotFoundError, OSError):
        return False
    joined = " ".join(lines).upper()
    return (
        "MODEL COVARIANCE" in joined
        or "AUTOREGRESSION" in joined
        or ("SMOOTHING" in joined and "MASK" in joined)
    )


def is_control_file(path: PathLike) -> bool:
    """Return True if *path* is a ModEM inversion control file.

    The control file uses key : value lines; the most reliable marker is
    the output-stem line or the search-step specification.
    """
    try:
        lines = _head(path, 10)
    except (FileNotFoundError, OSError):
        return False
    joined = " ".join(lines).upper()
    return (
        "MODEL AND DATA OUTPUT FILE NAME" in joined
        or "INITIAL SEARCH STEP IN MODEL UNITS" in joined
        or ("INITIAL DAMPING FACTOR LAMBDA" in joined and ":" in joined)
    )


def is_log_file(path: PathLike) -> bool:
    """Return True if *path* is a ModEM log file."""
    try:
        lines = _head(path, 15)
    except (FileNotFoundError, OSError):
        return False
    joined = " ".join(lines).upper()
    return (
        "NLCG ITERATION" in joined
        or "COMPLETED NLCG" in joined
        or ("DAMPING PARAMETER LAMBDA" in joined and "RMS" in joined)
    )


def detect_file_type(path: PathLike) -> str:
    """Return a ``ModEmFileType`` constant for *path*, or ``'unknown'``."""
    for fn, typ in [
        (is_log_file,        ModEmFileType.LOG),
        (is_data_file,       ModEmFileType.DATA),
        (is_model_3d_file,   ModEmFileType.MODEL_3D),
        (is_model_2d_file,   ModEmFileType.MODEL_2D),
        (is_covariance_file, ModEmFileType.COVARIANCE),
        (is_control_file,    ModEmFileType.CONTROL),
    ]:
        if fn(path):
            return typ
    return ModEmFileType.UNKNOWN
