# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""File-type detection for MARE2DEM input and output files."""

from __future__ import annotations

from enum import Enum, auto
from pathlib import Path

__all__ = [
    "Mare2DEMFileType",
    "detect_file_type",
    "is_emdata_file",
    "is_resistivity_file",
    "is_settings_file",
    "is_log_file",
    "is_response_file",
    "is_sensitivity_file",
]


class Mare2DEMFileType(Enum):
    """Enumeration of recognized MARE2DEM file types."""

    EMDATA = auto()       # .emdata              — observed data
    RESISTIVITY = auto()  # .resistivity          — mesh + model
    POLY = auto()         # .poly                 — Triangle polygon mesh
    SETTINGS = auto()     # .settings             — inversion settings
    LOG = auto()          # .log / .logfile       — run log
    RESPONSE = auto()     # *.resp / *_MARE2DEM.emdata — predicted data
    SENSITIVITY = auto()  # .sensitivity          — sensitivity file
    UNKNOWN = auto()


def detect_file_type(path: str | Path) -> Mare2DEMFileType:
    """Detect the MARE2DEM file type from the file extension and name.

    Parameters
    ----------
    path : path-like
        Path to the file to classify.

    Returns
    -------
    Mare2DEMFileType
        Detected file type, or :attr:`Mare2DEMFileType.UNKNOWN`.
    """
    p = Path(path)
    stem = p.stem.lower()
    suffix = p.suffix.lower()

    if suffix == ".emdata":
        if stem.endswith("_mare2dem"):
            return Mare2DEMFileType.RESPONSE
        return Mare2DEMFileType.EMDATA
    if suffix == ".resp":
        return Mare2DEMFileType.RESPONSE
    if suffix == ".resistivity":
        return Mare2DEMFileType.RESISTIVITY
    if suffix == ".poly":
        return Mare2DEMFileType.POLY
    if suffix == ".settings":
        return Mare2DEMFileType.SETTINGS
    if suffix in (".log", ".logfile"):
        return Mare2DEMFileType.LOG
    if suffix == ".sensitivity":
        return Mare2DEMFileType.SENSITIVITY
    return Mare2DEMFileType.UNKNOWN


def is_emdata_file(path: str | Path) -> bool:
    """Return ``True`` when *path* is a MARE2DEM observed-data file."""
    return detect_file_type(path) == Mare2DEMFileType.EMDATA


def is_resistivity_file(path: str | Path) -> bool:
    """Return ``True`` when *path* is a MARE2DEM resistivity model file."""
    return detect_file_type(path) == Mare2DEMFileType.RESISTIVITY


def is_settings_file(path: str | Path) -> bool:
    """Return ``True`` when *path* is a MARE2DEM settings file."""
    return detect_file_type(path) == Mare2DEMFileType.SETTINGS


def is_log_file(path: str | Path) -> bool:
    """Return ``True`` when *path* is a MARE2DEM log file."""
    return detect_file_type(path) == Mare2DEMFileType.LOG


def is_response_file(path: str | Path) -> bool:
    """Return ``True`` when *path* is a MARE2DEM predicted-response file."""
    return detect_file_type(path) == Mare2DEMFileType.RESPONSE


def is_sensitivity_file(path: str | Path) -> bool:
    """Return ``True`` when *path* is a MARE2DEM sensitivity file."""
    return detect_file_type(path) == Mare2DEMFileType.SENSITIVITY
