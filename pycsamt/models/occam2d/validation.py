# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""File-format validators for Occam2D files.

Each ``is_*`` function accepts a path and returns ``True`` if the file
matches the expected Occam2D format signature, ``False`` otherwise.
Raises ``ValueError`` when the path does not exist.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

PathLike = Union[str, Path]

__all__ = [
    "is_data_file",
    "is_mesh_file",
    "is_model_file",
    "is_startup_file",
    "is_iter_file",
    "is_response_file",
    "is_log_file",
    "OccamFileType",
    "detect_file_type",
]

# Magic strings that identify each file type (first non-blank line prefix)
_SIGNATURES: dict[str, str] = {
    "data":    "FORMAT:",       # OCCAM2MTDATA_1.0
    "model":   "FORMAT:",       # OCCAM2MTMOD_1.0  (checked with second token)
    "startup": "FORMAT:",       # OCCAMITER_FLEX
    "iter":    "FORMAT:",       # OCCAMITER_FLEX  + Iteration: > 0
    "mesh":    "MESH FILE",     # first line of Occam2DMesh
    "response":"FORMAT:",       # similar to data
    "log":     "OCCAM ITERATION",
}


class OccamFileType:
    DATA     = "data"
    MESH     = "mesh"
    MODEL    = "model"
    STARTUP  = "startup"
    ITER     = "iter"
    RESPONSE = "response"
    LOG      = "log"
    UNKNOWN  = "unknown"


def _first_lines(path: PathLike, n: int = 5) -> list[str]:
    p = Path(path)
    if not p.exists():
        raise ValueError(f"File not found: {p}")
    lines = []
    with p.open("r", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped:
                lines.append(stripped.upper())
            if len(lines) >= n:
                break
    return lines


def is_data_file(path: PathLike) -> bool:
    """Return True if *path* is an Occam2DMT data file."""
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    return any("OCCAM2MTDATA" in ln for ln in lines)


def is_mesh_file(path: PathLike) -> bool:
    """Return True if *path* is an Occam2DMesh file."""
    try:
        lines = _first_lines(path, 2)
    except ValueError:
        return False
    # Mesh files begin with two integers (nx_cols, nz_rows) or a comment
    if lines:
        tokens = lines[0].split()
        if len(tokens) >= 2 and all(t.lstrip("-").isdigit() for t in tokens[:2]):
            return True
    return False


def is_model_file(path: PathLike) -> bool:
    """Return True if *path* is an Occam2DModel file."""
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    return any("OCCAM2MTMOD" in ln for ln in lines)


def is_startup_file(path: PathLike) -> bool:
    """Return True if *path* is a Startup file (iteration == 0)."""
    try:
        lines = _first_lines(path, 15)
    except ValueError:
        return False
    fmt = any("OCCAMITER_FLEX" in ln for ln in lines)
    if not fmt:
        return False
    for ln in lines:
        if ln.startswith("ITERATION:") and ln.split(":")[1].strip() == "0":
            return True
    return False


def is_iter_file(path: PathLike) -> bool:
    """Return True if *path* is an Occam .iter file (iteration > 0)."""
    try:
        lines = _first_lines(path, 15)
    except ValueError:
        return False
    fmt = any("OCCAMITER_FLEX" in ln for ln in lines)
    if not fmt:
        return False
    for ln in lines:
        if ln.startswith("ITERATION:"):
            try:
                return int(ln.split(":")[1].strip()) > 0
            except ValueError:
                return False
    return False


def is_response_file(path: PathLike) -> bool:
    """Return True if *path* is an Occam response (.resp) file."""
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    return any("OCCAM2MTRESP" in ln or "RESPONSE" in ln for ln in lines)


def is_log_file(path: PathLike) -> bool:
    """Return True if *path* is an Occam log file."""
    try:
        lines = _first_lines(path, 5)
    except ValueError:
        return False
    return any("OCCAM ITERATION" in ln or "FMINOCC" in ln for ln in lines)


def detect_file_type(path: PathLike) -> str:
    """Return an ``OccamFileType`` constant for *path*, or ``'unknown'``."""
    for fn, typ in [
        (is_data_file,     OccamFileType.DATA),
        (is_model_file,    OccamFileType.MODEL),
        (is_startup_file,  OccamFileType.STARTUP),
        (is_iter_file,     OccamFileType.ITER),
        (is_mesh_file,     OccamFileType.MESH),
        (is_response_file, OccamFileType.RESPONSE),
        (is_log_file,      OccamFileType.LOG),
    ]:
        if fn(path):
            return typ
    return OccamFileType.UNKNOWN
