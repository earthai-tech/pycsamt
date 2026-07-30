# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""File-format validators for Occam2D files.

Each ``is_*`` function accepts a path and returns ``True`` if the file
matches the expected Occam2D format signature, ``False`` otherwise
(including when the path does not exist).

The validators are intentionally lightweight. They read only the first
few non-empty lines and look for stable Occam2D signatures such as
``OCCAM2MTDATA_1.0``, ``OCCAM2MTMOD_1.0``, ``OCCAMITER_FLEX``, mesh
control lines, or log iteration markers. They are useful when scanning
an inversion directory before choosing the appropriate parser.
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
    "data": "FORMAT:",  # OCCAM2MTDATA_1.0
    "model": "FORMAT:",  # OCCAM2MTMOD_1.0  (checked with second token)
    "startup": "FORMAT:",  # OCCAMITER_FLEX
    "iter": "FORMAT:",  # OCCAMITER_FLEX  + Iteration: > 0
    "mesh": "MESH FILE",  # first line of Occam2DMesh
    "response": "FORMAT:",  # similar to data
    "log": "OCCAM ITERATION",
}


class OccamFileType:
    """String constants returned by :func:`detect_file_type`.

    Attributes
    ----------
    DATA, MESH, MODEL, STARTUP, ITER, RESPONSE, LOG : str
        Known Occam2D file categories.
    UNKNOWN : str
        Returned when no validator recognizes a file.
    """

    DATA = "data"
    MESH = "mesh"
    MODEL = "model"
    STARTUP = "startup"
    ITER = "iter"
    RESPONSE = "response"
    LOG = "log"
    UNKNOWN = "unknown"


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
    """Return whether ``path`` is an Occam2D data file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when an early non-empty line contains the
        ``OCCAM2MTDATA`` format tag.
    """
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    return any("OCCAM2MTDATA" in ln for ln in lines)


def is_mesh_file(path: PathLike) -> bool:
    """Return whether ``path`` is an Occam2D mesh file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when the file looks like an ``Occam2DMesh``
        file, either by a ``MESH FILE`` comment or by an early
        numeric control line.
    """
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    for ln in lines:
        # Comment line always starts with MESH FILE
        if "MESH FILE" in ln:
            return True
        # Control line: first two tokens are integers (could be line 1 or 2)
        tokens = ln.split()
        if len(tokens) >= 2 and all(
            t.lstrip("-").isdigit() for t in tokens[:2]
        ):
            return True
    return False


def is_model_file(path: PathLike) -> bool:
    """Return whether ``path`` is an Occam2D model file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when an early non-empty line contains the
        ``OCCAM2MTMOD`` model-format tag.
    """
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    return any("OCCAM2MTMOD" in ln for ln in lines)


def is_startup_file(path: PathLike) -> bool:
    """Return whether ``path`` is an Occam2D startup file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when the file has the ``OCCAMITER_FLEX`` tag
        and an ``Iteration: 0`` header.
    """
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
    """Return whether ``path`` is an Occam2D iteration file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when the file has the ``OCCAMITER_FLEX`` tag
        and an ``Iteration`` value greater than zero.
    """
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
    """Return whether ``path`` is an Occam2D response file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when the early lines contain a response
        marker such as ``OCCAM2MTRESP`` or ``RESPONSE``.
    """
    try:
        lines = _first_lines(path, 3)
    except ValueError:
        return False
    return any("OCCAM2MTRESP" in ln or "RESPONSE" in ln for ln in lines)


def is_log_file(path: PathLike) -> bool:
    """Return whether ``path`` is an Occam2D log file.

    Parameters
    ----------
    path : path-like
        Candidate file path. Missing files return ``False``.

    Returns
    -------
    bool
        ``True`` when early lines contain common Occam log
        markers such as ``** ITERATION`` or ``FMINOCC``.
    """
    try:
        lines = _first_lines(path, 5)
    except ValueError:
        return False
    return any(
        "OCCAM ITERATION" in ln
        or "** ITERATION" in ln
        or "FMINOCC" in ln
        or "STARTING R.M.S." in ln
        for ln in lines
    )


def detect_file_type(path: PathLike) -> str:
    """Detect the Occam2D file type represented by ``path``.

    The detector applies the public validators in a fixed
    order and returns a string constant from
    :class:`OccamFileType`. Missing or unrecognized files
    return ``OccamFileType.UNKNOWN``.

    Parameters
    ----------
    path : path-like
        Candidate file path. The value may be a string,
        :class:`pathlib.Path`, or any object accepted by
        :class:`pathlib.Path`.

    Returns
    -------
    str
        One of ``"data"``, ``"mesh"``, ``"model"``,
        ``"startup"``, ``"iter"``, ``"response"``, ``"log"``,
        or ``"unknown"``.

    See Also
    --------
    is_data_file
        Recognizes Occam data files.
    is_startup_file
        Distinguishes startup files from non-zero iterations.
    is_iter_file
        Recognizes non-zero Occam iteration files.

    Examples
    --------
    >>> from pycsamt.models.occam2d.validation import (
    ...     detect_file_type,
    ... )
    >>> detect_file_type("occam_run/Occam2DMesh")
    'mesh'
    """
    for fn, typ in [
        (is_data_file, OccamFileType.DATA),
        (is_model_file, OccamFileType.MODEL),
        (is_startup_file, OccamFileType.STARTUP),
        (is_iter_file, OccamFileType.ITER),
        (is_mesh_file, OccamFileType.MESH),
        (is_response_file, OccamFileType.RESPONSE),
        (is_log_file, OccamFileType.LOG),
    ]:
        if fn(path):
            return typ
    return OccamFileType.UNKNOWN
