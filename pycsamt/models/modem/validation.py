# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""File-format validators for ModEM files.

This module provides lightweight predicates for recognizing
common ModEM input and output files before a full parser is
called. The checks are intentionally signature based: they
read only the first few non-empty lines and look for stable
tokens in the header. They are therefore fast and safe for
directory scans, but they do not replace structural parsing by
the dedicated data, model, covariance, control, or log classes.

Examples
--------
>>> from pycsamt.models.modem.validation import detect_file_type
>>> from pycsamt.models.modem.validation import ModEmFileType
>>> detect_file_type("missing.dat") == ModEmFileType.UNKNOWN
True

See Also
--------
pycsamt.models.modem.data.ModEmData
    Full reader and writer for ModEM data files.
pycsamt.models.modem.model2d.ModEmModel2D
    Full reader and writer for 2-D model files.
pycsamt.models.modem.model3d.ModEmModel3D
    Full reader and writer for 3-D model files.
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
_DATA_COMPONENT_TYPES = frozenset(
    {
        "FULL_IMPEDANCE",
        "OFF_DIAGONAL_IMPEDANCE",
        "DETERMINANT_IMPEDANCE",
        "TE_IMPEDANCE",
        "TM_IMPEDANCE",
        "FULL_VERTICAL_COMPONENTS",
        "PHASE_TENSOR",
    }
)


class ModEmFileType:
    """String constants returned by :func:`detect_file_type`."""

    DATA = "data"
    MODEL_2D = "model_2d"
    MODEL_3D = "model_3d"
    COVARIANCE = "covariance"
    CONTROL = "control"
    LOG = "log"
    UNKNOWN = "unknown"


def _head(path: PathLike, n_lines: int = 20) -> list[str]:
    """Return up to ``n_lines`` non-empty stripped file lines."""
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
    """Check whether a file looks like a ModEM data file.

    ModEM data files describe observed or predicted
    electromagnetic responses for a set of stations and
    periods. The validator uses header signatures rather than
    a full parse. A file is accepted when its leading lines
    contain either a known component declaration such as
    ``> Full_Impedance`` or a tabular header containing
    ``Period(s)``, ``Code``, and ``GG_Lat``.

    Parameters
    ----------
    path : path-like
        Candidate file path. The value may be a string or
        :class:`pathlib.Path`. Missing files and unreadable
        files are treated as non-matches.

    Returns
    -------
    bool
        ``True`` when the leading file content matches a
        recognized ModEM data-file signature, otherwise
        ``False``.

    Examples
    --------
    >>> from pycsamt.models.modem.validation import is_data_file
    >>> is_data_file("d0.dat")
    False

    See Also
    --------
    detect_file_type
        Return the most likely ModEM file category.
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
        marker = (
            "PERIOD(S)" in ln.upper()
            and "CODE" in ln.upper()
            and "GG_LAT" in ln.upper()
        )
        if marker:
            return True
    return False


def _parse_model_header(lines: list[str]):
    """Return integer count and encoding marker from model header."""
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
    """Check whether a file looks like a ModEM 2-D model.

    A 2-D ModEM model header begins with two integer cell
    counts, commonly :math:`N_x` and :math:`N_z`, followed by a
    resistivity encoding token such as ``LOGE``, ``LOG10``, or
    ``LINEAR``. This predicate checks only that leading header
    signature; it does not validate the full number of mesh or
    resistivity values.

    Parameters
    ----------
    path : path-like
        Candidate model file. Missing or unreadable files
        return ``False``.

    Returns
    -------
    bool
        ``True`` if the file has a 2-D ModEM model header,
        otherwise ``False``.

    Examples
    --------
    >>> from pycsamt.models.modem.validation import is_model_2d_file
    >>> is_model_2d_file("m0.rho")
    False
    """
    try:
        lines = _head(path, 5)
    except (FileNotFoundError, OSError):
        return False
    n_ints, loge = _parse_model_header(lines)
    return n_ints == 2 and loge


def is_model_3d_file(path: PathLike) -> bool:
    """Check whether a file looks like a ModEM 3-D model.

    A 3-D ModEM model header begins with at least three
    integer cell counts, commonly :math:`N_x`, :math:`N_y`,
    and :math:`N_z`, plus a resistivity encoding token. The
    validator is intentionally lightweight and reads only the
    first few non-empty lines.

    Parameters
    ----------
    path : path-like
        Candidate model file. Missing or unreadable files
        return ``False``.

    Returns
    -------
    bool
        ``True`` if the file has a 3-D ModEM model header,
        otherwise ``False``.

    Examples
    --------
    >>> from pycsamt.models.modem.validation import is_model_3d_file
    >>> is_model_3d_file("m0.ws")
    False
    """
    try:
        lines = _head(path, 5)
    except (FileNotFoundError, OSError):
        return False
    n_ints, loge = _parse_model_header(lines)
    return n_ints >= 3 and loge


def is_model_file(path: PathLike) -> bool:
    """Check whether a file looks like any ModEM model file.

    Parameters
    ----------
    path : path-like
        Candidate 2-D or 3-D model file.

    Returns
    -------
    bool
        ``True`` when either :func:`is_model_2d_file` or
        :func:`is_model_3d_file` accepts the file.
    """
    return is_model_2d_file(path) or is_model_3d_file(path)


def is_covariance_file(path: PathLike) -> bool:
    """Check whether a file looks like a ModEM covariance file.

    Covariance files describe regularization smoothing and
    active-cell masks for 3-D inversions. This validator looks
    for typical header words such as ``Model Covariance``,
    ``Autoregression``, or the combination of ``Smoothing``
    and ``Mask``.

    Parameters
    ----------
    path : path-like
        Candidate covariance file. Missing or unreadable files
        return ``False``.

    Returns
    -------
    bool
        ``True`` when the leading lines contain a covariance
        signature, otherwise ``False``.
    """
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
    r"""Check whether a file looks like a ModEM control file.

    ModEM control files are key-value text files that define
    nonlinear inversion settings. The most reliable signatures
    are the output-stem line, the initial search-step line, or
    the initial damping parameter :math:`\lambda` line.

    Parameters
    ----------
    path : path-like
        Candidate inversion-control file. Missing or unreadable
        files return ``False``.

    Returns
    -------
    bool
        ``True`` when the header contains recognized control
        parameters, otherwise ``False``.
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
    """Check whether a file looks like a ModEM run log.

    Run logs contain nonlinear conjugate-gradient progress
    records, including RMS values and damping parameters. The
    predicate checks for common strings written by ModEM logs,
    such as ``NLCG iteration``, ``Completed NLCG``, or the
    combination of ``Damping parameter lambda`` and ``RMS``.

    Parameters
    ----------
    path : path-like
        Candidate log file. Missing or unreadable files return
        ``False``.

    Returns
    -------
    bool
        ``True`` when the file has a recognizable ModEM log
        signature, otherwise ``False``.
    """
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
    """Detect the most likely ModEM file category.

    The detector applies the public validators in a stable
    order and returns the first matching
    :class:`ModEmFileType` constant. Log and data signatures
    are checked before model signatures, followed by
    covariance and control files. If no predicate matches, the
    file is reported as ``ModEmFileType.UNKNOWN``.

    Parameters
    ----------
    path : path-like
        Candidate ModEM file. Missing files are allowed and
        return ``ModEmFileType.UNKNOWN``.

    Returns
    -------
    str
        One of ``"data"``, ``"model_2d"``, ``"model_3d"``,
        ``"covariance"``, ``"control"``, ``"log"``, or
        ``"unknown"``.

    Examples
    --------
    >>> from pycsamt.models.modem.validation import detect_file_type
    >>> from pycsamt.models.modem.validation import ModEmFileType
    >>> detect_file_type("missing.dat") == ModEmFileType.UNKNOWN
    True

    See Also
    --------
    is_data_file
        Validate ModEM observed or predicted data files.
    is_model_file
        Validate either 2-D or 3-D ModEM model files.
    is_control_file
        Validate ModEM inversion-control files.
    """
    for fn, typ in [
        (is_log_file, ModEmFileType.LOG),
        (is_data_file, ModEmFileType.DATA),
        (is_model_3d_file, ModEmFileType.MODEL_3D),
        (is_model_2d_file, ModEmFileType.MODEL_2D),
        (is_covariance_file, ModEmFileType.COVARIANCE),
        (is_control_file, ModEmFileType.CONTROL),
    ]:
        if fn(path):
            return typ
    return ModEmFileType.UNKNOWN
