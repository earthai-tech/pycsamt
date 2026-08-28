# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Semantic-version parsing and compatibility checks for PCBH."""

from __future__ import annotations

import warnings

__all__ = ["parse_pcbh_version", "check_pcbh_version"]


def parse_pcbh_version(text: str) -> tuple[int, int, int]:
    """Parse a PCBH ``MAJOR.MINOR.PATCH`` version.

    Parameters
    ----------
    text : str
        Version text to parse.

    Returns
    -------
    tuple of int
        The major, minor, and patch components.

    Raises
    ------
    ValueError
        If *text* is not a non-negative semantic version.
    """
    if not isinstance(text, str):
        raise ValueError("pcbh_version must be a string")
    parts = text.split(".")
    if len(parts) != 3:
        raise ValueError(
            f"malformed pcbh_version {text!r}; expected MAJOR.MINOR.PATCH"
        )
    try:
        version = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise ValueError(
            f"malformed pcbh_version {text!r}; expected MAJOR.MINOR.PATCH"
        ) from exc
    if any(part < 0 for part in version):
        raise ValueError(
            f"malformed pcbh_version {text!r}; components must be >= 0"
        )
    return version  # type: ignore[return-value]


def check_pcbh_version(file_version: str, reader_version: str) -> None:
    """Enforce PCBH major/minor reader compatibility.

    Parameters
    ----------
    file_version, reader_version : str
        Semantic versions for the input document and reader.

    Raises
    ------
    ValueError
        If the major versions differ or either version is malformed.

    Warns
    -----
    UserWarning
        If the file uses a newer minor version than the reader.
    """
    current_major, current_minor, _ = parse_pcbh_version(reader_version)
    file_major, file_minor, _ = parse_pcbh_version(file_version)
    if file_major != current_major:
        raise ValueError(
            f"unsupported pcbh_version {file_version!r}: this reader "
            f"supports major version {current_major}"
        )
    if file_minor > current_minor:
        warnings.warn(
            f"PCBH file version {file_version!r} is newer than reader "
            f"version {reader_version!r}; unknown optional fields may "
            "be preserved but not interpreted",
            UserWarning,
            stacklevel=2,
        )
