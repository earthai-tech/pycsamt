# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared semver parsing/gating for PCSF's two on-disk encodings.

:mod:`pycsamt.format.io` (HDF5, ``.pcsf``) and
:mod:`pycsamt.format.text` (ASCII, ``.pcsm``) serialize the exact same
:class:`~pycsamt.format.schema.PCSFModel`, so they share one
versioning policy (``pycsamt/format/SPEC.md`` section 5) and one
implementation of it rather than two copies that could drift apart.
"""

from __future__ import annotations

import warnings

__all__ = ["parse_pcsf_version", "check_pcsf_version"]


def parse_pcsf_version(text: str) -> tuple[int, int, int]:
    parts = text.split(".")
    if len(parts) != 3:
        raise ValueError(
            f"malformed pcsf_version {text!r}; expected MAJOR.MINOR.PATCH"
        )
    try:
        major, minor, patch = (int(part) for part in parts)
    except ValueError as exc:
        raise ValueError(
            f"malformed pcsf_version {text!r}; expected MAJOR.MINOR.PATCH"
        ) from exc
    return major, minor, patch


def check_pcsf_version(file_version: str, reader_version: str) -> None:
    """Enforce the versioning policy in ``pycsamt/format/SPEC.md`` (S5):
    reject an unrecognised MAJOR version, warn on a newer MINOR one.
    """
    current_major, current_minor, _ = parse_pcsf_version(reader_version)
    file_major, file_minor, _ = parse_pcsf_version(file_version)
    if file_major != current_major:
        raise ValueError(
            f"unsupported pcsf_version {file_version!r}: this reader only "
            f"recognises major version {current_major} (reader "
            f"PCSF_VERSION={reader_version!r}); see pycsamt/format/SPEC.md "
            "section 5 for the versioning policy"
        )
    if file_minor > current_minor:
        warnings.warn(
            f"file was written with pcsf_version {file_version!r}, newer "
            f"than this reader's {reader_version!r}; fields added since "
            "then will be ignored",
            stacklevel=2,
        )
