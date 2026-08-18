# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Common native-I/O dispatcher for airborne EM deliveries.

No vendor reader is registered by default. Technology adapters already accept
decoded scientific arrays; native readers should enter this registry only when
their actual delivery schema has been verified from a representative file or
an authoritative format specification.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .base import AirborneEMDataset
from .registry import (
    AirborneFormatDetectionError,
    AirborneRegistryError,
    detect_airborne_format,
    get_airborne_format,
    get_airborne_technology,
    list_airborne_formats,
)

__all__ = [
    "AirborneIOError",
    "read_airborne",
    "write_airborne",
    "available_airborne_readers",
    "available_airborne_writers",
]


class AirborneIOError(RuntimeError):
    """Raised when no defensible native airborne I/O path is available."""


def _resolved_format(
    source: Any,
    *,
    format: str | None,
    technology: str | None,
):
    if format is not None:
        definition = get_airborne_format(format)
        if definition is None:
            raise AirborneIOError(f"unknown airborne format: {format!r}")
    else:
        try:
            detected = detect_airborne_format(source)
        except AirborneFormatDetectionError as exc:
            raise AirborneIOError(str(exc)) from exc
        if detected is None:
            tech_text = ""
            if technology is not None:
                tech_text = f" for technology {technology!r}"
            raise AirborneIOError(
                "no registered native airborne format recognized the source"
                f"{tech_text}; register a verified reader after obtaining a "
                "representative delivery file"
            )
        definition = get_airborne_format(detected)
        if definition is None:  # pragma: no cover - registry invariant
            raise AirborneIOError(f"detected unknown format: {detected!r}")

    if technology is not None:
        tech = get_airborne_technology(technology)
        if tech is None:
            raise AirborneIOError(f"unknown technology: {technology!r}")
        if definition.technology != tech.name:
            raise AirborneIOError(
                f"format {definition.name!r} belongs to "
                f"{definition.technology!r}, not {tech.name!r}"
            )
    return definition


def read_airborne(
    source: Any,
    *,
    format: str | None = None,
    technology: str | None = None,
    **kwargs: Any,
) -> AirborneEMDataset:
    """Read one verified native airborne delivery into the common dataset.

    Passing an existing :class:`AirborneEMDataset` is an intentional no-op
    when no explicit native format is requested.
    """
    if isinstance(source, AirborneEMDataset) and format is None:
        return source
    definition = _resolved_format(
        source,
        format=format,
        technology=technology,
    )
    if definition.reader is None:
        raise AirborneIOError(
            f"airborne format {definition.name!r} has no registered reader"
        )
    result = definition.reader(source, **kwargs)
    if not isinstance(result, AirborneEMDataset):
        raise AirborneIOError(
            f"reader {definition.name!r} did not return AirborneEMDataset"
        )
    return result


def write_airborne(
    dataset: AirborneEMDataset,
    target: Any,
    *,
    format: str | None = None,
    technology: str | None = None,
    **kwargs: Any,
) -> Any:
    """Write a dataset through a verified native airborne writer."""
    if not isinstance(dataset, AirborneEMDataset):
        raise TypeError("dataset must be an AirborneEMDataset")
    if format is None:
        if isinstance(target, (str, Path)):
            detected = detect_airborne_format(target)
        else:
            detected = None
        if detected is None:
            raise AirborneIOError(
                "output format is required because no registered native "
                "airborne format can be inferred"
            )
        format = detected
    definition = _resolved_format(
        target,
        format=format,
        technology=technology,
    )
    if definition.writer is None:
        raise AirborneIOError(
            f"airborne format {definition.name!r} has no registered writer"
        )
    return definition.writer(dataset, target, **kwargs)


def available_airborne_readers(
    *,
    technology: str | None = None,
) -> tuple[str, ...]:
    """Return native formats with a registered reader."""
    try:
        formats = list_airborne_formats(technology=technology)
    except AirborneRegistryError as exc:
        raise AirborneIOError(str(exc)) from exc
    return tuple(fmt.name for fmt in formats if fmt.readable)


def available_airborne_writers(
    *,
    technology: str | None = None,
) -> tuple[str, ...]:
    """Return native formats with a registered writer."""
    try:
        formats = list_airborne_formats(technology=technology)
    except AirborneRegistryError as exc:
        raise AirborneIOError(str(exc)) from exc
    return tuple(fmt.name for fmt in formats if fmt.writable)
