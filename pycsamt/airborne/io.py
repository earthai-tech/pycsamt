# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Common native-I/O dispatcher for airborne EM deliveries.

No vendor reader is registered by default. Technology adapters already accept
decoded scientific arrays; native readers should enter this registry only when
their actual delivery schema has been verified from a representative file or
an authoritative format specification.

Architecturally this module mirrors :mod:`pycsamt.io.formats` /
:mod:`pycsamt.io.transfer` (:func:`~pycsamt.io.formats.detect_tf_format`,
:func:`~pycsamt.io.transfer.read_transfer_function`): a
:class:`~pycsamt.airborne.registry.AirborneFormatDefinition` registry
resolved by detector-then-extension, dispatched through one stable
public entry point. The two are deliberately *not* the same registry,
because they operate on different units of work: ``pycsamt.io``
readers/writers exchange one site's transfer function (EDI/EMTF XML),
while this module exchanges a whole-survey
:class:`~pycsamt.airborne.base.AirborneEMDataset` (many flight lines,
each with many samples). :class:`AirborneIOError` is a ``RuntimeError``
rather than :class:`~pycsamt.io.formats.TransferFunctionFormatError`'s
``ValueError`` for the same reason it exists at all right now: every
failure currently reachable here is "no native reader/writer is
registered for this technology yet" -- a capability gap, not bad user
input -- because, per the project roadmap, no vendor has supplied a
delivery sample to validate a native decoder against (see
:mod:`pycsamt.airborne.mobilemt` for why that is permanent for
MobileMT specifically). Should a genuine format-detection-from-bad-
input failure mode be added later, reconsider this rather than
assuming that day is today.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .base import AirborneEMDataset
from .registry import (
    AirborneFormatDefinition,
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
    """Raised when no defensible native airborne I/O path is available.

    See the module docstring for why this is a ``RuntimeError`` rather
    than a ``ValueError``.
    """


def _resolved_format(
    source: Any,
    *,
    format: str | None,
    technology: str | None,
) -> AirborneFormatDefinition:
    """Resolve one :class:`AirborneFormatDefinition` for *source*.

    Parameters
    ----------
    source : Any
        Candidate to resolve a format for: an explicit *format* name
        skips inspecting this value entirely, otherwise it is passed
        to :func:`~pycsamt.airborne.registry.detect_airborne_format`.
    format : str, optional
        Explicit registered format name or alias. When ``None``,
        content/extension-based detection is used instead.
    technology : str, optional
        When given, cross-checked against the resolved format's own
        ``technology`` so a caller cannot silently read/write a
        MobileMT file through a ZTEM-scoped call, for example.

    Returns
    -------
    AirborneFormatDefinition
        The resolved format definition.

    Raises
    ------
    AirborneIOError
        If *format* is given but unregistered; if detection is
        ambiguous, finds nothing, or *format* is omitted; or if
        *technology* is given but does not own the resolved format.
    """
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

    Parameters
    ----------
    source : Any
        Delivery to read: typically a path, though the concrete type
        accepted depends on the registered reader. Passing an existing
        :class:`AirborneEMDataset` is an intentional no-op when *format*
        is not explicitly requested, so pipeline code can call this
        uniformly whether it already has a dataset or a raw delivery.
    format : str, optional
        Explicit registered format name or alias. When omitted,
        content/extension-based detection selects the format.
    technology : str, optional
        Restrict resolution to one technology's formats; see
        :func:`_resolved_format`.
    **kwargs
        Forwarded to the selected format's registered reader.

    Returns
    -------
    AirborneEMDataset
        The dataset produced by the resolved reader, or *source*
        itself when it already was one and *format* was omitted.

    Raises
    ------
    AirborneIOError
        If no reader can be resolved for *source* (see
        :func:`_resolved_format`), or if the resolved format has no
        registered reader, or if that reader does not return an
        :class:`AirborneEMDataset`.
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
    """Write a dataset through a verified native airborne writer.

    Parameters
    ----------
    dataset : AirborneEMDataset
        Dataset to serialize.
    target : Any
        Output destination. When *format* is omitted, *target* must be
        a path/string with an extension that resolves unambiguously
        via :func:`~pycsamt.airborne.registry.detect_airborne_format`.
    format : str, optional
        Explicit registered format name or alias.
    technology : str, optional
        Restrict resolution to one technology's formats; see
        :func:`_resolved_format`.
    **kwargs
        Forwarded to the selected format's registered writer.

    Returns
    -------
    Any
        Whatever the resolved writer returns; not constrained by this
        dispatcher.

    Raises
    ------
    TypeError
        If *dataset* is not an :class:`AirborneEMDataset`.
    AirborneIOError
        If *format* is omitted and cannot be inferred from *target*,
        if no writer can otherwise be resolved (see
        :func:`_resolved_format`), or if the resolved format has no
        registered writer.
    """
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
    """Return native formats with a registered reader.

    Parameters
    ----------
    technology : str, optional
        Restrict to one technology's formats; ``None`` lists across
        all registered technologies.

    Returns
    -------
    tuple of str
        Canonical format names currently readable. Empty until a
        native reader has been registered for at least one format;
        see the module docstring.

    Raises
    ------
    AirborneIOError
        If *technology* is given but not registered.
    """
    try:
        formats = list_airborne_formats(technology=technology)
    except AirborneRegistryError as exc:
        raise AirborneIOError(str(exc)) from exc
    return tuple(fmt.name for fmt in formats if fmt.readable)


def available_airborne_writers(
    *,
    technology: str | None = None,
) -> tuple[str, ...]:
    """Return native formats with a registered writer.

    Parameters
    ----------
    technology : str, optional
        Restrict to one technology's formats; ``None`` lists across
        all registered technologies.

    Returns
    -------
    tuple of str
        Canonical format names currently writable. Empty until a
        native writer has been registered for at least one format;
        see the module docstring.

    Raises
    ------
    AirborneIOError
        If *technology* is given but not registered.
    """
    try:
        formats = list_airborne_formats(technology=technology)
    except AirborneRegistryError as exc:
        raise AirborneIOError(str(exc)) from exc
    return tuple(fmt.name for fmt in formats if fmt.writable)
