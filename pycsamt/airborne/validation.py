# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Shared parameter validation and normalization for ``pycsamt.airborne``.

This module centralizes the boundary-validation logic that would otherwise
be reimplemented independently by :mod:`pycsamt.airborne.base`,
:mod:`pycsamt.airborne.navigation`, and each technology adapter
(:mod:`pycsamt.airborne.mobilemt`, :mod:`pycsamt.airborne.ztem`,
:mod:`pycsamt.airborne.afmag`). Every helper accepts an ``error_cls``
keyword so a technology adapter can keep raising its own public exception
type (for example ``MobileMTValidationError``) while sharing one
implementation, instead of every adapter re-declaring an equivalent
private helper under a different name.

Structural, type, and shape validation is the majority of this module.
Scientific invariants that are specific to one technology (for example
the MobileMT 3x2 admittance matrix shape) remain the responsibility of
the owning adapter module, consistent with the pyCSAMT convention that
public methods validate near the API boundary while private helpers
may assume already-validated internal state. A small number of helpers
(:func:`reference_station_mapping`, :func:`merge_remote_reference_processing`)
go slightly beyond pure validation into shared EMTF-metadata assembly,
because that assembly logic was independently duplicated across
technology adapters just as much as the shape checks were -- keeping
it out of this module would not have made the package simpler, only
moved the duplication somewhere less discoverable.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from types import MappingProxyType
from typing import Any

import numpy as np

__all__ = [
    "emtf_class",
    "clean_identifier",
    "normalize_optional_identifier",
    "normalize_positive_float",
    "normalize_count_range",
    "normalize_numeric_vector",
    "normalize_object_vector",
    "normalize_frequency",
    "resolve_frequency_or_periods",
    "resolve_line_frequency_grid",
    "normalize_frequency_range",
    "normalize_fixed_channels",
    "normalize_estimate_array",
    "normalize_sample_axis_array",
    "normalize_record_mask",
    "reference_station_mapping",
    "merge_remote_reference_processing",
]


def emtf_class():
    """Return :class:`pycsamt.emtf.EMTF`, imported lazily.

    Both :mod:`pycsamt.airborne.base` and :mod:`pycsamt.airborne.qc`
    need ``EMTF`` only for ``isinstance`` checks. Importing it lazily
    here, once, avoids a hard import-time dependency on
    :mod:`pycsamt.emtf` from either module while keeping the check
    itself in one place instead of duplicated per module.
    """
    from ..emtf.document import EMTF

    return EMTF


def clean_identifier(
    value: Any,
    *,
    name: str,
    error_cls: type[Exception] = ValueError,
) -> str:
    """Return a stripped, non-empty string identifier.

    Parameters
    ----------
    value : Any
        Candidate identifier. ``str(value)`` is used before stripping.
    name : str
        Parameter name used in the error message.
    error_cls : type, default ValueError
        Exception type raised when *value* strips to an empty string.

    Returns
    -------
    str
        The stripped identifier.

    Raises
    ------
    error_cls
        If *value* is empty after stripping.
    """
    text = str(value).strip()
    if not text:
        raise error_cls(f"{name} must be non-empty")
    return text


def normalize_optional_identifier(value: Any | None) -> str | None:
    """Return a stripped identifier, or ``None`` if absent/blank.

    Unlike :func:`clean_identifier`, an empty result is not an error:
    this is for genuinely optional identifiers (for example a
    reference station's ``station_id`` before it falls back to a
    :class:`~pycsamt.metadata.SiteMeta` name), where "not supplied"
    and "supplied but blank" should collapse to the same ``None``
    rather than one being valid and the other raising.
    """
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def normalize_numeric_vector(
    value: Any | None,
    *,
    name: str,
    size: int,
    error_cls: type[Exception] = ValueError,
) -> np.ndarray | None:
    """Return a validated 1-D float vector aligned to *size*, or ``None``.

    Parameters
    ----------
    value : array-like or None
        Candidate numeric vector. ``None`` passes through unchanged
        rather than being treated as a physical zero.
    name : str
        Parameter name used in error messages.
    size : int
        Required vector length.
    error_cls : type, default ValueError
        Exception type raised on shape/finiteness violations.

    Returns
    -------
    ndarray of shape (size,), or None
        ``NaN`` values are permitted: an individual missing sample is
        represented by ``nan``, never by a fabricated zero.

    Raises
    ------
    error_cls
        If *value* is not 1-D, its length does not match *size*, or it
        contains an infinite value.
    """
    if value is None:
        return None
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise error_cls(f"{name} must be a 1-D array")
    if arr.size != size:
        raise error_cls(
            f"{name} length must match sample_ids: {arr.size} != {size}"
        )
    if np.any(np.isinf(arr)):
        raise error_cls(f"{name} must not contain infinite values")
    return arr


def normalize_object_vector(
    value: Any | None,
    *,
    name: str,
    size: int,
    error_cls: type[Exception] = ValueError,
) -> tuple[Any, ...] | None:
    """Return a validated length-*size* tuple of opaque values, or ``None``.

    Used for sample-aligned fields, such as timestamps, whose element
    type is not itself numeric or geophysical.
    """
    if value is None:
        return None
    out = tuple(value)
    if len(out) != size:
        raise error_cls(
            f"{name} length must match sample_ids: {len(out)} != {size}"
        )
    return out


def normalize_frequency(
    value: Any,
    *,
    name: str = "frequency",
    error_cls: type[Exception] = ValueError,
) -> np.ndarray:
    """Return a validated 1-D, finite, strictly positive vector.

    Parameters
    ----------
    value : array-like
        Candidate frequency (Hz) or period (s) values. A scalar is
        promoted to a length-1 array.
    name : str, default "frequency"
        Parameter name used in error messages; pass ``"periods"`` when
        validating a period axis instead of a frequency axis.
    error_cls : type, default ValueError
        Exception type raised on shape/positivity violations.

    Returns
    -------
    ndarray of shape (n,)
        Finite, strictly positive values.

    Raises
    ------
    error_cls
        If *value* is not 1-D, is empty, or contains a non-finite or
        non-positive value.
    """
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise error_cls(f"{name} must be a 1-D array")
    if arr.size == 0:
        raise error_cls(f"{name} must be non-empty")
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise error_cls(f"{name} must contain finite positive values")
    return arr


def resolve_frequency_or_periods(
    *,
    frequency: Any | None,
    periods: Any | None,
    error_cls: type[Exception] = ValueError,
) -> tuple[np.ndarray, np.ndarray]:
    """Return validated ``(frequency, periods)`` from exactly one input.

    Exactly one of *frequency* or *periods* must be supplied; the other
    axis is derived as its reciprocal. This is the shared contract used
    by every airborne technology adapter that accepts either axis.

    Returns
    -------
    (ndarray, ndarray)
        ``(frequency, periods)``, each of shape ``(n,)``.

    Raises
    ------
    error_cls
        If both or neither of *frequency*/*periods* are supplied, or if
        the supplied axis fails :func:`normalize_frequency`.
    """
    if (frequency is None) == (periods is None):
        raise error_cls(
            "exactly one of frequency or periods must be supplied"
        )
    if frequency is not None:
        freq = normalize_frequency(frequency, error_cls=error_cls)
        return freq, 1.0 / freq
    period_arr = normalize_frequency(
        periods, name="periods", error_cls=error_cls
    )
    return 1.0 / period_arr, period_arr


def resolve_line_frequency_grid(
    frequency: Any,
    *,
    n_samples: int,
    n_frequency: int,
    error_cls: type[Exception] = ValueError,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Return ``(common_frequency, frequency_rows)`` for one flight line.

    Parameters
    ----------
    frequency : array-like
        Either one shared frequency vector of shape ``(n_frequency,)``
        or a per-sample grid of shape ``(n_samples, n_frequency)``.
    n_samples : int
        Number of navigation samples on the line.
    n_frequency : int
        Number of frequency samples expected from the response data.
    error_cls : type, default ValueError
        Exception type raised on shape violations.

    Returns
    -------
    (ndarray or None, ndarray or None)
        Exactly one of ``(common_frequency, None)`` or
        ``(None, frequency_rows)``.

    Raises
    ------
    error_cls
        If *frequency* is 1-D but its length does not match
        *n_frequency*, or if it is neither a valid ``(n_frequency,)``
        vector nor a ``(n_samples, n_frequency)`` matrix.

    Notes
    -----
    The shared 1-D case is fully validated here via
    :func:`normalize_frequency`. The per-sample 2-D case is only
    shape-checked: each row's own finiteness/positivity is validated
    lazily, once per attached sample, by the caller (typically via
    another :func:`normalize_frequency` call inside its per-sample
    loop). This means a row for a sample excluded by a line's
    ``record_mask`` is never required to be valid -- consistent with
    the rest of this module's missing-is-not-invalid stance, since a
    masked-out sample contributes no record for that row to belong to.
    """
    freq = np.asarray(frequency, dtype=float)
    if freq.ndim == 1:
        common = normalize_frequency(freq, error_cls=error_cls)
        if common.size != n_frequency:
            raise error_cls(
                "shared frequency length does not match response data: "
                f"{common.size} != {n_frequency}"
            )
        return common, None
    if freq.ndim == 2 and freq.shape == (n_samples, n_frequency):
        return None, freq
    raise error_cls(
        "frequency must have shape (nf,) or (n_samples, nf)"
    )


def normalize_positive_float(
    value: Any,
    *,
    name: str,
    error_cls: type[Exception] = ValueError,
) -> float:
    """Return *value* as a finite, strictly positive ``float``.

    Shared by every technology ``SystemSpec`` publishing a single
    positive descriptive rate/measurement (a sampling rate, an output
    rate, a coil angle, ...) as opposed to a range; see
    :func:`normalize_frequency_range` for the two-value case.
    """
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise error_cls(f"{name} must be positive")
    return number


def normalize_count_range(
    value: Any,
    *,
    name: str,
    error_cls: type[Exception] = ValueError,
) -> tuple[int, int]:
    """Return a validated ``(low, high)`` count with ``0 < low <= high``.

    Shared by every technology ``SystemSpec`` publishing a typical
    minimum/maximum item count (for example a typical processed
    frequency-window count). Unlike :func:`normalize_frequency_range`,
    equality (``low == high``) is valid here: a system with a fixed
    count still has "typical min == typical max".
    """
    try:
        low, high = (int(v) for v in value)
    except (TypeError, ValueError) as exc:
        raise error_cls(
            f"{name} must contain exactly two values"
        ) from exc
    if low <= 0 or high < low:
        raise error_cls(f"{name} must satisfy 0 < low <= high")
    return low, high


def normalize_frequency_range(
    value: Any,
    *,
    name: str,
    error_cls: type[Exception] = ValueError,
) -> tuple[float, float]:
    """Return a validated ``(low, high)`` band with ``0 < low < high``.

    Shared by every technology ``SystemSpec`` that publishes a nominal or
    practical frequency band as descriptive metadata.
    """
    try:
        low, high = (float(v) for v in value)
    except (TypeError, ValueError) as exc:
        raise error_cls(
            f"{name} must contain exactly two values"
        ) from exc
    if not np.isfinite(low) or not np.isfinite(high):
        raise error_cls(f"{name} must be finite")
    if low <= 0.0 or high <= low:
        raise error_cls(f"{name} must satisfy 0 < low < high")
    return low, high


def normalize_fixed_channels(
    value: Iterable[Any],
    *,
    expected: tuple[str, ...],
    name: str,
    error_cls: type[Exception] = ValueError,
) -> tuple[str, ...]:
    """Return *value* as a stripped tuple, requiring it to equal *expected*.

    Several technology contracts (for example MobileMT's ``Ex``/``Ey``
    admittance inputs) fix the channel layout by scientific definition
    rather than by user choice. Centralizing this "must equal" check
    keeps the behavior and message consistent across adapters.
    """
    channels = tuple(str(v).strip() for v in value)
    if channels != expected:
        raise error_cls(f"{name} must be {expected}, got {channels}")
    return channels


def normalize_estimate_array(
    value: Any,
    *,
    n_frequency: int,
    tail: tuple[int, int],
    name: str,
    error_cls: type[Exception] = ValueError,
) -> np.ndarray:
    """Return a validated ``(n_frequency, *tail)`` estimate array.

    A single ``tail``-shaped matrix is promoted to one frequency. This is
    the shared contract for ``VAR``/``INVSIGCOV``/``RESIDCOV`` payloads
    across technology adapters.
    """
    arr = np.asarray(value)
    if arr.ndim == 2 and arr.shape == tail:
        arr = arr[None, ...]
    if arr.ndim != 3 or arr.shape != (n_frequency, *tail):
        raise error_cls(
            f"{name} must have shape {(n_frequency, *tail)}, "
            f"got {arr.shape}"
        )
    if arr.dtype.kind not in "biufc":
        raise error_cls(f"{name} must be numeric")
    return arr


def normalize_sample_axis_array(
    value: Any,
    *,
    name: str,
    n_samples: int,
    expected: tuple[int, ...],
    error_cls: type[Exception] = ValueError,
) -> np.ndarray:
    """Return a validated array with a leading sample axis of *n_samples*.

    When ``n_samples == 1`` a caller may omit the leading axis; it is
    promoted automatically. This is the shared contract used to attach
    per-sample transfer-function arrays (admittance, tipper, tensor, ...)
    and per-sample statistical estimates to a flight line.
    """
    arr = np.asarray(value)
    if n_samples == 1 and arr.shape == expected[1:]:
        arr = arr[None, ...]
    if arr.shape != expected:
        raise error_cls(
            f"{name} must have shape {expected}, got {arr.shape}"
        )
    return arr


def normalize_record_mask(
    value: Any | None,
    *,
    n_samples: int,
    error_cls: type[Exception] = ValueError,
) -> np.ndarray:
    """Return a boolean record-coverage mask of length *n_samples*.

    ``None`` means every navigation sample has an attached EM record.
    """
    if value is None:
        return np.ones(n_samples, dtype=bool)
    mask = np.asarray(value, dtype=bool)
    if mask.ndim != 1 or mask.size != n_samples:
        raise error_cls(
            "record_mask must have one boolean per navigation sample"
        )
    return mask


def reference_station_mapping(
    reference_station: Any | None,
    *,
    channel_fields: tuple[str, ...],
) -> dict[str, Any] | None:
    """Return a duck-typed reference station as a plain ``EMTF.attrs`` map.

    Every technology's ``*ReferenceStation`` metadata class (MobileMT's,
    ZTEM's, AFMAG's) exposes ``preferred_id``, one or more channel-name
    tuples, an optional ``site``, and an optional ``attrs`` dict; this
    is the shared shape those adapters build into
    ``EMTF.attrs["<technology>"]["reference_station"]``.

    Parameters
    ----------
    reference_station : Any or None
        Duck-typed reference-station object exposing ``preferred_id``,
        ``site``, ``attrs``, and one attribute per name in
        *channel_fields*. ``None`` returns ``None``.
    channel_fields : tuple of str
        Attribute names to read off *reference_station* and include
        verbatim under the same key, for example ``("electric_channels",)``
        for MobileMT or ``("measured_channels", "transfer_input_channels")``
        for AirMt.

    Returns
    -------
    dict or None
        ``None`` when *reference_station* is ``None``; otherwise a
        mapping with ``"station_id"``, one entry per *channel_fields*
        name, and ``"site"``/``"attrs"`` when those are non-empty.
    """
    if reference_station is None:
        return None
    out: dict[str, Any] = {"station_id": reference_station.preferred_id}
    for field_name in channel_fields:
        out[field_name] = list(getattr(reference_station, field_name))
    if reference_station.site is not None:
        out["site"] = reference_station.site.to_dict(max_depth=2)
    if reference_station.attrs:
        out["attrs"] = dict(reference_station.attrs)
    return out


def merge_remote_reference_processing(
    reference_station: Any | None,
    processing: Any | None,
    *,
    reference_type: str,
    technology: str,
    extra: Mapping[str, Any] = MappingProxyType({}),
    error_cls: type[Exception] = ValueError,
):
    """Merge a duck-typed reference station into remote-reference metadata.

    ``reference_station`` is each technology's scientifically typed way
    to supply its fixed ground reference; a caller-supplied
    ``processing`` is the general
    :class:`~pycsamt.metadata.ProcessingMeta` way. When both are given,
    they must describe the same reference site -- this raises rather
    than silently letting one win. When only *reference_station* is
    given, a :class:`~pycsamt.metadata.ProcessingMeta` is synthesized
    around it so :attr:`~pycsamt.emtf.EMTF.processing` is always the
    one place downstream code (for example
    :func:`~pycsamt.airborne.qc.assess_airborne_qc`) looks for
    reference metadata, regardless of which technology built the EMTF.

    Parameters
    ----------
    reference_station : Any or None
        Duck-typed reference-station object exposing ``preferred_id``.
        ``None`` returns *processing* unchanged.
    processing : ProcessingMeta or None
        Caller-supplied processing metadata to merge into.
    reference_type : str
        ``RemoteReferenceMeta.reference_type`` to record, for example
        ``"fixed_ground_magnetic"``.
    technology : str
        Recorded as ``extra["technology"]`` on the synthesized
        :class:`~pycsamt.metadata.RemoteReferenceMeta`.
    extra : Mapping, optional
        Additional technology-specific ``RemoteReferenceMeta.extra``
        entries (for example measured/transfer channel names).
    error_cls : type, default ValueError
        Exception type raised on a genuine site conflict. The type
        check on *processing* itself always raises ``TypeError``,
        regardless of *error_cls*, since that failure is never a
        reference-station/processing conflict.

    Returns
    -------
    ProcessingMeta or None
        *processing*, either unchanged, merged with the synthesized
        remote reference, or newly created around it.

    Raises
    ------
    TypeError
        If *processing* is supplied and is not a
        :class:`~pycsamt.metadata.ProcessingMeta`.
    error_cls
        If *processing* already has a remote reference whose ``site``
        conflicts with *reference_station*'s.
    """
    from ..metadata import ProcessingMeta, RemoteReferenceMeta

    if processing is not None and not isinstance(processing, ProcessingMeta):
        raise TypeError("processing must be a ProcessingMeta or None")
    if reference_station is None:
        return processing

    remote = RemoteReferenceMeta(
        reference_type=reference_type,
        site=reference_station.preferred_id,
        extra={**dict(extra), "technology": technology},
    )
    if processing is None:
        return ProcessingMeta(remote_reference=remote)
    if processing.remote_reference is None:
        return ProcessingMeta(
            sign_convention=processing.sign_convention,
            processed_by=processing.processed_by,
            software=processing.software,
            remote_reference=remote,
            processing_tag=processing.processing_tag,
            run_list=(
                None
                if processing.run_list is None
                else list(processing.run_list)
            ),
            extra=dict(processing.extra),
        )

    existing = processing.remote_reference
    reference_id = reference_station.preferred_id
    if existing.site is not None and reference_id is not None:
        if str(existing.site) != str(reference_id):
            raise error_cls(
                "processing remote-reference site conflicts with the "
                "supplied reference station metadata"
            )
    return processing
