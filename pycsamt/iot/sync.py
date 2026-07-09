"""Clock synchronisation audit for IoT field acquisition.

Offset alone is not enough to certify field-grade timing. This module
adds clock *drift* (ppm), *jitter* (ms), GPS-lock accounting, and an
overall :class:`SyncQuality` grade, plus batch assessment and GPS-dropout
detection across a deployment.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..api.property import PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c

__all__ = [
    "SyncConfig",
    "SyncQuality",
    "SyncStatus",
    "ClockSynchronizer",
    "estimate_clock_offset_ms",
    "estimate_clock_drift_ppm",
    "estimate_clock_jitter_ms",
    "assess_sync_quality",
    "batch_assess_sync",
    "detect_gps_dropout",
    "sync_status_table",
]


class SyncQuality(str, Enum):
    """Overall synchronisation grade for a device."""

    EXCELLENT = "excellent"
    GOOD = "good"
    FAIR = "fair"
    POOR = "poor"
    UNKNOWN = "unknown"


@dataclass
class SyncConfig(PyCSAMTObject):
    """Clock-synchronisation tolerances."""

    tolerance_ms: float = 1.0
    reference: str = "gps"
    max_drift_ppm: Optional[float] = None
    max_jitter_ms: Optional[float] = None
    min_reference_points: int = 2

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise sync tolerances."""
        self.tolerance_ms = _c.as_positive(self.tolerance_ms, "tolerance_ms")
        self.reference = _c.as_nonempty_str(self.reference, "reference").lower()
        self.max_drift_ppm = _c.as_optional_positive(
            self.max_drift_ppm, "max_drift_ppm"
        )
        self.max_jitter_ms = _c.as_optional_positive(
            self.max_jitter_ms, "max_jitter_ms"
        )
        self.min_reference_points = int(self.min_reference_points)
        if self.min_reference_points < 1:
            raise ValueError("min_reference_points must be >= 1.")


@dataclass
class SyncStatus(PyCSAMTObject):
    """Clock-synchronisation status for one device."""

    device_id: str
    offset_ms: float
    within_tolerance: bool
    reference: str = "gps"
    drift_ppm: float = float("nan")
    jitter_ms: float = float("nan")
    gps_lock: Optional[bool] = None
    n_reference_points: int = 0
    quality: SyncQuality | str = SyncQuality.UNKNOWN

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise sync-status fields."""
        self.device_id = _c.as_nonempty_str(self.device_id, "device_id")
        self.offset_ms = float(self.offset_ms)
        self.within_tolerance = bool(self.within_tolerance)
        self.reference = str(self.reference or "gps").lower()
        self.drift_ppm = float(self.drift_ppm)
        self.jitter_ms = float(self.jitter_ms)
        if self.gps_lock is not None:
            self.gps_lock = _c.as_bool(self.gps_lock)
        self.n_reference_points = int(self.n_reference_points)
        if self.n_reference_points < 0:
            raise ValueError("n_reference_points must be >= 0.")
        self.quality = _c.normalise_enum(self.quality, SyncQuality, "quality")

    def as_dict(self) -> Dict[str, Any]:
        """Return a flat, serialisable status dictionary."""
        return dict(
            device_id=self.device_id,
            offset_ms=self.offset_ms,
            within_tolerance=self.within_tolerance,
            reference=self.reference,
            drift_ppm=self.drift_ppm,
            jitter_ms=self.jitter_ms,
            gps_lock=self.gps_lock,
            n_reference_points=self.n_reference_points,
            quality=self.quality.value if isinstance(self.quality, SyncQuality)
            else str(self.quality),
        )


def _paired(local: Any, reference: Any) -> Tuple[np.ndarray, np.ndarray]:
    loc = np.asarray(local, dtype=float).ravel()
    ref = np.asarray(reference, dtype=float).ravel()
    n = min(loc.size, ref.size)
    if n == 0:
        return np.empty(0), np.empty(0)
    loc, ref = loc[:n], ref[:n]
    finite = np.isfinite(loc) & np.isfinite(ref)
    return loc[finite], ref[finite]


def estimate_clock_offset_ms(
    local_timestamps: Any,
    reference_timestamps: Any,
) -> float:
    """Estimate median local-reference clock offset in milliseconds."""
    loc, ref = _paired(local_timestamps, reference_timestamps)
    if loc.size == 0:
        return float("nan")
    return float(np.median((loc - ref) * 1000.0))


def estimate_clock_drift_ppm(
    local_timestamps: Any,
    reference_timestamps: Any,
) -> float:
    """Estimate clock drift in parts-per-million.

    Fits the local-reference offset (in seconds) against reference time
    and returns the slope scaled to ppm. Requires at least two distinct
    reference points; otherwise returns ``nan``.
    """
    loc, ref = _paired(local_timestamps, reference_timestamps)
    if loc.size < 2 or np.ptp(ref) <= 0:
        return float("nan")
    offset_s = loc - ref
    slope = float(np.polyfit(ref, offset_s, 1)[0])  # seconds per second
    return slope * 1.0e6


def estimate_clock_jitter_ms(
    local_timestamps: Any,
    reference_timestamps: Any,
) -> float:
    """Estimate timing jitter as the std of drift-corrected offsets (ms)."""
    loc, ref = _paired(local_timestamps, reference_timestamps)
    if loc.size < 2:
        return float("nan")
    offset_ms = (loc - ref) * 1000.0
    if loc.size >= 3 and np.ptp(ref) > 0:
        slope, intercept = np.polyfit(ref, offset_ms, 1)
        residual = offset_ms - (slope * ref + intercept)
        return float(np.std(residual))
    return float(np.std(offset_ms))


def assess_sync_quality(
    *,
    offset_ms: float,
    drift_ppm: float = float("nan"),
    jitter_ms: float = float("nan"),
    gps_lock: Optional[bool] = None,
    tolerance_ms: float = 1.0,
) -> SyncQuality:
    """Grade synchronisation from offset, drift, jitter, and GPS lock."""
    if not np.isfinite(offset_ms):
        return SyncQuality.UNKNOWN
    tol = float(tolerance_ms)
    abs_offset = abs(offset_ms)
    abs_drift = abs(drift_ppm) if np.isfinite(drift_ppm) else 0.0
    jit = jitter_ms if np.isfinite(jitter_ms) else 0.0

    if gps_lock is False:
        # Free-running clock caps achievable quality regardless of offset.
        if abs_offset <= 5.0 * tol:
            return SyncQuality.FAIR
        return SyncQuality.POOR
    if abs_offset <= tol and abs_drift <= 2.0 and jit <= 0.5 * tol:
        return SyncQuality.EXCELLENT
    if abs_offset <= tol and abs_drift <= 10.0 and jit <= tol:
        return SyncQuality.GOOD
    if abs_offset <= 5.0 * tol and abs_drift <= 50.0:
        return SyncQuality.FAIR
    return SyncQuality.POOR


class ClockSynchronizer:
    """Evaluate device clock status against a reference."""

    def __init__(self, config: Optional[SyncConfig] = None) -> None:
        self.config = config or SyncConfig()
        self.config.validate()

    def assess(
        self,
        device_id: str,
        local_timestamps: Any,
        reference_timestamps: Any,
        *,
        gps_lock: Optional[bool] = None,
    ) -> SyncStatus:
        """Return offset, drift, jitter, and an overall quality grade."""
        loc, ref = _paired(local_timestamps, reference_timestamps)
        offset = estimate_clock_offset_ms(loc, ref)
        drift = estimate_clock_drift_ppm(loc, ref)
        jitter = estimate_clock_jitter_ms(loc, ref)
        within = bool(
            np.isfinite(offset)
            and abs(offset) <= self.config.tolerance_ms
            and (self.config.max_drift_ppm is None
                 or not np.isfinite(drift)
                 or abs(drift) <= self.config.max_drift_ppm)
            and (self.config.max_jitter_ms is None
                 or not np.isfinite(jitter)
                 or jitter <= self.config.max_jitter_ms)
        )
        quality = assess_sync_quality(
            offset_ms=offset,
            drift_ppm=drift,
            jitter_ms=jitter,
            gps_lock=gps_lock,
            tolerance_ms=self.config.tolerance_ms,
        )
        return SyncStatus(
            device_id=str(device_id),
            offset_ms=offset,
            within_tolerance=within,
            reference=self.config.reference,
            drift_ppm=drift,
            jitter_ms=jitter,
            gps_lock=gps_lock,
            n_reference_points=int(loc.size),
            quality=quality,
        )


def batch_assess_sync(
    references: Mapping[str, Any],
    *,
    config: Optional[SyncConfig] = None,
    api: bool | None = None,
) -> Any:
    """Assess many devices at once and return a status table.

    Parameters
    ----------
    references : mapping
        ``{device_id: spec}`` where *spec* is either a ``(local,
        reference)`` pair or a mapping with ``local`` / ``reference``
        (and optional ``gps_lock``) keys.
    config : SyncConfig, optional
        Shared tolerances.
    api : bool, optional
        Frame-wrapping override passed to :func:`maybe_wrap_frame`.
    """
    synchronizer = ClockSynchronizer(config)
    statuses: List[SyncStatus] = []
    for device_id, spec in dict(references).items():
        local, reference, gps_lock = _unpack_reference_spec(spec)
        statuses.append(
            synchronizer.assess(
                device_id, local, reference, gps_lock=gps_lock
            )
        )
    return sync_status_table(statuses, api=api)


def _unpack_reference_spec(
    spec: Any,
) -> Tuple[Any, Any, Optional[bool]]:
    if isinstance(spec, Mapping):
        local = spec.get("local", spec.get("local_timestamps"))
        reference = spec.get("reference", spec.get("reference_timestamps"))
        gps = spec.get("gps_lock")
        return local, reference, (None if gps is None else _c.as_bool(gps))
    if isinstance(spec, Sequence) and not isinstance(spec, (str, bytes)):
        items = list(spec)
        if len(items) < 2:
            raise ValueError(
                "reference spec must provide (local, reference) sequences."
            )
        gps = items[2] if len(items) > 2 else None
        return items[0], items[1], (None if gps is None else _c.as_bool(gps))
    raise TypeError(
        "reference spec must be a (local, reference) pair or a mapping."
    )


def detect_gps_dropout(
    gps_lock: Iterable[Any],
    timestamps: Optional[Iterable[Any]] = None,
    *,
    min_lock_fraction: float = 0.9,
) -> Dict[str, Any]:
    """Summarise GPS-lock dropouts across a sample sequence.

    Parameters
    ----------
    gps_lock : iterable
        Per-sample lock flags (bool/0-1). Non-finite entries are treated
        as unlocked.
    timestamps : iterable, optional
        Matching sample timestamps (seconds), used to report the longest
        dropout duration.
    min_lock_fraction : float
        Threshold below which ``ok`` is ``False``.

    Returns
    -------
    dict
        Keys: ``n_samples``, ``n_locked``, ``lock_fraction``,
        ``n_dropout_events``, ``longest_dropout_samples``,
        ``longest_dropout_s``, and ``ok``.
    """
    flags = []
    for value in gps_lock:
        try:
            flags.append(bool(_c.as_bool(value)))
        except ValueError:
            flags.append(False)
    n = len(flags)
    min_lock_fraction = _c.as_probability(min_lock_fraction, "min_lock_fraction")
    if n == 0:
        return dict(
            n_samples=0, n_locked=0, lock_fraction=float("nan"),
            n_dropout_events=0, longest_dropout_samples=0,
            longest_dropout_s=float("nan"), ok=False,
        )
    times = None
    if timestamps is not None:
        times = np.asarray(list(timestamps), dtype=float).ravel()
        if times.size != n:
            times = None

    n_locked = int(sum(flags))
    events = 0
    longest = 0
    longest_s = 0.0
    run_start: Optional[int] = None
    for i, locked in enumerate(flags):
        if not locked and run_start is None:
            run_start = i
            events += 1
        if (locked or i == n - 1) and run_start is not None:
            end = i if locked else i + 1
            run_len = end - run_start
            longest = max(longest, run_len)
            if times is not None and end - 1 < n:
                span = float(times[min(end, n - 1)] - times[run_start])
                longest_s = max(longest_s, abs(span))
            run_start = None
    lock_fraction = n_locked / n
    return dict(
        n_samples=n,
        n_locked=n_locked,
        lock_fraction=lock_fraction,
        n_dropout_events=events,
        longest_dropout_samples=longest,
        longest_dropout_s=(longest_s if times is not None else float("nan")),
        ok=bool(lock_fraction >= min_lock_fraction),
    )


def sync_status_table(
    statuses: SyncStatus | Iterable[SyncStatus],
    *,
    api: bool | None = None,
) -> Any:
    """Return one row per device from sync-status objects."""
    items = [statuses] if isinstance(statuses, SyncStatus) else list(statuses)
    df = pd.DataFrame.from_records([status.as_dict() for status in items])
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_sync_status",
        kind="iot.sync.status",
        source=items,
        description="Clock-synchronisation status by device.",
    )
