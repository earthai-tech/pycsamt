"""CSEM controlled-source electromagnetic edge quality control.

CSEM shares the controlled-source primitives in
:mod:`~pycsamt.iot.edge_csamt` (skin-depth field zones, transmitter
frequency comb, source stability), but its defining measurement is
different: a dipole source is recorded by a *receiver array* and the
signature data product is the response as a function of source-receiver
**offset** at each frequency -- magnitude-versus-offset (MVO) and
phase-versus-offset (PVO).

This module provides the offset-domain QC that is specific to CSEM:

* :func:`field_vs_offset` builds an MVO/PVO curve, flags where the signal
  drops below the noise floor (the detectability limit), and checks that
  the amplitude decays monotonically with offset -- a bump usually means a
  bad receiver, a geometry error, or genuine 3-D structure worth a second
  look;
* :func:`csem_edge_report` / :func:`csem_edge_table` collate one or more
  frequency responses.

Everything is numpy-only, consistent with the rest of the edge modules.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from ..api.property import PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c

__all__ = [
    "OffsetResponse",
    "field_vs_offset",
    "csem_edge_report",
    "csem_edge_table",
]


@dataclass
class OffsetResponse(PyCSAMTObject):
    """Result of :func:`field_vs_offset` -- one MVO/PVO curve."""

    frequency_hz: float | None
    offsets_m: list[float] = field(default_factory=list)
    amplitudes: list[float] = field(default_factory=list)
    phases_deg: list[float] | None = None
    above_noise: list[bool] = field(default_factory=list)
    max_detectable_offset_m: float | None = None
    monotonic_decay: bool = False
    dynamic_range_db: float = float("nan")

    @property
    def n_offsets(self) -> int:
        return len(self.offsets_m)

    @property
    def n_detectable(self) -> int:
        return int(sum(1 for ok in self.above_noise if ok))

    def as_dict(self) -> dict[str, Any]:
        return dict(
            frequency_hz=self.frequency_hz,
            n_offsets=self.n_offsets,
            n_detectable=self.n_detectable,
            max_detectable_offset_m=self.max_detectable_offset_m,
            monotonic_decay=self.monotonic_decay,
            dynamic_range_db=self.dynamic_range_db,
            offsets_m=list(self.offsets_m),
            amplitudes=list(self.amplitudes),
            phases_deg=(
                list(self.phases_deg) if self.phases_deg is not None else None
            ),
            above_noise=list(self.above_noise),
        )


def field_vs_offset(
    offsets_m: Any,
    amplitudes: Any,
    *,
    phases_deg: Any = None,
    noise_floor: float | None = None,
    frequency_hz: float | None = None,
    monotonic_tol: float = 0.05,
) -> OffsetResponse:
    """Build a magnitude/phase-versus-offset curve and QC it.

    Parameters
    ----------
    offsets_m : array-like of float
        Source-receiver offsets in metres (need not be sorted).
    amplitudes : array-like of float
        Field amplitude at each offset (e.g. normalised E-field). Must be
        non-negative and the same length as *offsets_m*.
    phases_deg : array-like of float, optional
        Phase at each offset in degrees, carried through for PVO.
    noise_floor : float, optional
        Amplitude below which a reading is not detectable. When given, the
        largest offset above the floor is reported as the detectability
        limit and only detectable points drive the decay/dynamic-range QC.
    frequency_hz : float, optional
        Frequency this curve was measured at, recorded for reference.
    monotonic_tol : float
        Fractional tolerance when checking that amplitude never increases
        with offset (``a[i+1] <= a[i] * (1 + tol)``).

    Returns
    -------
    OffsetResponse
    """
    off = np.asarray(offsets_m, dtype=float).ravel()
    amp = np.asarray(amplitudes, dtype=float).ravel()
    if off.shape != amp.shape:
        raise ValueError("offsets_m and amplitudes must have equal length.")
    ph = None
    if phases_deg is not None:
        ph = np.asarray(phases_deg, dtype=float).ravel()
        if ph.shape != off.shape:
            raise ValueError("phases_deg must match the length of offsets_m.")
    monotonic_tol = _c.as_nonnegative(monotonic_tol, "monotonic_tol")
    if frequency_hz is not None:
        frequency_hz = _c.as_positive(frequency_hz, "frequency_hz")

    result = OffsetResponse(frequency_hz=frequency_hz)
    if off.size == 0:
        return result

    finite = np.isfinite(off) & np.isfinite(amp)
    if ph is not None:
        finite &= np.isfinite(ph)
    off, amp = off[finite], amp[finite]
    if ph is not None:
        ph = ph[finite]
    if off.size == 0:
        return result
    if np.any(amp < 0):
        raise ValueError("amplitudes must be non-negative.")

    order = np.argsort(off)
    off, amp = off[order], amp[order]
    if ph is not None:
        ph = ph[order]

    if noise_floor is not None:
        floor = _c.as_positive(noise_floor, "noise_floor")
        above = amp > floor
    else:
        above = np.ones(off.shape, dtype=bool)

    result.offsets_m = [float(v) for v in off]
    result.amplitudes = [float(v) for v in amp]
    result.phases_deg = None if ph is None else [float(v) for v in ph]
    result.above_noise = [bool(v) for v in above]

    detectable_off = off[above]
    result.max_detectable_offset_m = (
        float(detectable_off.max()) if detectable_off.size else None
    )

    det_amp = amp[above]
    if det_amp.size >= 2:
        # Monotonic decay: amplitude should not grow as offset increases.
        increases = det_amp[1:] > det_amp[:-1] * (1.0 + monotonic_tol)
        result.monotonic_decay = not bool(np.any(increases))
        lo = float(np.min(det_amp))
        hi = float(np.max(det_amp))
        result.dynamic_range_db = (
            float(20.0 * np.log10(hi / lo)) if lo > 0 else float("inf")
        )
    elif det_amp.size == 1:
        result.monotonic_decay = True
        result.dynamic_range_db = 0.0
    return result


def csem_edge_report(
    offsets_m: Any,
    amplitudes: Any,
    *,
    phases_deg: Any = None,
    noise_floor: float | None = None,
    frequency_hz: float | None = None,
) -> dict[str, Any]:
    """Run the CSEM offset-domain diagnostics for one frequency."""
    response = field_vs_offset(
        offsets_m,
        amplitudes,
        phases_deg=phases_deg,
        noise_floor=noise_floor,
        frequency_hz=frequency_hz,
    )
    return dict(offset_response=response.as_dict())


def csem_edge_table(
    reports: dict[str, dict[str, Any]] | Iterable[tuple[str, dict[str, Any]]],
    *,
    api: bool | None = None,
) -> Any:
    """Flatten one or more :func:`csem_edge_report` results into a table.

    Accepts a ``{label: report}`` mapping (or ``(label, report)`` pairs);
    the label is typically a frequency or receiver-line identifier.
    """
    items = (
        list(reports.items()) if isinstance(reports, dict) else list(reports)
    )
    rows: list[dict[str, Any]] = []
    for label, report in items:
        resp = report.get("offset_response", {})
        rows.append(
            dict(
                label=str(label),
                frequency_hz=resp.get("frequency_hz"),
                n_offsets=resp.get("n_offsets"),
                n_detectable=resp.get("n_detectable"),
                max_detectable_offset_m=resp.get("max_detectable_offset_m"),
                monotonic_decay=resp.get("monotonic_decay"),
                dynamic_range_db=resp.get("dynamic_range_db"),
            )
        )
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_csem_edge_table",
        kind="iot.edge.csem",
        source=items,
        description="CSEM magnitude/phase-versus-offset edge QC by frequency.",
    )
