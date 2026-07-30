"""CSAMT / controlled-source specific edge quality-control metrics.

:mod:`~pycsamt.iot.edge_amt` covers natural-source AMT/MT concerns
(powerline harmonics, passive electrode contact, continuous spectra).
Controlled-source audio-magnetotellurics adds a transmitter and therefore
a different set of field diagnostics:

* **field-zone classification** — the single most important CSAMT QC:
  whether a frequency is recorded in the far field (where the plane-wave
  MT assumption holds) or has slipped into the transition/near field,
  where apparent resistivities roll off and a near-field correction is
  required. Driven by the skin depth versus the transmitter-receiver
  offset;
* **transmitter-frequency comb detection** — CSAMT transmits a discrete
  set of frequencies rather than a natural continuum, so QC checks that
  energy is actually present, and resolvable, at each expected line;
* **source-signal stability** — transmitter current/voltage steadiness
  and on/off keying, which bound the quality of every sounding.

Like :mod:`~pycsamt.iot.edge_amt`, everything here is numpy-only (the
spectral helpers are shared with that module) so the subpackage imports
without SciPy.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np
import pandas as pd

from ..api.property import PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c
from .edge_amt import _positive_sample_rate, _prep_signal, _welch_psd

__all__ = [
    "FieldZone",
    "FieldZoneCoverage",
    "SourceLine",
    "TransmitterComb",
    "SourceStability",
    "skin_depth_m",
    "classify_field_zones",
    "detect_transmitter_frequencies",
    "assess_source_stability",
    "csamt_edge_report",
    "csamt_edge_table",
]

# 503.29 = sqrt(2 / (mu0 * (2*pi))) * ... the standard EM skin-depth constant
# so that delta[m] = 503.29 * sqrt(rho[ohm.m] / f[Hz]).
_SKIN_DEPTH_CONST = 503.29


class FieldZone(str, Enum):
    """Wave-propagation regime of a CSAMT measurement at one frequency."""

    NEAR = "near"
    TRANSITION = "transition"
    FAR = "far"


# ---------------------------------------------------------------------------
# skin depth and field zones
# ---------------------------------------------------------------------------
def skin_depth_m(resistivity: Any, freq: Any) -> np.ndarray:
    r"""Electromagnetic skin depth in metres.

    :math:`\delta = 503.29\,\sqrt{\rho / f}` with :math:`\rho` in
    :math:`\Omega\cdot m` and :math:`f` in Hz. Scalars and arrays
    broadcast together; non-positive or non-finite inputs yield ``nan``.
    """
    rho = np.asarray(resistivity, dtype=float)
    f = np.asarray(freq, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = _SKIN_DEPTH_CONST * np.sqrt(rho / f)
    delta = np.where(np.isfinite(delta) & (rho > 0) & (f > 0), delta, np.nan)
    return delta


@dataclass
class FieldZoneCoverage(PyCSAMTObject):
    """Result of :func:`classify_field_zones`."""

    offset_m: float
    freq_hz: list[float] = field(default_factory=list)
    skin_depth_m: list[float] = field(default_factory=list)
    offset_ratio: list[float] = field(default_factory=list)
    zones: list[str] = field(default_factory=list)
    n_near: int = 0
    n_transition: int = 0
    n_far: int = 0

    @property
    def far_fraction(self) -> float:
        """Fraction of frequencies recorded in the far field."""
        total = len(self.zones)
        return (self.n_far / total) if total else float("nan")

    @property
    def all_far_field(self) -> bool:
        """True when every frequency satisfies the far-field assumption."""
        return bool(self.zones) and self.n_far == len(self.zones)

    @property
    def correction_recommended(self) -> bool:
        """True when any frequency needs a near-field correction."""
        return (self.n_near + self.n_transition) > 0

    def first_far_field_hz(self) -> float | None:
        """Lowest frequency that reaches the far field, if any."""
        far = [
            f
            for f, z in zip(self.freq_hz, self.zones)
            if z == FieldZone.FAR.value
        ]
        return min(far) if far else None

    def as_dict(self) -> dict[str, Any]:
        return dict(
            offset_m=self.offset_m,
            n_near=self.n_near,
            n_transition=self.n_transition,
            n_far=self.n_far,
            far_fraction=self.far_fraction,
            all_far_field=self.all_far_field,
            correction_recommended=self.correction_recommended,
            first_far_field_hz=self.first_far_field_hz(),
            freq_hz=list(self.freq_hz),
            skin_depth_m=list(self.skin_depth_m),
            offset_ratio=list(self.offset_ratio),
            zones=list(self.zones),
        )


def classify_field_zones(
    freq: Any,
    resistivity: Any,
    offset_m: float,
    *,
    near_ratio: float = 1.0,
    far_ratio: float = 3.0,
) -> FieldZoneCoverage:
    r"""Classify each frequency as near, transition, or far field.

    The controlling quantity is the transmitter-receiver offset expressed
    in skin depths, :math:`r / \delta`. A larger ratio means the receiver
    is electrically farther from the source and the plane-wave assumption
    is safer.

    Parameters
    ----------
    freq : array-like of float
        Frequencies in Hz.
    resistivity : float or array-like
        (Apparent) resistivity in :math:`\Omega\cdot m`. A scalar is used
        as a half-space value for every frequency; an array must match
        *freq*.
    offset_m : float
        Transmitter-receiver separation :math:`r` in metres.
    near_ratio : float
        ``r / delta`` at or below which a frequency is *near field*.
    far_ratio : float
        ``r / delta`` at or above which a frequency is *far field*.
        Values in between are the *transition zone*.

    Returns
    -------
    FieldZoneCoverage
    """
    f = np.asarray(freq, dtype=float).ravel()
    offset_m = _c.as_positive(offset_m, "offset_m")
    near_ratio = _c.as_positive(near_ratio, "near_ratio")
    far_ratio = _c.as_positive(far_ratio, "far_ratio")
    if far_ratio < near_ratio:
        raise ValueError("far_ratio must be >= near_ratio.")

    rho = np.asarray(resistivity, dtype=float).ravel()
    if rho.size == 1:
        rho = np.full(f.shape, float(rho.item()))
    elif rho.shape != f.shape:
        raise ValueError(
            "resistivity must be a scalar or match the length of freq."
        )

    delta = skin_depth_m(rho, f)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(delta > 0, offset_m / delta, np.nan)

    coverage = FieldZoneCoverage(offset_m=offset_m)
    for fi, di, ri in zip(f, delta, ratio):
        if not np.isfinite(ri):
            zone = FieldZone.NEAR  # undefined skin depth -> treat cautiously
        elif ri >= far_ratio:
            zone = FieldZone.FAR
        elif ri <= near_ratio:
            zone = FieldZone.NEAR
        else:
            zone = FieldZone.TRANSITION
        coverage.freq_hz.append(float(fi))
        coverage.skin_depth_m.append(
            float(di) if np.isfinite(di) else float("nan")
        )
        coverage.offset_ratio.append(
            float(ri) if np.isfinite(ri) else float("nan")
        )
        coverage.zones.append(zone.value)
        if zone is FieldZone.FAR:
            coverage.n_far += 1
        elif zone is FieldZone.TRANSITION:
            coverage.n_transition += 1
        else:
            coverage.n_near += 1
    return coverage


# ---------------------------------------------------------------------------
# transmitter frequency comb
# ---------------------------------------------------------------------------
@dataclass
class SourceLine(PyCSAMTObject):
    """Detection result for one expected transmitter frequency."""

    frequency_hz: float
    snr_db: float
    power_ratio: float
    detected: bool

    def as_dict(self) -> dict[str, Any]:
        return dict(
            frequency_hz=self.frequency_hz,
            snr_db=self.snr_db,
            power_ratio=self.power_ratio,
            detected=self.detected,
        )


@dataclass
class TransmitterComb(PyCSAMTObject):
    """Result of :func:`detect_transmitter_frequencies`."""

    lines: list[SourceLine] = field(default_factory=list)

    @property
    def n_expected(self) -> int:
        return len(self.lines)

    @property
    def n_detected(self) -> int:
        return sum(1 for line in self.lines if line.detected)

    @property
    def detection_fraction(self) -> float:
        return (
            (self.n_detected / self.n_expected) if self.lines else float("nan")
        )

    def missing(self) -> list[float]:
        """Expected frequencies whose energy was not resolvable."""
        return [ln.frequency_hz for ln in self.lines if not ln.detected]

    def as_dict(self) -> dict[str, Any]:
        return dict(
            n_expected=self.n_expected,
            n_detected=self.n_detected,
            detection_fraction=self.detection_fraction,
            missing_hz=self.missing(),
            lines=[ln.as_dict() for ln in self.lines],
        )


def detect_transmitter_frequencies(
    data: Any,
    sample_rate: float,
    tx_frequencies: Iterable[float],
    *,
    bandwidth_hz: float = 1.0,
    snr_threshold_db: float = 6.0,
) -> TransmitterComb:
    """Check that recorded energy is present at each transmitted frequency.

    For every expected transmitter frequency the in-band power is compared
    with a local out-of-band noise estimate to form a per-line SNR. A line
    counts as detected when its SNR meets ``snr_threshold_db``; lines above
    Nyquist are reported as undetected.

    Parameters
    ----------
    data : array-like
        Single-channel time series.
    sample_rate : float
        Sampling frequency in Hz.
    tx_frequencies : iterable of float
        Expected transmitter frequencies in Hz.
    bandwidth_hz : float
        Half-width of the signal band around each line.
    snr_threshold_db : float
        Minimum in-band/out-of-band SNR for a line to count as detected.

    Returns
    -------
    TransmitterComb
    """
    fs = _positive_sample_rate(sample_rate)
    bandwidth_hz = _c.as_positive(bandwidth_hz, "bandwidth_hz")
    lines_hz = [_c.as_positive(v, "tx_frequency") for v in tx_frequencies]
    comb = TransmitterComb()
    x = _prep_signal(data)
    freqs, psd = _welch_psd(x, fs)
    nyquist = fs / 2.0
    if freqs.size == 0:
        for fc in lines_hz:
            comb.lines.append(
                SourceLine(fc, float("nan"), float("nan"), False)
            )
        return comb

    df = float(freqs[1] - freqs[0]) if freqs.size >= 2 else float("nan")
    # A robust background floor: the median PSD. The comb lines are a small
    # minority of bins, so they do not bias the median, and neighbouring
    # transmitter lines never contaminate each other's noise estimate.
    positive = psd[psd > 0]
    floor = float(np.median(positive)) if positive.size else 0.0
    for fc in lines_hz:
        if fc >= nyquist or floor <= 0:
            comb.lines.append(SourceLine(fc, float("nan"), 0.0, False))
            continue
        bw = max(bandwidth_hz, 1.5 * df) if np.isfinite(df) else bandwidth_hz
        mask = (freqs >= fc - bw) & (freqs <= fc + bw)
        peak = float(np.max(psd[mask])) if mask.any() else 0.0
        ratio = float(peak / floor) if floor > 0 else 0.0
        snr_db = 10.0 * np.log10(ratio) if ratio > 0 else float("-inf")
        detected = np.isfinite(snr_db) and snr_db >= snr_threshold_db
        comb.lines.append(
            SourceLine(
                frequency_hz=float(fc),
                snr_db=float(snr_db),
                power_ratio=ratio,
                detected=bool(detected),
            )
        )
    return comb


# ---------------------------------------------------------------------------
# source signal stability
# ---------------------------------------------------------------------------
@dataclass
class SourceStability(PyCSAMTObject):
    """Result of :func:`assess_source_stability`."""

    current_mean_a: float
    current_cv: float
    current_drift_a: float
    on_fraction: float
    voltage_mean_v: float | None = None
    stable: bool = False
    flags: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return dict(
            current_mean_a=self.current_mean_a,
            current_cv=self.current_cv,
            current_drift_a=self.current_drift_a,
            on_fraction=self.on_fraction,
            voltage_mean_v=self.voltage_mean_v,
            stable=self.stable,
            flags=list(self.flags),
        )


def assess_source_stability(
    tx_current: Any,
    tx_voltage: Any = None,
    *,
    max_cv: float = 0.05,
    on_threshold_frac: float = 0.5,
) -> SourceStability:
    """Assess transmitter current (and optional voltage) steadiness.

    The transmitter current sets the signal level of every sounding, so
    its steadiness bounds data quality. This reports the coefficient of
    variation and linear drift of the *on-state* current, together with
    the on-fraction of an on/off-keyed source.

    Parameters
    ----------
    tx_current : array-like
        Transmitter current samples (A).
    tx_voltage : array-like, optional
        Transmitter voltage samples (V).
    max_cv : float
        Maximum on-state current coefficient of variation for a stable
        source.
    on_threshold_frac : float
        Fraction of the peak current above which the source counts as
        "on"; on-state statistics use only those samples.

    Returns
    -------
    SourceStability
    """
    max_cv = _c.as_positive(max_cv, "max_cv")
    on_threshold_frac = _c.as_probability(
        on_threshold_frac, "on_threshold_frac"
    )
    x = _prep_signal(tx_current)
    flags: list[str] = []
    if x.size < 3:
        return SourceStability(
            current_mean_a=float("nan"),
            current_cv=float("nan"),
            current_drift_a=float("nan"),
            on_fraction=float("nan"),
            stable=False,
            flags=["insufficient_samples"],
        )

    peak = float(np.max(np.abs(x)))
    on_mask = (
        np.abs(x) >= on_threshold_frac * peak
        if peak > 0
        else np.zeros(x.shape, dtype=bool)
    )
    on_fraction = float(np.mean(on_mask))
    on = x[on_mask]
    if on.size < 3:
        on = x  # source never keyed off; use the whole record

    mean = float(np.mean(on))
    std = float(np.std(on))
    cv = float(std / abs(mean)) if mean != 0 else float("inf")
    t = np.arange(on.size, dtype=float)
    slope = float(np.polyfit(t, on, 1)[0]) if np.ptp(t) > 0 else 0.0
    drift = abs(slope) * on.size

    if not np.isfinite(cv) or cv > max_cv:
        flags.append("current_unstable")
    if mean == 0:
        flags.append("no_current")

    voltage_mean = None
    if tx_voltage is not None:
        v = _prep_signal(tx_voltage)
        if v.size:
            voltage_mean = float(np.mean(v))

    stable = len(flags) == 0
    return SourceStability(
        current_mean_a=mean,
        current_cv=cv,
        current_drift_a=drift,
        on_fraction=on_fraction,
        voltage_mean_v=voltage_mean,
        stable=stable,
        flags=flags,
    )


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------
def csamt_edge_report(
    data: Any,
    sample_rate: float,
    *,
    tx_frequencies: Iterable[float],
    offset_m: float,
    resistivity: Any,
    tx_current: Any = None,
    tx_voltage: Any = None,
    near_ratio: float = 1.0,
    far_ratio: float = 3.0,
) -> dict[str, Any]:
    """Run the core CSAMT edge diagnostics on one channel and collate them."""
    lines = list(tx_frequencies)
    comb = detect_transmitter_frequencies(data, sample_rate, lines)
    zones = classify_field_zones(
        lines,
        resistivity,
        offset_m,
        near_ratio=near_ratio,
        far_ratio=far_ratio,
    )
    report: dict[str, Any] = dict(
        transmitter=comb.as_dict(),
        field_zones=zones.as_dict(),
    )
    if tx_current is not None:
        report["source_stability"] = assess_source_stability(
            tx_current, tx_voltage
        ).as_dict()
    return report


def csamt_edge_table(
    reports: dict[str, dict[str, Any]] | Iterable[tuple[str, dict[str, Any]]],
    *,
    api: bool | None = None,
) -> Any:
    """Flatten one or more :func:`csamt_edge_report` results into a table."""
    items = (
        list(reports.items()) if isinstance(reports, dict) else list(reports)
    )
    rows: list[dict[str, Any]] = []
    for channel, report in items:
        comb = report.get("transmitter", {})
        zones = report.get("field_zones", {})
        source = report.get("source_stability", {})
        rows.append(
            dict(
                channel=str(channel).lower(),
                n_tx_expected=comb.get("n_expected"),
                n_tx_detected=comb.get("n_detected"),
                tx_detection_fraction=comb.get("detection_fraction"),
                far_fraction=zones.get("far_fraction"),
                all_far_field=zones.get("all_far_field"),
                correction_recommended=zones.get("correction_recommended"),
                first_far_field_hz=zones.get("first_far_field_hz"),
                source_stable=source.get("stable"),
                current_cv=source.get("current_cv"),
            )
        )
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_csamt_edge_table",
        kind="iot.edge.csamt",
        source=items,
        description="CSAMT controlled-source edge quality-control metrics by channel.",
    )
