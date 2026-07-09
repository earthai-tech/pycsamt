from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Iterable, List, Mapping, Optional

import numpy as np
import pandas as pd

from ..api.property import PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c
from .core import DeviceConfig, PacketKind, TelemetryPacket


class EdgeDecision(str, Enum):
    """Acceptance state assigned by edge-side quality control.

    ``ACCEPT`` and ``REJECT`` are the hard decisions. ``WARNING`` marks a
    window that is usable but marginal, ``REPEAT`` requests re-occupation,
    and ``UNKNOWN`` is used when QC could not be evaluated.
    """

    ACCEPT = "accept"
    REJECT = "reject"
    WARNING = "warning"
    REPEAT = "repeat"
    UNKNOWN = "unknown"


def _as_probability(value: Any, name: str) -> float:
    out = float(value)
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1.")
    return out


def _as_positive_int(value: Any, name: str) -> int:
    out = int(value)
    if out < 1:
        raise ValueError(f"{name} must be >= 1.")
    return out


def _as_bool(value: Any) -> bool:
    # Parse string booleans explicitly so "false"/"0"/"off" are not truthy.
    return _c.as_bool(value)


def _finite_rms(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(finite ** 2)))


def _finite_stat(values: np.ndarray, fn: str) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    if fn == "mean":
        return float(np.mean(finite))
    if fn == "std":
        return float(np.std(finite))
    if fn == "min":
        return float(np.min(finite))
    if fn == "max":
        return float(np.max(finite))
    raise ValueError(f"Unknown statistic {fn!r}.")


def _robust_spike_fraction(
    values: np.ndarray,
    *,
    threshold: Optional[float],
) -> float:
    if threshold is None:
        return 0.0
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 3:
        return 0.0
    centre = float(np.median(finite))
    mad = float(np.median(np.abs(finite - centre)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.std(finite))
    if not np.isfinite(scale) or scale <= 0:
        return 0.0
    spikes = np.abs(finite - centre) > float(threshold) * scale
    return float(np.mean(spikes))


def _channel_names(
    n_channel: int,
    provided: Optional[Iterable[str]],
) -> List[str]:
    if provided is None:
        return [f"ch{i}" for i in range(n_channel)]
    names = [str(name).strip().lower() for name in provided]
    if len(names) != n_channel:
        raise ValueError(
            "channel_names length must match the number of channels."
        )
    if any(not name for name in names):
        raise ValueError("channel_names cannot contain empty labels.")
    return names


@dataclass
class EdgeProcessingConfig(PyCSAMTObject):
    """Configuration for lightweight field-side processing."""

    decimation: int = 1
    finite_threshold: float = 0.9
    sample_axis: int = 0
    channel_names: Optional[List[str]] = None
    compute_rms: bool = True
    compute_coverage: bool = True
    compute_spikes: bool = True
    spike_threshold: Optional[float] = 6.0
    max_spike_fraction: float = 0.05
    warn_finite_threshold: Optional[float] = None
    warn_spike_fraction: Optional[float] = None
    retain_payload_samples: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise edge-processing options."""
        self.decimation = _as_positive_int(self.decimation, "decimation")
        self.finite_threshold = _as_probability(
            self.finite_threshold,
            "finite_threshold",
        )
        self.sample_axis = int(self.sample_axis)
        if self.channel_names is not None:
            self.channel_names = [
                str(name).strip().lower()
                for name in self.channel_names
            ]
            if any(not name for name in self.channel_names):
                raise ValueError("channel_names cannot contain empty labels.")
        self.compute_rms = _as_bool(self.compute_rms)
        self.compute_coverage = _as_bool(self.compute_coverage)
        self.compute_spikes = _as_bool(self.compute_spikes)
        if self.spike_threshold is not None:
            self.spike_threshold = float(self.spike_threshold)
            if self.spike_threshold <= 0:
                raise ValueError("spike_threshold must be positive.")
        self.max_spike_fraction = _as_probability(
            self.max_spike_fraction,
            "max_spike_fraction",
        )
        if self.warn_finite_threshold is not None:
            self.warn_finite_threshold = _as_probability(
                self.warn_finite_threshold,
                "warn_finite_threshold",
            )
            if self.warn_finite_threshold < self.finite_threshold:
                raise ValueError(
                    "warn_finite_threshold must be >= finite_threshold."
                )
        if self.warn_spike_fraction is not None:
            self.warn_spike_fraction = _as_probability(
                self.warn_spike_fraction,
                "warn_spike_fraction",
            )
            if self.warn_spike_fraction > self.max_spike_fraction:
                raise ValueError(
                    "warn_spike_fraction must be <= max_spike_fraction."
                )
        self.retain_payload_samples = _as_bool(self.retain_payload_samples)
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})


@dataclass
class EdgeChannelSummary(PyCSAMTObject):
    """Quality-control summary for one edge data channel."""

    channel: str
    n_sample: int
    finite_coverage: float
    rms: float
    mean: float
    std: float
    min: float
    max: float
    spike_fraction: float
    accepted: bool
    reasons: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise the channel summary."""
        self.channel = str(self.channel).strip().lower()
        if not self.channel:
            raise ValueError("channel cannot be empty.")
        self.n_sample = int(self.n_sample)
        if self.n_sample < 0:
            raise ValueError("n_sample must be >= 0.")
        self.finite_coverage = _as_probability(
            self.finite_coverage,
            "finite_coverage",
        )
        self.rms = float(self.rms)
        self.mean = float(self.mean)
        self.std = float(self.std)
        self.min = float(self.min)
        self.max = float(self.max)
        self.spike_fraction = _as_probability(
            self.spike_fraction,
            "spike_fraction",
        )
        self.accepted = bool(self.accepted)
        self.reasons = [str(reason) for reason in list(self.reasons or [])]

    def as_dict(self) -> Dict[str, Any]:
        """Return a serialisable channel summary."""
        return dict(
            channel=self.channel,
            n_sample=self.n_sample,
            finite_coverage=self.finite_coverage,
            rms=self.rms,
            mean=self.mean,
            std=self.std,
            min=self.min,
            max=self.max,
            spike_fraction=self.spike_fraction,
            accepted=self.accepted,
            reasons=";".join(self.reasons),
        )


@dataclass
class EdgeProcessingResult(PyCSAMTObject):
    """Summary returned by :class:`EdgeProcessor`."""

    data: np.ndarray
    metrics: Dict[str, Any]
    accepted: bool
    reasons: List[str] = field(default_factory=list)
    channels: List[EdgeChannelSummary] = field(default_factory=list)
    decision_override: Optional[EdgeDecision | str] = None

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise the processing result."""
        self.data = np.asarray(self.data)
        self.metrics = dict(self.metrics or {})
        self.accepted = bool(self.accepted)
        self.reasons = [str(reason) for reason in list(self.reasons or [])]
        self.channels = [
            ch if isinstance(ch, EdgeChannelSummary)
            else EdgeChannelSummary(**dict(ch))
            for ch in list(self.channels or [])
        ]
        if self.decision_override is not None:
            self.decision_override = _c.normalise_enum(
                self.decision_override, EdgeDecision, "decision_override"
            )

    @property
    def decision(self) -> EdgeDecision:
        """Return the compact QC decision.

        Falls back to a plain accept/reject derived from ``accepted``
        unless an explicit ``decision_override`` (e.g. ``WARNING``) was
        assigned by the processor.
        """
        if self.decision_override is not None:
            return self.decision_override
        return EdgeDecision.ACCEPT if self.accepted else EdgeDecision.REJECT

    def payload(self, *, include_data: Optional[bool] = None) -> Dict[str, Any]:
        """Return a serialisable QC payload for telemetry."""
        keep_data = bool(include_data)
        out = dict(
            metrics=dict(self.metrics),
            accepted=bool(self.accepted),
            decision=self.decision.value,
            reasons=list(self.reasons),
            channels=[ch.as_dict() for ch in self.channels],
        )
        if keep_data:
            out["data"] = np.asarray(self.data).tolist()
        return out

    def to_packet(
        self,
        device: DeviceConfig,
        *,
        timestamp: float,
        survey_id: Optional[str] = None,
        qos: int = 0,
        retained: bool = False,
        include_data: Optional[bool] = None,
    ) -> TelemetryPacket:
        """Encode this edge result as a QC telemetry packet."""
        return TelemetryPacket.from_device(
            device,
            timestamp=timestamp,
            payload=self.payload(include_data=include_data),
            kind=PacketKind.QC,
            survey_id=survey_id,
            qos=qos,
            retained=retained,
        )


class EdgeProcessor(PyCSAMTObject):
    """Small edge-processing block for telemetry payload reduction."""

    def __init__(self, config: Optional[EdgeProcessingConfig] = None) -> None:
        self.config = config or EdgeProcessingConfig()
        self.config.validate()

    def process(
        self,
        data: Any,
        *,
        channel_names: Optional[Iterable[str]] = None,
    ) -> EdgeProcessingResult:
        """Decimate numeric data and compute simple quality metrics."""
        self.config.validate()
        arr = np.asarray(data, dtype=float)
        if arr.ndim == 0:
            arr = arr.reshape(1)
        sample_axis = self._normalise_sample_axis(arr.ndim)
        indices = np.arange(0, arr.shape[sample_axis], self.config.decimation)
        reduced = np.take(arr, indices, axis=sample_axis).copy()
        finite = np.isfinite(reduced)
        coverage = float(np.mean(finite)) if reduced.size else 0.0
        summaries = self._summarise_channels(
            reduced,
            sample_axis=sample_axis,
            channel_names=channel_names,
        )
        max_spike = max(
            (summary.spike_fraction for summary in summaries),
            default=0.0,
        )
        reasons = self._decision_reasons(
            reduced,
            finite_coverage=coverage,
            max_spike_fraction=max_spike,
        )
        metrics: Dict[str, Any] = dict(
            original_samples=int(arr.shape[sample_axis]),
            emitted_samples=int(reduced.shape[sample_axis]),
            finite_coverage=coverage,
            n_channel=len(summaries),
            spike_fraction_max=max_spike,
        )
        if self.config.compute_rms:
            metrics["rms"] = _finite_rms(reduced)
        accepted = len(reasons) == 0
        warnings = self._warning_reasons(
            finite_coverage=coverage,
            max_spike_fraction=max_spike,
        ) if accepted else []
        if accepted and warnings:
            decision = EdgeDecision.WARNING
        elif accepted:
            decision = EdgeDecision.ACCEPT
        else:
            decision = EdgeDecision.REJECT
        metrics["accepted"] = bool(accepted)
        metrics["decision"] = decision.value
        metrics["reasons"] = ";".join(reasons)
        if warnings:
            metrics["warnings"] = ";".join(warnings)
        return EdgeProcessingResult(
            data=reduced,
            metrics=metrics,
            accepted=bool(accepted),
            reasons=reasons,
            channels=summaries,
            decision_override=(
                decision if decision is EdgeDecision.WARNING else None
            ),
        )

    def _warning_reasons(
        self,
        *,
        finite_coverage: float,
        max_spike_fraction: float,
    ) -> List[str]:
        """Return marginal-quality reasons for an otherwise accepted result."""
        warnings: List[str] = []
        warn_cov = self.config.warn_finite_threshold
        warn_spike = self.config.warn_spike_fraction
        if warn_cov is not None and finite_coverage < warn_cov:
            warnings.append("finite_coverage_marginal")
        if warn_spike is not None and max_spike_fraction > warn_spike:
            warnings.append("spike_fraction_marginal")
        return warnings

    def _normalise_sample_axis(self, ndim: int) -> int:
        axis = int(self.config.sample_axis)
        if axis < 0:
            axis += ndim
        if axis < 0 or axis >= ndim:
            raise ValueError("sample_axis is outside the data dimensions.")
        return axis

    def _summarise_channels(
        self,
        data: np.ndarray,
        *,
        sample_axis: int,
        channel_names: Optional[Iterable[str]],
    ) -> List[EdgeChannelSummary]:
        sample_first = np.moveaxis(np.asarray(data, dtype=float), sample_axis, 0)
        if sample_first.ndim == 1:
            matrix = sample_first.reshape(sample_first.shape[0], 1)
        else:
            matrix = sample_first.reshape(sample_first.shape[0], -1)
        names = _channel_names(
            matrix.shape[1],
            channel_names or self.config.channel_names,
        )
        summaries: List[EdgeChannelSummary] = []
        for idx, name in enumerate(names):
            values = matrix[:, idx]
            finite = np.isfinite(values)
            coverage = float(np.mean(finite)) if values.size else 0.0
            spike_fraction = (
                _robust_spike_fraction(
                    values,
                    threshold=self.config.spike_threshold,
                )
                if self.config.compute_spikes else 0.0
            )
            reasons: List[str] = []
            if self.config.compute_coverage and (
                coverage < self.config.finite_threshold
            ):
                reasons.append("finite_coverage_below_threshold")
            if self.config.compute_spikes and (
                spike_fraction > self.config.max_spike_fraction
            ):
                reasons.append("spike_fraction_above_threshold")
            summaries.append(
                EdgeChannelSummary(
                    channel=name,
                    n_sample=int(values.size),
                    finite_coverage=coverage,
                    rms=_finite_rms(values),
                    mean=_finite_stat(values, "mean"),
                    std=_finite_stat(values, "std"),
                    min=_finite_stat(values, "min"),
                    max=_finite_stat(values, "max"),
                    spike_fraction=spike_fraction,
                    accepted=len(reasons) == 0,
                    reasons=reasons,
                )
            )
        return summaries

    def _decision_reasons(
        self,
        data: np.ndarray,
        *,
        finite_coverage: float,
        max_spike_fraction: float,
    ) -> List[str]:
        reasons: List[str] = []
        if data.size == 0:
            reasons.append("empty_data")
        if self.config.compute_coverage and (
            finite_coverage < self.config.finite_threshold
        ):
            reasons.append("finite_coverage_below_threshold")
        if self.config.compute_spikes and (
            max_spike_fraction > self.config.max_spike_fraction
        ):
            reasons.append("spike_fraction_above_threshold")
        return reasons


def edge_summary_table(
    result: EdgeProcessingResult | Iterable[EdgeProcessingResult],
    *,
    api: bool | None = None,
) -> Any:
    """Return one row per channel from edge-processing results."""
    results = (
        [result] if isinstance(result, EdgeProcessingResult)
        else list(result)
    )
    rows: List[Mapping[str, Any]] = []
    for idx, item in enumerate(results):
        item.validate()
        for channel in item.channels:
            row = channel.as_dict()
            row.update(
                result_index=idx,
                decision=item.decision.value,
                result_accepted=item.accepted,
                result_reasons=";".join(item.reasons),
            )
            rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_edge_summary",
        kind="iot.edge",
        source=results,
        description="Edge processing quality-control summaries by channel.",
    )


__all__ = [
    "EdgeDecision",
    "EdgeChannelSummary",
    "EdgeProcessingConfig",
    "EdgeProcessingResult",
    "EdgeProcessor",
    "edge_summary_table",
]
