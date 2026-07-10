from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import (
    Any,
)

import numpy as np
import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..api.view import maybe_wrap_frame
from .core import PacketKind, TelemetryPacket


class EMMethod(str, Enum):
    """Electromagnetic survey methods recognised by IoT monitoring."""

    AMT = "amt"
    MT = "mt"
    CSAMT = "csamt"
    TDEM = "tdem"
    TEM = "tem"
    UNKNOWN = "unknown"


class MonitoringLevel(str, Enum):
    """Overall status level for a monitored telemetry stream."""

    OK = "ok"
    WARNING = "warning"
    CRITICAL = "critical"
    NO_DATA = "no_data"


def _normalise_enum(value: Any, enum_cls: type[Enum], name: str) -> Enum:
    if isinstance(value, enum_cls):
        return value
    text = str(value).strip().lower()
    for item in enum_cls:
        if text in {item.value, item.name.lower()}:
            return item
    allowed = ", ".join(item.value for item in enum_cls)
    raise ValueError(f"{name} must be one of: {allowed}.")


def _as_probability(value: Any, name: str) -> float:
    out = float(value)
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1.")
    return out


def _as_optional_positive(value: Any, name: str) -> float | None:
    if value is None:
        return None
    out = float(value)
    if out <= 0:
        raise ValueError(f"{name} must be positive.")
    return out


def _payload(packet: TelemetryPacket | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(packet, TelemetryPacket):
        return dict(packet.payload or {})
    payload = dict(packet).get("payload", {})
    return dict(payload or {})


def _packet_dict(packet: TelemetryPacket | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(packet, TelemetryPacket):
        return packet.as_dict()
    out = dict(packet)
    out.setdefault("payload", dict(out.get("payload") or {}))
    out.setdefault("kind", str(out.get("kind", "")) or PacketKind.DATA.value)
    out.setdefault("qos", int(out.get("qos", 0)))
    out.setdefault("retained", bool(out.get("retained", False)))
    return out


def _payload_value(
    payload: Mapping[str, Any],
    keys: Iterable[str],
    default: Any = None,
) -> Any:
    for key in keys:
        if key in payload and payload[key] is not None:
            return payload[key]
    return default


def _as_float_or_nan(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _safe_mean(values: Iterable[Any]) -> float:
    arr = np.asarray([_as_float_or_nan(v) for v in values], dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr)) if arr.size else float("nan")


def _safe_max(values: Iterable[Any]) -> float:
    arr = np.asarray([_as_float_or_nan(v) for v in values], dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.max(arr)) if arr.size else float("nan")


def _safe_min(values: Iterable[Any]) -> float:
    arr = np.asarray([_as_float_or_nan(v) for v in values], dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.min(arr)) if arr.size else float("nan")


def _method_from_payload(payload: Mapping[str, Any]) -> EMMethod:
    value = _payload_value(payload, ("method", "survey_method", "em_method"))
    if value is None:
        return EMMethod.UNKNOWN
    try:
        return _normalise_enum(value, EMMethod, "method")  # type: ignore[return-value]
    except ValueError:
        return EMMethod.UNKNOWN


@dataclass
class MonitoringConfig(PyCSAMTObject, MetadataMixin):
    """Thresholds used to monitor AMT/MT/CSAMT telemetry streams."""

    method: EMMethod | str = EMMethod.UNKNOWN
    expected_interval_s: float | None = None
    max_latency_s: float | None = 30.0
    max_gap_s: float | None = None
    min_packet_success_rate: float = 0.95
    min_edge_acceptance_rate: float = 0.85
    min_battery_v: float | None = None
    max_clock_offset_ms: float | None = None
    frequency_band_hz: tuple[float, float] | None = None
    required_channels: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise monitoring thresholds."""
        self.method = _normalise_enum(self.method, EMMethod, "method")
        self.expected_interval_s = _as_optional_positive(
            self.expected_interval_s,
            "expected_interval_s",
        )
        self.max_latency_s = _as_optional_positive(
            self.max_latency_s,
            "max_latency_s",
        )
        self.max_gap_s = _as_optional_positive(self.max_gap_s, "max_gap_s")
        self.min_packet_success_rate = _as_probability(
            self.min_packet_success_rate,
            "min_packet_success_rate",
        )
        self.min_edge_acceptance_rate = _as_probability(
            self.min_edge_acceptance_rate,
            "min_edge_acceptance_rate",
        )
        self.min_battery_v = _as_optional_positive(
            self.min_battery_v,
            "min_battery_v",
        )
        self.max_clock_offset_ms = _as_optional_positive(
            self.max_clock_offset_ms,
            "max_clock_offset_ms",
        )
        if self.frequency_band_hz is not None:
            lo, hi = [float(v) for v in self.frequency_band_hz]
            if lo <= 0 or hi <= 0 or lo > hi:
                raise ValueError(
                    "frequency_band_hz must be positive and ordered."
                )
            self.frequency_band_hz = (lo, hi)
        self.required_channels = [
            str(ch).strip().lower()
            for ch in list(self.required_channels or [])
        ]
        if any(not ch for ch in self.required_channels):
            raise ValueError("required_channels cannot contain empty labels.")
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})


@dataclass
class MonitoringStatus(PyCSAMTObject):
    """Status returned by :class:`TelemetryMonitor`."""

    level: MonitoringLevel | str
    n_packet: int
    packet_success_rate: float
    edge_acceptance_rate: float
    mean_latency_s: float
    max_gap_s: float
    battery_min_v: float
    clock_offset_max_ms: float
    methods: list[str] = field(default_factory=list)
    stations: list[str] = field(default_factory=list)
    channels: list[str] = field(default_factory=list)
    frequency_min_hz: float = float("nan")
    frequency_max_hz: float = float("nan")
    issues: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise status fields."""
        self.level = _normalise_enum(self.level, MonitoringLevel, "level")
        self.n_packet = int(self.n_packet)
        if self.n_packet < 0:
            raise ValueError("n_packet must be >= 0.")
        self.packet_success_rate = _as_probability(
            self.packet_success_rate,
            "packet_success_rate",
        )
        self.edge_acceptance_rate = _as_probability(
            self.edge_acceptance_rate,
            "edge_acceptance_rate",
        )
        self.mean_latency_s = float(self.mean_latency_s)
        self.max_gap_s = float(self.max_gap_s)
        self.battery_min_v = float(self.battery_min_v)
        self.clock_offset_max_ms = float(self.clock_offset_max_ms)
        self.methods = sorted({str(v).lower() for v in self.methods if v})
        self.stations = sorted({str(v) for v in self.stations if v})
        self.channels = sorted({str(v).lower() for v in self.channels if v})
        self.frequency_min_hz = float(self.frequency_min_hz)
        self.frequency_max_hz = float(self.frequency_max_hz)
        self.issues = [str(issue) for issue in list(self.issues or [])]

    @property
    def ok(self) -> bool:
        """Return whether the monitored stream is acceptable."""
        return self.level == MonitoringLevel.OK

    def as_dict(self) -> dict[str, Any]:
        """Return a flat serialisable status dictionary."""
        return dict(
            level=self.level.value if isinstance(self.level, MonitoringLevel)
            else str(self.level),
            n_packet=self.n_packet,
            packet_success_rate=self.packet_success_rate,
            edge_acceptance_rate=self.edge_acceptance_rate,
            mean_latency_s=self.mean_latency_s,
            max_gap_s=self.max_gap_s,
            battery_min_v=self.battery_min_v,
            clock_offset_max_ms=self.clock_offset_max_ms,
            methods=";".join(self.methods),
            stations=";".join(self.stations),
            channels=";".join(self.channels),
            frequency_min_hz=self.frequency_min_hz,
            frequency_max_hz=self.frequency_max_hz,
            issues=";".join(self.issues),
        )


class TelemetryMonitor(PyCSAMTObject):
    """Monitor field telemetry from AMT, MT, CSAMT, and related surveys."""

    def __init__(self, config: MonitoringConfig | None = None) -> None:
        self.config = config or MonitoringConfig()
        self.config.validate()

    def assess(
        self,
        packets: Iterable[TelemetryPacket | Mapping[str, Any]],
        *,
        now: float | None = None,
    ) -> MonitoringStatus:
        """Assess a telemetry packet stream against configured thresholds."""
        self.config.validate()
        rows = [_packet_dict(packet) for packet in packets]
        if not rows:
            return MonitoringStatus(
                level=MonitoringLevel.NO_DATA,
                n_packet=0,
                packet_success_rate=0.0,
                edge_acceptance_rate=0.0,
                mean_latency_s=float("nan"),
                max_gap_s=float("nan"),
                battery_min_v=float("nan"),
                clock_offset_max_ms=float("nan"),
                issues=["no_packets"],
            )
        enriched = [self._enrich_row(row, now=now) for row in rows]
        issues = self._issues(enriched)
        level = self._level(issues)
        return MonitoringStatus(
            level=level,
            n_packet=len(enriched),
            packet_success_rate=self._packet_success_rate(enriched),
            edge_acceptance_rate=self._edge_acceptance_rate(enriched),
            mean_latency_s=_safe_mean(row.get("latency_s")
                                      for row in enriched),
            max_gap_s=self._max_gap(enriched),
            battery_min_v=_safe_min(row.get("battery_v")
                                    for row in enriched),
            clock_offset_max_ms=_safe_max(
                abs(_as_float_or_nan(row.get("clock_offset_ms")))
                for row in enriched
            ),
            methods=[
                str(row.get("method", EMMethod.UNKNOWN.value))
                for row in enriched
            ],
            stations=[
                str(row.get("station"))
                for row in enriched
                if row.get("station") is not None
            ],
            channels=[
                ch for row in enriched
                for ch in list(row.get("channels") or [])
            ],
            frequency_min_hz=_safe_min(row.get("frequency_min_hz")
                                       for row in enriched),
            frequency_max_hz=_safe_max(row.get("frequency_max_hz")
                                       for row in enriched),
            issues=issues,
        )

    def table(
        self,
        packets: Iterable[TelemetryPacket | Mapping[str, Any]],
        *,
        now: float | None = None,
        api: bool | None = None,
    ) -> Any:
        """Return enriched packet rows used by the monitor."""
        rows = [
            self._enrich_row(_packet_dict(packet), now=now)
            for packet in packets
        ]
        df = pd.DataFrame.from_records(rows)
        return maybe_wrap_frame(
            df,
            api=api,
            name="iot_monitoring_packets",
            kind="iot.monitoring.packets",
            source=packets,
            description="Telemetry packets enriched with EM survey metadata.",
        )

    def _enrich_row(
        self,
        row: Mapping[str, Any],
        *,
        now: float | None,
    ) -> dict[str, Any]:
        out = dict(row)
        payload = dict(out.get("payload") or {})
        method = _method_from_payload(payload)
        station = _payload_value(payload, ("station", "site", "station_id"))
        channels = _payload_value(payload, ("channels", "channel"), [])
        if isinstance(channels, str):
            channels = [channels]
        channels = [str(ch).strip().lower() for ch in list(channels or [])]
        freq = _payload_value(
            payload,
            ("frequency_hz", "freq_hz", "frequency"),
        )
        band = _payload_value(
            payload,
            ("frequency_band_hz", "band_hz", "freq_band_hz"),
        )
        fmin = fmax = _as_float_or_nan(freq)
        if band is not None:
            try:
                f0, f1 = list(band)[:2]
                fmin, fmax = float(f0), float(f1)
            except Exception:
                pass
        latency = _payload_value(payload, ("latency_s", "latency"))
        if latency is None and now is not None:
            latency = float(now) - float(out.get("timestamp"))
        out.update(
            method=method.value,
            station=station,
            channels=channels,
            frequency_min_hz=fmin,
            frequency_max_hz=fmax,
            latency_s=_as_float_or_nan(latency),
            battery_v=_as_float_or_nan(
                _payload_value(payload, ("battery_v", "battery_voltage_v"))
            ),
            clock_offset_ms=_as_float_or_nan(
                _payload_value(payload, ("clock_offset_ms", "offset_ms"))
            ),
            ack_ok=bool(_payload_value(payload, ("ack_ok", "ok"), True)),
            edge_accepted=self._edge_accepted(payload),
        )
        return out

    def _edge_accepted(self, payload: Mapping[str, Any]) -> bool | None:
        value = _payload_value(
            payload,
            ("accepted", "edge_accepted", "qc_accepted"),
        )
        decision = _payload_value(payload, ("decision", "edge_decision"))
        if value is not None:
            return bool(value)
        if decision is not None:
            return str(decision).strip().lower() in {"accept", "ok", "pass"}
        return None

    def _issues(self, rows: list[Mapping[str, Any]]) -> list[str]:
        issues: list[str] = []
        success = self._packet_success_rate(rows)
        edge = self._edge_acceptance_rate(rows)
        max_gap = self._max_gap(rows)
        mean_latency = _safe_mean(row.get("latency_s") for row in rows)
        battery_min = _safe_min(row.get("battery_v") for row in rows)
        clock_max = _safe_max(
            abs(_as_float_or_nan(row.get("clock_offset_ms")))
            for row in rows
        )
        if success < self.config.min_packet_success_rate:
            issues.append("packet_success_rate_below_threshold")
        if edge < self.config.min_edge_acceptance_rate:
            issues.append("edge_acceptance_rate_below_threshold")
        if self.config.max_gap_s is not None and np.isfinite(max_gap):
            if max_gap > self.config.max_gap_s:
                issues.append("packet_gap_above_threshold")
        if self.config.expected_interval_s is not None and np.isfinite(max_gap):
            if max_gap > 2.5 * self.config.expected_interval_s:
                issues.append("packet_gap_exceeds_expected_interval")
        if self.config.max_latency_s is not None and np.isfinite(mean_latency):
            if mean_latency > self.config.max_latency_s:
                issues.append("latency_above_threshold")
        if self.config.min_battery_v is not None and np.isfinite(battery_min):
            if battery_min < self.config.min_battery_v:
                issues.append("battery_below_threshold")
        if (
            self.config.max_clock_offset_ms is not None
            and np.isfinite(clock_max)
            and clock_max > self.config.max_clock_offset_ms
        ):
            issues.append("clock_offset_above_threshold")
        if self.config.method != EMMethod.UNKNOWN:
            methods = {str(row.get("method")) for row in rows}
            if self.config.method.value not in methods:
                issues.append("method_mismatch")
        channels = {
            ch for row in rows for ch in list(row.get("channels") or [])
        }
        missing = set(self.config.required_channels) - channels
        if missing:
            issues.append("required_channels_missing")
        if self.config.frequency_band_hz is not None:
            lo, hi = self.config.frequency_band_hz
            seen = [
                (
                    _as_float_or_nan(row.get("frequency_min_hz")),
                    _as_float_or_nan(row.get("frequency_max_hz")),
                )
                for row in rows
            ]
            in_band = [
                f0 >= lo and f1 <= hi
                for f0, f1 in seen
                if np.isfinite(f0) and np.isfinite(f1)
            ]
            if in_band and not all(in_band):
                issues.append("frequency_outside_configured_band")
        return sorted(dict.fromkeys(issues))

    def _level(self, issues: list[str]) -> MonitoringLevel:
        if not issues:
            return MonitoringLevel.OK
        critical = {
            "packet_success_rate_below_threshold",
            "edge_acceptance_rate_below_threshold",
            "battery_below_threshold",
            "clock_offset_above_threshold",
            "method_mismatch",
            "required_channels_missing",
        }
        if any(issue in critical for issue in issues):
            return MonitoringLevel.CRITICAL
        return MonitoringLevel.WARNING

    def _packet_success_rate(self, rows: list[Mapping[str, Any]]) -> float:
        if not rows:
            return 0.0
        return float(np.mean([bool(row.get("ack_ok", True)) for row in rows]))

    def _edge_acceptance_rate(self, rows: list[Mapping[str, Any]]) -> float:
        values = [
            row.get("edge_accepted")
            for row in rows
            if row.get("edge_accepted") is not None
        ]
        if not values:
            return 1.0
        return float(np.mean([bool(value) for value in values]))

    def _max_gap(self, rows: list[Mapping[str, Any]]) -> float:
        times = sorted(
            _as_float_or_nan(row.get("timestamp"))
            for row in rows
            if np.isfinite(_as_float_or_nan(row.get("timestamp")))
        )
        if len(times) < 2:
            return 0.0
        return float(np.max(np.diff(times)))


def packet_table(
    packets: Iterable[TelemetryPacket | Mapping[str, Any]],
    *,
    api: bool | None = None,
) -> Any:
    """Return telemetry packets as a pyCSAMT table."""
    rows: list[dict] = []
    for packet in packets:
        row = _packet_dict(packet)
        row["payload_keys"] = ";".join(sorted(row["payload"].keys()))
        rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_packet_table",
        kind="iot.telemetry.packets",
        source=packets,
        description="Recorded IoT telemetry packets.",
    )


def telemetry_summary(
    packets: Iterable[TelemetryPacket | Mapping[str, Any]],
    *,
    api: bool | None = None,
) -> Any:
    """Summarise telemetry packet counts by device and topic."""
    rows = [_packet_dict(p) for p in packets]
    if not rows:
        df = pd.DataFrame(columns=["device_id", "topic", "n_packet"])
    else:
        df0 = pd.DataFrame.from_records(rows)
        df = (
            df0.groupby(["device_id", "topic"], dropna=False)
            .size()
            .reset_index(name="n_packet")
        )
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_telemetry_summary",
        kind="iot.telemetry.summary",
        source=packets,
        description="Telemetry packet counts by device and topic.",
    )


def monitoring_status_table(
    statuses: MonitoringStatus | Iterable[MonitoringStatus],
    *,
    api: bool | None = None,
) -> Any:
    """Return monitoring statuses as a pyCSAMT table."""
    rows = (
        [statuses] if isinstance(statuses, MonitoringStatus)
        else list(statuses)
    )
    df = pd.DataFrame.from_records([status.as_dict() for status in rows])
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_monitoring_status",
        kind="iot.monitoring.status",
        source=statuses,
        description="IoT monitoring status for EM field telemetry.",
    )


def assess_telemetry(
    packets: Iterable[TelemetryPacket | Mapping[str, Any]],
    *,
    config: MonitoringConfig | None = None,
    now: float | None = None,
) -> MonitoringStatus:
    """Convenience wrapper around :class:`TelemetryMonitor`."""
    return TelemetryMonitor(config).assess(packets, now=now)


__all__ = [
    "EMMethod",
    "MonitoringConfig",
    "MonitoringLevel",
    "MonitoringStatus",
    "TelemetryMonitor",
    "assess_telemetry",
    "monitoring_status_table",
    "packet_table",
    "telemetry_summary",
]
