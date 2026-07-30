"""Formal telemetry payload schemas.

Raw IoT payloads are free-form dictionaries, which is convenient but
fragile: the same quantity may arrive as ``battery_v``, ``battery_voltage``,
or ``voltage``. These schemas define a canonical field set per packet
kind, absorb common key aliases, and validate ranges, so downstream
monitoring and provenance code can rely on stable keys.

Each schema is a light dataclass exposing:

* ``from_mapping(payload)`` — tolerant parser that pulls known aliases
  into canonical fields and preserves anything else under ``extra``.
* ``as_dict(drop_none=...)`` — flat, serialisable payload with canonical
  keys merged with ``extra``.

Use :func:`validate_payload` / :func:`parse_payload` to route a payload
through the schema registered for a given :class:`~pycsamt.iot.core.PacketKind`.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from ..api.property import PyCSAMTObject
from . import _common as _c
from .core import PacketKind

__all__ = [
    "EventSeverity",
    "TelemetryPayload",
    "HealthPayload",
    "QCPayload",
    "PowerPayload",
    "SourcePayload",
    "SyncPayload",
    "EventPayload",
    "AcquisitionPayload",
    "PAYLOAD_SCHEMAS",
    "schema_for",
    "parse_payload",
    "validate_payload",
]


class EventSeverity(str, Enum):
    """Severity levels for event telemetry."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


def _first(mapping: Mapping[str, Any], aliases: tuple[str, ...]) -> Any:
    """Return the first present, non-null value among *aliases*."""
    for key in aliases:
        if key in mapping and mapping[key] is not None:
            return mapping[key]
    return None


def _known_keys(*alias_groups: tuple[str, ...]) -> set[str]:
    return {key for group in alias_groups for key in group}


@dataclass
class TelemetryPayload(PyCSAMTObject):
    """Base class for canonical telemetry payloads."""

    #: Canonical packet kind the schema describes. Overridden per subclass.
    kind: PacketKind = PacketKind.DATA
    station: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    # ---- shared field alias groups -------------------------------------
    _STATION_ALIASES = ("station", "site", "station_id", "station_name")

    def _base_validate(self) -> None:
        self.station = _c.as_optional_str(self.station, "station")
        if not isinstance(self.extra, dict):
            self.extra = dict(self.extra or {})

    @classmethod
    def _consume_station(cls, data: dict[str, Any]) -> str | None:
        return _c.as_optional_str(
            _first(data, cls._STATION_ALIASES), "station"
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        """Return a flat payload dictionary with canonical keys."""
        raise NotImplementedError  # pragma: no cover - overridden

    def _finish(
        self,
        canonical: dict[str, Any],
        *,
        drop_none: bool,
    ) -> dict[str, Any]:
        out = dict(self.extra)
        for key, value in canonical.items():
            if drop_none and value is None:
                continue
            out[key] = value
        return out


@dataclass
class HealthPayload(TelemetryPayload):
    """Device-health telemetry (battery, temperature, link quality)."""

    kind: PacketKind = PacketKind.HEALTH
    battery_v: float | None = None
    temperature_c: float | None = None
    uptime_s: float | None = None
    free_storage_mb: float | None = None
    rssi_dbm: float | None = None
    firmware: str | None = None

    _BATTERY = (
        "battery_v",
        "battery_voltage",
        "battery_voltage_v",
        "voltage",
        "batt_v",
    )
    _TEMPERATURE = ("temperature_c", "temperature", "temp_c", "temp")
    _UPTIME = ("uptime_s", "uptime", "uptime_seconds")
    _STORAGE = ("free_storage_mb", "storage_mb", "free_mb")
    _RSSI = ("rssi_dbm", "rssi", "signal_dbm")
    _FIRMWARE = ("firmware", "firmware_version", "fw")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        self.battery_v = _c.as_optional_finite_float(
            self.battery_v, "battery_v"
        )
        self.temperature_c = _c.as_optional_finite_float(
            self.temperature_c, "temperature_c"
        )
        if self.uptime_s is not None:
            self.uptime_s = _c.as_nonnegative(self.uptime_s, "uptime_s")
        self.free_storage_mb = _c.as_optional_finite_float(
            self.free_storage_mb, "free_storage_mb"
        )
        self.rssi_dbm = _c.as_optional_finite_float(self.rssi_dbm, "rssi_dbm")
        self.firmware = _c.as_optional_str(self.firmware, "firmware")

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> HealthPayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._BATTERY,
            cls._TEMPERATURE,
            cls._UPTIME,
            cls._STORAGE,
            cls._RSSI,
            cls._FIRMWARE,
        )
        return cls(
            station=cls._consume_station(data),
            battery_v=_first(data, cls._BATTERY),
            temperature_c=_first(data, cls._TEMPERATURE),
            uptime_s=_first(data, cls._UPTIME),
            free_storage_mb=_first(data, cls._STORAGE),
            rssi_dbm=_first(data, cls._RSSI),
            firmware=_first(data, cls._FIRMWARE),
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        return self._finish(
            dict(
                station=self.station,
                battery_v=self.battery_v,
                temperature_c=self.temperature_c,
                uptime_s=self.uptime_s,
                free_storage_mb=self.free_storage_mb,
                rssi_dbm=self.rssi_dbm,
                firmware=self.firmware,
            ),
            drop_none=drop_none,
        )


@dataclass
class QCPayload(TelemetryPayload):
    """Edge quality-control telemetry for one acquisition window."""

    kind: PacketKind = PacketKind.QC
    accepted: bool | None = None
    decision: str | None = None
    snr_db: float | None = None
    finite_coverage: float | None = None
    spike_fraction: float | None = None
    rms: float | None = None
    method: str | None = None
    channels: list[str] = field(default_factory=list)
    frequency_band_hz: tuple[float, float] | None = None
    reasons: list[str] = field(default_factory=list)

    _ACCEPTED = ("accepted", "edge_accepted", "qc_accepted", "ok")
    _DECISION = ("decision", "edge_decision", "qc_decision")
    _SNR = ("snr_db", "snr", "channel_snr_db")
    _COVERAGE = ("finite_coverage", "coverage")
    _SPIKE = ("spike_fraction", "spikes", "spike_frac")
    _RMS = ("rms", "rms_value")
    _METHOD = ("method", "survey_method", "em_method")
    _CHANNELS = ("channels", "channel")
    _BAND = ("frequency_band_hz", "band_hz", "freq_band_hz")
    _REASONS = ("reasons", "qc_reasons")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        if self.accepted is not None:
            self.accepted = _c.as_bool(self.accepted)
        self.decision = _c.as_optional_str(self.decision, "decision")
        if self.decision is not None:
            self.decision = self.decision.lower()
        self.snr_db = _c.as_optional_finite_float(self.snr_db, "snr_db")
        if self.finite_coverage is not None:
            self.finite_coverage = _c.as_probability(
                self.finite_coverage, "finite_coverage"
            )
        if self.spike_fraction is not None:
            self.spike_fraction = _c.as_probability(
                self.spike_fraction, "spike_fraction"
            )
        self.rms = _c.as_optional_finite_float(self.rms, "rms")
        self.method = _c.as_optional_str(self.method, "method")
        if self.method is not None:
            self.method = self.method.lower()
        self.channels = _c.as_channel_list(self.channels)
        self.frequency_band_hz = _normalise_band(self.frequency_band_hz)
        self.reasons = _as_str_list(self.reasons)

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> QCPayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._ACCEPTED,
            cls._DECISION,
            cls._SNR,
            cls._COVERAGE,
            cls._SPIKE,
            cls._RMS,
            cls._METHOD,
            cls._CHANNELS,
            cls._BAND,
            cls._REASONS,
        )
        return cls(
            station=cls._consume_station(data),
            accepted=_first(data, cls._ACCEPTED),
            decision=_first(data, cls._DECISION),
            snr_db=_first(data, cls._SNR),
            finite_coverage=_first(data, cls._COVERAGE),
            spike_fraction=_first(data, cls._SPIKE),
            rms=_first(data, cls._RMS),
            method=_first(data, cls._METHOD),
            channels=_first(data, cls._CHANNELS) or [],
            frequency_band_hz=_first(data, cls._BAND),
            reasons=_first(data, cls._REASONS) or [],
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        return self._finish(
            dict(
                station=self.station,
                accepted=self.accepted,
                decision=self.decision,
                snr_db=self.snr_db,
                finite_coverage=self.finite_coverage,
                spike_fraction=self.spike_fraction,
                rms=self.rms,
                method=self.method,
                channels=list(self.channels),
                frequency_band_hz=(
                    list(self.frequency_band_hz)
                    if self.frequency_band_hz is not None
                    else None
                ),
                reasons=list(self.reasons),
            ),
            drop_none=drop_none,
        )


@dataclass
class PowerPayload(TelemetryPayload):
    """Energy-budget telemetry for a field node."""

    kind: PacketKind = PacketKind.POWER
    battery_v: float | None = None
    state: str | None = None
    runtime_days: float | None = None
    net_wh_per_day: float | None = None
    solar_w: float | None = None
    load_w: float | None = None

    _BATTERY = ("battery_v", "battery_voltage", "voltage")
    _STATE = ("state", "power_state")
    _RUNTIME = ("runtime_days", "autonomy_days", "runtime")
    _NET = ("net_wh_per_day", "net_energy_wh_per_day")
    _SOLAR = ("solar_w", "solar_power_w", "harvest_w")
    _LOAD = ("load_w", "load_power_w", "average_power_w")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        self.battery_v = _c.as_optional_finite_float(
            self.battery_v, "battery_v"
        )
        self.state = _c.as_optional_str(self.state, "state")
        if self.state is not None:
            self.state = self.state.lower()
        self.runtime_days = _c.as_optional_finite_float(
            self.runtime_days, "runtime_days"
        )
        self.net_wh_per_day = _c.as_optional_finite_float(
            self.net_wh_per_day, "net_wh_per_day"
        )
        self.solar_w = _c.as_optional_finite_float(self.solar_w, "solar_w")
        self.load_w = _c.as_optional_finite_float(self.load_w, "load_w")

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> PowerPayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._BATTERY,
            cls._STATE,
            cls._RUNTIME,
            cls._NET,
            cls._SOLAR,
            cls._LOAD,
        )
        return cls(
            station=cls._consume_station(data),
            battery_v=_first(data, cls._BATTERY),
            state=_first(data, cls._STATE),
            runtime_days=_first(data, cls._RUNTIME),
            net_wh_per_day=_first(data, cls._NET),
            solar_w=_first(data, cls._SOLAR),
            load_w=_first(data, cls._LOAD),
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        return self._finish(
            dict(
                station=self.station,
                battery_v=self.battery_v,
                state=self.state,
                runtime_days=self.runtime_days,
                net_wh_per_day=self.net_wh_per_day,
                solar_w=self.solar_w,
                load_w=self.load_w,
            ),
            drop_none=drop_none,
        )


@dataclass
class SourcePayload(TelemetryPayload):
    """Transmitter telemetry for a controlled-source (CSAMT/TDEM) survey.

    Reported by a transmitter node so the receiver-side QC and provenance
    know the state of the source that produced each sounding: the injected
    current and voltage, the frequency currently being transmitted, the
    grounded-dipole geometry, and the keyed on/off state.
    """

    kind: PacketKind = PacketKind.SOURCE
    source_id: str | None = None
    tx_current_a: float | None = None
    tx_voltage_v: float | None = None
    tx_frequency_hz: float | None = None
    tx_power_w: float | None = None
    dipole_length_m: float | None = None
    duty_cycle: float | None = None
    on: bool | None = None
    offset_m: float | None = None
    azimuth_deg: float | None = None

    _SOURCE_ID = ("source_id", "tx_id", "transmitter_id")
    _CURRENT = ("tx_current_a", "current_a", "tx_current", "current")
    _VOLTAGE = ("tx_voltage_v", "voltage_v", "tx_voltage")
    _FREQ = ("tx_frequency_hz", "frequency_hz", "tx_freq_hz", "frequency")
    _POWER = ("tx_power_w", "power_w", "tx_power")
    _DIPOLE = ("dipole_length_m", "tx_dipole_length_m", "ab_length_m", "ab_m")
    _DUTY = ("duty_cycle", "duty")
    _ON = ("on", "keyed_on", "tx_on", "transmitting")
    _OFFSET = ("offset_m", "tx_rx_offset_m", "tx_rx_offset", "rx_offset_m")
    _AZIMUTH = ("azimuth_deg", "tx_azimuth_deg", "bearing_deg")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        self.source_id = _c.as_optional_str(self.source_id, "source_id")
        self.tx_current_a = _c.as_optional_finite_float(
            self.tx_current_a, "tx_current_a"
        )
        self.tx_voltage_v = _c.as_optional_finite_float(
            self.tx_voltage_v, "tx_voltage_v"
        )
        self.tx_frequency_hz = _c.as_optional_positive(
            self.tx_frequency_hz, "tx_frequency_hz"
        )
        self.tx_power_w = _c.as_optional_finite_float(
            self.tx_power_w, "tx_power_w"
        )
        self.dipole_length_m = _c.as_optional_positive(
            self.dipole_length_m, "dipole_length_m"
        )
        if self.duty_cycle is not None:
            self.duty_cycle = _c.as_probability(self.duty_cycle, "duty_cycle")
        if self.on is not None:
            self.on = _c.as_bool(self.on)
        self.offset_m = _c.as_optional_positive(self.offset_m, "offset_m")
        self.azimuth_deg = _c.as_optional_finite_float(
            self.azimuth_deg, "azimuth_deg"
        )

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> SourcePayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._SOURCE_ID,
            cls._CURRENT,
            cls._VOLTAGE,
            cls._FREQ,
            cls._POWER,
            cls._DIPOLE,
            cls._DUTY,
            cls._ON,
            cls._OFFSET,
            cls._AZIMUTH,
        )
        return cls(
            station=cls._consume_station(data),
            source_id=_first(data, cls._SOURCE_ID),
            tx_current_a=_first(data, cls._CURRENT),
            tx_voltage_v=_first(data, cls._VOLTAGE),
            tx_frequency_hz=_first(data, cls._FREQ),
            tx_power_w=_first(data, cls._POWER),
            dipole_length_m=_first(data, cls._DIPOLE),
            duty_cycle=_first(data, cls._DUTY),
            on=_first(data, cls._ON),
            offset_m=_first(data, cls._OFFSET),
            azimuth_deg=_first(data, cls._AZIMUTH),
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        return self._finish(
            dict(
                station=self.station,
                source_id=self.source_id,
                tx_current_a=self.tx_current_a,
                tx_voltage_v=self.tx_voltage_v,
                tx_frequency_hz=self.tx_frequency_hz,
                tx_power_w=self.tx_power_w,
                dipole_length_m=self.dipole_length_m,
                duty_cycle=self.duty_cycle,
                on=self.on,
                offset_m=self.offset_m,
                azimuth_deg=self.azimuth_deg,
            ),
            drop_none=drop_none,
        )


@dataclass
class SyncPayload(TelemetryPayload):
    """Clock-synchronisation telemetry for a field node."""

    kind: PacketKind = PacketKind.SYNC
    offset_ms: float | None = None
    drift_ppm: float | None = None
    jitter_ms: float | None = None
    gps_lock: bool | None = None
    n_reference_points: int | None = None
    reference: str | None = None
    # Controlled-source (CSAMT/CSEM) transmitter-receiver timing lock.
    tx_locked: bool | None = None
    tx_sync_offset_ms: float | None = None
    tx_id: str | None = None

    _OFFSET = ("offset_ms", "clock_offset_ms")
    _DRIFT = ("drift_ppm", "clock_drift_ppm")
    _JITTER = ("jitter_ms", "clock_jitter_ms")
    _GPS = ("gps_lock", "gps_locked", "has_gps")
    _NREF = ("n_reference_points", "n_reference", "n_ref")
    _REFERENCE = ("reference", "clock_reference", "ref")
    _TX_LOCKED = ("tx_locked", "transmitter_locked", "source_locked")
    _TX_OFFSET = ("tx_sync_offset_ms", "tx_offset_ms", "source_offset_ms")
    _TX_ID = ("tx_id", "transmitter_id", "source_id")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        self.offset_ms = _c.as_optional_finite_float(
            self.offset_ms, "offset_ms"
        )
        self.drift_ppm = _c.as_optional_finite_float(
            self.drift_ppm, "drift_ppm"
        )
        if self.jitter_ms is not None:
            self.jitter_ms = _c.as_nonnegative(self.jitter_ms, "jitter_ms")
        if self.gps_lock is not None:
            self.gps_lock = _c.as_bool(self.gps_lock)
        if self.n_reference_points is not None:
            self.n_reference_points = int(self.n_reference_points)
            if self.n_reference_points < 0:
                raise ValueError("n_reference_points must be >= 0.")
        self.reference = _c.as_optional_str(self.reference, "reference")
        if self.tx_locked is not None:
            self.tx_locked = _c.as_bool(self.tx_locked)
        self.tx_sync_offset_ms = _c.as_optional_finite_float(
            self.tx_sync_offset_ms, "tx_sync_offset_ms"
        )
        self.tx_id = _c.as_optional_str(self.tx_id, "tx_id")

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> SyncPayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._OFFSET,
            cls._DRIFT,
            cls._JITTER,
            cls._GPS,
            cls._NREF,
            cls._REFERENCE,
            cls._TX_LOCKED,
            cls._TX_OFFSET,
            cls._TX_ID,
        )
        return cls(
            station=cls._consume_station(data),
            offset_ms=_first(data, cls._OFFSET),
            drift_ppm=_first(data, cls._DRIFT),
            jitter_ms=_first(data, cls._JITTER),
            gps_lock=_first(data, cls._GPS),
            n_reference_points=_first(data, cls._NREF),
            reference=_first(data, cls._REFERENCE),
            tx_locked=_first(data, cls._TX_LOCKED),
            tx_sync_offset_ms=_first(data, cls._TX_OFFSET),
            tx_id=_first(data, cls._TX_ID),
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        return self._finish(
            dict(
                station=self.station,
                offset_ms=self.offset_ms,
                drift_ppm=self.drift_ppm,
                jitter_ms=self.jitter_ms,
                gps_lock=self.gps_lock,
                n_reference_points=self.n_reference_points,
                reference=self.reference,
                tx_locked=self.tx_locked,
                tx_sync_offset_ms=self.tx_sync_offset_ms,
                tx_id=self.tx_id,
            ),
            drop_none=drop_none,
        )


@dataclass
class EventPayload(TelemetryPayload):
    """Discrete field event (state change, alarm, operator note)."""

    kind: PacketKind = PacketKind.EVENT
    event: str | None = None
    severity: EventSeverity | str = EventSeverity.INFO
    message: str | None = None
    code: str | None = None

    _EVENT = ("event", "event_type", "name")
    _SEVERITY = ("severity", "level")
    _MESSAGE = ("message", "detail", "description")
    _CODE = ("code", "event_code")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        self.event = _c.as_optional_str(self.event, "event")
        self.severity = _c.normalise_enum(
            self.severity, EventSeverity, "severity"
        )
        self.message = _c.as_optional_str(self.message, "message")
        self.code = _c.as_optional_str(self.code, "code")

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> EventPayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._EVENT,
            cls._SEVERITY,
            cls._MESSAGE,
            cls._CODE,
        )
        return cls(
            station=cls._consume_station(data),
            event=_first(data, cls._EVENT),
            severity=_first(data, cls._SEVERITY) or EventSeverity.INFO,
            message=_first(data, cls._MESSAGE),
            code=_first(data, cls._CODE),
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        severity = (
            self.severity.value
            if isinstance(self.severity, EventSeverity)
            else str(self.severity)
        )
        return self._finish(
            dict(
                station=self.station,
                event=self.event,
                severity=severity,
                message=self.message,
                code=self.code,
            ),
            drop_none=drop_none,
        )


@dataclass
class AcquisitionPayload(TelemetryPayload):
    """Metadata describing one raw acquisition record/window."""

    kind: PacketKind = PacketKind.DATA
    method: str | None = None
    channels: list[str] = field(default_factory=list)
    sample_rate_hz: float | None = None
    frequency_hz: float | None = None
    frequency_band_hz: tuple[float, float] | None = None
    n_samples: int | None = None
    gain: float | None = None
    duration_s: float | None = None

    _METHOD = ("method", "survey_method", "em_method")
    _CHANNELS = ("channels", "channel")
    _RATE = ("sample_rate_hz", "sample_rate", "fs_hz", "fs")
    _FREQ = ("frequency_hz", "freq_hz", "frequency")
    _BAND = ("frequency_band_hz", "band_hz", "freq_band_hz")
    _NSAMPLES = ("n_samples", "nsamples", "n_sample")
    _GAIN = ("gain", "gain_db")
    _DURATION = ("duration_s", "duration", "window_s")

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self._base_validate()
        self.method = _c.as_optional_str(self.method, "method")
        if self.method is not None:
            self.method = self.method.lower()
        self.channels = _c.as_channel_list(self.channels)
        self.sample_rate_hz = _c.as_optional_positive(
            self.sample_rate_hz, "sample_rate_hz"
        )
        self.frequency_hz = _c.as_optional_positive(
            self.frequency_hz, "frequency_hz"
        )
        self.frequency_band_hz = _normalise_band(self.frequency_band_hz)
        if self.n_samples is not None:
            self.n_samples = int(self.n_samples)
            if self.n_samples < 0:
                raise ValueError("n_samples must be >= 0.")
        self.gain = _c.as_optional_finite_float(self.gain, "gain")
        self.duration_s = _c.as_optional_positive(
            self.duration_s, "duration_s"
        )

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> AcquisitionPayload:
        data = dict(payload or {})
        known = _known_keys(
            cls._STATION_ALIASES,
            cls._METHOD,
            cls._CHANNELS,
            cls._RATE,
            cls._FREQ,
            cls._BAND,
            cls._NSAMPLES,
            cls._GAIN,
            cls._DURATION,
        )
        return cls(
            station=cls._consume_station(data),
            method=_first(data, cls._METHOD),
            channels=_first(data, cls._CHANNELS) or [],
            sample_rate_hz=_first(data, cls._RATE),
            frequency_hz=_first(data, cls._FREQ),
            frequency_band_hz=_first(data, cls._BAND),
            n_samples=_first(data, cls._NSAMPLES),
            gain=_first(data, cls._GAIN),
            duration_s=_first(data, cls._DURATION),
            extra={k: v for k, v in data.items() if k not in known},
        )

    def as_dict(self, *, drop_none: bool = False) -> dict[str, Any]:
        return self._finish(
            dict(
                station=self.station,
                method=self.method,
                channels=list(self.channels),
                sample_rate_hz=self.sample_rate_hz,
                frequency_hz=self.frequency_hz,
                frequency_band_hz=(
                    list(self.frequency_band_hz)
                    if self.frequency_band_hz is not None
                    else None
                ),
                n_samples=self.n_samples,
                gain=self.gain,
                duration_s=self.duration_s,
            ),
            drop_none=drop_none,
        )


def _normalise_band(value: Any) -> tuple[float, float] | None:
    if value is None:
        return None
    try:
        lo, hi = list(value)[:2]
    except Exception as exc:  # noqa: BLE001 - report a clear message
        raise ValueError(
            "frequency_band_hz must be a (low, high) pair."
        ) from exc
    lo = _c.as_positive(lo, "frequency_band_hz[0]")
    hi = _c.as_positive(hi, "frequency_band_hz[1]")
    if lo > hi:
        raise ValueError("frequency_band_hz must be ordered (low <= high).")
    return (lo, hi)


def _as_str_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        # Support ";"-joined reason strings emitted by edge tables.
        parts = [p.strip() for p in value.split(";")]
        return [p for p in parts if p]
    return [str(v) for v in list(value)]


#: Map each packet kind to its canonical payload schema.
PAYLOAD_SCHEMAS: dict[PacketKind, type[TelemetryPayload]] = {
    PacketKind.HEALTH: HealthPayload,
    PacketKind.QC: QCPayload,
    PacketKind.POWER: PowerPayload,
    PacketKind.SYNC: SyncPayload,
    PacketKind.EVENT: EventPayload,
    PacketKind.DATA: AcquisitionPayload,
    PacketKind.SOURCE: SourcePayload,
}


def schema_for(kind: PacketKind | str) -> type[TelemetryPayload]:
    """Return the payload schema registered for *kind*."""
    pkt_kind = _c.normalise_enum(kind, PacketKind, "kind")
    try:
        return PAYLOAD_SCHEMAS[pkt_kind]
    except KeyError as exc:  # pragma: no cover - all kinds are mapped
        raise ValueError(f"No payload schema for kind {pkt_kind!r}.") from exc


def parse_payload(
    kind: PacketKind | str,
    payload: Mapping[str, Any],
) -> TelemetryPayload:
    """Parse *payload* into the schema registered for *kind*."""
    return schema_for(kind).from_mapping(payload)


def validate_payload(
    kind: PacketKind | str,
    payload: Mapping[str, Any],
    *,
    drop_none: bool = True,
) -> dict[str, Any]:
    """Return a canonicalised payload dictionary for *kind*.

    Aliased keys are folded into canonical names and values are range
    checked. Unknown keys are preserved. Raises :class:`ValueError` when
    a field fails validation.
    """
    return parse_payload(kind, payload).as_dict(drop_none=drop_none)
