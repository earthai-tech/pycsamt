from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import (
    Any,
)

import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c


class IoTCapability(str, Enum):
    """IoT capability flags used in deployment reports."""

    TELEMETRY = "telemetry"
    EDGE_PROCESSING = "edge_processing"
    SYNCHRONIZATION = "synchronization"
    ENERGY_MANAGEMENT = "energy_management"
    HEALTH_MONITORING = "health_monitoring"


class DeviceRole(str, Enum):
    """Role of an IoT device in a field deployment."""

    SENSOR_NODE = "sensor_node"
    RECORDER = "recorder"
    GATEWAY = "gateway"
    BASE_STATION = "base_station"
    REMOTE_REFERENCE = "remote_reference"


class PacketKind(str, Enum):
    """Telemetry packet categories used by pyCSAMT IoT workflows."""

    DATA = "data"
    QC = "qc"
    HEALTH = "health"
    SYNC = "sync"
    POWER = "power"
    EVENT = "event"


def _as_nonempty_str(value: Any, name: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} cannot be empty.")
    return text


def _normalise_enum(value: Any, enum_cls: type[Enum], name: str) -> Enum:
    if isinstance(value, enum_cls):
        return value
    text = str(value).strip().lower()
    for item in enum_cls:
        if item.value == text or item.name.lower() == text:
            return item
    allowed = ", ".join(item.value for item in enum_cls)
    raise ValueError(f"{name} must be one of: {allowed}.")


def _normalise_capabilities(
    values: Iterable[IoTCapability | str],
) -> list[IoTCapability]:
    out: list[IoTCapability] = []
    for value in values:
        cap = _normalise_enum(value, IoTCapability, "capability")
        if cap not in out:
            out.append(cap)
    return out


@dataclass
class DeviceConfig(PyCSAMTObject, MetadataMixin):
    """Configuration for one field IoT node or recorder gateway."""

    device_id: str
    station: str | None = None
    protocol: str = "mqtt"
    sample_rate_hz: float | None = None
    channels: list[str] = field(default_factory=list)
    role: DeviceRole | str = DeviceRole.SENSOR_NODE
    enabled: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise device configuration."""
        self.device_id = _as_nonempty_str(self.device_id, "device_id")
        if self.station is not None:
            self.station = _as_nonempty_str(self.station, "station")
        self.protocol = _as_nonempty_str(self.protocol, "protocol").lower()
        if self.sample_rate_hz is not None:
            self.sample_rate_hz = float(self.sample_rate_hz)
            if self.sample_rate_hz <= 0:
                raise ValueError("sample_rate_hz must be positive.")
        self.channels = [
            _as_nonempty_str(ch, "channel").lower()
            for ch in list(self.channels or [])
        ]
        self.role = _normalise_enum(self.role, DeviceRole, "role")
        self.enabled = _c.as_bool(self.enabled)
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})

    @property
    def n_channels(self) -> int:
        """Number of declared acquisition channels."""
        return len(self.channels)

    def topic(
        self,
        kind: PacketKind | str,
        *,
        survey_id: str | None = None,
    ) -> str:
        """Return a canonical telemetry topic for this device."""
        pkt_kind = _normalise_enum(kind, PacketKind, "kind")
        prefix = f"pycsamt/{survey_id}" if survey_id else "pycsamt"
        station = self.station or "unassigned"
        return f"{prefix}/{station}/{self.device_id}/{pkt_kind.value}"

    def as_dict(self) -> dict[str, Any]:
        """Return a serialisable representation."""
        return dict(
            device_id=self.device_id,
            station=self.station,
            protocol=self.protocol,
            sample_rate_hz=self.sample_rate_hz,
            channels=list(self.channels),
            role=self.role.value if isinstance(self.role, DeviceRole)
            else str(self.role),
            enabled=bool(self.enabled),
            metadata=dict(self.metadata),
        )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> DeviceConfig:
        """Build a device config from a mapping."""
        return cls(**dict(data))


@dataclass
class TelemetryPacket(PyCSAMTObject):
    """A timestamped message emitted by an IoT-enabled field node."""

    device_id: str
    timestamp: float
    topic: str
    payload: dict[str, Any]
    qos: int = 0
    kind: PacketKind | str = PacketKind.DATA
    retained: bool = False

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise packet fields."""
        self.device_id = _as_nonempty_str(self.device_id, "device_id")
        # Reject NaN/inf/negative time: a common sign of an un-synced clock.
        self.timestamp = _c.as_timestamp(self.timestamp, "timestamp")
        self.topic = _as_nonempty_str(self.topic, "topic")
        self.payload = dict(self.payload or {})
        self.qos = _c.as_qos(self.qos, "qos")
        self.kind = _normalise_enum(self.kind, PacketKind, "kind")
        self.retained = _c.as_bool(self.retained)

    def as_dict(self) -> dict[str, Any]:
        """Return a flat packet dictionary suitable for logging."""
        return dict(
            device_id=self.device_id,
            timestamp=float(self.timestamp),
            topic=self.topic,
            qos=int(self.qos),
            kind=self.kind.value if isinstance(self.kind, PacketKind)
            else str(self.kind),
            retained=bool(self.retained),
            payload=dict(self.payload),
        )

    @classmethod
    def from_device(
        cls,
        device: DeviceConfig,
        *,
        timestamp: float,
        payload: Mapping[str, Any],
        kind: PacketKind | str = PacketKind.DATA,
        survey_id: str | None = None,
        qos: int = 0,
        retained: bool = False,
    ) -> TelemetryPacket:
        """Create a packet using a device's canonical topic."""
        return cls(
            device_id=device.device_id,
            timestamp=timestamp,
            topic=device.topic(kind, survey_id=survey_id),
            payload=dict(payload),
            qos=qos,
            kind=kind,
            retained=retained,
        )


@dataclass
class DeploymentConfig(PyCSAMTObject, MetadataMixin):
    """Survey-level IoT deployment metadata."""

    survey_id: str
    devices: list[DeviceConfig] = field(default_factory=list)
    capabilities: list[IoTCapability] = field(default_factory=list)
    notes: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise deployment fields."""
        self.survey_id = _as_nonempty_str(self.survey_id, "survey_id")
        self.devices = [
            device if isinstance(device, DeviceConfig)
            else DeviceConfig.from_mapping(device)
            for device in list(self.devices or [])
        ]
        ids = [device.device_id for device in self.devices]
        if len(ids) != len(set(ids)):
            raise ValueError("device_id values must be unique.")
        self.capabilities = _normalise_capabilities(self.capabilities)
        self.notes = str(self.notes or "")
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})

    @property
    def n_devices(self) -> int:
        """Number of devices in the deployment."""
        return len(self.devices)

    def add_device(self, device: DeviceConfig | Mapping[str, Any]) -> None:
        """Append a device and revalidate uniqueness."""
        self.devices.append(
            device if isinstance(device, DeviceConfig)
            else DeviceConfig.from_mapping(device)
        )
        self.validate()

    def has_capability(self, capability: IoTCapability | str) -> bool:
        """Return whether the deployment declares *capability*."""
        cap = _normalise_enum(capability, IoTCapability, "capability")
        return cap in self.capabilities

    def device(self, device_id: str) -> DeviceConfig:
        """Return a device by id."""
        key = _as_nonempty_str(device_id, "device_id")
        for device in self.devices:
            if device.device_id == key:
                return device
        raise KeyError(f"Unknown device_id {key!r}.")

    def as_dict(self) -> dict[str, Any]:
        """Return a serialisable deployment dictionary."""
        return dict(
            survey_id=self.survey_id,
            devices=[device.as_dict() for device in self.devices],
            capabilities=[cap.value for cap in self.capabilities],
            notes=self.notes,
            metadata=dict(self.metadata),
        )


def deployment_report(
    deployment: DeploymentConfig,
    *,
    api: bool | None = None,
) -> Any:
    """Return one row per device describing declared IoT capabilities."""
    deployment.validate()
    caps = [c.value if isinstance(c, IoTCapability) else str(c)
            for c in deployment.capabilities]
    rows = []
    for device in deployment.devices:
        row = device.as_dict()
        row.update(
            survey_id=deployment.survey_id,
            capabilities=";".join(caps),
            n_capabilities=len(caps),
            n_channels=device.n_channels,
            notes=deployment.notes,
        )
        rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_deployment_report",
        kind="iot.deployment",
        source=deployment,
        description="IoT deployment devices and declared capabilities.",
    )


__all__ = [
    "IoTCapability",
    "DeviceRole",
    "PacketKind",
    "DeviceConfig",
    "TelemetryPacket",
    "DeploymentConfig",
    "deployment_report",
]
