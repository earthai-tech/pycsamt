"""Field acquisition session — the hub of the IoT subpackage.

A :class:`FieldSession` connects devices, stations, telemetry packets,
edge QC, and synchronisation, and exports the result into pyCSAMT-facing
structures: station descriptors (:meth:`to_sites`), a processing-pipeline
input (:meth:`to_pipeline_input`), and a reproducible provenance manifest
(:meth:`export_manifest`).

Example
-------
>>> from pycsamt.iot import FieldSession, DeviceConfig
>>> session = FieldSession(
...     "SSL2026",
...     devices=[DeviceConfig("node-1", station="S01", channels=["ex", "hy"])],
... )
>>> _ = session.add_packet(
...     {
...         "device_id": "node-1",
...         "timestamp": 100.0,
...         "topic": "t",
...         "kind": "qc",
...         "payload": {
...             "station": "S01",
...             "accepted": True,
...             "channels": ["ex", "hy"],
...             "frequency_band_hz": [1.0, 1000.0],
...         },
...     }
... )
>>> status = session.assess()
>>> sites = session.to_sites()
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import (
    Any,
)

from ..api.property import MetadataMixin, PyCSAMTObject
from . import _common as _c
from .core import DeviceConfig, PacketKind, TelemetryPacket
from .monitoring import (
    MonitoringConfig,
    MonitoringStatus,
    TelemetryMonitor,
)
from .provenance import (
    AcquisitionManifest,
    ProvenanceRecord,
    log_qc_decision,
)
from .schemas import parse_payload
from .station import StationConfig, station_table
from .sync import SyncConfig

__all__ = ["FieldSession"]


@dataclass
class _StationAggregate:
    """Rolling per-station summary derived from telemetry packets."""

    station_id: str
    channels: list[str] = field(default_factory=list)
    method: str | None = None
    n_packets: int = 0
    n_accept: int = 0
    n_reject: int = 0
    accepted_band_hz: tuple | None = None
    battery_v: float | None = None
    sample_rate_hz: float | None = None
    first_ts: float | None = None
    last_ts: float | None = None
    last_decision: str | None = None
    qc_decisions: list[dict[str, Any]] = field(default_factory=list)

    @property
    def acceptance_rate(self) -> float | None:
        total = self.n_accept + self.n_reject
        return (self.n_accept / total) if total else None


class FieldSession(PyCSAMTObject, MetadataMixin):
    """Stateful container for one IoT-enabled field acquisition session."""

    __repr_fields__ = ("survey_id", "n_devices", "n_stations", "n_packets")

    def __init__(
        self,
        survey_id: str,
        *,
        devices: Iterable[DeviceConfig | Mapping[str, Any]] | None = None,
        stations: Iterable[StationConfig | Mapping[str, Any]] | None = None,
        monitoring_config: MonitoringConfig | None = None,
        sync_config: SyncConfig | None = None,
        operator: str | None = None,
        method: str | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> None:
        self.survey_id = _c.as_nonempty_str(survey_id, "survey_id")
        self.devices: dict[str, DeviceConfig] = {}
        self.stations: dict[str, StationConfig] = {}
        self.packets: list[TelemetryPacket] = []
        self.monitoring_config = monitoring_config or MonitoringConfig()
        self.sync_config = sync_config or SyncConfig()
        self.operator = _c.as_optional_str(operator, "operator")
        self.method = _c.as_optional_str(method, "method")
        self.metadata = dict(metadata or {})
        for device in devices or []:
            self.add_device(device)
        for station in stations or []:
            self.add_station(station)

    # -- properties --------------------------------------------------------
    @property
    def n_devices(self) -> int:
        return len(self.devices)

    @property
    def n_stations(self) -> int:
        return len(self.stations)

    @property
    def n_packets(self) -> int:
        return len(self.packets)

    def validate(self) -> None:
        """Validate the session (required by :class:`PyCSAMTObject`)."""
        self.survey_id = _c.as_nonempty_str(self.survey_id, "survey_id")

    # -- registration ------------------------------------------------------
    def add_device(
        self, device: DeviceConfig | Mapping[str, Any]
    ) -> DeviceConfig:
        """Register a device, auto-creating its station when named."""
        cfg = (
            device
            if isinstance(device, DeviceConfig)
            else DeviceConfig.from_mapping(device)
        )
        self.devices[cfg.device_id] = cfg
        if cfg.station:
            station = self.stations.get(cfg.station)
            if station is None:
                station = StationConfig(
                    station_id=cfg.station,
                    channels=list(cfg.channels),
                    device_ids=[cfg.device_id],
                )
                self.stations[station.station_id] = station
            else:
                station.attach_device(cfg.device_id)
                for channel in cfg.channels:
                    if channel not in station.channels:
                        station.channels.append(channel)
        return cfg

    def add_station(
        self, station: StationConfig | Mapping[str, Any]
    ) -> StationConfig:
        """Register (or merge into) a station."""
        cfg = (
            station
            if isinstance(station, StationConfig)
            else StationConfig.from_mapping(station)
        )
        existing = self.stations.get(cfg.station_id)
        if existing is None:
            self.stations[cfg.station_id] = cfg
            return cfg
        # Merge device links and channels into the existing station.
        for device_id in cfg.device_ids:
            existing.attach_device(device_id)
        for channel in cfg.channels:
            if channel not in existing.channels:
                existing.channels.append(channel)
        # Explicit geospatial/orientation values override an auto-created
        # placeholder created earlier from a device's ``station`` field.
        for attr in (
            "lat",
            "lon",
            "elevation",
            "profile",
            "position_m",
            "dipole_length_m",
            "ex_azimuth_deg",
            "ey_azimuth_deg",
            "operator",
        ):
            value = getattr(cfg, attr)
            if value is not None:
                setattr(existing, attr, value)
        if cfg.notes:
            existing.notes = cfg.notes
        existing.validate()
        return existing

    def add_packet(
        self, packet: TelemetryPacket | Mapping[str, Any]
    ) -> TelemetryPacket:
        """Append one telemetry packet, registering unknown devices."""
        pkt = (
            packet
            if isinstance(packet, TelemetryPacket)
            else TelemetryPacket(**dict(packet))
        )
        if pkt.device_id not in self.devices:
            self.devices[pkt.device_id] = DeviceConfig(device_id=pkt.device_id)
        self.packets.append(pkt)
        # Register the station this telemetry belongs to so the session
        # inventory (``n_stations``) reflects everything that reported in.
        station_id = self._resolve_station_id(pkt)
        if station_id is not None:
            station = self.stations.get(station_id)
            if station is None:
                station = StationConfig(station_id=station_id)
                self.stations[station_id] = station
            station.attach_device(pkt.device_id)
        return pkt

    def add_packets(
        self, packets: Iterable[TelemetryPacket | Mapping[str, Any]]
    ) -> None:
        """Append many telemetry packets."""
        for packet in packets:
            self.add_packet(packet)

    # -- assessment --------------------------------------------------------
    def assess(self, *, now: float | None = None) -> MonitoringStatus:
        """Assess telemetry quality with the session monitoring config."""
        return TelemetryMonitor(self.monitoring_config).assess(
            self.packets, now=now
        )

    def station_table(self, *, api: bool | None = None) -> Any:
        """Return the registered stations as a pyCSAMT table."""
        return station_table(list(self.stations.values()), api=api)

    # -- pyCSAMT-facing exports -------------------------------------------
    def to_sites(self) -> list[StationConfig]:
        """Return acquisition-level station descriptors.

        These are :class:`StationConfig` objects enriched from telemetry
        (channels and method observed in packets). They are acquisition
        descriptors, not EDI-backed :class:`pycsamt.site.base.Site`
        objects, which require processed impedance data.

        Once impedance is available, :meth:`to_edifiles` and
        :meth:`to_sites_collection` promote these descriptors into
        EDI-backed sites ready for the processing pipeline.
        """
        aggregates = self._station_aggregates()
        sites: list[StationConfig] = []
        for station_id in self._ordered_station_ids():
            station = self.stations.get(station_id) or StationConfig(
                station_id=station_id
            )
            agg = aggregates.get(station_id)
            if agg is not None:
                for channel in agg.channels:
                    if channel not in station.channels:
                        station.channels.append(channel)
            sites.append(station)
        return sites

    def to_pipeline_input(self) -> dict[str, Any]:
        """Return a serialisable hand-off for the processing pipeline.

        Contains, per station, the location, channels, observed method,
        accepted frequency band, and QC acceptance — everything the
        downstream pyCSAMT processing flow needs to know from acquisition.
        """
        aggregates = self._station_aggregates()
        stations_out: list[dict[str, Any]] = []
        for station_id in self._ordered_station_ids():
            station = self.stations.get(station_id)
            agg = aggregates.get(station_id)
            channels = list(station.channels) if station else []
            if agg is not None:
                for channel in agg.channels:
                    if channel not in channels:
                        channels.append(channel)
            stations_out.append(
                dict(
                    station_id=station_id,
                    coords=(station.coords if station else None),
                    profile=(station.profile if station else None),
                    position_m=(station.position_m if station else None),
                    channels=channels,
                    method=(agg.method if agg else None),
                    accepted_band_hz=(
                        list(agg.accepted_band_hz)
                        if agg and agg.accepted_band_hz
                        else None
                    ),
                    acceptance_rate=(agg.acceptance_rate if agg else None),
                    n_packets=(agg.n_packets if agg else 0),
                    last_decision=(agg.last_decision if agg else None),
                )
            )
        return dict(
            survey_id=self.survey_id,
            method=self.method,
            operator=self.operator,
            n_stations=len(stations_out),
            stations=stations_out,
        )

    def to_edifiles(
        self,
        impedance: Mapping[str, Any],
        freq: Any = None,
        *,
        method: str | None = None,
        savepath: str | None = None,
        write: bool = False,
        acqby: str = "pycsamt.iot",
    ) -> dict[str, Any]:
        """Promote per-station impedance into ``EDIFile`` objects.

        A thin convenience over
        :func:`pycsamt.iot.bridge.field_session_to_edifiles`: each EDI is
        enriched with the geometry this session recorded for the station.
        See that function for the accepted *impedance* mapping forms and
        the meaning of *write*/*savepath*.
        """
        from .bridge import field_session_to_edifiles

        return field_session_to_edifiles(
            self,
            impedance,
            freq,
            savepath=savepath,
            method=method or self.method,
            write=write,
            acqby=acqby,
        )

    def to_sites_collection(
        self,
        impedance: Mapping[str, Any],
        freq: Any = None,
        *,
        method: str | None = None,
    ) -> Any:
        """Return an EDI-backed :class:`pycsamt.site.base.Sites` collection.

        This is the processing hand-off: given per-station impedance, it
        builds one EDI per station and wraps them into a ``Sites``
        collection that can be fed straight to
        :meth:`pycsamt.pipeline.Pipeline.run`. Requires the optional
        geospatial stack (``pyproj``).
        """
        edifiles = self.to_edifiles(
            impedance, freq, method=method, write=False
        )
        try:
            from ..site.base import Sites
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "to_sites_collection requires the pyCSAMT geospatial stack "
                "(e.g. pyproj). Use to_edifiles for the raw EDIFile objects, "
                "which have no geospatial dependency."
            ) from exc
        return Sites(list(edifiles.values()))

    def to_manifest(self) -> AcquisitionManifest:
        """Build a reproducible :class:`AcquisitionManifest`."""
        aggregates = self._station_aggregates()
        records: list[ProvenanceRecord] = []
        all_qc: list[dict[str, Any]] = []
        method = self.method
        for station_id in self._ordered_station_ids():
            station = self.stations.get(station_id)
            agg = aggregates.get(station_id)
            if agg is not None and agg.method and method is None:
                method = agg.method
            record = ProvenanceRecord(
                station_id=station_id,
                operator=self.operator,
                lat=(station.lat if station else None),
                lon=(station.lon if station else None),
                elevation=(station.elevation if station else None),
                ex_azimuth_deg=(station.ex_azimuth_deg if station else None),
                ey_azimuth_deg=(station.ey_azimuth_deg if station else None),
                occupation_start=(agg.first_ts if agg else None),
                occupation_end=(agg.last_ts if agg else None),
                sample_rate_hz=(agg.sample_rate_hz if agg else None),
                battery_status=(
                    f"{agg.battery_v:.2f} V"
                    if agg and agg.battery_v is not None
                    else None
                ),
                accepted_band_hz=(agg.accepted_band_hz if agg else None),
                qc_decisions=(list(agg.qc_decisions) if agg else []),
            )
            records.append(record)
            if agg is not None:
                all_qc.extend(agg.qc_decisions)
        return AcquisitionManifest(
            survey_id=self.survey_id,
            method=method,
            operator=self.operator,
            records=records,
            devices=[d.as_dict() for d in self.devices.values()],
            qc_decisions=all_qc,
            metadata=dict(self.metadata),
        )

    def export_manifest(self, path: str, *, indent: int = 2) -> str:
        """Write the acquisition manifest to *path* and return its path."""
        return self.to_manifest().write(path, indent=indent)

    # -- serialisation -----------------------------------------------------
    def to_dict(self) -> dict[str, Any]:  # type: ignore[override]
        """Return a serialisable representation of the whole session."""
        return dict(
            survey_id=self.survey_id,
            operator=self.operator,
            method=self.method,
            devices=[d.as_dict() for d in self.devices.values()],
            stations=[s.as_dict() for s in self.stations.values()],
            packets=[p.as_dict() for p in self.packets],
            metadata=dict(self.metadata),
        )

    def to_json(self, *, indent: int = 2) -> str:
        """Return the session as a JSON string."""
        return json.dumps(self.to_dict(), indent=indent, default=str)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FieldSession:
        """Reconstruct a session from :meth:`to_dict` output."""
        data = dict(data)
        session = cls(
            survey_id=data["survey_id"],
            operator=data.get("operator"),
            method=data.get("method"),
            metadata=data.get("metadata"),
        )
        for device in data.get("devices", []):
            session.add_device(device)
        for station in data.get("stations", []):
            session.add_station(station)
        for packet in data.get("packets", []):
            session.add_packet(packet)
        return session

    @classmethod
    def from_json(cls, text: str) -> FieldSession:
        """Reconstruct a session from a JSON string."""
        return cls.from_dict(json.loads(text))

    def write_json(self, path: str, *, indent: int = 2) -> str:
        """Write the session JSON to *path* and return its path."""
        import os

        parent = os.path.dirname(os.path.abspath(path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(self.to_json(indent=indent))
        return os.path.abspath(path)

    @classmethod
    def read_json(cls, path: str) -> FieldSession:
        """Load a session from a JSON file."""
        with open(path, encoding="utf-8") as handle:
            return cls.from_json(handle.read())

    # -- internals ---------------------------------------------------------
    def _resolve_station_id(self, packet: TelemetryPacket) -> str | None:
        payload = packet.payload or {}
        for key in ("station", "site", "station_id", "station_name"):
            value = payload.get(key)
            if value:
                return str(value)
        device = self.devices.get(packet.device_id)
        if device is not None and device.station:
            return device.station
        return None

    def _ordered_station_ids(self) -> list[str]:
        ids = list(self.stations.keys())
        seen = set(ids)
        for packet in self.packets:
            sid = self._resolve_station_id(packet)
            if sid and sid not in seen:
                seen.add(sid)
                ids.append(sid)

        def _key(sid: str):
            station = self.stations.get(sid)
            profile = station.profile if station else None
            position = (
                station.position_m
                if station and station.position_m is not None
                else float("inf")
            )
            return (profile or "", position, sid)

        return sorted(ids, key=_key)

    def _station_aggregates(self) -> dict[str, _StationAggregate]:
        aggregates: dict[str, _StationAggregate] = {}
        for packet in self.packets:
            sid = self._resolve_station_id(packet)
            if sid is None:
                continue
            agg = aggregates.setdefault(sid, _StationAggregate(station_id=sid))
            agg.n_packets += 1
            ts = float(packet.timestamp)
            agg.first_ts = (
                ts if agg.first_ts is None else min(agg.first_ts, ts)
            )
            agg.last_ts = ts if agg.last_ts is None else max(agg.last_ts, ts)
            parsed = parse_payload(packet.kind, packet.payload)
            self._fold_payload(agg, packet.kind, parsed)
        return aggregates

    def _fold_payload(
        self,
        agg: _StationAggregate,
        kind: PacketKind,
        parsed: Any,
    ) -> None:
        for channel in getattr(parsed, "channels", []) or []:
            if channel not in agg.channels:
                agg.channels.append(channel)
        method = getattr(parsed, "method", None)
        if method and agg.method is None:
            agg.method = method
        rate = getattr(parsed, "sample_rate_hz", None)
        if rate is not None:
            agg.sample_rate_hz = rate
        battery = getattr(parsed, "battery_v", None)
        if battery is not None:
            agg.battery_v = battery

        if kind is PacketKind.QC:
            decision = getattr(parsed, "decision", None)
            accepted = getattr(parsed, "accepted", None)
            if accepted is None and decision is not None:
                accepted = decision in {"accept", "ok", "pass", "warning"}
            if accepted is True:
                agg.n_accept += 1
            elif accepted is False:
                agg.n_reject += 1
            if decision is not None:
                agg.last_decision = decision
            band = getattr(parsed, "frequency_band_hz", None)
            if band is not None and accepted is not False:
                agg.accepted_band_hz = tuple(band)
            agg.qc_decisions.append(
                log_qc_decision(
                    station=agg.station_id,
                    decision=(
                        decision or ("accept" if accepted else "reject")
                    ),
                    reasons=getattr(parsed, "reasons", None),
                    operator=self.operator,
                )
            )
