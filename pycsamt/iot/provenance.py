"""Acquisition provenance and reproducibility for IoT field sessions.

For publication and reviewer credibility, every IoT field session should
be able to emit a reproducible audit trail: what was recorded, where, by
whom, with which instrument, which windows were rejected, and the exact
bytes of the raw files. This module builds that trail as plain JSON-able
structures.

The routines are dependency-free (standard library only) and do not
import :mod:`pycsamt.iot.session`, so they can be reused independently.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from ..api.property import PyCSAMTObject
from . import _common as _c

__all__ = [
    "ProvenanceRecord",
    "AcquisitionManifest",
    "hash_raw_file",
    "hash_bytes",
    "hash_mapping",
    "log_qc_decision",
    "build_acquisition_manifest",
    "export_acquisition_manifest",
    "export_station_audit",
    "export_reproducibility_bundle",
]


def _tool_version() -> str:
    try:
        from importlib.metadata import version

        return version("pycsamt")
    except Exception:  # noqa: BLE001 - version metadata may be absent
        try:
            import pycsamt

            return getattr(pycsamt, "__version__", "unknown")
        except Exception:  # noqa: BLE001
            return "unknown"


def _utc_iso(timestamp: Optional[float] = None) -> str:
    if timestamp is None:
        dt = datetime.now(timezone.utc)
    else:
        dt = datetime.fromtimestamp(float(timestamp), tz=timezone.utc)
    return dt.isoformat()


def hash_bytes(data: bytes, *, algo: str = "sha256") -> str:
    """Return the hex digest of *data* using *algo*."""
    hasher = hashlib.new(algo)
    hasher.update(data)
    return hasher.hexdigest()


def hash_mapping(mapping: Mapping[str, Any], *, algo: str = "sha256") -> str:
    """Return a stable hash of *mapping* (key-sorted JSON, UTF-8)."""
    encoded = json.dumps(
        dict(mapping), sort_keys=True, default=str, separators=(",", ":")
    ).encode("utf-8")
    return hash_bytes(encoded, algo=algo)


def hash_raw_file(
    path: str,
    *,
    algo: str = "sha256",
    chunk_size: int = 1 << 20,
) -> Dict[str, Any]:
    """Hash a raw acquisition file and return its integrity record.

    Returns a mapping with ``path``, ``bytes``, ``algo``, ``digest``, and
    ``modified_utc``. Raises :class:`FileNotFoundError` when the file is
    missing.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"raw file not found: {path}")
    hasher = hashlib.new(algo)
    size = 0
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            hasher.update(block)
            size += len(block)
    return dict(
        path=os.path.abspath(path),
        name=os.path.basename(path),
        bytes=size,
        algo=algo,
        digest=hasher.hexdigest(),
        modified_utc=_utc_iso(os.path.getmtime(path)),
    )


def log_qc_decision(
    station: str,
    decision: str,
    *,
    channel: Optional[str] = None,
    reasons: Optional[Iterable[str]] = None,
    operator: Optional[str] = None,
    timestamp: Optional[float] = None,
    window: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """Return a normalised QC-decision record for the audit trail."""
    return dict(
        station=_c.as_nonempty_str(station, "station"),
        channel=(str(channel).lower() if channel is not None else None),
        decision=_c.as_nonempty_str(decision, "decision").lower(),
        reasons=list(reasons or []),
        operator=_c.as_optional_str(operator, "operator"),
        logged_utc=_utc_iso(timestamp),
        window=(list(window) if window is not None else None),
    )


@dataclass
class ProvenanceRecord(PyCSAMTObject):
    """Per-station occupation provenance for one field session."""

    station_id: str
    instrument_serial: Optional[str] = None
    firmware: Optional[str] = None
    operator: Optional[str] = None
    lat: Optional[float] = None
    lon: Optional[float] = None
    elevation: Optional[float] = None
    ex_azimuth_deg: Optional[float] = None
    ey_azimuth_deg: Optional[float] = None
    occupation_start: Optional[float] = None
    occupation_end: Optional[float] = None
    sample_rate_hz: Optional[float] = None
    gps_quality: Optional[str] = None
    battery_status: Optional[str] = None
    accepted_band_hz: Optional[Tuple[float, float]] = None
    field_notes: str = ""
    qc_decisions: List[Dict[str, Any]] = field(default_factory=list)
    rejected_windows: List[Any] = field(default_factory=list)
    processing_steps: List[str] = field(default_factory=list)
    raw_files: List[Dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.station_id = _c.as_nonempty_str(self.station_id, "station_id")
        self.lat = _c.as_latitude(self.lat)
        self.lon = _c.as_longitude(self.lon)
        self.elevation = _c.as_elevation(self.elevation)
        if self.sample_rate_hz is not None:
            self.sample_rate_hz = _c.as_positive(
                self.sample_rate_hz, "sample_rate_hz"
            )
        if self.occupation_start is not None:
            self.occupation_start = _c.as_timestamp(
                self.occupation_start, "occupation_start"
            )
        if self.occupation_end is not None:
            self.occupation_end = _c.as_timestamp(
                self.occupation_end, "occupation_end"
            )

    @property
    def occupation_seconds(self) -> Optional[float]:
        """Occupation duration in seconds, if both bounds are known."""
        if self.occupation_start is None or self.occupation_end is None:
            return None
        return float(self.occupation_end - self.occupation_start)

    def add_raw_file(self, path: str, *, algo: str = "sha256") -> Dict[str, Any]:
        """Hash *path* and attach the integrity record."""
        record = hash_raw_file(path, algo=algo)
        self.raw_files.append(record)
        return record

    def as_dict(self) -> Dict[str, Any]:
        return dict(
            station_id=self.station_id,
            instrument_serial=self.instrument_serial,
            firmware=self.firmware,
            operator=self.operator,
            lat=self.lat,
            lon=self.lon,
            elevation=self.elevation,
            ex_azimuth_deg=self.ex_azimuth_deg,
            ey_azimuth_deg=self.ey_azimuth_deg,
            occupation_start=self.occupation_start,
            occupation_end=self.occupation_end,
            occupation_seconds=self.occupation_seconds,
            sample_rate_hz=self.sample_rate_hz,
            gps_quality=self.gps_quality,
            battery_status=self.battery_status,
            accepted_band_hz=(
                list(self.accepted_band_hz)
                if self.accepted_band_hz is not None else None
            ),
            field_notes=self.field_notes,
            qc_decisions=list(self.qc_decisions),
            rejected_windows=list(self.rejected_windows),
            processing_steps=list(self.processing_steps),
            raw_files=list(self.raw_files),
        )


@dataclass
class AcquisitionManifest(PyCSAMTObject):
    """Reproducible manifest describing one IoT acquisition session."""

    survey_id: str
    created_utc: str = field(default_factory=_utc_iso)
    tool: str = "pycsamt.iot"
    tool_version: str = field(default_factory=_tool_version)
    method: Optional[str] = None
    operator: Optional[str] = None
    records: List[ProvenanceRecord] = field(default_factory=list)
    devices: List[Dict[str, Any]] = field(default_factory=list)
    files: List[Dict[str, Any]] = field(default_factory=list)
    qc_decisions: List[Dict[str, Any]] = field(default_factory=list)
    processing_steps: List[str] = field(default_factory=list)
    environment: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.survey_id = _c.as_nonempty_str(self.survey_id, "survey_id")
        self.records = [
            r if isinstance(r, ProvenanceRecord)
            else ProvenanceRecord(**dict(r))
            for r in list(self.records or [])
        ]
        if not self.environment:
            self.environment = _default_environment()

    def add_record(self, record: ProvenanceRecord | Mapping[str, Any]) -> None:
        self.records.append(
            record if isinstance(record, ProvenanceRecord)
            else ProvenanceRecord(**dict(record))
        )

    def add_file(self, path: str, *, algo: str = "sha256") -> Dict[str, Any]:
        record = hash_raw_file(path, algo=algo)
        self.files.append(record)
        return record

    def add_qc_decision(self, **kwargs: Any) -> Dict[str, Any]:
        decision = log_qc_decision(**kwargs)
        self.qc_decisions.append(decision)
        return decision

    def add_processing_step(self, step: str) -> None:
        self.processing_steps.append(_c.as_nonempty_str(step, "step"))

    def as_dict(self) -> Dict[str, Any]:
        payload = dict(
            survey_id=self.survey_id,
            created_utc=self.created_utc,
            tool=self.tool,
            tool_version=self.tool_version,
            method=self.method,
            operator=self.operator,
            environment=dict(self.environment),
            devices=list(self.devices),
            records=[r.as_dict() for r in self.records],
            files=list(self.files),
            qc_decisions=list(self.qc_decisions),
            processing_steps=list(self.processing_steps),
            metadata=dict(self.metadata),
        )
        payload["content_hash"] = hash_mapping(payload)
        return payload

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.as_dict(), indent=indent, default=str)

    def write(self, path: str, *, indent: int = 2) -> str:
        parent = os.path.dirname(os.path.abspath(path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(self.to_json(indent=indent))
        return os.path.abspath(path)


def _default_environment() -> Dict[str, Any]:
    return dict(
        python=platform.python_version(),
        platform=platform.platform(),
        machine=platform.machine(),
    )


def build_acquisition_manifest(
    survey_id: str,
    *,
    records: Optional[Iterable[ProvenanceRecord | Mapping[str, Any]]] = None,
    devices: Optional[Iterable[Mapping[str, Any]]] = None,
    qc_decisions: Optional[Iterable[Mapping[str, Any]]] = None,
    processing_steps: Optional[Iterable[str]] = None,
    method: Optional[str] = None,
    operator: Optional[str] = None,
    metadata: Optional[Mapping[str, Any]] = None,
) -> AcquisitionManifest:
    """Assemble an :class:`AcquisitionManifest` from session components."""
    manifest = AcquisitionManifest(
        survey_id=survey_id,
        method=method,
        operator=operator,
        records=[
            r if isinstance(r, ProvenanceRecord) else ProvenanceRecord(**dict(r))
            for r in list(records or [])
        ],
        devices=[dict(d) for d in list(devices or [])],
        qc_decisions=[dict(q) for q in list(qc_decisions or [])],
        processing_steps=list(processing_steps or []),
        metadata=dict(metadata or {}),
    )
    return manifest


def export_acquisition_manifest(
    survey_id: str,
    path: str,
    *,
    indent: int = 2,
    **kwargs: Any,
) -> str:
    """Build an acquisition manifest and write it to *path*.

    Thin convenience over :func:`build_acquisition_manifest` +
    :meth:`AcquisitionManifest.write`; ``kwargs`` are forwarded to the
    builder (``records``, ``devices``, ``qc_decisions``, ``method``, ...).
    Returns the written path.
    """
    manifest = build_acquisition_manifest(survey_id, **kwargs)
    return manifest.write(path, indent=indent)


def export_station_audit(
    record: ProvenanceRecord | Mapping[str, Any],
    path: str,
    *,
    indent: int = 2,
) -> str:
    """Write a single-station provenance audit to *path* (JSON)."""
    prov = (
        record if isinstance(record, ProvenanceRecord)
        else ProvenanceRecord(**dict(record))
    )
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    audit = prov.as_dict()
    audit["audit_hash"] = hash_mapping(audit)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=indent, default=str)
    return os.path.abspath(path)


def export_reproducibility_bundle(
    manifest: AcquisitionManifest,
    out_dir: str,
    *,
    include_raw: bool = False,
    zip_bundle: bool = False,
) -> Dict[str, Any]:
    """Write a manifest and per-station audits into *out_dir*.

    Parameters
    ----------
    manifest : AcquisitionManifest
        The session manifest to export.
    out_dir : str
        Destination directory (created if needed).
    include_raw : bool
        When ``True``, copy referenced raw files into ``out_dir/raw`` if
        they exist on disk.
    zip_bundle : bool
        When ``True``, also produce ``<out_dir>.zip`` of the bundle.

    Returns
    -------
    dict
        Paths written: ``manifest``, ``audits``, optional ``raw`` and
        ``zip``.
    """
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, "acquisition_manifest.json")
    manifest.write(manifest_path)

    audits_dir = os.path.join(out_dir, "audits")
    os.makedirs(audits_dir, exist_ok=True)
    audit_paths: List[str] = []
    for record in manifest.records:
        audit_path = os.path.join(audits_dir, f"{record.station_id}.json")
        export_station_audit(record, audit_path)
        audit_paths.append(audit_path)

    result: Dict[str, Any] = dict(
        manifest=os.path.abspath(manifest_path),
        audits=audit_paths,
    )

    if include_raw:
        import shutil

        raw_dir = os.path.join(out_dir, "raw")
        os.makedirs(raw_dir, exist_ok=True)
        copied: List[str] = []
        seen: List[Dict[str, Any]] = list(manifest.files)
        for record in manifest.records:
            seen.extend(record.raw_files)
        for entry in seen:
            src = entry.get("path")
            if src and os.path.isfile(src):
                dst = os.path.join(raw_dir, os.path.basename(src))
                shutil.copy2(src, dst)
                copied.append(dst)
        result["raw"] = copied

    if zip_bundle:
        zip_path = os.path.abspath(out_dir.rstrip(os.sep)) + ".zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(out_dir):
                for name in files:
                    full = os.path.join(root, name)
                    zf.write(full, os.path.relpath(full, out_dir))
        result["zip"] = zip_path

    return result
