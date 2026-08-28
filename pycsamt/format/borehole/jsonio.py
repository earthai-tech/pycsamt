# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Canonical UTF-8 JSON reader and writer for PCBH 0.1."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from collections.abc import Mapping, Sequence
from os import PathLike
from pathlib import Path
from typing import Any

from ._version import check_pcbh_version
from .schema import (
    PCBH_VERSION,
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    StructureObservation,
    SurveyStation,
    Trajectory,
    UnitSystem,
    VocabularyEntry,
)

__all__ = [
    "PCBH_SCHEMA_URI",
    "DEFAULT_MAX_BYTES",
    "DEFAULT_MAX_BOREHOLES",
    "DEFAULT_MAX_INTERVALS",
    "DEFAULT_MAX_NESTING",
    "pcbh_to_dict",
    "pcbh_from_dict",
    "read_pcbh",
    "write_pcbh",
]

PCBH_SCHEMA_URI = "https://pycsamt.org/schemas/pcbh/0.1/schema.json"
DEFAULT_MAX_BYTES = 16 * 1024 * 1024
DEFAULT_MAX_BOREHOLES = 10_000
DEFAULT_MAX_INTERVALS = 1_000_000
DEFAULT_MAX_NESTING = 32

logger = logging.getLogger(__name__)


def _positive_limit(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{path} must be a JSON object")
    return value


def _sequence(value: Any, path: str) -> Sequence[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{path} must be a JSON array")
    return value


def _reject_unknown(
    value: Mapping[str, Any], allowed: set[str], path: str
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(
            f"{path} contains unsupported field(s): {', '.join(unknown)}"
        )


def _required(value: Mapping[str, Any], name: str, path: str) -> Any:
    if name not in value:
        raise ValueError(f"{path}.{name} is required")
    return value[name]


def _without_none(value: dict[str, Any]) -> dict[str, Any]:
    return {key: item for key, item in value.items() if item is not None}


def _vocabulary_to_dict(entry: VocabularyEntry) -> dict[str, Any]:
    return _without_none(
        {
            "code": entry.code,
            "name": entry.name,
            "color": entry.color,
            "description": entry.description,
            "external_ids": dict(entry.external_ids),
            "properties": dict(entry.properties),
        }
    )


def _station_to_dict(station: SurveyStation) -> dict[str, Any]:
    return {
        "md": station.md,
        "azimuth_deg": station.azimuth_deg,
        "inclination_deg": station.inclination_deg,
    }


def _interval_to_dict(interval: LogInterval) -> dict[str, Any]:
    return _without_none(
        {
            "from_md": interval.from_md,
            "to_md": interval.to_md,
            "code": interval.code,
            "label": interval.label,
            "description": interval.description,
            "resistivity_ohm_m": interval.resistivity_ohm_m,
            "data_nature": interval.data_nature,
            "confidence": interval.confidence,
            "properties": dict(interval.properties),
        }
    )


def _structure_to_dict(
    structure: StructureObservation,
) -> dict[str, Any]:
    return _without_none(
        {
            "kind": structure.kind,
            "at_md": structure.at_md,
            "from_md": structure.from_md,
            "to_md": structure.to_md,
            "orientation_representation": (
                structure.orientation_representation
            ),
            "strike_deg": structure.strike_deg,
            "dip_deg": structure.dip_deg,
            "dip_direction_deg": structure.dip_direction_deg,
            "trend_deg": structure.trend_deg,
            "plunge_deg": structure.plunge_deg,
            "alpha_deg": structure.alpha_deg,
            "beta_deg": structure.beta_deg,
            "aperture_m": structure.aperture_m,
            "fill": structure.fill,
            "data_nature": structure.data_nature,
            "confidence": structure.confidence,
        }
    )


def _borehole_to_dict(borehole: PCBHBorehole) -> dict[str, Any]:
    collar = _without_none(
        {
            "x": borehole.collar.x,
            "y": borehole.collar.y,
            "z": borehole.collar.z,
            "longitude": borehole.collar.longitude,
            "latitude": borehole.collar.latitude,
            "position_uncertainty_m": (borehole.collar.position_uncertainty_m),
        }
    )
    trajectory = {
        "method": borehole.trajectory.method,
        "north_reference": borehole.trajectory.north_reference,
        "desurvey_method": borehole.trajectory.desurvey_method,
        "stations": [
            _station_to_dict(station)
            for station in borehole.trajectory.stations
        ],
    }
    result = _without_none(
        {
            "id": borehole.id,
            "name": borehole.name,
            "kind": borehole.kind,
            "status": borehole.status,
            "aliases": list(borehole.aliases),
            "collar": collar,
            "total_depth_md": borehole.total_depth_md,
            "diameter": borehole.diameter,
            "trajectory": trajectory,
            "interval_logs": {
                family: [_interval_to_dict(item) for item in intervals]
                for family, intervals in borehole.interval_logs.items()
            },
            "structures": [
                _structure_to_dict(item) for item in borehole.structures
            ],
            "metadata": dict(borehole.metadata),
            "extensions": dict(borehole.extensions),
        }
    )
    return result


def pcbh_to_dict(
    document: PCBHDocument, *, validate: bool = True
) -> dict[str, Any]:
    """Convert a PCBH document to its canonical JSON-compatible mapping.

    Parameters
    ----------
    document : PCBHDocument
        In-memory document to convert.
    validate : bool, default True
        Run semantic validation before conversion.

    Returns
    -------
    dict
        Canonically ordered JSON-compatible mapping.

    Raises
    ------
    TypeError
        If *document* is not a :class:`PCBHDocument`.
    PCBHValidationError
        If semantic validation fails.
    ValueError
        If metadata or extension values are not finite JSON values.
    """
    if not isinstance(document, PCBHDocument):
        raise TypeError("document must be a PCBHDocument")
    if validate:
        document.validate()
    result = {
        "$schema": PCBH_SCHEMA_URI,
        "pcbh_version": document.pcbh_version,
        "document_id": document.document_id,
        "title": document.title,
        "description": document.description,
        "created_at": document.created_at,
        "created_by": document.created_by,
        "crs": {
            "horizontal": document.crs.horizontal,
            "vertical": document.crs.vertical,
            "axis_order": document.crs.axis_order,
            "coordinate_unit": document.crs.coordinate_unit,
        },
        "units": {
            "depth": document.units.depth,
            "diameter": document.units.diameter,
            "angle": document.units.angle,
            "resistivity": document.units.resistivity,
        },
        "conventions": {
            "depth_reference": document.depth_reference,
            "z_positive": document.z_positive,
            "azimuth_direction": document.azimuth_direction,
            "inclination_reference": document.inclination_reference,
        },
        "dictionaries": {
            "lithologies": [
                _vocabulary_to_dict(entry) for entry in document.lithologies
            ],
            "formations": [
                _vocabulary_to_dict(entry) for entry in document.formations
            ],
        },
        "boreholes": [
            _borehole_to_dict(borehole) for borehole in document.boreholes
        ],
        "metadata": dict(document.metadata),
        "extensions": dict(document.extensions),
    }
    _check_json_value(result, "$", max_nesting=DEFAULT_MAX_NESTING)
    return result


def _vocabulary_from_dict(value: Any, path: str) -> VocabularyEntry:
    item = _mapping(value, path)
    _reject_unknown(
        item,
        {
            "code",
            "name",
            "color",
            "description",
            "external_ids",
            "properties",
        },
        path,
    )
    return VocabularyEntry(
        code=_required(item, "code", path),
        name=_required(item, "name", path),
        color=item.get("color"),
        description=item.get("description", ""),
        external_ids=dict(_mapping(item.get("external_ids", {}), path)),
        properties=dict(_mapping(item.get("properties", {}), path)),
    )


def _collar_from_dict(value: Any, path: str) -> Collar:
    item = _mapping(value, path)
    _reject_unknown(
        item,
        {
            "x",
            "y",
            "z",
            "longitude",
            "latitude",
            "position_uncertainty_m",
        },
        path,
    )
    return Collar(
        x=_required(item, "x", path),
        y=_required(item, "y", path),
        z=_required(item, "z", path),
        longitude=item.get("longitude"),
        latitude=item.get("latitude"),
        position_uncertainty_m=item.get("position_uncertainty_m"),
    )


def _trajectory_from_dict(value: Any, path: str) -> Trajectory:
    item = _mapping(value, path)
    _reject_unknown(
        item,
        {"method", "north_reference", "desurvey_method", "stations"},
        path,
    )
    stations = []
    station_values = _sequence(item.get("stations", []), f"{path}.stations")
    for index, station in enumerate(station_values):
        station_path = f"{path}.stations[{index}]"
        station_item = _mapping(station, station_path)
        _reject_unknown(
            station_item,
            {"md", "azimuth_deg", "inclination_deg"},
            station_path,
        )
        stations.append(
            SurveyStation(
                md=_required(station_item, "md", station_path),
                azimuth_deg=_required(
                    station_item, "azimuth_deg", station_path
                ),
                inclination_deg=_required(
                    station_item, "inclination_deg", station_path
                ),
            )
        )
    return Trajectory(
        method=_required(item, "method", path),
        north_reference=_required(item, "north_reference", path),
        desurvey_method=_required(item, "desurvey_method", path),
        stations=stations,
    )


def _interval_from_dict(value: Any, path: str) -> LogInterval:
    item = _mapping(value, path)
    _reject_unknown(
        item,
        {
            "from_md",
            "to_md",
            "code",
            "label",
            "description",
            "resistivity_ohm_m",
            "data_nature",
            "confidence",
            "properties",
        },
        path,
    )
    return LogInterval(
        from_md=_required(item, "from_md", path),
        to_md=_required(item, "to_md", path),
        code=item.get("code"),
        label=item.get("label"),
        description=item.get("description", ""),
        resistivity_ohm_m=item.get("resistivity_ohm_m"),
        data_nature=item.get("data_nature", "unknown"),
        confidence=item.get("confidence"),
        properties=dict(_mapping(item.get("properties", {}), path)),
    )


def _structure_from_dict(value: Any, path: str) -> StructureObservation:
    item = _mapping(value, path)
    fields = {
        "kind",
        "at_md",
        "from_md",
        "to_md",
        "orientation_representation",
        "strike_deg",
        "dip_deg",
        "dip_direction_deg",
        "trend_deg",
        "plunge_deg",
        "alpha_deg",
        "beta_deg",
        "aperture_m",
        "fill",
        "data_nature",
        "confidence",
    }
    _reject_unknown(item, fields, path)
    return StructureObservation(
        kind=_required(item, "kind", path),
        at_md=item.get("at_md"),
        from_md=item.get("from_md"),
        to_md=item.get("to_md"),
        orientation_representation=item.get(
            "orientation_representation", "none"
        ),
        strike_deg=item.get("strike_deg"),
        dip_deg=item.get("dip_deg"),
        dip_direction_deg=item.get("dip_direction_deg"),
        trend_deg=item.get("trend_deg"),
        plunge_deg=item.get("plunge_deg"),
        alpha_deg=item.get("alpha_deg"),
        beta_deg=item.get("beta_deg"),
        aperture_m=item.get("aperture_m"),
        fill=item.get("fill"),
        data_nature=item.get("data_nature", "unknown"),
        confidence=item.get("confidence"),
    )


def _borehole_from_dict(value: Any, path: str) -> PCBHBorehole:
    item = _mapping(value, path)
    _reject_unknown(
        item,
        {
            "id",
            "name",
            "kind",
            "status",
            "aliases",
            "collar",
            "total_depth_md",
            "diameter",
            "trajectory",
            "interval_logs",
            "structures",
            "metadata",
            "extensions",
        },
        path,
    )
    logs_raw = _mapping(item.get("interval_logs", {}), f"{path}.interval_logs")
    interval_logs = {
        str(family): [
            _interval_from_dict(
                interval, f"{path}.interval_logs.{family}[{index}]"
            )
            for index, interval in enumerate(
                _sequence(intervals, f"{path}.interval_logs.{family}")
            )
        ]
        for family, intervals in logs_raw.items()
    }
    return PCBHBorehole(
        id=_required(item, "id", path),
        name=_required(item, "name", path),
        kind=_required(item, "kind", path),
        status=_required(item, "status", path),
        aliases=list(_sequence(item.get("aliases", []), f"{path}.aliases")),
        collar=_collar_from_dict(
            _required(item, "collar", path), f"{path}.collar"
        ),
        total_depth_md=_required(item, "total_depth_md", path),
        diameter=item.get("diameter"),
        trajectory=_trajectory_from_dict(
            _required(item, "trajectory", path), f"{path}.trajectory"
        ),
        interval_logs=interval_logs,
        structures=[
            _structure_from_dict(structure, f"{path}.structures[{index}]")
            for index, structure in enumerate(
                _sequence(item.get("structures", []), f"{path}.structures")
            )
        ],
        metadata=dict(_mapping(item.get("metadata", {}), f"{path}.metadata")),
        extensions=dict(
            _mapping(item.get("extensions", {}), f"{path}.extensions")
        ),
    )


def pcbh_from_dict(
    value: Mapping[str, Any],
    *,
    validate: bool = True,
    max_boreholes: int = DEFAULT_MAX_BOREHOLES,
    max_intervals: int = DEFAULT_MAX_INTERVALS,
    max_nesting: int = DEFAULT_MAX_NESTING,
) -> PCBHDocument:
    """Build a PCBH document from a decoded JSON mapping.

    Parameters
    ----------
    value : mapping
        Decoded PCBH root object.
    validate : bool, default True
        Run semantic validation before returning.
    max_boreholes, max_intervals, max_nesting : int
        Positive resource limits applied before object construction.

    Returns
    -------
    PCBHDocument
        Parsed in-memory document.

    Raises
    ------
    ValueError
        If structure, version, or resource limits are invalid.
    PCBHValidationError
        If semantic validation fails.
    """
    max_boreholes = _positive_limit("max_boreholes", max_boreholes)
    max_intervals = _positive_limit("max_intervals", max_intervals)
    max_nesting = _positive_limit("max_nesting", max_nesting)
    root = _mapping(value, "$")
    _check_json_value(root, "$", max_nesting=max_nesting)
    _reject_unknown(
        root,
        {
            "$schema",
            "pcbh_version",
            "document_id",
            "title",
            "description",
            "created_at",
            "created_by",
            "crs",
            "units",
            "conventions",
            "dictionaries",
            "boreholes",
            "metadata",
            "extensions",
        },
        "$",
    )
    version = _required(root, "pcbh_version", "$")
    check_pcbh_version(version, PCBH_VERSION)
    boreholes_raw = _sequence(_required(root, "boreholes", "$"), "$.boreholes")
    if len(boreholes_raw) > max_boreholes:
        raise ValueError(
            f"PCBH contains {len(boreholes_raw)} boreholes; limit is "
            f"{max_boreholes}"
        )
    interval_count = _count_intervals(boreholes_raw)
    if interval_count > max_intervals:
        raise ValueError(
            f"PCBH contains {interval_count} intervals; limit is "
            f"{max_intervals}"
        )
    crs_raw = _mapping(_required(root, "crs", "$"), "$.crs")
    _reject_unknown(
        crs_raw,
        {"horizontal", "vertical", "axis_order", "coordinate_unit"},
        "$.crs",
    )
    units_raw = _mapping(_required(root, "units", "$"), "$.units")
    _reject_unknown(
        units_raw,
        {"depth", "diameter", "angle", "resistivity"},
        "$.units",
    )
    conventions = _mapping(
        _required(root, "conventions", "$"), "$.conventions"
    )
    _reject_unknown(
        conventions,
        {
            "depth_reference",
            "z_positive",
            "azimuth_direction",
            "inclination_reference",
        },
        "$.conventions",
    )
    dictionaries = _mapping(root.get("dictionaries", {}), "$.dictionaries")
    _reject_unknown(
        dictionaries, {"lithologies", "formations"}, "$.dictionaries"
    )
    document = PCBHDocument(
        pcbh_version=version,
        document_id=_required(root, "document_id", "$"),
        title=root.get("title", ""),
        description=root.get("description", ""),
        created_at=_required(root, "created_at", "$"),
        created_by=_required(root, "created_by", "$"),
        crs=CoordinateReferenceSystem(
            horizontal=_required(crs_raw, "horizontal", "$.crs"),
            vertical=_required(crs_raw, "vertical", "$.crs"),
            axis_order=_required(crs_raw, "axis_order", "$.crs"),
            coordinate_unit=_required(crs_raw, "coordinate_unit", "$.crs"),
        ),
        units=UnitSystem(
            depth=_required(units_raw, "depth", "$.units"),
            diameter=_required(units_raw, "diameter", "$.units"),
            angle=_required(units_raw, "angle", "$.units"),
            resistivity=_required(units_raw, "resistivity", "$.units"),
        ),
        depth_reference=_required(
            conventions, "depth_reference", "$.conventions"
        ),
        z_positive=_required(conventions, "z_positive", "$.conventions"),
        azimuth_direction=_required(
            conventions, "azimuth_direction", "$.conventions"
        ),
        inclination_reference=_required(
            conventions, "inclination_reference", "$.conventions"
        ),
        lithologies=[
            _vocabulary_from_dict(item, f"$.dictionaries.lithologies[{index}]")
            for index, item in enumerate(
                _sequence(
                    dictionaries.get("lithologies", []),
                    "$.dictionaries.lithologies",
                )
            )
        ],
        formations=[
            _vocabulary_from_dict(item, f"$.dictionaries.formations[{index}]")
            for index, item in enumerate(
                _sequence(
                    dictionaries.get("formations", []),
                    "$.dictionaries.formations",
                )
            )
        ],
        boreholes=[
            _borehole_from_dict(item, f"$.boreholes[{index}]")
            for index, item in enumerate(boreholes_raw)
        ],
        metadata=dict(_mapping(root.get("metadata", {}), "$.metadata")),
        extensions=dict(_mapping(root.get("extensions", {}), "$.extensions")),
    )
    if validate:
        document.validate()
    return document


def _count_intervals(boreholes: Sequence[Any]) -> int:
    count = 0
    for index, value in enumerate(boreholes):
        hole = _mapping(value, f"$.boreholes[{index}]")
        logs = _mapping(
            hole.get("interval_logs", {}),
            f"$.boreholes[{index}].interval_logs",
        )
        for family, intervals in logs.items():
            count += len(
                _sequence(
                    intervals,
                    f"$.boreholes[{index}].interval_logs.{family}",
                )
            )
    return count


def _check_json_value(value: Any, path: str, *, max_nesting: int) -> None:
    if max_nesting < 0:
        raise ValueError(f"JSON nesting exceeds configured limit at {path}")
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError(f"non-finite JSON number at {path}")
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(f"JSON object key at {path} must be a string")
            _check_json_value(
                item,
                f"{path}.{key}",
                max_nesting=max_nesting - 1,
            )
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _check_json_value(
                item,
                f"{path}[{index}]",
                max_nesting=max_nesting - 1,
            )
        return
    raise ValueError(
        f"value at {path} is not JSON serializable: {type(value).__name__}"
    )


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key {key!r}")
        result[key] = value
    return result


def read_pcbh(
    path: str | PathLike[str],
    *,
    validate: bool = True,
    max_bytes: int = DEFAULT_MAX_BYTES,
    max_boreholes: int = DEFAULT_MAX_BOREHOLES,
    max_intervals: int = DEFAULT_MAX_INTERVALS,
    max_nesting: int = DEFAULT_MAX_NESTING,
) -> PCBHDocument:
    """Read a canonical PCBH JSON file.

    Parameters
    ----------
    path : path-like
        Input `.pcbh.json` file.
    validate : bool, default True
        Run semantic validation before returning.
    max_bytes, max_boreholes, max_intervals, max_nesting : int
        Positive limits for untrusted input.

    Returns
    -------
    PCBHDocument
        Parsed and optionally validated document.

    Raises
    ------
    OSError
        If the file cannot be read.
    UnicodeError
        If the file is not UTF-8.
    ValueError
        If JSON, structure, version, or limits are invalid.
    PCBHValidationError
        If semantic validation fails.
    """
    max_bytes = _positive_limit("max_bytes", max_bytes)
    source = Path(path)
    size = source.stat().st_size
    if size > max_bytes:
        raise ValueError(
            f"PCBH file is {size} bytes; max_bytes is {max_bytes}"
        )
    raw = source.read_bytes()
    if len(raw) > max_bytes:
        raise ValueError(
            f"PCBH file is {len(raw)} bytes; max_bytes is {max_bytes}"
        )
    try:
        value = json.loads(
            raw.decode("utf-8"), object_pairs_hook=_unique_object
        )
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"{source} is not valid PCBH JSON: {exc.msg} at line "
            f"{exc.lineno}, column {exc.colno}"
        ) from exc
    document = pcbh_from_dict(
        value,
        validate=validate,
        max_boreholes=max_boreholes,
        max_intervals=max_intervals,
        max_nesting=max_nesting,
    )
    logger.debug("Read PCBH document %s from %s", document.document_id, source)
    return document


def write_pcbh(
    document: PCBHDocument,
    path: str | PathLike[str],
    *,
    validate: bool = True,
    indent: int = 2,
) -> Path:
    """Atomically write a canonical PCBH JSON file.

    Parameters
    ----------
    document : PCBHDocument
        Document to serialize.
    path : path-like
        Destination path. Parent directories are created when needed.
    validate : bool, default True
        Run semantic validation before writing.
    indent : int, default 2
        Positive JSON indentation width.

    Returns
    -------
    pathlib.Path
        Destination path.

    Raises
    ------
    TypeError
        If *document* is not a PCBH document.
    ValueError
        If validation, JSON values, or *indent* are invalid.
    OSError
        If the atomic write fails.

    Notes
    -----
    The temporary file is created beside the destination so
    :func:`os.replace` remains an atomic same-filesystem operation.
    """
    if isinstance(indent, bool) or not isinstance(indent, int) or indent <= 0:
        raise ValueError("indent must be a positive integer")
    value = pcbh_to_dict(document, validate=validate)
    text = (
        json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            indent=indent,
            sort_keys=False,
        )
        + "\n"
    )
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
    logger.debug(
        "Wrote PCBH document %s to %s", document.document_id, destination
    )
    return destination
