# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Manifest-driven relational CSV import and export for PCBH."""

from __future__ import annotations

import csv
import hashlib
import math
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from ...api.property import PyCSAMTObject
from .mapping import ImportReport, PCBHCSVImportError
from .schema import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    PCBHValidationError,
    StructureObservation,
    SurveyStation,
    Trajectory,
    UnitSystem,
    VocabularyEntry,
)

__all__ = [
    "RELATIONAL_VERSION",
    "RelationalManifest",
    "boreholes_from_csv_directory",
    "write_csv_directory",
]

RELATIONAL_VERSION = "0.1.0"
_TABLES = (
    "collars",
    "surveys",
    "lithology",
    "structures",
    "samples",
    "assays",
)


@dataclass(frozen=True, repr=False)
class RelationalManifest(PyCSAMTObject):
    """Validated description of a PCBH relational CSV directory."""

    crs_horizontal: str
    tables: dict[str, str]
    depth_unit: str = "m"
    resistivity_unit: str = "ohm.m"
    version: str = RELATIONAL_VERSION

    @classmethod
    def read(cls, path: str | Path) -> RelationalManifest:
        """Read a manifest with ``yaml.safe_load`` and validate it."""
        source = Path(path)
        data = yaml.safe_load(source.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("relational manifest must be a mapping")
        units = data.get("units", {})
        result = cls(
            version=str(data.get("version", "")),
            crs_horizontal=str(data.get("crs", {}).get("horizontal", "")),
            depth_unit=str(units.get("depth", "m")),
            resistivity_unit=str(units.get("resistivity", "ohm.m")),
            tables=dict(data.get("tables", {})),
        )
        result.validate()
        return result

    def validate(self) -> None:
        """Validate version, units, table names, and safe relative paths."""
        if self.version != RELATIONAL_VERSION:
            raise ValueError(
                f"manifest version must be {RELATIONAL_VERSION!r}"
            )
        if not self.crs_horizontal.strip():
            raise ValueError("manifest CRS horizontal value is required")
        if self.depth_unit != "m" or self.resistivity_unit != "ohm.m":
            raise ValueError(
                "relational import currently requires m and ohm.m"
            )
        if "collars" not in self.tables:
            raise ValueError("manifest tables must include collars")
        for name, value in self.tables.items():
            if name not in _TABLES:
                raise ValueError(f"unknown relational table {name!r}")
            path = Path(value)
            if path.is_absolute() or ".." in path.parts or path.name != value:
                raise ValueError(
                    f"table path must be a safe filename: {value!r}"
                )


def boreholes_from_csv_directory(
    directory: str | Path,
    *,
    manifest: str = "import.yaml",
    strict: bool = True,
    document_id: str | None = None,
    created_by: str = "pycsamt relational CSV importer",
) -> tuple[PCBHDocument, ImportReport]:
    """Import joined collar, survey, log, structure, sample, assay tables."""
    root = Path(directory).resolve()
    manifest_path = (root / manifest).resolve()
    if manifest_path.parent != root:
        raise ValueError("manifest must be inside the project directory")
    config = RelationalManifest.read(manifest_path)
    report = ImportReport(
        source=str(root),
        source_sha256=_sha256(manifest_path),
        delimiter=",",
        strict=strict,
    )
    tables: dict[str, list[dict[str, str]]] = {}
    for name, filename in config.tables.items():
        path = (root / filename).resolve()
        if path.parent != root:
            raise ValueError("table path escapes the project directory")
        if not path.exists():
            report.add(
                "error", "csv.table_missing", f"missing table {filename!r}"
            )
            continue
        tables[name] = _read_table(path)
        report.source_files[filename] = _sha256(path)
        report.rows_read += len(tables[name])
    if report.errors:
        raise PCBHCSVImportError(report)

    holes: dict[str, PCBHBorehole] = {}
    for row_number, row in enumerate(tables["collars"], start=2):
        try:
            hole_id = _required(row, "borehole_id")
            if hole_id in holes:
                raise ValueError(f"duplicate collar for {hole_id!r}")
            holes[hole_id] = PCBHBorehole(
                id=hole_id,
                name=row.get("name", "").strip() or hole_id,
                kind=row.get("kind", "").strip() or "unknown",
                status=row.get("status", "").strip() or "unknown",
                collar=Collar(
                    _float(row, "x"), _float(row, "y"), _float(row, "z")
                ),
                total_depth_md=_float(row, "total_depth_md"),
                diameter=_optional_float(row, "diameter"),
            )
            report.rows_accepted += 1
        except (TypeError, ValueError) as error:
            _reject(report, "collars", row_number, error)

    surveys: dict[str, list[SurveyStation]] = defaultdict(list)
    for row_number, row in enumerate(tables.get("surveys", []), start=2):
        hole_id = row.get("borehole_id", "").strip()
        if hole_id not in holes:
            _reject(report, "surveys", row_number, KeyError(hole_id))
            continue
        try:
            surveys[hole_id].append(
                SurveyStation(
                    _float(row, "md"),
                    _float(row, "azimuth_deg"),
                    _float(row, "inclination_deg"),
                )
            )
            report.rows_accepted += 1
        except (TypeError, ValueError) as error:
            _reject(report, "surveys", row_number, error)
    for hole_id, stations in surveys.items():
        north = next(
            (
                row.get("north_reference", "").strip()
                for row in tables.get("surveys", [])
                if row.get("borehole_id", "").strip() == hole_id
                and row.get("north_reference", "").strip()
            ),
            "unknown",
        )
        holes[hole_id].trajectory = Trajectory(
            method="survey",
            north_reference=north,
            stations=sorted(stations, key=lambda item: item.md),
        )

    vocab: dict[str, VocabularyEntry] = {}
    for row_number, row in enumerate(tables.get("lithology", []), start=2):
        hole = _joined_hole(holes, row, "lithology", row_number, report)
        if hole is None:
            continue
        try:
            label = row.get("label", "").strip()
            code = row.get("code", "").strip() or _code(label)
            interval = LogInterval(
                _float(row, "from_md"),
                _float(row, "to_md"),
                code=code,
                label=label or code,
                description=row.get("description", "").strip(),
                resistivity_ohm_m=_optional_float(row, "resistivity_ohm_m"),
                data_nature=row.get("data_nature", "").strip() or "observed",
                confidence=_optional_float(row, "confidence"),
            )
            hole.interval_logs.setdefault("lithology", []).append(interval)
            existing = vocab.get(code)
            if existing and existing.name != interval.label:
                raise ValueError(f"conflicting lithology code {code!r}")
            vocab.setdefault(
                code, VocabularyEntry(code, interval.label or code)
            )
            report.rows_accepted += 1
        except (TypeError, ValueError) as error:
            _reject(report, "lithology", row_number, error)

    for row_number, row in enumerate(tables.get("structures", []), start=2):
        hole = _joined_hole(holes, row, "structures", row_number, report)
        if hole is None:
            continue
        try:
            hole.structures.append(_structure(row))
            report.rows_accepted += 1
        except (TypeError, ValueError) as error:
            _reject(report, "structures", row_number, error)

    samples = _joined_records(
        tables.get("samples", []), holes, "sample_id", report
    )
    sample_ids = {item["sample_id"] for item in samples}
    assays = []
    assay_keys: set[tuple[str, str]] = set()
    for index, row in enumerate(tables.get("assays", []), start=2):
        sample_id = row.get("sample_id", "").strip()
        analyte = row.get("analyte", "").strip()
        key = (sample_id, analyte.casefold())
        if sample_id not in sample_ids:
            _reject(report, "assays", index, KeyError(row.get("sample_id")))
        elif not analyte or key in assay_keys:
            _reject(
                report,
                "assays",
                index,
                ValueError("missing or duplicate sample/analyte result"),
            )
        else:
            assays.append(_clean_record(row))
            assay_keys.add(key)
            report.rows_accepted += 1
    if strict and report.errors:
        raise PCBHCSVImportError(report)
    if not holes:
        report.add("error", "csv.no_valid_rows", "no valid collars remain")
        raise PCBHCSVImportError(report)
    for hole in holes.values():
        for intervals in hole.interval_logs.values():
            intervals.sort(key=lambda item: item.from_md)
    extensions = {}
    if samples:
        extensions["pcbh:samples"] = samples
    if assays:
        extensions["pcbh:assays"] = assays
    document = PCBHDocument(
        document_id=document_id or f"csv-directory:{root.name}",
        created_at=datetime.now(timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        created_by=created_by,
        crs=CoordinateReferenceSystem(config.crs_horizontal),
        units=UnitSystem(),
        boreholes=list(holes.values()),
        lithologies=list(vocab.values()),
        extensions=extensions,
    )
    try:
        document.validate()
    except PCBHValidationError as error:
        report.add("error", "csv.document_invalid", str(error))
        raise PCBHCSVImportError(report) from error
    return document, report


def write_csv_directory(document: PCBHDocument, directory: str | Path) -> Path:
    """Export supported PCBH content as a relational CSV directory."""
    document.validate()
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    tables = {"collars": "collars.csv"}
    survey_rows = _survey_rows(document)
    lithology_rows = _lithology_rows(document)
    structure_rows = _structure_rows(document)
    if survey_rows:
        tables["surveys"] = "surveys.csv"
    if lithology_rows:
        tables["lithology"] = "lithology.csv"
    if structure_rows:
        tables["structures"] = "structures.csv"
    if document.extensions.get("pcbh:samples"):
        tables["samples"] = "samples.csv"
    if document.extensions.get("pcbh:assays"):
        tables["assays"] = "assays.csv"
    _write_rows(
        root / tables["collars"],
        [_collar_row(hole) for hole in document.boreholes],
    )
    if survey_rows:
        _write_rows(root / tables["surveys"], survey_rows)
    if lithology_rows:
        _write_rows(root / tables["lithology"], lithology_rows)
    if structure_rows:
        _write_rows(root / tables["structures"], structure_rows)
    for extension, table in (
        ("pcbh:samples", "samples"),
        ("pcbh:assays", "assays"),
    ):
        if table in tables:
            _write_rows(root / tables[table], document.extensions[extension])
    manifest = {
        "version": RELATIONAL_VERSION,
        "crs": {"horizontal": document.crs.horizontal},
        "units": {
            "depth": document.units.depth,
            "resistivity": document.units.resistivity,
        },
        "tables": tables,
    }
    (root / "import.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )
    return root


def _read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f"CSV table has no header: {path}")
        return [
            dict(row)
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]


def _required(row: dict[str, str], key: str) -> str:
    value = row.get(key, "").strip()
    if not value:
        raise ValueError(f"{key} is required")
    return value


def _float(row: dict[str, str], key: str) -> float:
    value = float(_required(row, key))
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def _optional_float(row: dict[str, str], key: str) -> float | None:
    return None if not row.get(key, "").strip() else _float(row, key)


def _reject(
    report: ImportReport, table: str, row: int, error: Exception
) -> None:
    report.rows_rejected += 1
    report.add("error", "csv.join_or_value", f"{table}: {error}", row=row)


def _joined_hole(holes, row, table, number, report):
    hole_id = row.get("borehole_id", "").strip()
    if hole_id not in holes:
        _reject(report, table, number, KeyError(hole_id))
        return None
    return holes[hole_id]


def _structure(row: dict[str, str]) -> StructureObservation:
    return StructureObservation(
        kind=_required(row, "kind"),
        at_md=_optional_float(row, "at_md"),
        from_md=_optional_float(row, "from_md"),
        to_md=_optional_float(row, "to_md"),
        orientation_representation=row.get(
            "orientation_representation", ""
        ).strip()
        or "none",
        strike_deg=_optional_float(row, "strike_deg"),
        dip_deg=_optional_float(row, "dip_deg"),
        dip_direction_deg=_optional_float(row, "dip_direction_deg"),
        trend_deg=_optional_float(row, "trend_deg"),
        plunge_deg=_optional_float(row, "plunge_deg"),
        alpha_deg=_optional_float(row, "alpha_deg"),
        beta_deg=_optional_float(row, "beta_deg"),
        aperture_m=_optional_float(row, "aperture_m"),
        fill=row.get("fill", "").strip() or None,
        data_nature=row.get("data_nature", "").strip() or "observed",
        confidence=_optional_float(row, "confidence"),
    )


def _joined_records(rows, holes, identifier, report):
    output, seen = [], set()
    for number, row in enumerate(rows, start=2):
        hole = _joined_hole(holes, row, "samples", number, report)
        record_id = row.get(identifier, "").strip()
        if hole is None:
            continue
        if not record_id or record_id in seen:
            _reject(
                report,
                "samples",
                number,
                ValueError("missing or duplicate sample_id"),
            )
            continue
        try:
            start, end = _float(row, "from_md"), _float(row, "to_md")
            if not 0 <= start < end <= hole.total_depth_md:
                raise ValueError("sample bounds exceed borehole")
            record = _clean_record(row)
            record["from_md"], record["to_md"] = start, end
            output.append(record)
            seen.add(record_id)
            report.rows_accepted += 1
        except ValueError as error:
            _reject(report, "samples", number, error)
    return output


def _clean_record(row):
    return {
        key: value.strip()
        for key, value in row.items()
        if value and value.strip()
    }


def _code(label: str) -> str:
    import re

    return re.sub(r"[^A-Z0-9]+", "_", label.upper()).strip("_") or "UNKNOWN"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    headers = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def _collar_row(hole):
    return {
        "borehole_id": hole.id,
        "name": hole.name,
        "x": hole.collar.x,
        "y": hole.collar.y,
        "z": hole.collar.z,
        "total_depth_md": hole.total_depth_md,
        "kind": hole.kind,
        "status": hole.status,
        "diameter": hole.diameter if hole.diameter is not None else "",
    }


def _survey_rows(document):
    return [
        {
            "borehole_id": hole.id,
            "md": station.md,
            "azimuth_deg": station.azimuth_deg,
            "inclination_deg": station.inclination_deg,
            "north_reference": hole.trajectory.north_reference,
        }
        for hole in document.boreholes
        for station in hole.trajectory.stations
    ]


def _lithology_rows(document):
    return [
        {
            "borehole_id": hole.id,
            "from_md": item.from_md,
            "to_md": item.to_md,
            "code": item.code or "",
            "label": item.label or "",
            "description": item.description,
            "resistivity_ohm_m": (
                item.resistivity_ohm_m
                if item.resistivity_ohm_m is not None
                else ""
            ),
            "data_nature": item.data_nature,
            "confidence": item.confidence
            if item.confidence is not None
            else "",
        }
        for hole in document.boreholes
        for item in hole.interval_logs.get("lithology", [])
    ]


def _structure_rows(document):
    fields = (
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
    )
    return [
        {
            "borehole_id": hole.id,
            **{
                name: getattr(item, name)
                if getattr(item, name) is not None
                else ""
                for name in fields
            },
        }
        for hole in document.boreholes
        for item in hole.structures
    ]
