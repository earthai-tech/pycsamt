# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Application-neutral draft, mapping, and validation helpers for PCBH."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from ...api.property import PyCSAMTObject
from .jsonio import pcbh_from_dict, pcbh_to_dict
from .mapping import (
    CANONICAL_CSV_FIELDS,
    ImportReport,
    resolve_csv_columns,
)
from .schema import PCBHDocument, PCBHValidationError, ValidationIssue

__all__ = [
    "CSVMappingProfile",
    "BuilderDiagnostic",
    "BuilderValidation",
    "new_builder_draft",
    "document_from_builder",
    "validate_builder",
    "preview_csv_mapping",
    "document_to_builder",
]

_TABLE_EXTENSIONS = ("water", "construction", "samples", "assays")


@dataclass(frozen=True, repr=False)
class CSVMappingProfile(PyCSAMTObject):
    """Reusable source-column mapping and constant values."""

    name: str
    columns: dict[str, str] = field(default_factory=dict)
    constants: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("mapping profile name must be non-empty")
        unknown = set(self.columns) - set(CANONICAL_CSV_FIELDS)
        if unknown:
            raise ValueError(
                f"unknown canonical mapping fields: {sorted(unknown)}"
            )

    def to_dict(self) -> dict[str, Any]:
        self.validate()
        return {
            "name": self.name,
            "columns": dict(self.columns),
            "constants": dict(self.constants),
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> CSVMappingProfile:
        profile = cls(
            str(value.get("name", "")),
            dict(value.get("columns", {})),
            dict(value.get("constants", {})),
        )
        profile.validate()
        return profile


@dataclass(frozen=True, repr=False)
class BuilderDiagnostic(PyCSAMTObject):
    """Validation issue augmented with editor/table navigation."""

    severity: str
    code: str
    message: str
    path: str
    editor: str
    row: int | None = None


@dataclass(frozen=True, repr=False)
class BuilderValidation(PyCSAMTObject):
    """Result of constructing and validating a builder draft."""

    document: PCBHDocument | None
    diagnostics: tuple[BuilderDiagnostic, ...]

    @property
    def ok(self) -> bool:
        return self.document is not None and not any(
            item.severity == "error" for item in self.diagnostics
        )


def new_builder_draft() -> dict[str, Any]:
    """Return a JSON-safe empty draft suitable for local browser storage."""
    return {
        "project": {
            "document_id": "pcbh:new-project",
            "created_by": "pyCSAMT Borehole Builder",
            "crs_horizontal": "EPSG:4326",
            "crs_vertical": "unknown",
            "coordinate_unit": "m",
            "depth_unit": "m",
            "diameter_unit": "m",
        },
        "boreholes": [],
        "surveys": [],
        "intervals": [],
        "structures": [],
        "water": [],
        "construction": [],
        "samples": [],
        "assays": [],
        "mapping_profiles": {},
    }


def document_from_builder(draft: dict[str, Any]) -> PCBHDocument:
    """Construct the canonical PCBH object from editable flat tables."""
    if not isinstance(draft, dict):
        raise TypeError("builder draft must be a mapping")
    project = dict(draft.get("project") or {})
    borehole_rows = list(draft.get("boreholes") or [])
    surveys = _group_rows(draft.get("surveys"), "borehole_id")
    intervals = _group_rows(draft.get("intervals"), "borehole_id")
    structures = _group_rows(draft.get("structures"), "borehole_id")
    extension_groups = {
        name: _group_rows(draft.get(name), "borehole_id")
        for name in _TABLE_EXTENSIONS
    }
    vocabulary: dict[str, dict[str, Any]] = {}
    boreholes = []
    for row in borehole_rows:
        hole_id = _required(row, "id")
        hole_intervals: dict[str, list[dict[str, Any]]] = {}
        for item in intervals.get(hole_id, []):
            family = str(item.get("family") or "lithology")
            interval = {
                "from_md": _number(item.get("from_md")),
                "to_md": _number(item.get("to_md")),
                "code": _optional_text(item.get("code")),
                "label": _optional_text(item.get("label")),
                "description": _optional_text(item.get("description")),
                "data_nature": str(item.get("data_nature") or "unknown"),
            }
            hole_intervals.setdefault(family, []).append(interval)
            if family == "lithology" and interval["code"]:
                vocabulary.setdefault(
                    interval["code"],
                    {
                        "code": interval["code"],
                        "name": interval["label"] or interval["code"],
                        "color": _optional_text(item.get("color")),
                    },
                )
        survey_rows = surveys.get(hole_id, [])
        trajectory = {
            "method": "survey" if survey_rows else "vertical",
            "north_reference": str(row.get("north_reference") or "unknown"),
            "desurvey_method": "minimum_curvature",
            "stations": [
                {
                    "md": _number(item.get("md")),
                    "azimuth_deg": _number(item.get("azimuth_deg")),
                    "inclination_deg": _number(item.get("inclination_deg")),
                }
                for item in survey_rows
            ],
        }
        extensions = {}
        for name, groups in extension_groups.items():
            rows = groups.get(hole_id, [])
            if rows:
                extensions[f"pcbh:{name}"] = [
                    {
                        key: value
                        for key, value in item.items()
                        if key != "borehole_id"
                    }
                    for item in rows
                ]
        boreholes.append(
            {
                "id": hole_id,
                "name": str(row.get("name") or hole_id),
                "kind": str(row.get("kind") or "unknown"),
                "status": str(row.get("status") or "unknown"),
                "collar": {
                    "x": _number(row.get("x")),
                    "y": _number(row.get("y")),
                    "z": _number(row.get("z")),
                },
                "total_depth_md": _number(row.get("total_depth_md")),
                "diameter": _optional_number(row.get("diameter")),
                "trajectory": trajectory,
                "interval_logs": hole_intervals,
                "structures": [
                    _structure(item) for item in structures.get(hole_id, [])
                ],
                "extensions": extensions,
            }
        )
    now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    payload = {
        "pcbh_version": "0.1.0",
        "document_id": str(project.get("document_id") or "pcbh:new-project"),
        "created_at": str(project.get("created_at") or now),
        "created_by": str(
            project.get("created_by") or "pyCSAMT Borehole Builder"
        ),
        "crs": {
            "horizontal": str(project.get("crs_horizontal") or ""),
            "vertical": str(project.get("crs_vertical") or "unknown"),
            "axis_order": "xy",
            "coordinate_unit": str(project.get("coordinate_unit") or "m"),
        },
        "units": {
            "depth": str(project.get("depth_unit") or "m"),
            "diameter": str(project.get("diameter_unit") or "m"),
            "angle": "deg",
            "resistivity": "ohm.m",
        },
        "conventions": {
            "depth_reference": "collar",
            "z_positive": "up",
            "azimuth_direction": "clockwise",
            "inclination_reference": "vertical_down",
        },
        "dictionaries": {
            "lithologies": list(vocabulary.values()),
            "formations": [],
        },
        "boreholes": boreholes,
    }
    return pcbh_from_dict(payload)


def validate_builder(draft: dict[str, Any]) -> BuilderValidation:
    """Validate a draft and return diagnostics mapped to editor rows."""
    try:
        document = document_from_builder(draft)
    except PCBHValidationError as error:
        return BuilderValidation(
            None,
            tuple(_diagnostic(item) for item in error.issues),
        )
    except (TypeError, ValueError, KeyError) as error:
        return BuilderValidation(
            None,
            (
                BuilderDiagnostic(
                    "error",
                    "builder.construction",
                    str(error),
                    "$",
                    "project",
                ),
            ),
        )
    return BuilderValidation(document, ())


def preview_csv_mapping(
    data: bytes,
    *,
    profile: CSVMappingProfile | None = None,
    max_rows: int = 20,
) -> dict[str, Any]:
    """Return bounded rows, inferred mapping, and diagnostics before import."""
    if not isinstance(data, bytes):
        raise TypeError("CSV preview data must be bytes")
    if len(data) > 10 * 1024 * 1024:
        raise ValueError("CSV preview exceeds 10485760 bytes")
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError("CSV preview must be UTF-8") from error
    dialect = csv.Sniffer().sniff(text[:65536], delimiters=",;\t|")
    reader = csv.DictReader(io.StringIO(text), dialect=dialect)
    headers = list(reader.fieldnames or [])
    if profile is not None:
        profile.validate()
    report = ImportReport(
        "browser-preview",
        "",
        dialect.delimiter,
        False,
    )
    mapping = resolve_csv_columns(
        headers,
        explicit=profile.columns if profile else None,
        constants=profile.constants if profile else {},
        report=report,
    )
    rows = []
    for index, row in enumerate(reader):
        if index >= max_rows:
            break
        rows.append(dict(row))
    return {
        "headers": headers,
        "rows": rows,
        "mapping": mapping,
        "delimiter": dialect.delimiter,
        "issues": [item.__dict__ for item in report.issues],
    }


def document_to_builder(document: PCBHDocument) -> dict[str, Any]:
    """Flatten a valid document for editing and browser-session recovery."""
    document.validate()
    canonical = pcbh_to_dict(document)
    draft = new_builder_draft()
    draft["project"] = {
        "document_id": document.document_id,
        "created_at": document.created_at,
        "created_by": document.created_by,
        "crs_horizontal": document.crs.horizontal,
        "crs_vertical": document.crs.vertical,
        "coordinate_unit": document.crs.coordinate_unit,
        "depth_unit": document.units.depth,
        "diameter_unit": document.units.diameter,
    }
    colors = {entry.code: entry.color for entry in document.lithologies}
    for hole, raw in zip(document.boreholes, canonical["boreholes"]):
        draft["boreholes"].append(
            {
                "id": hole.id,
                "name": hole.name,
                "kind": hole.kind,
                "status": hole.status,
                "x": hole.collar.x,
                "y": hole.collar.y,
                "z": hole.collar.z,
                "total_depth_md": hole.total_depth_md,
                "diameter": hole.diameter,
                "north_reference": hole.trajectory.north_reference,
            }
        )
        for station in raw["trajectory"]["stations"]:
            draft["surveys"].append({"borehole_id": hole.id, **station})
        for family, intervals in raw["interval_logs"].items():
            for interval in intervals:
                if family == "lithology" and interval.get("code"):
                    interval["color"] = colors.get(interval["code"])
                draft["intervals"].append(
                    {"borehole_id": hole.id, "family": family, **interval}
                )
        for item in raw["structures"]:
            draft["structures"].append({"borehole_id": hole.id, **item})
        for name in _TABLE_EXTENSIONS:
            for item in hole.extensions.get(f"pcbh:{name}", []):
                draft[name].append({"borehole_id": hole.id, **item})
    return draft


def _group_rows(value, key):
    result: dict[str, list[dict[str, Any]]] = {}
    for row in list(value or []):
        name = str(row.get(key) or "")
        result.setdefault(name, []).append(dict(row))
    return result


def _required(row, key):
    value = str(row.get(key) or "").strip()
    if not value:
        raise ValueError(f"borehole {key} is required")
    return value


def _number(value):
    if value is None or value == "":
        raise ValueError("required numeric value is missing")
    return float(value)


def _optional_number(value):
    return None if value is None or value == "" else float(value)


def _optional_text(value):
    return None if value is None or value == "" else str(value)


def _structure(item):
    return {
        "kind": str(item.get("kind") or "unknown"),
        "at_md": _optional_number(item.get("at_md")),
        "from_md": _optional_number(item.get("from_md")),
        "to_md": _optional_number(item.get("to_md")),
        "orientation_representation": str(
            item.get("orientation_representation") or "none"
        ),
        "dip_deg": _optional_number(item.get("dip_deg")),
        "dip_direction_deg": _optional_number(item.get("dip_direction_deg")),
    }


def _diagnostic(issue: ValidationIssue) -> BuilderDiagnostic:
    path = issue.path
    editor = "project"
    row = None
    if ".boreholes[" in path:
        editor = "boreholes"
        try:
            row = int(path.split(".boreholes[", 1)[1].split("]", 1)[0])
        except ValueError:
            row = None
        if ".trajectory.stations[" in path:
            editor = "surveys"
        elif ".interval_logs." in path:
            editor = "intervals"
        elif ".structures[" in path:
            editor = "structures"
    return BuilderDiagnostic(
        issue.severity,
        issue.code,
        issue.message,
        issue.path,
        editor,
        row,
    )
