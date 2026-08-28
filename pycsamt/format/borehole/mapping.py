# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Column mapping and reporting primitives for PCBH tabular imports."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from ...api.property import PyCSAMTObject

__all__ = [
    "CANONICAL_CSV_FIELDS",
    "REQUIRED_CSV_FIELDS",
    "ImportIssue",
    "ImportReport",
    "PCBHCSVImportError",
    "resolve_csv_columns",
]

CANONICAL_CSV_FIELDS = (
    "borehole.id",
    "borehole.name",
    "borehole.kind",
    "borehole.status",
    "borehole.total_depth_md",
    "collar.x",
    "collar.y",
    "collar.z",
    "crs.horizontal",
    "interval.from_md",
    "interval.to_md",
    "interval.code",
    "interval.lithology",
    "interval.description",
    "interval.resistivity_ohm_m",
    "interval.data_nature",
)

REQUIRED_CSV_FIELDS = (
    "borehole.id",
    "collar.x",
    "collar.y",
    "collar.z",
    "interval.from_md",
    "interval.to_md",
    "interval.lithology",
)

_ALIASES = {
    "borehole.id": (
        "borehole_id",
        "holeid",
        "hole_id",
        "well",
        "well_id",
        "bhid",
    ),
    "borehole.name": ("borehole_name", "hole_name", "well_name"),
    "borehole.kind": ("kind", "borehole_type", "hole_type"),
    "borehole.status": ("status", "hole_status"),
    "borehole.total_depth_md": ("total_depth_md", "total_depth", "td"),
    "collar.x": ("x", "easting", "east", "longitude", "lon"),
    "collar.y": ("y", "northing", "north", "latitude", "lat"),
    "collar.z": ("z", "rl", "elevation", "elev", "collar_elevation"),
    "crs.horizontal": ("crs", "horizontal_crs", "epsg"),
    "interval.from_md": ("from_md", "from", "top", "depth_from"),
    "interval.to_md": ("to_md", "to", "bottom", "depth_to"),
    "interval.code": ("lithology_code", "lith_code", "rock_code", "code"),
    "interval.lithology": ("lithology", "lith", "rock", "formation"),
    "interval.description": ("description", "lithology_description"),
    "interval.resistivity_ohm_m": (
        "resistivity_ohm_m",
        "resistivity",
        "tres",
        "rho",
    ),
    "interval.data_nature": ("data_nature", "nature"),
}


@dataclass(frozen=True, repr=False)
class ImportIssue(PyCSAMTObject):
    """One localized CSV import diagnostic."""

    severity: str
    code: str
    message: str
    row: int | None = None
    column: str | None = None

    def __post_init__(self) -> None:
        if self.severity not in {"error", "warning", "info"}:
            raise ValueError("severity must be error, warning, or info")


@dataclass(repr=False)
class ImportReport(PyCSAMTObject):
    """Structured provenance and diagnostics for one CSV import."""

    source: str
    source_sha256: str
    delimiter: str
    strict: bool
    source_files: dict[str, str] = field(default_factory=dict)
    rows_read: int = 0
    rows_accepted: int = 0
    rows_skipped: int = 0
    rows_rejected: int = 0
    column_mapping: dict[str, str] = field(default_factory=dict)
    constants: dict[str, Any] = field(default_factory=dict)
    inferred_values: list[str] = field(default_factory=list)
    unit_conversions: list[str] = field(default_factory=list)
    conflict_resolutions: list[str] = field(default_factory=list)
    unresolved_labels: list[str] = field(default_factory=list)
    issues: list[ImportIssue] = field(default_factory=list)

    @property
    def errors(self) -> tuple[ImportIssue, ...]:
        """Return error diagnostics in source order."""
        return tuple(item for item in self.issues if item.severity == "error")

    @property
    def warnings(self) -> tuple[ImportIssue, ...]:
        """Return warning diagnostics in source order."""
        return tuple(
            item for item in self.issues if item.severity == "warning"
        )

    @property
    def ok(self) -> bool:
        """Whether the import recorded no errors."""
        return not self.errors

    def add(
        self,
        severity: str,
        code: str,
        message: str,
        *,
        row: int | None = None,
        column: str | None = None,
    ) -> None:
        """Append one localized diagnostic."""
        self.issues.append(
            ImportIssue(severity, code, message, row=row, column=column)
        )


class PCBHCSVImportError(ValueError):
    """Raised when strict CSV import records one or more errors."""

    def __init__(self, report: ImportReport):
        self.report = report
        preview = "; ".join(item.message for item in report.errors[:3])
        super().__init__(preview or "PCBH CSV import failed")


def resolve_csv_columns(
    headers: list[str],
    *,
    explicit: dict[str, str] | None,
    constants: dict[str, Any],
    report: ImportReport,
) -> dict[str, str]:
    """Resolve canonical fields from explicit mappings and aliases."""
    normalized: dict[str, list[str]] = {}
    for header in headers:
        normalized.setdefault(_normalize(header), []).append(header)
    duplicates = [values for values in normalized.values() if len(values) > 1]
    for values in duplicates:
        report.add(
            "error",
            "csv.duplicate_header",
            f"headers normalize to the same name: {values!r}",
        )

    result: dict[str, str] = {}
    for canonical, source in (explicit or {}).items():
        if canonical not in CANONICAL_CSV_FIELDS:
            report.add(
                "error",
                "csv.mapping_field",
                f"unknown canonical field {canonical!r}",
            )
            continue
        if source not in headers:
            report.add(
                "error",
                "csv.mapping_column",
                f"mapped source column {source!r} is absent",
                column=source,
            )
            continue
        if canonical in constants:
            report.add(
                "error",
                "csv.mapping_constant_conflict",
                f"{canonical!r} is supplied by both a column and a constant",
            )
            continue
        result[canonical] = source

    explicit_sources: dict[str, list[str]] = {}
    for canonical, source in result.items():
        explicit_sources.setdefault(source, []).append(canonical)
    for source, canonicals in explicit_sources.items():
        if len(canonicals) > 1:
            report.add(
                "error",
                "csv.mapping_source_reused",
                f"source column {source!r} maps to multiple fields "
                f"{canonicals!r}",
                column=source,
            )

    used_sources = set(result.values())
    for canonical in CANONICAL_CSV_FIELDS:
        if canonical in result or canonical in constants:
            continue
        matches: list[str] = []
        candidates = {_normalize(canonical.rsplit(".", 1)[-1])}
        candidates.update(_normalize(item) for item in _ALIASES[canonical])
        for candidate in candidates:
            matches.extend(normalized.get(candidate, []))
        matches = list(dict.fromkeys(matches))
        if len(matches) > 1:
            report.add(
                "error",
                "csv.mapping_ambiguous",
                f"{canonical!r} matches multiple columns {matches!r}; "
                "provide an explicit mapping",
            )
        elif matches and matches[0] not in used_sources:
            result[canonical] = matches[0]
            used_sources.add(matches[0])
            if matches[0] != canonical:
                report.add(
                    "warning",
                    "csv.alias_mapping",
                    f"mapped {matches[0]!r} to {canonical!r}",
                    column=matches[0],
                )

    for canonical in REQUIRED_CSV_FIELDS:
        if canonical not in result and canonical not in constants:
            report.add(
                "error",
                "csv.required_column",
                f"no column or constant supplies {canonical!r}",
            )
    if "crs.horizontal" not in result and "crs.horizontal" not in constants:
        report.add(
            "error",
            "csv.required_crs",
            "supply a CRS column or constants['crs.horizontal']",
        )
    report.column_mapping = dict(result)
    return result


def _normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().casefold()).strip("_")
