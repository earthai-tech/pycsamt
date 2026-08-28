# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Combined interval-CSV importer for PCBH 0.1."""

from __future__ import annotations

import csv
import hashlib
import io
import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .mapping import (
    CANONICAL_CSV_FIELDS,
    ImportReport,
    PCBHCSVImportError,
    resolve_csv_columns,
)
from .schema import (
    BOREHOLE_KINDS,
    BOREHOLE_STATUSES,
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    PCBHValidationError,
    Trajectory,
    UnitSystem,
    VocabularyEntry,
)

__all__ = [
    "DEFAULT_CSV_MAX_BYTES",
    "DEFAULT_CSV_MAX_ROWS",
    "boreholes_from_csv",
]

DEFAULT_CSV_MAX_BYTES = 10 * 1024 * 1024
DEFAULT_CSV_MAX_ROWS = 250_000
_DELIMITERS = (",", ";", "\t", "|")
_MISSING = {"", "na", "n/a", "nan", "none", "null"}
_CONSTANT_FIELDS = set(CANONICAL_CSV_FIELDS) | {
    "units.depth",
    "units.resistivity",
}


@dataclass
class _HoleRows:
    borehole_id: str
    name: str
    x: float
    y: float
    z: float
    crs: str
    kind: str
    status: str
    total_depth_md: float | None
    intervals: list[LogInterval] = field(default_factory=list)


def boreholes_from_csv(
    path: str | Path,
    *,
    columns: dict[str, str] | None = None,
    constants: dict[str, Any] | None = None,
    strict: bool = True,
    delimiter: str | None = None,
    document_id: str | None = None,
    created_by: str = "pycsamt CSV importer",
    max_bytes: int = DEFAULT_CSV_MAX_BYTES,
    max_rows: int = DEFAULT_CSV_MAX_ROWS,
) -> tuple[PCBHDocument, ImportReport]:
    """Import a combined collar-and-interval CSV as a PCBH document.

    Parameters
    ----------
    path : path-like
        UTF-8 CSV containing repeated collar fields and interval rows.
    columns : dict, optional
        Explicit ``{canonical_field: source_header}`` mapping. Canonical names
        use dotted paths such as ``borehole.id`` and ``interval.from_md``.
    constants : dict, optional
        Constant canonical values, commonly ``crs.horizontal``. Constants
        take precedence over mapped row values.
    strict : bool, default=True
        Raise :class:`PCBHCSVImportError` if any error is recorded. In
        permissive mode, invalid rows are rejected and valid rows returned.
    delimiter : {',', ';', '\\t', '|'}, optional
        Explicit delimiter. If omitted, detection is restricted to this set.
    document_id : str, optional
        PCBH document identifier. Defaults to ``csv:<file stem>``.
    created_by : str, default='pycsamt CSV importer'
        Provenance name stored on the document.
    max_bytes : int, default=10485760
        Maximum source size in bytes.
    max_rows : int, default=250000
        Maximum number of data rows.

    Returns
    -------
    document : PCBHDocument
        Valid document containing every accepted row.
    report : ImportReport
        Source checksum, mappings, inferences, counts, and diagnostics.

    Raises
    ------
    PCBHCSVImportError
        If the file structure is unusable, no valid boreholes remain, or
        strict mode records an error. The exception exposes ``report``.
    ValueError
        If a resource limit or delimiter parameter is invalid.

    Notes
    -----
    Missing tokens are normalized before conversion and never stringified.
    Repeated collar and total-depth values must agree within a borehole.
    """
    if columns is not None and not isinstance(columns, dict):
        raise TypeError("columns must be a dict or None")
    if constants is not None and not isinstance(constants, dict):
        raise TypeError("constants must be a dict or None")
    if not isinstance(strict, bool):
        raise TypeError("strict must be a boolean")
    if not isinstance(created_by, str) or not created_by.strip():
        raise ValueError("created_by must be a non-empty string")
    source = Path(path)
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int):
        raise TypeError("max_bytes must be an integer")
    if max_bytes <= 0:
        raise ValueError("max_bytes must be greater than zero")
    if isinstance(max_rows, bool) or not isinstance(max_rows, int):
        raise TypeError("max_rows must be an integer")
    if max_rows <= 0:
        raise ValueError("max_rows must be greater than zero")
    raw = source.read_bytes()
    if len(raw) > max_bytes:
        raise ValueError(
            f"CSV input is {len(raw)} bytes; limit is {max_bytes} bytes"
        )
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError("PCBH CSV input must be UTF-8") from error
    chosen_delimiter = _delimiter(text, delimiter)
    supplied_constants = dict(constants or {})
    report = ImportReport(
        source=str(source),
        source_sha256=hashlib.sha256(raw).hexdigest(),
        delimiter=chosen_delimiter,
        strict=strict,
        constants=supplied_constants,
    )
    for key in supplied_constants:
        if key not in _CONSTANT_FIELDS:
            report.add(
                "error",
                "csv.constant_field",
                f"unknown constant field {key!r}",
            )

    reader = csv.reader(io.StringIO(text), delimiter=chosen_delimiter)
    try:
        headers = next(reader)
    except StopIteration:
        report.add("error", "csv.empty", "CSV file is empty")
        raise PCBHCSVImportError(report)
    headers = [header.strip() for header in headers]
    if len(headers) < 2 or any(not header for header in headers):
        report.add(
            "error",
            "csv.header",
            "CSV requires a non-empty header row with at least two columns",
            row=1,
        )
    mapping = resolve_csv_columns(
        headers,
        explicit=columns,
        constants=supplied_constants,
        report=report,
    )
    if report.errors:
        raise PCBHCSVImportError(report)

    groups: dict[str, _HoleRows] = {}
    vocabulary: dict[str, VocabularyEntry] = {}
    label_codes: dict[str, str] = {}
    global_crs: str | None = None
    for row_number, values in enumerate(reader, start=2):
        report.rows_read += 1
        if report.rows_read > max_rows:
            report.add(
                "error",
                "csv.row_limit",
                f"CSV exceeds the {max_rows}-row limit",
                row=row_number,
            )
            break
        if not values or all(_missing(value) for value in values):
            report.rows_skipped += 1
            continue
        errors_before = len(report.errors)
        if len(values) != len(headers):
            report.add(
                "error",
                "csv.row_width",
                f"expected {len(headers)} columns, found {len(values)}",
                row=row_number,
            )
            report.rows_rejected += 1
            continue
        row = dict(zip(headers, values))
        parsed = _parse_row(
            row,
            row_number,
            mapping,
            supplied_constants,
            report,
        )
        if parsed is None:
            report.rows_rejected += 1
            continue
        hole_id = parsed["borehole.id"]
        group = groups.get(hole_id)
        if group is None:
            group = _new_group(parsed, report)
            groups[hole_id] = group
        elif not _consistent_group(group, parsed, row_number, report):
            report.rows_rejected += 1
            continue
        if global_crs is not None and parsed["crs.horizontal"] != global_crs:
            report.add(
                "error",
                "csv.crs_conflict",
                f"CRS {parsed['crs.horizontal']!r} conflicts with "
                f"document CRS {global_crs!r}",
                row=row_number,
                column=mapping.get("crs.horizontal"),
            )
            report.rows_rejected += 1
            if not group.intervals:
                groups.pop(hole_id, None)
            continue

        interval = parsed["interval"]
        if _overlaps(interval, group.intervals):
            report.add(
                "error",
                "csv.interval_overlap",
                f"interval [{interval.from_md}, {interval.to_md}) overlaps "
                f"another interval for {hole_id!r}",
                row=row_number,
            )
            report.rows_rejected += 1
            continue
        if (
            group.total_depth_md is not None
            and interval.to_md > group.total_depth_md
        ):
            report.add(
                "error",
                "csv.interval_beyond_td",
                "interval extends beyond total_depth_md",
                row=row_number,
                column=mapping.get("interval.to_md"),
            )
            report.rows_rejected += 1
            continue
        if len(report.errors) > errors_before:
            report.rows_rejected += 1
            continue
        if not _register_lithology(
            interval,
            vocabulary,
            label_codes,
            row_number,
            report,
        ):
            report.rows_rejected += 1
            continue
        if global_crs is None:
            global_crs = parsed["crs.horizontal"]
        group.intervals.append(interval)
        report.rows_accepted += 1

    if strict and report.errors:
        raise PCBHCSVImportError(report)
    groups = {key: value for key, value in groups.items() if value.intervals}
    if not groups or global_crs is None:
        report.add(
            "error",
            "csv.no_valid_rows",
            "CSV import produced no valid boreholes",
        )
        raise PCBHCSVImportError(report)

    boreholes = [_build_borehole(group, report) for group in groups.values()]
    units = UnitSystem(
        depth=str(supplied_constants.get("units.depth", "m")),
        resistivity=str(supplied_constants.get("units.resistivity", "ohm.m")),
    )
    if units.depth != "m" or units.resistivity != "ohm.m":
        report.add(
            "error",
            "csv.units_unsupported",
            "combined CSV currently requires metres and ohm metres; "
            "convert source units before import",
        )
        raise PCBHCSVImportError(report)
    document = PCBHDocument(
        document_id=document_id or f"csv:{source.stem}",
        created_at=datetime.now(timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        created_by=created_by,
        crs=CoordinateReferenceSystem(global_crs),
        units=units,
        boreholes=boreholes,
        lithologies=list(vocabulary.values()),
        metadata={
            "csv_source": str(source),
            "csv_source_sha256": report.source_sha256,
        },
    )
    try:
        document.validate()
    except PCBHValidationError as error:
        for issue in error.issues:
            report.add(
                "error",
                f"csv.document.{issue.code}",
                f"{issue.path}: {issue.message}",
            )
        raise PCBHCSVImportError(report) from error
    return document, report


def _delimiter(text: str, supplied: str | None) -> str:
    if supplied is not None:
        if supplied not in _DELIMITERS:
            raise ValueError(f"delimiter must be one of {_DELIMITERS!r}")
        return supplied
    sample = text[:65536]
    try:
        return (
            csv.Sniffer()
            .sniff(sample, delimiters="".join(_DELIMITERS))
            .delimiter
        )
    except csv.Error as error:
        raise ValueError(
            "could not detect CSV delimiter; pass delimiter explicitly"
        ) from error


def _parse_row(
    row: dict[str, str],
    row_number: int,
    mapping: dict[str, str],
    constants: dict[str, Any],
    report: ImportReport,
) -> dict[str, Any] | None:
    def value(canonical: str) -> Any:
        raw = constants.get(canonical)
        if canonical not in constants:
            source = mapping.get(canonical)
            raw = row.get(source, "") if source else None
        return None if _missing(raw) else raw

    hole_id = value("borehole.id")
    lithology = value("interval.lithology")
    if hole_id is None:
        report.add(
            "error", "csv.missing_id", "borehole ID is missing", row=row_number
        )
    elif not str(hole_id).strip():
        report.add(
            "error", "csv.missing_id", "borehole ID is blank", row=row_number
        )
    if lithology is None:
        report.add(
            "error",
            "csv.missing_lithology",
            "lithology is missing",
            row=row_number,
            column=mapping.get("interval.lithology"),
        )
    elif not str(lithology).strip():
        report.add(
            "error",
            "csv.missing_lithology",
            "lithology is blank",
            row=row_number,
            column=mapping.get("interval.lithology"),
        )
    numeric: dict[str, float | None] = {}
    for canonical in (
        "collar.x",
        "collar.y",
        "collar.z",
        "borehole.total_depth_md",
        "interval.from_md",
        "interval.to_md",
        "interval.resistivity_ohm_m",
    ):
        raw = value(canonical)
        required = canonical in {
            "collar.x",
            "collar.y",
            "collar.z",
            "interval.from_md",
            "interval.to_md",
        }
        numeric[canonical] = _number(
            raw,
            canonical,
            row_number,
            mapping.get(canonical),
            report,
            required=required,
        )
    crs = value("crs.horizontal")
    if crs is None:
        report.add(
            "error", "csv.missing_crs", "CRS is missing", row=row_number
        )
    elif not str(crs).strip():
        report.add("error", "csv.missing_crs", "CRS is blank", row=row_number)
    if report.errors and report.errors[-1].row == row_number:
        return None
    from_md = numeric["interval.from_md"]
    to_md = numeric["interval.to_md"]
    resistivity = numeric["interval.resistivity_ohm_m"]
    if from_md < 0 or to_md <= from_md:
        report.add(
            "error",
            "csv.interval_bounds",
            "interval requires 0 <= from_md < to_md",
            row=row_number,
        )
        return None
    if resistivity is not None and resistivity <= 0:
        report.add(
            "error",
            "csv.resistivity",
            "resistivity must be greater than zero",
            row=row_number,
        )
        return None
    total_depth = numeric["borehole.total_depth_md"]
    if total_depth is not None and total_depth <= 0:
        report.add(
            "error",
            "csv.total_depth",
            "total_depth_md must be greater than zero",
            row=row_number,
        )
        return None
    name_value = value("borehole.name")
    kind_value = value("borehole.kind")
    status_value = value("borehole.status")
    inferred = {
        canonical
        for canonical, raw in (
            ("borehole.name", name_value),
            ("borehole.kind", kind_value),
            ("borehole.status", status_value),
        )
        if raw is None
    }
    kind = str(kind_value or "unknown").strip()
    status = str(status_value or "unknown").strip()
    if not _controlled(kind, BOREHOLE_KINDS):
        report.add(
            "error",
            "csv.borehole_kind",
            f"unsupported borehole kind {kind!r}",
            row=row_number,
        )
        return None
    if not _controlled(status, BOREHOLE_STATUSES):
        report.add(
            "error",
            "csv.borehole_status",
            f"unsupported borehole status {status!r}",
            row=row_number,
        )
        return None
    data_nature = str(value("interval.data_nature") or "observed")
    code_value = value("interval.code")
    interval = LogInterval(
        from_md=from_md,
        to_md=to_md,
        code=str(code_value).strip() if code_value else None,
        label=str(lithology).strip(),
        description=str(value("interval.description") or "").strip(),
        resistivity_ohm_m=resistivity,
        data_nature=data_nature,
    )
    interval_errors = [
        issue
        for issue in interval.collect_issues()
        if issue.severity == "error"
    ]
    if interval_errors:
        for issue in interval_errors:
            report.add(
                "error",
                f"csv.{issue.code}",
                issue.message,
                row=row_number,
            )
        return None
    return {
        "borehole.id": str(hole_id).strip(),
        "borehole.name": str(name_value or hole_id).strip(),
        "borehole.kind": kind,
        "borehole.status": status,
        "borehole.total_depth_md": total_depth,
        "collar.x": numeric["collar.x"],
        "collar.y": numeric["collar.y"],
        "collar.z": numeric["collar.z"],
        "crs.horizontal": str(crs).strip(),
        "interval": interval,
        "_inferred": inferred,
    }


def _number(
    value: Any,
    canonical: str,
    row: int,
    column: str | None,
    report: ImportReport,
    *,
    required: bool,
) -> float | None:
    if value is None:
        if required:
            report.add(
                "error",
                "csv.missing_number",
                f"{canonical} is missing",
                row=row,
                column=column,
            )
        return None
    if isinstance(value, bool):
        number = math.nan
    else:
        try:
            number = float(value)
        except (TypeError, ValueError):
            number = math.nan
    if not math.isfinite(number):
        report.add(
            "error",
            "csv.invalid_number",
            f"{canonical} must be a finite number, got {value!r}",
            row=row,
            column=column,
        )
        return None
    return number


def _new_group(parsed: dict[str, Any], report: ImportReport) -> _HoleRows:
    hole_id = parsed["borehole.id"]
    for field_name, default in (
        ("borehole.name", hole_id),
        ("borehole.kind", "unknown"),
        ("borehole.status", "unknown"),
    ):
        if field_name in parsed["_inferred"]:
            text = f"{field_name}={default!r} for {hole_id!r}"
            if text not in report.inferred_values:
                report.inferred_values.append(text)
    return _HoleRows(
        borehole_id=hole_id,
        name=parsed["borehole.name"],
        x=parsed["collar.x"],
        y=parsed["collar.y"],
        z=parsed["collar.z"],
        crs=parsed["crs.horizontal"],
        kind=parsed["borehole.kind"],
        status=parsed["borehole.status"],
        total_depth_md=parsed["borehole.total_depth_md"],
    )


def _consistent_group(
    group: _HoleRows,
    parsed: dict[str, Any],
    row: int,
    report: ImportReport,
) -> bool:
    comparisons = (
        ("collar.x", group.x),
        ("collar.y", group.y),
        ("collar.z", group.z),
        ("crs.horizontal", group.crs),
        ("borehole.name", group.name),
        ("borehole.kind", group.kind),
        ("borehole.status", group.status),
        ("borehole.total_depth_md", group.total_depth_md),
    )
    for canonical, expected in comparisons:
        actual = parsed[canonical]
        if canonical == "borehole.total_depth_md" and actual is None:
            continue
        if expected is None and canonical == "borehole.total_depth_md":
            group.total_depth_md = actual
            continue
        equal = (
            math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-9)
            if isinstance(expected, float)
            else actual == expected
        )
        if not equal:
            message = (
                f"{canonical}={actual!r} conflicts with first value "
                f"{expected!r} for {group.borehole_id!r}"
            )
            report.add("error", "csv.collar_conflict", message, row=row)
            report.conflict_resolutions.append(
                message + "; first value retained and row rejected"
            )
            return False
    return True


def _overlaps(interval: LogInterval, existing: list[LogInterval]) -> bool:
    return any(
        interval.from_md < other.to_md and interval.to_md > other.from_md
        for other in existing
    )


def _register_lithology(
    interval: LogInterval,
    vocabulary: dict[str, VocabularyEntry],
    label_codes: dict[str, str],
    row: int,
    report: ImportReport,
) -> bool:
    label = interval.label or "unknown"
    folded = label.casefold()
    code = interval.code or label_codes.get(folded)
    if code is None:
        base = re.sub(r"[^A-Z0-9]+", "_", label.upper()).strip("_")
        base = base or "UNKNOWN"
        code = base
        suffix = 2
        while (
            code in vocabulary and vocabulary[code].name.casefold() != folded
        ):
            code = f"{base}_{suffix}"
            suffix += 1
        report.inferred_values.append(
            f"interval.code={code!r} for lithology {label!r}"
        )
    current = vocabulary.get(code)
    if current is not None and current.name.casefold() != folded:
        report.add(
            "error",
            "csv.lithology_code_conflict",
            f"lithology code {code!r} identifies both {current.name!r} "
            f"and {label!r}",
            row=row,
        )
        return False
    interval.code = code
    label_codes[folded] = code
    vocabulary.setdefault(code, VocabularyEntry(code=code, name=label))
    return True


def _build_borehole(group: _HoleRows, report: ImportReport) -> PCBHBorehole:
    intervals = sorted(group.intervals, key=lambda item: item.from_md)
    total_depth = group.total_depth_md
    if total_depth is None:
        total_depth = max(item.to_md for item in intervals)
        report.inferred_values.append(
            f"borehole.total_depth_md={total_depth!r} for "
            f"{group.borehole_id!r} from deepest interval"
        )
    return PCBHBorehole(
        id=group.borehole_id,
        name=group.name,
        kind=group.kind,
        status=group.status,
        collar=Collar(group.x, group.y, group.z),
        total_depth_md=total_depth,
        trajectory=Trajectory(method="vertical"),
        interval_logs={"lithology": intervals},
    )


def _missing(value: Any) -> bool:
    return value is None or (
        isinstance(value, str) and value.strip().casefold() in _MISSING
    )


def _controlled(value: str, allowed: tuple[str, ...]) -> bool:
    left, separator, right = value.partition(":")
    return value in allowed or bool(separator and left and right)
