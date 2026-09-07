# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Spreadsheet (``.xlsx``) importer for PCBH 0.1.

Field borehole logs almost never arrive in the tidy one-interval-per-row
shape that :func:`pycsamt.format.borehole.boreholes_from_csv` expects: a
typical workbook has a multi-row merged header, a *Rock name* column, a
*From*/*To* depth pair, bilingual labels, and **no collar coordinates**
at all.  This module turns such a sheet into a canonical
:class:`~pycsamt.format.borehole.schema.PCBHDocument` by

1. flattening one chosen sheet to a header row + value rows,
2. resolving canonical PCBH fields from that header (the same alias
   machinery as the CSV importer, plus spreadsheet-specific hints),
3. injecting an external collar position (per hole, or one shared
   collar) supplied out-of-band, and
4. delegating to
   :func:`pycsamt.format.borehole.csvio.assemble_interval_document`.

No pandas dependency: ``openpyxl`` is imported lazily, read-only.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ...api.property import PyCSAMTObject
from .csvio import (
    DEFAULT_CSV_MAX_BYTES,
    DEFAULT_CSV_MAX_ROWS,
    assemble_interval_document,
)
from .mapping import (
    CANONICAL_CSV_FIELDS,
    ImportReport,
    PCBHCSVImportError,
    resolve_csv_columns,
)
from .schema import PCBHDocument

__all__ = [
    "SheetOutline",
    "WorkbookOutline",
    "inspect_workbook",
    "boreholes_from_xlsx",
    "sheet_rows",
]

_PREVIEW_ROWS = 15
_CONSTANT_FIELDS = set(CANONICAL_CSV_FIELDS) | {
    "units.depth",
    "units.resistivity",
}
_MISSING = {"", "na", "n/a", "nan", "none", "null"}


# ---------------------------------------------------------------------------
# workbook inspection
# ---------------------------------------------------------------------------


@dataclass(repr=False)
class SheetOutline(PyCSAMTObject):
    """Lightweight description of one worksheet."""

    name: str
    n_rows: int
    n_cols: int
    preview: list[list[str]] = field(default_factory=list)
    header_candidates: list[int] = field(default_factory=list)


@dataclass(repr=False)
class WorkbookOutline(PyCSAMTObject):
    """All sheets in a workbook, previews only — no bulk data read."""

    source: str
    sheets: list[SheetOutline] = field(default_factory=list)

    def sheet(self, key: str | int) -> SheetOutline:
        """Return one sheet by name or positional index."""
        if isinstance(key, int):
            return self.sheets[key]
        for outline in self.sheets:
            if outline.name == key:
                return outline
        raise KeyError(f"no sheet named {key!r}")


def inspect_workbook(source: str | Path | bytes) -> WorkbookOutline:
    """Return sheet names, shapes, previews, and header-row guesses.

    Parameters
    ----------
    source : path-like or bytes
        An ``.xlsx`` workbook. Bytes are accepted so an application can
        pass an uploaded payload without touching disk.

    Returns
    -------
    WorkbookOutline
        One :class:`SheetOutline` per worksheet.
    """
    raw, name = _read_bytes(source, DEFAULT_CSV_MAX_BYTES)
    workbook = _load_workbook(raw)
    try:
        outlines: list[SheetOutline] = []
        for worksheet in workbook.worksheets:
            rows = _trimmed_rows(worksheet)
            preview = [
                [_cell_text(cell) for cell in row]
                for row in rows[:_PREVIEW_ROWS]
            ]
            outlines.append(
                SheetOutline(
                    name=worksheet.title,
                    n_rows=len(rows),
                    n_cols=max((len(row) for row in rows), default=0),
                    preview=preview,
                    header_candidates=_header_candidates(rows),
                )
            )
    finally:
        workbook.close()
    return WorkbookOutline(source=name, sheets=outlines)


def sheet_rows(
    source: str | Path | bytes,
    *,
    sheet: str | int | None = None,
    header_row: int | None = None,
) -> tuple[str, list[str], list[list[str]]]:
    """Flatten one worksheet to ``(sheet_name, headers, data_rows)``.

    A thin public wrapper around the workbook flattening used by
    :func:`boreholes_from_xlsx`, for callers that resolve their own
    canonical fields (e.g. :mod:`pycsamt.format.pointset`). ``header_row``
    is 0-based; auto-detected when omitted.
    """
    raw, _ = _read_bytes(source, DEFAULT_CSV_MAX_BYTES)
    workbook = _load_workbook(raw)
    try:
        worksheet = _resolve_sheet(workbook, sheet)
        name = worksheet.title
        rows = _trimmed_rows(worksheet)
    finally:
        workbook.close()
    if not rows:
        return name, [], []
    if header_row is None:
        header_row = next(
            (i for i, r in enumerate(rows) if _looks_like_header(r)), 0
        )
    header_row = max(0, min(int(header_row), len(rows) - 1))
    headers = _header_names(rows[header_row])
    data = [
        _normalise_row(r, len(headers))
        for r in rows[header_row + 1 :]
    ]
    return name, headers, data


# ---------------------------------------------------------------------------
# import
# ---------------------------------------------------------------------------


def boreholes_from_xlsx(
    source: str | Path | bytes,
    *,
    sheet: str | int | None = None,
    header_row: int | None = None,
    columns: dict[str, Any] | None = None,
    constants: dict[str, Any] | None = None,
    collars: dict[str, Any] | str | Path | None = None,
    default_borehole_id: str | None = None,
    strict: bool = True,
    document_id: str | None = None,
    created_by: str = "pycsamt XLSX importer",
    max_bytes: int = DEFAULT_CSV_MAX_BYTES,
    max_rows: int = DEFAULT_CSV_MAX_ROWS,
) -> tuple[PCBHDocument, ImportReport]:
    """Import one worksheet of a borehole workbook as a PCBH document.

    Parameters
    ----------
    source : path-like or bytes
        The ``.xlsx`` workbook.
    sheet : str or int, optional
        Worksheet name or 0-based index. Defaults to the first sheet.
    header_row : int, optional
        0-based row index of the header. Auto-detected when omitted
        (the first row with at least two text cells).
    columns : dict, optional
        Explicit ``{canonical_field: source}`` mapping. ``source`` may be
        a header string, a 0-based column index, or an Excel column
        letter. Canonical names use dotted paths such as ``interval.from_md``.
    constants : dict, optional
        Constant canonical values (commonly ``borehole.id`` and
        ``crs.horizontal`` for a single-hole sheet).
    collars : dict or path-like, optional
        External collar positions, keyed by borehole id
        (``{"ZK2203": {"x": .., "y": .., "z": .., "crs": "EPSG:32648"}}``),
        or a single ``{"x": .., "y": .., "z": ..}`` applied to every hole,
        or a small CSV with ``id,x,y,z[,crs]`` columns. Spreadsheets rarely
        carry coordinates, so this is how they enter the document.
    strict : bool, default=True
        Raise :class:`PCBHCSVImportError` on any recorded error. In
        permissive mode a missing collar/CRS becomes a flagged placeholder
        and invalid rows are dropped.
    document_id : str, optional
        Defaults to ``xlsx:<file stem>``.
    created_by : str, default='pycsamt XLSX importer'
        Provenance name stored on the document.
    max_bytes, max_rows : int
        Resource limits, shared with the CSV importer.

    Returns
    -------
    document : PCBHDocument
    report : ImportReport
        Includes the resolved sheet, header row, column mapping, injected
        collars, and every rejected row.
    """
    if columns is not None and not isinstance(columns, dict):
        raise TypeError("columns must be a dict or None")
    if constants is not None and not isinstance(constants, dict):
        raise TypeError("constants must be a dict or None")
    if not isinstance(strict, bool):
        raise TypeError("strict must be a boolean")

    raw, name = _read_bytes(source, max_bytes)
    supplied_constants = dict(constants or {})
    report = ImportReport(
        source=name,
        source_sha256=hashlib.sha256(raw).hexdigest(),
        delimiter="xlsx",
        strict=strict,
        constants=supplied_constants,
    )
    for key in supplied_constants:
        if key not in _CONSTANT_FIELDS:
            report.add(
                "error",
                "xlsx.constant_field",
                f"unknown constant field {key!r}",
            )

    workbook = _load_workbook(raw)
    try:
        worksheet = _resolve_sheet(workbook, sheet)
        report.source_files["sheet"] = worksheet.title
        rows = _trimmed_rows(worksheet)
    finally:
        workbook.close()

    if not rows:
        report.add("error", "xlsx.empty", "worksheet has no rows")
        raise PCBHCSVImportError(report)

    header_index = _resolve_header_row(rows, header_row, report)
    headers = _header_names(rows[header_index])
    report.inferred_values.append(
        f"header row = {header_index + 1} (1-based)"
    )
    if len(headers) < 2:
        report.add(
            "error",
            "xlsx.header",
            "header row needs at least two columns",
            row=header_index + 1,
        )

    data_rows = [
        _normalise_row(row, len(headers))
        for row in rows[header_index + 1 :]
    ]

    explicit = _resolve_column_tokens(columns, headers, report)
    # Preflight the alias resolver on a throwaway report so collar
    # injection can key on the borehole-id column even when it was
    # matched by alias rather than named explicitly.
    scratch = ImportReport(
        source=name, source_sha256="", delimiter="xlsx", strict=False
    )
    preflight = resolve_csv_columns(
        headers,
        explicit=explicit,
        constants=supplied_constants,
        report=scratch,
    )
    if (
        default_borehole_id
        and "borehole.id" not in preflight
        and "borehole.id" not in explicit
        and "borehole.id" not in supplied_constants
    ):
        supplied_constants["borehole.id"] = str(default_borehole_id)
        report.inferred_values.append(
            f"borehole.id = {default_borehole_id!r} (no id column found)"
        )
    headers, data_rows, explicit = _inject_collars(
        headers,
        data_rows,
        explicit,
        collars=collars,
        constants=supplied_constants,
        id_source=explicit.get("borehole.id") or preflight.get("borehole.id"),
        report=report,
        strict=strict,
    )
    if strict:
        _guard_required_present(explicit, supplied_constants, report)

    mapping = resolve_csv_columns(
        headers,
        explicit=explicit,
        constants=supplied_constants,
        report=report,
    )
    if report.errors:
        raise PCBHCSVImportError(report)

    return assemble_interval_document(
        headers,
        iter(data_rows),
        mapping=mapping,
        constants=supplied_constants,
        report=report,
        strict=strict,
        document_id=document_id or f"xlsx:{Path(name).stem}",
        created_by=created_by,
        source_str=name,
        max_rows=max_rows,
        row_offset=header_index + 2,
        source_kind="xlsx",
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _read_bytes(source: Any, max_bytes: int) -> tuple[bytes, str]:
    if isinstance(source, (bytes, bytearray)):
        raw = bytes(source)
        name = "workbook.xlsx"
    elif hasattr(source, "read"):
        raw = source.read()
        name = getattr(source, "name", "workbook.xlsx")
    else:
        path = Path(source)
        raw = path.read_bytes()
        name = str(path)
    if len(raw) > max_bytes:
        raise ValueError(
            f"workbook is {len(raw)} bytes; limit is {max_bytes} bytes"
        )
    return raw, name


def _load_workbook(raw: bytes):
    try:
        import openpyxl
    except ModuleNotFoundError as error:  # pragma: no cover - packaging
        raise ModuleNotFoundError(
            "reading .xlsx workbooks requires 'openpyxl'"
        ) from error
    try:
        return openpyxl.load_workbook(
            io.BytesIO(raw), read_only=True, data_only=True
        )
    except Exception as error:  # noqa: BLE001 - surface a clean message
        raise ValueError(f"could not open workbook: {error}") from error


def _resolve_sheet(workbook, sheet):
    names = workbook.sheetnames
    if sheet is None:
        return workbook[names[0]]
    if isinstance(sheet, int):
        if not 0 <= sheet < len(names):
            raise ValueError(
                f"sheet index {sheet} out of range (0..{len(names) - 1})"
            )
        return workbook[names[sheet]]
    if sheet not in names:
        raise ValueError(f"no sheet named {sheet!r}; have {names!r}")
    return workbook[sheet]


def _trimmed_rows(worksheet) -> list[list[Any]]:
    rows = [list(row) for row in worksheet.iter_rows(values_only=True)]
    while rows and all(_is_blank(cell) for cell in rows[-1]):
        rows.pop()
    return rows


def _is_blank(cell: Any) -> bool:
    return cell is None or (isinstance(cell, str) and not cell.strip())


def _cell_text(cell: Any) -> str:
    if cell is None:
        return ""
    if isinstance(cell, str):
        return cell.strip()
    if isinstance(cell, bool):
        return "true" if cell else "false"
    if isinstance(cell, float):
        if cell.is_integer():
            return str(int(cell))
        return repr(cell)
    if isinstance(cell, (_dt.datetime, _dt.date)):
        return cell.isoformat()
    return str(cell).strip()


def _normalise_row(row: list[Any], width: int) -> list[str]:
    values = [_cell_text(cell) for cell in row]
    if len(values) < width:
        values.extend([""] * (width - len(values)))
    elif len(values) > width:
        values = values[:width]
    return values


def _header_names(row: list[Any]) -> list[str]:
    names: list[str] = []
    seen: dict[str, int] = {}
    for index, cell in enumerate(row, start=1):
        text = _cell_text(cell) or f"column_{index}"
        if text in seen:
            seen[text] += 1
            text = f"{text}_{seen[text]}"
        else:
            seen[text] = 1
        names.append(text)
    return names


def _numeric_count(row: list[Any]) -> int:
    count = 0
    for cell in row:
        if isinstance(cell, (int, float)) and not isinstance(cell, bool):
            count += 1
        elif isinstance(cell, str):
            try:
                float(cell.strip())
                count += 1
            except ValueError:
                pass
    return count


def _text_count(row: list[Any]) -> int:
    return sum(
        1
        for cell in row
        if isinstance(cell, str)
        and cell.strip()
        and not _is_number_text(cell)
    )


def _is_number_text(cell: str) -> bool:
    try:
        float(cell.strip())
        return True
    except ValueError:
        return False


def _looks_like_header(row: list[Any]) -> bool:
    return _text_count(row) >= 2 and _numeric_count(row) == 0


def _header_candidates(rows: list[list[Any]]) -> list[int]:
    candidates: list[int] = []
    for index, row in enumerate(rows[:_PREVIEW_ROWS]):
        if not _looks_like_header(row):
            continue
        following = rows[index + 1 : index + 6]
        if any(_numeric_count(later) >= 2 for later in following):
            candidates.append(index)
    return candidates


def _resolve_header_row(
    rows: list[list[Any]], header_row: int | None, report: ImportReport
) -> int:
    if header_row is not None:
        if not 0 <= header_row < len(rows):
            report.add(
                "error",
                "xlsx.header_row",
                f"header_row {header_row} is outside 0..{len(rows) - 1}",
            )
            raise PCBHCSVImportError(report)
        return header_row
    for index, row in enumerate(rows):
        if _looks_like_header(row):
            return index
    return 0


def _resolve_column_tokens(
    columns: dict[str, Any] | None,
    headers: list[str],
    report: ImportReport,
) -> dict[str, str]:
    resolved: dict[str, str] = {}
    for canonical, token in (columns or {}).items():
        header = _token_to_header(token, headers)
        if header is None:
            report.add(
                "error",
                "xlsx.column_token",
                f"cannot resolve column {token!r} for {canonical!r}",
            )
            continue
        resolved[canonical] = header
    return resolved


def _token_to_header(token: Any, headers: list[str]) -> str | None:
    if isinstance(token, int) and not isinstance(token, bool):
        return headers[token] if 0 <= token < len(headers) else None
    text = str(token).strip()
    if text in headers:
        return text
    if text.isalpha():
        index = 0
        for char in text.upper():
            index = index * 26 + (ord(char) - ord("A") + 1)
        index -= 1
        return headers[index] if 0 <= index < len(headers) else None
    return None


def _normalise_collars(
    collars: dict[str, Any] | str | Path | None,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any] | None]:
    """Return ``(per_id, shared)`` collar tables."""
    if collars is None:
        return {}, None
    if isinstance(collars, (str, Path)):
        import csv as _csv

        per_id: dict[str, dict[str, Any]] = {}
        with open(collars, newline="", encoding="utf-8-sig") as handle:
            for row in _csv.DictReader(handle):
                key = str(
                    row.get("id") or row.get("borehole_id") or ""
                ).strip()
                if not key:
                    continue
                per_id[key] = {
                    "x": row.get("x"),
                    "y": row.get("y"),
                    "z": row.get("z"),
                    "crs": row.get("crs"),
                }
        return per_id, None
    if not isinstance(collars, dict):
        raise TypeError("collars must be a dict, path, or None")
    lowered = {str(k).lower() for k in collars}
    if {"x", "y", "z"} <= lowered:
        return {}, {str(k).lower(): v for k, v in collars.items()}
    return {str(k): dict(v) for k, v in collars.items()}, None


def _inject_collars(
    headers: list[str],
    data_rows: list[list[str]],
    explicit: dict[str, str],
    *,
    collars: dict[str, Any] | str | Path | None,
    constants: dict[str, Any],
    id_source: str | None,
    report: ImportReport,
    strict: bool,
) -> tuple[list[str], list[list[str]], dict[str, str]]:
    per_id, shared = _normalise_collars(collars)
    have_collar = per_id or shared
    id_index = headers.index(id_source) if id_source in headers else None
    id_constant = str(constants.get("borehole.id") or "").strip()

    synth = {
        "collar.x": "__pcbh_collar_x",
        "collar.y": "__pcbh_collar_y",
        "collar.z": "__pcbh_collar_z",
    }
    crs_synth = "__pcbh_crs"
    add_crs = have_collar and (
        shared and shared.get("crs")
        or any(entry.get("crs") for entry in per_id.values())
    )

    if not have_collar:
        if strict:
            return headers, data_rows, explicit
        # permissive: placeholder collar + CRS so a preview still renders
        for canonical, value in (
            ("collar.x", 0.0),
            ("collar.y", 0.0),
            ("collar.z", 0.0),
        ):
            if canonical not in explicit and canonical not in constants:
                constants[canonical] = value
        if (
            "crs.horizontal" not in explicit
            and "crs.horizontal" not in constants
        ):
            constants["crs.horizontal"] = "LOCAL:unknown"
        report.add(
            "warning",
            "xlsx.placeholder_collar",
            "no collar coordinates supplied; using a (0, 0, 0) "
            "LOCAL:unknown placeholder — set real collars before use",
        )
        return headers, data_rows, explicit

    new_headers = list(headers) + list(synth.values())
    if add_crs:
        new_headers.append(crs_synth)
    injected: set[str] = set()
    for row in data_rows:
        key = ""
        if id_index is not None and id_index < len(row):
            key = row[id_index].strip()
        key = key or id_constant
        entry = per_id.get(key) if per_id else shared
        if entry is None and shared is not None:
            entry = shared
        if entry:
            injected.add(key or "<all>")
            row.append(_cell_or_blank(entry.get("x")))
            row.append(_cell_or_blank(entry.get("y")))
            row.append(_cell_or_blank(entry.get("z")))
            if add_crs:
                row.append(_cell_or_blank(entry.get("crs")))
        else:
            row.extend(["", "", ""])
            if add_crs:
                row.append("")

    explicit = dict(explicit)
    explicit.update(synth)
    if add_crs and "crs.horizontal" not in constants:
        explicit["crs.horizontal"] = crs_synth
    report.inferred_values.append(
        f"injected external collars for: {sorted(injected)}"
    )
    return new_headers, data_rows, explicit


def _cell_or_blank(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _guard_required_present(
    explicit: dict[str, str],
    constants: dict[str, Any],
    report: ImportReport,
) -> None:
    for canonical in ("collar.x", "collar.y", "collar.z", "crs.horizontal"):
        if canonical in explicit or canonical in constants:
            continue
        report.add(
            "error",
            "xlsx.missing_collar",
            f"{canonical!r} has no column, constant, or collar entry; "
            "pass collars= or constants=, or use strict=False",
        )
    if report.errors:
        raise PCBHCSVImportError(report)
