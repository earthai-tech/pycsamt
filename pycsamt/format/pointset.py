# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PCPT - pyCSAMT Points: a tiny format for targets and points of interest.

PCPT is deliberately small. It carries a flat list of named positions -
drill targets, planned collars, samples, anomalies, generic markers -
each with a location and optional depth window, so a map or a 3-D scene
can be annotated *independently* of a full :mod:`pycsamt.format.borehole`
document. No trajectory, no interval logs, no vocabularies.

Canonical encoding is UTF-8 JSON, canonical extension ``.pcpt.json``.
CSV and spreadsheet importers resolve the same canonical fields through
a small alias table.
"""

from __future__ import annotations

import csv
import io
import json
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from math import isfinite
from pathlib import Path
from typing import Any

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = [
    "PCPT_VERSION",
    "POINT_KINDS",
    "Point",
    "PointSet",
    "PointSetValidationError",
    "read_points",
    "write_points",
    "points_from_csv",
    "points_from_xlsx",
    "point_set_to_dict",
    "point_set_from_dict",
]

PCPT_VERSION = "0.1.0"
POINT_KINDS = (
    "target",
    "planned_borehole",
    "sample",
    "anomaly",
    "poi",
    "other",
)
_MISSING = {"", "na", "n/a", "nan", "none", "null"}

_FIELDS = (
    "id", "name", "x", "y", "longitude", "latitude", "z",
    "depth_top", "depth_bottom", "kind", "symbol", "color",
    "note", "category",
)
_ALIASES = {
    "id": ("id", "point_id", "target_id", "name_id", "code"),
    "name": ("name", "label", "title", "target", "target_name"),
    "x": ("x", "easting", "east", "utm_e", "utmx"),
    "y": ("y", "northing", "north", "utm_n", "utmy"),
    "longitude": ("longitude", "lon", "long", "lng"),
    "latitude": ("latitude", "lat"),
    "z": ("z", "elevation", "elev", "rl", "altitude"),
    "depth_top": ("depth_top", "depth_from", "from", "top", "from_m"),
    "depth_bottom": ("depth_bottom", "depth_to", "to", "bottom", "to_m"),
    "kind": ("kind", "type", "category_kind"),
    "symbol": ("symbol", "marker"),
    "color": ("color", "colour", "hex"),
    "note": ("note", "notes", "comment", "remark", "description"),
    "category": ("category", "group", "class", "domain"),
}


class PointSetValidationError(ValueError):
    """Raised when a PCPT document fails validation."""


@dataclass(repr=False)
class Point(PyCSAMTObject):
    """One named position with an optional depth window."""

    id: str
    name: str = ""
    x: float | None = None
    y: float | None = None
    longitude: float | None = None
    latitude: float | None = None
    z: float | None = None
    depth_top: float | None = None
    depth_bottom: float | None = None
    kind: str = "poi"
    symbol: str = "circle"
    color: str | None = None
    note: str = ""
    category: str = ""
    properties: dict[str, Any] = field(default_factory=dict)

    def issues(self) -> list[str]:
        out: list[str] = []
        if not str(self.id).strip():
            out.append("point id must be non-empty")
        has_xy = self.x is not None and self.y is not None
        has_ll = self.longitude is not None and self.latitude is not None
        if not has_xy and not has_ll:
            out.append(f"point {self.id!r} has neither x/y nor lon/lat")
        for name in ("x", "y", "longitude", "latitude", "z",
                     "depth_top", "depth_bottom"):
            value = getattr(self, name)
            if value is not None and not _finite(value):
                out.append(f"point {self.id!r}: {name} must be finite")
        if (
            _finite(self.longitude)
            and not -180.0 <= float(self.longitude) <= 180.0
        ):
            out.append(f"point {self.id!r}: longitude out of range")
        if (
            _finite(self.latitude)
            and not -90.0 <= float(self.latitude) <= 90.0
        ):
            out.append(f"point {self.id!r}: latitude out of range")
        if (
            _finite(self.depth_top)
            and _finite(self.depth_bottom)
            and float(self.depth_bottom) <= float(self.depth_top)
        ):
            out.append(
                f"point {self.id!r}: depth_bottom must exceed depth_top"
            )
        return out


@dataclass(repr=False)
class PointSet(PyCSAMTObject, MetadataMixin):
    """A validated flat list of points sharing one CRS."""

    document_id: str
    created_at: str
    created_by: str
    crs: str
    points: list[Point]
    pcpt_version: str = PCPT_VERSION
    title: str = ""
    description: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def issues(self) -> list[str]:
        out: list[str] = []
        if not str(self.document_id).strip():
            out.append("document_id must be non-empty")
        if not str(self.crs).strip():
            out.append("crs must be non-empty")
        try:
            parts = tuple(int(p) for p in self.pcpt_version.split("."))
        except (AttributeError, ValueError):
            parts = ()
        if len(parts) != 3 or parts[0] != 0:
            out.append(f"pcpt_version must be compatible with {PCPT_VERSION}")
        if not self.points:
            out.append("point set must contain at least one point")
        seen: set[str] = set()
        for point in self.points:
            if not isinstance(point, Point):
                out.append("every entry must be a Point")
                continue
            out.extend(point.issues())
            if point.id in seen:
                out.append(f"duplicate point id {point.id!r}")
            seen.add(point.id)
        return out

    def validate(self) -> None:
        problems = self.issues()
        if problems:
            raise PointSetValidationError("; ".join(problems[:6]))

    def lonlat(self) -> dict[str, tuple[float, float]]:
        """Return ``{id: (lat, lon)}`` for every locatable point."""
        result: dict[str, tuple[float, float]] = {}
        transformer = None
        horizontal = str(self.crs or "")
        want_transform = any(
            p.longitude is None or p.latitude is None for p in self.points
        )
        if (
            want_transform
            and horizontal
            and not horizontal.upper().startswith("LOCAL:")
        ):
            try:
                from pyproj import Transformer

                transformer = Transformer.from_crs(
                    horizontal, "EPSG:4326", always_xy=True
                )
            except Exception:  # noqa: BLE001
                transformer = None
        for point in self.points:
            if point.latitude is not None and point.longitude is not None:
                result[point.id] = (
                    float(point.latitude),
                    float(point.longitude),
                )
            elif (
                transformer is not None
                and point.x is not None
                and point.y is not None
            ):
                try:
                    lon, lat = transformer.transform(
                        float(point.x), float(point.y)
                    )
                except Exception:  # noqa: BLE001
                    continue
                if isfinite(lon) and isfinite(lat):
                    result[point.id] = (float(lat), float(lon))
        return result


# ---------------------------------------------------------------------------
# JSON I/O
# ---------------------------------------------------------------------------


def point_set_to_dict(point_set: PointSet) -> dict[str, Any]:
    """Return a canonical JSON-safe dict for *point_set*."""
    point_set.validate()
    return {
        "pcpt_version": point_set.pcpt_version,
        "document_id": point_set.document_id,
        "created_at": point_set.created_at,
        "created_by": point_set.created_by,
        "title": point_set.title,
        "description": point_set.description,
        "crs": point_set.crs,
        "metadata": dict(point_set.metadata),
        "points": [
            {
                key: value
                for key, value in {
                    "id": point.id,
                    "name": point.name,
                    "x": point.x,
                    "y": point.y,
                    "longitude": point.longitude,
                    "latitude": point.latitude,
                    "z": point.z,
                    "depth_top": point.depth_top,
                    "depth_bottom": point.depth_bottom,
                    "kind": point.kind,
                    "symbol": point.symbol,
                    "color": point.color,
                    "note": point.note,
                    "category": point.category,
                    "properties": dict(point.properties),
                }.items()
                if value not in (None, "", {})
                or key in ("id",)
            }
            for point in point_set.points
        ],
    }


def point_set_from_dict(payload: dict[str, Any]) -> PointSet:
    """Build and validate a :class:`PointSet` from a dict."""
    if not isinstance(payload, dict):
        raise PointSetValidationError("PCPT payload must be a JSON object")
    points = [
        Point(
            id=str(item.get("id", "")),
            name=str(item.get("name", "")),
            x=_num(item.get("x")),
            y=_num(item.get("y")),
            longitude=_num(item.get("longitude")),
            latitude=_num(item.get("latitude")),
            z=_num(item.get("z")),
            depth_top=_num(item.get("depth_top")),
            depth_bottom=_num(item.get("depth_bottom")),
            kind=str(item.get("kind", "poi")) or "poi",
            symbol=str(item.get("symbol", "circle")) or "circle",
            color=item.get("color") or None,
            note=str(item.get("note", "")),
            category=str(item.get("category", "")),
            properties=dict(item.get("properties", {}) or {}),
        )
        for item in payload.get("points", [])
    ]
    point_set = PointSet(
        document_id=str(payload.get("document_id", "pcpt:points")),
        created_at=str(payload.get("created_at") or _now()),
        created_by=str(payload.get("created_by", "pycsamt")),
        crs=str(payload.get("crs", "EPSG:4326")),
        points=points,
        pcpt_version=str(payload.get("pcpt_version", PCPT_VERSION)),
        title=str(payload.get("title", "")),
        description=str(payload.get("description", "")),
        metadata=dict(payload.get("metadata", {}) or {}),
    )
    point_set.validate()
    return point_set


def read_points(path: str | Path) -> PointSet:
    """Read a ``.pcpt.json`` file."""
    raw = Path(path).read_text(encoding="utf-8")
    return point_set_from_dict(json.loads(raw))


def write_points(point_set: PointSet, path: str | Path) -> Path:
    """Write *point_set* as canonical ``.pcpt.json``."""
    target = Path(path)
    payload = json.dumps(
        point_set_to_dict(point_set), ensure_ascii=False, indent=2
    )
    target.write_text(payload + "\n", encoding="utf-8")
    return target


# ---------------------------------------------------------------------------
# tabular import
# ---------------------------------------------------------------------------


def points_from_csv(
    path: str | Path,
    *,
    columns: dict[str, str] | None = None,
    crs: str | None = None,
    document_id: str | None = None,
) -> PointSet:
    """Import a points CSV (one point per row)."""
    text = Path(path).read_bytes().decode("utf-8-sig")
    reader = csv.reader(io.StringIO(text))
    try:
        headers = [h.strip() for h in next(reader)]
    except StopIteration as error:
        raise PointSetValidationError("CSV is empty") from error
    rows = [list(values) for values in reader]
    return _points_from_table(
        headers, rows, columns=columns, crs=crs,
        document_id=document_id or f"pcpt:{Path(path).stem}",
    )


def points_from_xlsx(
    source: str | Path | bytes,
    *,
    sheet: str | int | None = None,
    header_row: int | None = None,
    columns: dict[str, str] | None = None,
    crs: str | None = None,
    document_id: str | None = None,
) -> PointSet:
    """Import points from one worksheet of an ``.xlsx`` workbook."""
    from .borehole.xlsxio import sheet_rows

    name, headers, rows = sheet_rows(
        source, sheet=sheet, header_row=header_row
    )
    return _points_from_table(
        headers, rows, columns=columns, crs=crs,
        document_id=document_id or f"pcpt:{name}",
    )


def _points_from_table(
    headers: list[str],
    rows: list[list[str]],
    *,
    columns: dict[str, str] | None,
    crs: str | None,
    document_id: str,
) -> PointSet:
    mapping = _resolve_columns(headers, columns)
    points: list[Point] = []
    for index, values in enumerate(rows, start=2):
        record = dict(zip(headers, values))

        def get(field: str, _record: dict = record) -> Any:
            return _clean(_record.get(mapping.get(field)))

        pid = get("id") or get("name") or f"P{index - 1}"
        points.append(
            Point(
                id=str(pid),
                name=str(get("name") or ""),
                x=_num(get("x")),
                y=_num(get("y")),
                longitude=_num(get("longitude")),
                latitude=_num(get("latitude")),
                z=_num(get("z")),
                depth_top=_num(get("depth_top")),
                depth_bottom=_num(get("depth_bottom")),
                kind=str(get("kind") or "poi"),
                symbol=str(get("symbol") or "circle"),
                color=get("color") or None,
                note=str(get("note") or ""),
                category=str(get("category") or ""),
            )
        )
    points = [p for p in points if p.issues() == [] or _placeable(p)]
    if not points:
        raise PointSetValidationError(
            "no usable points — need id + (x,y) or (lon,lat) per row"
        )
    inferred_crs = crs or (
        "EPSG:4326"
        if any(p.longitude is not None for p in points)
        else "LOCAL:points-grid"
    )
    point_set = PointSet(
        document_id=document_id,
        created_at=_now(),
        created_by="pycsamt points importer",
        crs=inferred_crs,
        points=points,
    )
    point_set.validate()
    return point_set


def _resolve_columns(
    headers: list[str], explicit: dict[str, str] | None
) -> dict[str, str]:
    normalized = {_norm(h): h for h in headers}
    result = dict(explicit or {})
    for field_name, aliases in _ALIASES.items():
        if field_name in result:
            continue
        for alias in (field_name, *aliases):
            match = normalized.get(_norm(alias))
            if match is not None:
                result[field_name] = match
                break
    return result


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _finite(value: Any) -> bool:
    try:
        return isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _num(value: Any) -> float | None:
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if isfinite(number) else None


def _clean(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, str) and value.strip().casefold() in _MISSING:
        return None
    return value


def _placeable(point: Point) -> bool:
    return (
        (point.x is not None and point.y is not None)
        or (point.longitude is not None and point.latitude is not None)
    ) and bool(str(point.id).strip())


def _norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().casefold()).strip("_")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
