# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Station editing for the Edit ▸ Stations dialogs (Qt-free).

The Station Editor shows one table per concern -- coordinates, names,
lines, header metadata -- built by :func:`station_table`.  Whatever the
user changes is reduced by :func:`table_changes` to the cells that really
changed, and applied in one validated, audited batch through
:func:`pycsamt.site.metadata.update_metadata_all` (staged on copies;
coordinate ranges and duplicate names are checked before anything is
committed).  Lines are not stored in the sites but in the app's station
table, so they travel separately (``lines`` mappings).

Projected coordinates (easting / northing in any EPSG) are converted with
pyproj, so a table exported from a GPS or a survey report can be pasted as
it is.
"""

from __future__ import annotations

import io
import math
import re
from typing import Any

import pandas as pd

__all__ = [
    "HEAD_FIELDS",
    "apply_changes",
    "check_names",
    "detect_lines",
    "parse_pasted",
    "project_table",
    "rename_preview",
    "station_table",
    "table_changes",
    "unproject_table",
]

# EDI HEAD fields offered for editing (label, key)
HEAD_FIELDS = (
    ("Acquired by", "acqby"), ("Acq. date", "acqdate"),
    ("End date", "enddate"), ("Project", "project"), ("Survey", "survey"),
    ("Prospect", "prospect"), ("Location", "loc"), ("Country", "country"),
    ("State", "state"), ("County", "county"), ("Datum", "datum"),
    ("Declination", "declination"),
)
COORD_COLUMNS = ("lat", "lon", "elev")


def _num(v) -> float:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return math.nan
    return f


def _head_value(site, key: str) -> str:
    try:
        head = site.edi.get_section("head")
        v = getattr(head, key, None)
    except Exception:
        return ""
    if v is None:
        return ""
    s = str(v).strip().strip('"')
    return "" if s.lower() in ("none", "nan") else s


def station_table(sites, lines: dict | None = None) -> pd.DataFrame:
    """One row per station: name, line, lat, lon, elev and HEAD fields."""
    lines = lines or {}
    rows = []
    for site in sites:
        name = str(getattr(site, "name", ""))
        try:
            lat, lon, elev = site.coords
        except Exception:
            lat = lon = elev = math.nan
        row = {"station": name, "line": str(lines.get(name, "")),
               "lat": _num(lat), "lon": _num(lon), "elev": _num(elev)}
        for _label, key in HEAD_FIELDS:
            row[key] = _head_value(site, key)
        rows.append(row)
    return pd.DataFrame(rows)


def _same(a, b) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        fa, fb = _num(a), _num(b)
        if math.isnan(fa) and math.isnan(fb):
            return True
        return (not math.isnan(fa) and not math.isnan(fb)
                and abs(fa - fb) <= 1e-9 * max(1.0, abs(fa)))
    return str(a) == str(b)


def table_changes(before: pd.DataFrame, after: pd.DataFrame) -> dict:
    """``{station: {column: new_value}}`` for every changed cell.

    Rows are matched by position (the station column itself may be the
    edited one: a rename).  Column ``station`` changes become ``name``.
    """
    out: dict[str, dict] = {}
    for i in range(min(len(before), len(after))):
        old = before.iloc[i]
        new = after.iloc[i]
        diff = {}
        for col in after.columns:
            if col not in before.columns:
                continue
            if not _same(old[col], new[col]):
                diff["name" if col == "station" else col] = new[col]
        if diff:
            out[str(old["station"])] = diff
    return out


def apply_changes(sites, changes: dict):
    """Apply :func:`table_changes` output to *sites* (validated batch).

    Returns ``(new_sites, n_changed, lines)``: the updated sites (the input
    when only lines changed), how many stations changed, and the line
    changes ``{final_station_name: line}``.
    """
    from pycsamt.site.metadata import update_metadata_all

    rows, lines = [], {}
    for station, diff in changes.items():
        final = str(diff.get("name", station))
        if "line" in diff:
            lines[final] = str(diff["line"]).strip()
        spec = {"station": station}
        for key, value in diff.items():
            if key == "line":
                continue
            if key in COORD_COLUMNS:
                v = _num(value)
                if math.isnan(v):
                    raise ValueError(f"{station}: {key} must be a number "
                                     f"(got {value!r})")
                spec[key] = v
            elif key == "name":
                spec["name"] = str(value).strip()
            else:
                spec[f"head.{key}"] = str(value)
        if len(spec) > 1:
            rows.append(spec)
    if not rows:
        return sites, len(lines), lines
    frame = pd.DataFrame(rows)
    new = update_metadata_all(sites, frame, missing="raise",
                              on_error="raise")
    return new, len(changes), lines


# ── names ────────────────────────────────────────────────────────────────
def rename_preview(names, *, prefix: str = "", suffix: str = "",
                   find: str = "", replace: str = "", regex: bool = False,
                   case: str = "keep", pad: int = 0) -> list[str]:
    """New names for *names* under the rename rules.

    *pad* zero-pads the last number in each name to that many digits
    (``L1-3`` -> ``L1-003`` with ``pad=3``).  *case* is ``keep``,
    ``upper`` or ``lower``.
    """
    out = []
    for n in names:
        s = str(n)
        if find:
            s = re.sub(find, replace, s) if regex else s.replace(find,
                                                                  replace)
        if pad and pad > 0:
            m = list(re.finditer(r"\d+", s))
            if m:
                last = m[-1]
                s = s[:last.start()] + last.group().zfill(pad) + s[last.end():]
        if case == "upper":
            s = s.upper()
        elif case == "lower":
            s = s.lower()
        out.append(f"{prefix}{s}{suffix}")
    return out


def check_names(names) -> list[str]:
    """Problems with a list of final names (empty, duplicated)."""
    problems = []
    seen: dict[str, int] = {}
    for n in names:
        key = str(n).strip()
        if not key:
            problems.append("empty name")
        seen[key.casefold()] = seen.get(key.casefold(), 0) + 1
    dup = sorted(k for k, c in seen.items() if c > 1 and k)
    if dup:
        problems.append("duplicate names: " + ", ".join(dup[:6])
                        + (" …" if len(dup) > 6 else ""))
    return problems


# ── lines ────────────────────────────────────────────────────────────────
def detect_lines(names) -> dict[str, str]:
    """``{station: line}`` from the station names (``L18-001`` -> L18)."""
    from pycsamt.site.lines import detect_lines_from_station_ids

    groups = detect_lines_from_station_ids(list(map(str, names)))
    return {st: str(line) for line, members in groups.items()
            for st in members}


# ── pasting & projections ────────────────────────────────────────────────
def parse_pasted(text: str) -> list[list[str]]:
    """Cells of a block copied from a spreadsheet (tabs) or a CSV."""
    text = text.strip("\r\n")
    if not text:
        return []
    sep = "\t" if "\t" in text else ("," if "," in text else None)
    rows = []
    for line in io.StringIO(text):
        line = line.rstrip("\r\n")
        cells = line.split(sep) if sep else line.split()
        rows.append([c.strip() for c in cells])
    return rows


def _transformer(epsg: int, *, inverse: bool):
    try:
        from pyproj import Transformer
    except ImportError as exc:  # pragma: no cover - optional
        raise ImportError("projected coordinates need pyproj") from exc
    if inverse:
        return Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326",
                                    always_xy=True)
    return Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)


def project_table(table: pd.DataFrame, epsg: int) -> pd.DataFrame:
    """Add ``easting``/``northing`` in *epsg* from ``lat``/``lon``."""
    t = _transformer(int(epsg), inverse=False)
    out = table.copy()
    e, n = t.transform(out["lon"].astype(float).to_numpy(),
                       out["lat"].astype(float).to_numpy())
    out["easting"], out["northing"] = e, n
    return out


def unproject_table(table: pd.DataFrame, epsg: int) -> pd.DataFrame:
    """``lat``/``lon`` from ``easting``/``northing`` in *epsg*."""
    t = _transformer(int(epsg), inverse=True)
    out = table.copy()
    lon, lat = t.transform(out["easting"].astype(float).to_numpy(),
                           out["northing"].astype(float).to_numpy())
    out["lat"], out["lon"] = lat, lon
    return out


def utm_epsg(lat: float, lon: float) -> int:
    """The WGS84 UTM zone EPSG code for a point."""
    zone = int((lon + 180.0) // 6.0) + 1
    return (32600 if lat >= 0 else 32700) + max(1, min(zone, 60))


def carry_lines(old_lines: dict, changes: dict) -> dict:
    """The station -> line map after *changes*: renamed stations keep
    their line under the new name, edited lines take the new value."""
    out = dict(old_lines)
    for station, diff in changes.items():
        final = str(diff.get("name", station)).strip()
        line = out.pop(station, None) if final != station else \
            out.get(station)
        if "line" in diff:
            line = str(diff["line"]).strip()
        if line is not None:
            out[final] = line
    return out


def _sites_names(sites) -> list[str]:
    return [str(getattr(s, "name", "")) for s in sites]


def lines_from_frame(frame: Any) -> dict[str, str]:
    """``{ID: Line}`` from the app's station table."""
    if frame is None or "ID" not in frame or "Line" not in frame:
        return {}
    return dict(zip(frame["ID"].astype(str), frame["Line"].astype(str)))
