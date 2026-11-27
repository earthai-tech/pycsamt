# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Station coordinate tables: reading, CRS checks and EDI embedding.

Field crews often deliver station positions as a bare GPS table
(no header, a UTF-8 byte-order mark, trailing empty columns, a
projected CRS that is not stated) while the EDI files carry no
location at all. This module covers that path:

1. :func:`read_station_table` reads such tables into a standard
   frame (``station``, ``northing``, ``easting``, ``elev`` or
   ``lat``/``lon``), inferring columns only when unambiguous.
2. :func:`compare_crs` shows where a table lands under candidate
   projected CRS, so the right EPSG can be chosen knowingly (UTM
   and Gauss-Kruger on the same central meridian differ by ~1 km).
3. :func:`to_wgs84` converts the table to WGS84 latitude/longitude.
4. :func:`embed_station_coordinates` writes the positions into EDI
   headers (``>HEAD`` and ``>=DEFINEMEAS`` reference), matching
   stations by **name** (file stem or DATAID, extension and case
   insensitive) rather than by file order.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "read_station_table",
    "to_wgs84",
    "compare_crs",
    "station_key",
    "embed_station_coordinates",
]

_PROJ_COLS = ("station", "northing", "easting", "elev")


def station_key(name: Any) -> str:
    """Normalised station key used for name matching.

    Upper case, surrounding whitespace removed, and a trailing file
    or processing extension (``.edi``, ``.atimpd`` ...) dropped.

    Examples
    --------
    >>> from pycsamt.gis.stations import station_key
    >>> station_key(" l1-450.atimpd "), station_key("L1-450")
    ('L1-450', 'L1-450')
    """
    s = str(name).strip()
    s = re.sub(r"\.[A-Za-z][A-Za-z0-9]{0,7}$", "", s)
    return s.upper()


def read_station_table(
    path: str | Path,
    *,
    columns: Sequence[str] | None = None,
    header: str | bool = "auto",
    delimiter: str | None = None,
    sheet_name: int | str = 0,
) -> pd.DataFrame:
    """Read a station coordinate table (CSV, TXT, XLS/XLSX).

    Parameters
    ----------
    path : str or Path
        Table file.
    columns : sequence of str, optional
        Names for the columns, in file order (e.g.
        ``("station", "northing", "easting", "elev")``). Overrides
        any header; columns beyond the given names are dropped.
    header : {"auto", True, False}, default "auto"
        Whether the first row is a header. ``"auto"`` treats it as
        a header when a column that is numeric below row 1 holds
        text in row 1.
    delimiter : str, optional
        Field separator for text files (default: comma, or any
        whitespace when the file has no commas).
    sheet_name : int or str, default 0
        Sheet for Excel files.

    Returns
    -------
    pandas.DataFrame
        One row per station; ``station`` as stripped strings and
        coordinate columns as float. Fully empty columns and rows are
        dropped. Without *columns* and without a header, a table of
        one text column plus three numeric ones is labelled
        ``station, northing, easting, elev`` when unambiguous: the
        two columns with magnitudes above 1e4 are taken as the
        projected pair (the larger median is ``northing``, valid in
        the northern hemisphere with false eastings below 1e7) and
        the remaining one as ``elev``.

    Raises
    ------
    ValueError
        If the layout cannot be inferred unambiguously (pass
        *columns*), if *columns* has the wrong length, or if a
        ``station`` column is missing or duplicated.

    Examples
    --------
    >>> import pathlib, tempfile
    >>> from pycsamt.gis.stations import read_station_table
    >>> p = pathlib.Path(tempfile.mkdtemp()) / "l1.csv"
    >>> _ = p.write_text("\\ufeffL1-0,2830561.5,792692.4,373.7,,\\n"
    ...                  "L1-30,2830543.8,792716.7,366.8,,\\n",
    ...                  encoding="utf-8")
    >>> read_station_table(p).columns.tolist()
    ['station', 'northing', 'easting', 'elev']
    """
    p = Path(path)
    if p.suffix.lower() in (".xls", ".xlsx"):
        raw = pd.read_excel(p, header=None, sheet_name=sheet_name,
                            dtype=object)
    else:
        text = p.read_text(encoding="utf-8-sig", errors="replace")
        sep = delimiter or ("," if "," in text else r"\s+")
        raw = pd.read_csv(p, header=None, sep=sep, dtype=str,
                          encoding="utf-8-sig", skip_blank_lines=True,
                          engine="python")
    raw = raw.replace(r"^\s*$", np.nan, regex=True)
    raw = raw.dropna(axis=1, how="all").dropna(axis=0, how="all")
    raw = raw.reset_index(drop=True)
    if raw.empty:
        raise ValueError(f"{p}: no data rows.")

    def _num(col: pd.Series) -> pd.Series:
        return pd.to_numeric(col, errors="coerce")

    if header == "auto":
        body = raw.iloc[1:]
        has_header = bool(len(body)) and any(
            _num(body[c]).notna().all() and pd.isna(_num(
                raw[c].iloc[:1])).all()
            for c in raw.columns)
    else:
        has_header = bool(header)
    if has_header:
        names = [str(v).strip() for v in raw.iloc[0]]
        raw = raw.iloc[1:].reset_index(drop=True)
    else:
        names = None
    if columns is not None:
        if len(columns) > raw.shape[1]:
            raise ValueError(
                f"{len(columns)} column names for {raw.shape[1]} "
                "columns.")
        raw = raw.iloc[:, : len(columns)]
        names = list(columns)
    elif names is None:
        names = _infer_projected_layout(raw, _num)
        raw = raw.iloc[:, : len(names)]
    raw.columns = names
    if list(raw.columns).count("station") != 1:
        raise ValueError(
            "exactly one 'station' column is required (pass "
            "columns=...).")
    out = raw.copy()
    out["station"] = out["station"].astype(str).str.strip()
    for c in out.columns:
        if c != "station":
            out[c] = _num(out[c]).astype(float)
    return out.reset_index(drop=True)


def _infer_projected_layout(raw: pd.DataFrame, num) -> list[str]:
    numeric = [c for c in raw.columns if num(raw[c]).notna().all()]
    text = [c for c in raw.columns if c not in numeric]
    if len(text) != 1 or text[0] != raw.columns[0] or len(numeric) != 3:
        raise ValueError(
            "cannot infer the table layout (expected one station "
            "column followed by three numeric columns); pass "
            "columns=(...).")
    med = {c: float(np.nanmedian(np.abs(num(raw[c])))) for c in numeric}
    big = [c for c in numeric if med[c] > 1e4]
    if len(big) != 2:
        raise ValueError(
            "cannot tell the projected coordinate pair from the "
            "elevation; pass columns=(...).")
    north, east = sorted(big, key=lambda c: -med[c])
    elev = next(c for c in numeric if c not in big)
    order = {text[0]: "station", north: "northing", east: "easting",
             elev: "elev"}
    return [order[c] for c in raw.columns[:4]]


def to_wgs84(
    table: pd.DataFrame,
    epsg: int,
    *,
    x: str = "easting",
    y: str = "northing",
) -> pd.DataFrame:
    """Add WGS84 ``lat``/``lon`` columns to a projected table.

    Parameters
    ----------
    table : pandas.DataFrame
        Station table with projected coordinates.
    epsg : int
        EPSG code of the projected CRS (e.g. 4508 for CGCS2000 /
        Gauss-Kruger CM 111E, 32649 for WGS 84 / UTM 49N).
    x, y : str
        Easting and northing column names.

    Returns
    -------
    pandas.DataFrame
        A copy with ``lat`` and ``lon`` (decimal degrees).

    Raises
    ------
    ValueError
        If *epsg* is unknown to PROJ.
    """
    from pyproj import CRS, Transformer
    from pyproj.exceptions import CRSError

    try:
        CRS.from_epsg(int(epsg))
    except CRSError as exc:
        raise ValueError(f"EPSG:{epsg} is not a valid CRS.") from exc
    tr = Transformer.from_crs(int(epsg), 4326, always_xy=True)
    lon, lat = tr.transform(table[x].to_numpy(float),
                            table[y].to_numpy(float))
    out = table.copy()
    out["lat"], out["lon"] = lat, lon
    return out


def compare_crs(
    table: pd.DataFrame,
    candidates: Sequence[int],
    *,
    x: str = "easting",
    y: str = "northing",
) -> pd.DataFrame:
    """Where the stations land under each candidate CRS.

    Parameters
    ----------
    table : pandas.DataFrame
        Projected station table.
    candidates : sequence of int
        EPSG codes to compare; the first is the reference.
    x, y : str
        Easting and northing column names.

    Returns
    -------
    pandas.DataFrame
        One row per EPSG: CRS name, centroid ``lat``/``lon`` and the
        horizontal ``offset_m`` of the centroid from the reference
        candidate's (great-circle, m).

    Notes
    -----
    Projected coordinates do not identify their CRS. The usual
    ambiguity is UTM (scale factor 0.9996) versus Gauss-Kruger
    (scale 1) on the same central meridian, which differ by up to a
    kilometre; datum differences (CGCS2000, Xian 1980, Beijing 1954)
    are tens of metres. Confirm the CRS with the survey crew.

    Examples
    --------
    >>> import pandas as pd
    >>> from pycsamt.gis.stations import compare_crs
    >>> t = pd.DataFrame({"station": ["A"], "northing": [2830561.5],
    ...                   "easting": [792692.4]})
    >>> c = compare_crs(t, [32649, 4508])
    >>> round(float(c.offset_m.iloc[1]), -2)
    1100.0
    """
    from pyproj import CRS, Geod

    geod = Geod(ellps="WGS84")
    rows = []
    for epsg in candidates:
        ll = to_wgs84(table, epsg, x=x, y=y)
        rows.append({"epsg": int(epsg),
                     "name": CRS.from_epsg(int(epsg)).name,
                     "lat": float(ll["lat"].mean()),
                     "lon": float(ll["lon"].mean())})
    out = pd.DataFrame(rows)
    ref = out.iloc[0]
    out["offset_m"] = [
        float(geod.inv(ref.lon, ref.lat, r.lon, r.lat)[2])
        for r in out.itertuples()]
    return out


def _edi_keys(ed: Any, i: int) -> list[str]:
    """Candidate station keys of one EDI/Site (file stem first)."""
    edi = getattr(ed, "edi", ed)
    keys = []
    path = getattr(edi, "path", None)
    if path is not None:
        keys.append(station_key(Path(path).stem))
    for attr in ("station", "name"):
        v = getattr(edi, attr, None)
        if isinstance(v, str) and v:
            keys.append(station_key(v))
    head = edi.get_section("head") if hasattr(edi, "get_section") \
        else None
    v = getattr(head, "dataid", None)
    if isinstance(v, str) and v:
        keys.append(station_key(v))
    return list(dict.fromkeys(k for k in keys if k))


def embed_station_coordinates(
    sites: Any,
    table: pd.DataFrame,
    *,
    match: str = "name",
    set_reference: bool = True,
    note: str | None = None,
    require_all: bool = True,
) -> pd.DataFrame:
    """Write station positions into EDI headers.

    Parameters
    ----------
    sites : Sites or sequence of EDIFile
        Stations to update **in place**.
    table : pandas.DataFrame
        Must hold ``station``, ``lat``, ``lon`` and optionally
        ``elev`` (see :func:`read_station_table`, :func:`to_wgs84`).
    match : {"name", "order"}, default "name"
        ``"name"``: match each EDI by its file stem, then station /
        DATAID, against ``table.station`` with :func:`station_key`
        (extension- and case-insensitive). ``"order"``: row *i* goes
        to EDI *i* (only safe when both are in the same order).
    set_reference : bool, default True
        Also set ``>=DEFINEMEAS`` ``REFLAT/REFLONG/REFELEV``, which
        some programs read instead of ``>HEAD``.
    note : str, optional
        Free-text line added to ``>INFO`` (e.g. the source CRS). It
        is written without ``=`` signs being parsed as fields.
    require_all : bool, default True
        Raise when an EDI has no row in *table*.

    Returns
    -------
    pandas.DataFrame
        One row per EDI: ``key``, matched ``station``, ``lat``,
        ``lon``, ``elev`` and ``matched``.

    Raises
    ------
    ValueError
        If required columns are missing, *table* has duplicate
        station keys, *match* is invalid, the counts differ with
        ``match="order"``, or an EDI is unmatched while
        *require_all* is True.
    """
    from pycsamt.emtools._core import _iter_items
    from pycsamt.site.utils import set_coords

    missing = {"station", "lat", "lon"} - set(table.columns)
    if missing:
        raise ValueError(f"table lacks columns {sorted(missing)}.")
    if match not in ("name", "order"):
        raise ValueError("match must be 'name' or 'order'.")
    keys = [station_key(s) for s in table["station"]]
    if len(set(keys)) != len(keys):
        raise ValueError("table has duplicate station names.")
    rows = table.assign(_key=keys).set_index("_key")
    items = list(_iter_items(sites))
    if match == "order" and len(items) != len(table):
        raise ValueError(
            f"{len(items)} stations but {len(table)} table rows.")
    report = []
    for i, ed in enumerate(items):
        edi = getattr(ed, "edi", ed)
        cand = _edi_keys(ed, i)
        if match == "order":
            row = table.iloc[i]
        else:
            hit = next((k for k in cand if k in rows.index), None)
            if hit is None:
                if require_all:
                    raise ValueError(
                        f"no table row for station {cand or [i]}.")
                report.append({"key": cand[0] if cand else str(i),
                               "station": None, "lat": np.nan,
                               "lon": np.nan, "elev": np.nan,
                               "matched": False})
                continue
            row = rows.loc[hit]
        elev = float(row["elev"]) if "elev" in row.index and \
            np.isfinite(row["elev"]) else None
        set_coords(edi, lat=float(row["lat"]), lon=float(row["lon"]),
                   elev=elev, inplace=True)
        if set_reference and hasattr(edi, "get_section"):
            dm = edi.get_section("definemeas")
            if dm is not None:
                dm.reflat, dm.reflong = float(row["lat"]), float(
                    row["lon"])
                if elev is not None:
                    dm.refelev = elev
        if note and hasattr(edi, "get_section"):
            info = edi.get_section("info")
            if info is not None:
                info.info_text = [f"    {note}"] + list(
                    getattr(info, "info_text", []) or [])
        report.append({"key": cand[0] if cand else str(i),
                       "station": str(row["station"]),
                       "lat": float(row["lat"]),
                       "lon": float(row["lon"]),
                       "elev": np.nan if elev is None else elev,
                       "matched": True})
    return pd.DataFrame(report)
