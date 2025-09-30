# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence, Tuple, Iterable
import math
import re
import copy

import numpy as np

from ..constants import _EARTH_R
from ..core.config import get_config
from ..gis.utils import (
    assert_lat_value,
    assert_lon_value,
    assert_elevation_value,
)
from ..gis.config import HAS_GDAL  # noqa: F401
from ..seg.heads import Head  # type: ignore

# why not use the _get_head, _ensure_head  from pycsamt/seg/utils.py ? 

from .utils import _get_head, _ensure_head 

__all__ = [
    "Coord",
    "parse_lat",
    "parse_lon",
    "parse_elev",
    "ensure_head_coords",
    "apply_topography",
    "project",
    "distance",
    "bearing",
    "chainage_along",
]


@dataclass
class Coord:
    lat: float
    lon: float
    elev: float = 0.0


_HEMI_LAT = {"N": 1.0, "S": -1.0}
_HEMI_LON = {"E": 1.0, "W": -1.0}

_DMS_RX = re.compile(
    r"""
    ^\s*
    (?P<deg>[-+]?\d+(?:\.\d+)?)
    (?:[°:\s]+(?P<min>\d+(?:\.\d+)?))?
    (?:[':\s]+(?P<sec>\d+(?:\.\d+)?))?
    \s*(?P<hemi>[NSEW])?
    \s*$
    """,
    re.VERBOSE | re.IGNORECASE,
)


def _empty_value() -> float:
    try:
        return float(get_config().empty)
    except Exception:  # pragma: no cover
        return float("nan")


def _parse_angle(x: Any, hemi_map: dict[str, float]) -> float:
    if x is None:
        return math.nan
    if isinstance(x, (int, float, np.floating)):
        return float(x)
    s = str(x).strip()
    m = _DMS_RX.match(s)
    if m:
        d = float(m.group("deg"))
        mi = float(m.group("min")) if m.group("min") else 0.0
        se = float(m.group("sec")) if m.group("sec") else 0.0
        sign = 1.0
        h = m.group("hemi")
        if h:
            sign = hemi_map.get(h.upper(), 1.0)
        elif d < 0:
            sign = -1.0
        val = sign * (abs(d) + mi / 60.0 + se / 3600.0)
        return float(val)
    hemi = s[-1].upper() if s and s[-1].isalpha() else ""
    try:
        v = float(s[:-1]) if hemi else float(s)
    except Exception:
        return math.nan
    if hemi:
        v *= hemi_map.get(hemi, 1.0)
    return float(v)


def parse_lat(x: Any) -> float:
    return float(_parse_angle(x, _HEMI_LAT))


def parse_lon(x: Any) -> float:
    return float(_parse_angle(x, _HEMI_LON))


def parse_elev(x: Any) -> float:
    try:
        return float(x)
    except Exception:
        return math.nan

# why not use the _get_head, _ensure_head  from pycsamt/seg/utils.py ? 
# to apply DRY? 
def _get_head(ed: Any) -> Any:
    try:
        h = ed.get_section("head")  # type: ignore
        if h is not None:
            return h
    except Exception:
        pass
    h = getattr(ed, "Head", None)
    if h is not None:
        return h
    return None


def _ensure_head(ed: Any) -> Any:
    h = _get_head(ed)
    if h is not None:
        return h
    try:
        h = Head()
    except Exception:
        class _H:
            pass
        h = _H()
    try:
        ed.set_section("head", h)  # type: ignore
        return h
    except Exception:
        pass
    try:
        setattr(ed, "Head", h)
    except Exception:
        pass
    return h


def _get_station(ed: Any) -> str:
    keys = ("station", "dataid", "sitename", "name")
    h = _get_head(ed)
    if h is not None:
        for k in keys:
            v = getattr(h, k, None)
            if v:
                return str(v)
    for k in keys:
        v2 = getattr(ed, k, None)
        if v2:
            return str(v2)
    return ""


def ensure_head_coords(
    ed: Any,
    *,
    lat: Any | None = None,
    lon: Any | None = None,
    elev: Any | None = None,
    empty: float | None = None,
) -> Any:
    emp = _empty_value() if empty is None else float(empty)
    h = _ensure_head(ed)

    la0 = getattr(h, "lat", None)
    lo0 = getattr(h, "lon", getattr(h, "long", None))
    ev0 = getattr(h, "elev", getattr(h, "elevation", None))

    la = parse_lat(lat) if lat is not None else parse_lat(la0)
    lo = parse_lon(lon) if lon is not None else parse_lon(lo0)
    ev = parse_elev(elev) if elev is not None else parse_elev(ev0)

    la = assert_lat_value(la) if math.isfinite(la) else emp
    lo = assert_lon_value(lo) if math.isfinite(lo) else emp
    ev = assert_elevation_value(ev) if math.isfinite(ev) else emp

    try:
        setattr(h, "lat", float(la))
        setattr(h, "lon", float(lo))
        setattr(h, "elev", float(ev))
    except Exception:
        pass
    try:
        setattr(h, "long", float(lo))
    except Exception:
        pass
    return h


def _lower_cols(df: Any) -> Iterable[str]:
    try:
        return [str(c).strip() for c in df.columns]
    except Exception:
        return []


def _pick_col(df: Any, names: Sequence[str]) -> str | None:
    cols = _lower_cols(df)
    for n in names:
        for c in cols:
            if c.lower() == n.lower():
                return c
    return None


def _match_row(df: Any, sid: Any) -> Any:
    if df is None:
        return None
    try:
        if len(df) == 0:
            return None
    except Exception:
        return None
    idc = _pick_col(
        df, ("station", "site", "dataid", "id", "name"),
    )
    if not idc:
        return None
    try:
        row = df[df[idc] == sid]
        if getattr(row, "empty", True):
            row = df[df[idc] == str(sid)]
        if getattr(row, "empty", True):
            return None
        return row.iloc[0]
    except Exception:
        return None


def _set_coords_from_row(
    ed: Any,
    row: Any,
    *,
    empty: float,
) -> None:
    if row is None:
        return
    idx = row.index if hasattr(row, "index") else row
    latc = _pick_col(idx, ("latitude", "lat", "LAT"))
    lonc = _pick_col(idx, ("longitude", "lon", "long", "LON", "LONG"))
    elvc = _pick_col(idx, ("elevation", "elev", "alt", "ALT"))
    la = parse_lat(row[latc]) if latc else math.nan  # type: ignore
    lo = parse_lon(row[lonc]) if lonc else math.nan  # type: ignore
    ev = parse_elev(row[elvc]) if elvc else math.nan  # type: ignore
    ensure_head_coords(ed, lat=la, lon=lo, elev=ev, empty=empty)


def apply_topography(
    ed_or_sites: Any,
    frame: Any,
    *,
    empty: float | None = None,
    inplace: bool = True,
) -> Any:
    emp = _empty_value() if empty is None else float(empty)

    def _apply_one(ed: Any) -> Any:
        sid = _get_station(ed)
        row = _match_row(frame, sid)
        if row is not None:
            _ensure_head(ed)
            _set_coords_from_row(ed, row, empty=emp)
        return ed

    try:
        items = getattr(ed_or_sites, "_items", None)
        if items is not None:
            if not inplace:
                new = copy.deepcopy(ed_or_sites)
                it = getattr(new, "_items", [])
                for s in it:
                    _apply_one(getattr(s, "edi", s))
                return new
            for s in items:
                _apply_one(getattr(s, "edi", s))
            return ed_or_sites
    except Exception:
        pass

    try:
        if isinstance(ed_or_sites, (list, tuple)):
            vec = (
                [copy.deepcopy(x) for x in ed_or_sites]
                if not inplace else list(ed_or_sites)
            )
            for ed in vec:
                _apply_one(ed)
            return vec
    except Exception:
        pass

    return _apply_one(ed_or_sites)


def project(
    pts: Sequence[Tuple[float, float]] | Tuple[float, float],
    *,
    crs_from: Any,
    crs_to: Any,
) -> Tuple[np.ndarray, np.ndarray]:
    try:
        from pyproj import Transformer  # type: ignore
        use_pyproj = True
    except Exception:
        use_pyproj = False

    if not use_pyproj and not HAS_GDAL:
        raise RuntimeError(
            "project() requires pyproj or GDAL."
        )

    if isinstance(pts, tuple) and len(pts) == 2:
        xs, ys = [pts[0]], [pts[1]]
    else:
        xs, ys = zip(*pts)  # type: ignore[misc]

    xs_a = np.asarray(xs, float)
    ys_a = np.asarray(ys, float)

    if use_pyproj:
        tr = Transformer.from_crs(
            crs_from, crs_to, always_xy=True,
        )
        X, Y = tr.transform(xs_a, ys_a)
        return np.asarray(X, float), np.asarray(Y, float)

    # GDAL fallback
    try:
        from osgeo import osr  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError("GDAL not available.") from e

    def _to_srs(crs: Any) -> Any:
        srs = osr.SpatialReference()
        if isinstance(crs, int):
            srs.ImportFromEPSG(int(crs))
            return srs
        if isinstance(crs, str):
            cs = crs.strip()
            if cs.lower().startswith("epsg:"):
                code = int(cs.split(":")[1])
                srs.ImportFromEPSG(code)
            else:
                ok = srs.ImportFromWkt(cs)
                if ok != 0:
                    srs = osr.SpatialReference()
                    srs.ImportFromProj4(cs)
            return srs
        if hasattr(crs, "ExportToWkt"):
            return crs
        raise TypeError("Unsupported CRS spec.")

    s_from = _to_srs(crs_from)
    s_to = _to_srs(crs_to)

    try:
        s_from.SetAxisMappingStrategy(
            osr.OAMS_TRADITIONAL_GIS_ORDER
        )
        s_to.SetAxisMappingStrategy(
            osr.OAMS_TRADITIONAL_GIS_ORDER
        )
    except Exception:
        pass

    ct = osr.CoordinateTransformation(s_from, s_to)
    out_x: list[float] = []
    out_y: list[float] = []
    for x, y in zip(xs_a.tolist(), ys_a.tolist()):
        try:
            X, Y, _ = ct.TransformPoint(float(x), float(y))
        except Exception:
            X, Y = float("nan"), float("nan")
        out_x.append(float(X))
        out_y.append(float(Y))
    return np.asarray(out_x, float), np.asarray(out_y, float)


def _rad(x: float) -> float:
    return math.radians(float(x))

# my question is that , is it possible to compute the distance between 
# degree when EARTH_R is in metters? or we can use 1deg=111km? why not convert in utm first then 
# compute the distance ? or what is the best approach or both approach 
# is possible? if both is possible then give a kwarg parameter mode, 
# to do that ... 
# similar for bearing too . 

def distance(
    a: Coord | Tuple[float, float],
    b: Coord | Tuple[float, float],
) -> float:
    la1, lo1 = (a.lat, a.lon) if isinstance(a, Coord) else a
    la2, lo2 = (b.lat, b.lon) if isinstance(b, Coord) else b
    la1 = float(assert_lat_value(la1))
    lo1 = float(assert_lon_value(lo1))
    la2 = float(assert_lat_value(la2))
    lo2 = float(assert_lon_value(lo2))
    d1 = _rad(la2 - la1) # 
    d2 = _rad(lo2 - lo1)
    A = (
        math.sin(d1 / 2) ** 2
        + math.cos(_rad(la1)) * math.cos(_rad(la2))
        * math.sin(d2 / 2) ** 2
    )
    return 2.0 * _EARTH_R * math.asin(min(1.0, math.sqrt(A)))


def bearing(
    a: Coord | Tuple[float, float],
    b: Coord | Tuple[float, float],
) -> float:
    la1, lo1 = (a.lat, a.lon) if isinstance(a, Coord) else a
    la2, lo2 = (b.lat, b.lon) if isinstance(b, Coord) else b
    la1 = _rad(assert_lat_value(la1))  # type: ignore
    la2 = _rad(assert_lat_value(la2))  # type: ignore
    dlon = _rad(assert_lon_value(lo2) - assert_lon_value(lo1))
    y = math.sin(dlon) * math.cos(la2)
    x = (
        math.cos(la1) * math.sin(la2)
        - math.sin(la1) * math.cos(la2) * math.cos(dlon)
    )
    brg = math.degrees(math.atan2(y, x))
    return (brg + 360.0) % 360.0


def chainage_along(
    origin: Coord | Tuple[float, float],
    azimuth: float,
    pts: Sequence[Coord | Tuple[float, float]]
    | Coord
    | Tuple[float, float],
) -> np.ndarray | float:
    if isinstance(origin, Coord):
        la0, lo0 = origin.lat, origin.lon
    else:
        la0, lo0 = origin
    la0 = float(assert_lat_value(la0))
    lo0 = float(assert_lon_value(lo0))
    az = math.radians(float(azimuth))
    mperdeg = 111_000.0

    def _one(p: Coord | Tuple[float, float]) -> float:
        if isinstance(p, Coord):
            la, lo = p.lat, p.lon
        else:
            la, lo = p
        la = float(assert_lat_value(la))
        lo = float(assert_lon_value(lo))
        dy = (la - la0) * mperdeg
        dx = (lo - lo0) * mperdeg * math.cos(_rad(la0))
        return dx * math.cos(az) + dy * math.sin(az)

    if isinstance(pts, (tuple, Coord)):
        return float(_one(pts))
    out = [_one(p) for p in pts]  # type: ignore[arg-type]
    return np.asarray(out, float)
