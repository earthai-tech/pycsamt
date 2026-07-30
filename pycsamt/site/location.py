# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import copy
import math
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..constants import _EARTH_R
from ..gis.config import HAS_GDAL
from ..gis.utils import (
    assert_elevation_value,
    assert_lat_value,
    assert_lon_value,
)
from .utils import _ensure_head, _get_head

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


@dataclass
class Coord:
    r"""
    Lightweight geographic coordinate container.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees. North is positive.
    lon : float
        Longitude in decimal degrees. East is positive.
    elev : float, optional
        Elevation in meters. Defaults to 0.0.

    Notes
    -----
    This class is a convenience container used across the site
    tools. Values are not validated at construction time. Use
    ``pycsamt.gis.utils.assert_lat_value`` and related helpers
    if you need strict range checking.

    Examples
    --------
    >>> from pycsamt.site.location import Coord
    >>> c = Coord(10.0, 20.0, 100.0)
    >>> (c.lat, c.lon, c.elev)
    (10.0, 20.0, 100.0)
    """

    lat: float
    lon: float
    elev: float = 0.0


def parse_lat(x: Any) -> float:
    r"""
    Parse a latitude value from decimal or DMS text.

    Accepts floats, ints, or strings like
    ``"10"``, ``"-3.2"``, ``"12 30 00 N"``, ``"12:30:00N"``,
    ``"12d30mE"`` (hemisphere letter required for disambiguation
    when not using sign).

    Parameters
    ----------
    x : any
        Numeric or string representation of a latitude value.

    Returns
    -------
    float
        Parsed latitude in decimal degrees. Returns ``nan`` on
        parsing failure.

    Notes
    -----
    Hemisphere letters are interpreted as ``N=+1`` and ``S=-1``.
    Signed decimal values override hemisphere if both are given.

    Examples
    --------
    >>> from pycsamt.site.location import parse_lat
    >>> parse_lat("12 30 0 N")
    12.5
    >>> parse_lat("-7.25")
    -7.25

    See Also
    --------
    parse_lon, parse_elev
    """

    return float(_parse_angle(x, _HEMI_LAT))


def parse_lon(x: Any) -> float:
    r"""
    Parse a longitude value from decimal or DMS text.

    Accepts floats, ints, or strings like
    ``"20"``, ``"3.5W"``, ``"12 30 00 E"``, or
    ``"12:30:00W"``.

    Parameters
    ----------
    x : any
        Numeric or string representation of a longitude value.

    Returns
    -------
    float
        Parsed longitude in decimal degrees. Returns ``nan`` on
        parsing failure.

    Notes
    -----
    Hemisphere letters are interpreted as ``E=+1`` and ``W=-1``.
    Signed decimal values override hemisphere if both are given.

    Examples
    --------
    >>> from pycsamt.site.location import parse_lon
    >>> parse_lon("3.5W")
    -3.5

    See Also
    --------
    parse_lat, parse_elev
    """

    return float(_parse_angle(x, _HEMI_LON))


def parse_elev(x: Any) -> float:
    r"""
    Parse elevation from a numeric or string value.

    Parameters
    ----------
    x : any
        Numeric or string representation of elevation in meters.

    Returns
    -------
    float
        Elevation in meters. Returns ``nan`` on parsing failure.

    Examples
    --------
    >>> from pycsamt.site.location import parse_elev
    >>> parse_elev("123.4")
    123.4
    """
    try:
        return float(x)
    except Exception:
        return math.nan


def ensure_head_coords(
    ed: Any,
    *,
    lat: Any | None = None,
    lon: Any | None = None,
    elev: Any | None = None,
    empty: float | None = None,
) -> Any:
    r"""
    Create or update HEAD coordinates on an EDI-like object.

    This function guarantees that the EDI "head" section exists
    and that the ``lat``, ``lon`` and ``long`` aliases, and
    ``elev`` fields are present and numeric. Inputs are parsed
    with :func:`parse_lat`, :func:`parse_lon` and
    :func:`parse_elev`, validated with
    ``pycsamt.gis.utils.assert_*``, and then written back. When
    a value is missing or not finite, a default empty sentinel
    is used (0.0 by default).

    Parameters
    ----------
    ed : any
        EDI-like object providing ``get_section('head')`` and
        optionally ``set_section('head', head)``. A fallback
        attribute-based path is used when the API is not
        available.
    lat, lon, elev : any, optional
        Optional overrides. If omitted, the current values are
        read from the head and re-validated. If still missing,
        the ``empty`` value is used.
    empty : float, optional
        Empty sentinel for lat, lon and elev. Default is 0.0.

    Returns
    -------
    any
        The head section object that now carries ``lat``, ``lon``
        (and alias ``long``) and ``elev``.

    Notes
    -----
    This routine writes both ``lon`` and ``long`` to maximize
    compatibility with various EDI headers.

    Examples
    --------
    >>> from pycsamt.site.location import ensure_head_coords
    >>> class H:
    ...     pass
    >>> class EDI:
    ...     def __init__(self):
    ...         self._h = H()
    ...
    ...     def get_section(self, k):
    ...         return self._h if k == "head" else None
    ...
    ...     def set_section(self, k, v):
    ...         self._h = v
    >>> ed = EDI()
    >>> h = ensure_head_coords(ed, lat="12N", lon="3E", elev="10")
    >>> (h.lat, h.lon, h.elev)
    (12.0, 3.0, 10.0)

    See Also
    --------
    parse_lat, parse_lon, parse_elev
    pycsamt.site.utils.set_coords
    pycsamt.site.utils.get_coords
    """

    # for lat/lon/elev we use 0.0 as empty
    emp = 0.0 if empty is None else float(empty)
    h = _ensure_head(ed)

    la0 = getattr(h, "lat", None)
    lo0 = getattr(h, "long", getattr(h, "lon", None))
    ev0 = getattr(h, "elev", getattr(h, "elevation", None))

    # parse inputs; if None, keep existing; if still None → emp
    la = (
        parse_lat(lat)
        if lat is not None
        else (parse_lat(la0) if la0 is not None else emp)
    )
    lo = (
        parse_lon(lon)
        if lon is not None
        else (parse_lon(lo0) if lo0 is not None else emp)
    )
    ev = (
        parse_elev(elev)
        if elev is not None
        else (parse_elev(ev0) if ev0 is not None else emp)
    )

    la = assert_lat_value(la) if math.isfinite(la) else emp
    lo = assert_lon_value(lo) if math.isfinite(lo) else emp
    ev = assert_elevation_value(ev) if math.isfinite(ev) else emp

    # write both 'lon' and 'long' (and survive slots)
    try:
        h.lat = float(la)
    except Exception:
        pass

    wrote_lon = False
    try:
        h.lon = float(lo)
        wrote_lon = True
    except Exception:
        pass
    try:
        h.long = float(lo)
        wrote_lon = True or wrote_lon
    except Exception:
        pass
    try:
        h.elev = float(ev)
    except Exception:
        pass

    # after trying to set h.lon and h.long
    if not wrote_lon:
        try:
            nh = type("Head", (), {})()
            for k in (
                "lat",
                "long",
                "elev",
                "dataid",
                "station",
                "name",
                "sitename",
            ):
                if hasattr(h, k):
                    setattr(nh, k, getattr(h, k))
            nh.lon = float(lo)
            nh.long = float(lo)
            ed.set_section("head", nh)  # preferred API
            h = nh
        except Exception:
            # if set_section isn't available, also assign attribute
            # so that get_section can still discover it.
            try:
                ed.Head = nh  # fallback path
            except Exception:
                pass

    return h


def apply_topography(
    ed_or_sites: Any,
    frame: Any,
    *,
    empty: float | None = None,
    inplace: bool = True,
) -> Any:
    r"""
    Update site coordinates from a tabular frame by station id.

    Matches rows in ``frame`` against the EDI "station" or
    related identifiers and writes the associated ``latitude``,
    ``longitude`` and ``elevation`` into the EDI head section.

    Parameters
    ----------
    ed_or_sites : any or iterable
        A single EDI-like object, an iterable of EDI-like
        objects, or a container with a private ``._items`` list
        of site objects having an ``.edi`` attribute.
    frame : pandas.DataFrame or dict-like
        Table with columns naming a station id and coordinates.
        Station id columns are matched in order among
        ``["station","site","dataid","id","name"]``. Coordinate
        columns are searched among
        ``["latitude","lat"]``, ``["longitude","lon","long"]``,
        and ``["elevation","elev","alt"]``.
    empty : float, optional
        Empty sentinel applied when values are missing. Default
        is 0.0.
    inplace : bool, optional
        If ``True`` (default), update the provided objects. If
        ``False``, return a deep-copied updated structure.

    Returns
    -------
    any
        The updated object(s). For lists or containers, the
        return type mirrors the input.

    Notes
    -----
    Matching is case-insensitive and robust to whitespace. For
    containers, the function duck-types a ``._items`` attribute
    and attempts to update the underlying ``.edi`` objects.

    Examples
    --------
    >>> import pandas as pd
    >>> from pycsamt.site.location import apply_topography
    >>> # ed is an EDI-like object with a valid head and dataid
    >>> df = pd.DataFrame(
    ...     {
    ...         "station": ["S1"],
    ...         "latitude": [10.0],
    ...         "longitude": [2.0],
    ...         "elevation": [100.0],
    ...     }
    ... )
    >>> _ = apply_topography(ed, df, inplace=True)  # doctest: +SKIP
    """

    emp = 0.0 if empty is None else float(empty)

    def _apply_one(ed: Any) -> Any:
        sid = _get_station(ed)
        row = _match_row(frame, sid)
        if row is not None:
            _ensure_head(ed)
            _set_coords_from_row(ed, row, empty=emp)
        return ed

    # sites container (duck-typed)
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

    # list/tuple of EDIFile
    try:
        if isinstance(ed_or_sites, (list, tuple)):
            vec = (
                [copy.deepcopy(x) for x in ed_or_sites]
                if not inplace
                else list(ed_or_sites)
            )
            for ed in vec:
                _apply_one(ed)
            return vec
    except Exception:
        pass

    # single EDI
    return _apply_one(ed_or_sites)


def project(
    pts: Sequence[tuple[float, float]] | tuple[float, float],
    *,
    crs_from: Any,
    crs_to: Any,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Project points from one CRS to another using pyproj or GDAL.

    Parameters
    ----------
    pts : sequence of (float, float) or (float, float)
        Either a single coordinate pair or a sequence of pairs.
        Coordinates are interpreted according to ``crs_from``.
    crs_from : any
        Source CRS. For pyproj this can be an EPSG code, PROJ
        string, or CRS object. For GDAL, EPSG code, WKT, PROJ4,
        or a SpatialReference.
    crs_to : any
        Target CRS in the same form as ``crs_from``.

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Arrays ``(X, Y)`` holding transformed coordinates.

    Raises
    ------
    RuntimeError
        If neither pyproj nor GDAL are available.
    TypeError
        If CRS specification cannot be parsed.

    Notes
    -----
    When both pyproj and GDAL are available, pyproj is used. The
    axis order is forced to the traditional GIS order
    ``(x, y)``.

    Examples
    --------
    >>> from pycsamt.site.location import project
    >>> X, Y = project(
    ...     [(0.0, 0.0)], crs_from="EPSG:4326", crs_to="EPSG:3857"
    ... )  # doctest: +SKIP

    See Also
    --------
    pycsamt.site.location.distance
    pycsamt.site.location.bearing
    """

    try:
        from pyproj import Transformer  # type: ignore

        use_pyproj = True
    except Exception:
        use_pyproj = False

    if not use_pyproj and not HAS_GDAL:
        raise RuntimeError("project() requires pyproj or GDAL.")

    if isinstance(pts, tuple) and len(pts) == 2:
        xs, ys = [pts[0]], [pts[1]]
    else:
        xs, ys = zip(*pts)  # type: ignore[misc]
    xs_a = np.asarray(xs, float)
    ys_a = np.asarray(ys, float)

    if use_pyproj:
        tr = Transformer.from_crs(crs_from, crs_to, always_xy=True)
        X, Y = tr.transform(xs_a, ys_a)
        return np.asarray(X, float), np.asarray(Y, float)

    # GDAL fallback
    from osgeo import osr  # type: ignore

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
        s_from.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        s_to.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
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
    return (np.asarray(out_x, float), np.asarray(out_y, float))


def distance(
    a: Coord | tuple[float, float],
    b: Coord | tuple[float, float],
    *,
    mode: str = "geodetic",
    crs_to: Any | None = None,
) -> float:
    r"""
    Compute distance between two geographic points.

    Supports three modes:
    ``geodetic`` uses a haversine approximation on a spherical
    Earth of radius ``_EARTH_R``.
    ``flat`` uses a local small-extent planar approximation
    (about the mid-latitude).
    ``utm`` projects both points to UTM (auto zone unless
    ``crs_to`` is given) and computes the Euclidean distance.

    Parameters
    ----------
    a, b : Coord or (float, float)
        Points given as ``Coord`` or ``(lat, lon)``.
    mode : {'geodetic', 'flat', 'utm'}, optional
        Distance model. Default is ``'geodetic'``.
    crs_to : any, optional
        Target CRS for ``'utm'`` mode. If omitted, the UTM zone
        is inferred from the mean coordinate.

    Returns
    -------
    float
        Distance in meters.

    Raises
    ------
    ValueError
        If ``mode`` is not one of the supported options.

    Notes
    -----
    Geodetic mode implements the haversine formula:

    .. math::

       d = 2 R \arcsin\left( \sqrt{ \sin^2\Delta\phi/2 +
           \cos\phi_1 \cos\phi_2 \sin^2\Delta\lambda/2 } \right)

    where angles are in radians.

    Examples
    --------
    >>> from pycsamt.site.location import distance, Coord
    >>> distance(Coord(0, 0, 0), Coord(0, 1, 0), mode="geodetic")
    111000.0  # doctest: +ELLIPSIS

    See Also
    --------
    pycsamt.site.location.bearing
    pycsamt.site.location.project
    """

    la1, lo1 = (a.lat, a.lon) if isinstance(a, Coord) else a
    la2, lo2 = (b.lat, b.lon) if isinstance(b, Coord) else b
    la1 = float(assert_lat_value(la1))
    lo1 = float(assert_lon_value(lo1))
    la2 = float(assert_lat_value(la2))
    lo2 = float(assert_lon_value(lo2))

    m = mode.lower()
    if m == "geodetic":
        dlat = _rad(la2 - la1)
        dlon = _rad(lo2 - lo1)
        A = (
            math.sin(dlat / 2) ** 2
            + math.cos(_rad(la1))
            * math.cos(_rad(la2))
            * math.sin(dlon / 2) ** 2
        )
        return 2.0 * _EARTH_R * math.asin(min(1.0, math.sqrt(A)))
    if m == "flat":
        dx, dy = _flat_offsets_m(la1, lo1, la2, lo2)
        return math.hypot(dx, dy)
    if m == "utm":
        epsg = (
            _infer_utm_epsg((la1 + la2) * 0.5, (lo1 + lo2) * 0.5)
            if crs_to is None
            else crs_to
        )
        X, Y = project(
            [(lo1, la1), (lo2, la2)], crs_from="EPSG:4326", crs_to=epsg
        )
        return float(math.hypot(float(X[1] - X[0]), float(Y[1] - Y[0])))
    raise ValueError("mode must be 'geodetic','flat','utm'")


def bearing(
    a: Coord | tuple[float, float],
    b: Coord | tuple[float, float],
    *,
    mode: str = "geodetic",
    crs_to: Any | None = None,
) -> float:
    r"""
    Compute the azimuth from point ``a`` to point ``b``.

    Supports the same three modes as :func:`distance`:
    ``geodetic``, ``flat``, and ``utm``. The azimuth is expressed
    in degrees with 0 deg pointing to north and 90 deg to east.

    Parameters
    ----------
    a, b : Coord or (float, float)
        Points given as ``Coord`` or ``(lat, lon)``.
    mode : {'geodetic', 'flat', 'utm'}, optional
        Bearing model. Default is ``'geodetic'``.
    crs_to : any, optional
        Target CRS for ``'utm'`` mode. If omitted, the UTM zone
        is inferred.

    Returns
    -------
    float
        Azimuth in degrees in the range ``[0, 360)``.

    Raises
    ------
    ValueError
        If ``mode`` is not supported.

    Notes
    -----
    Geodetic mode uses the spherical forward azimuth:

    .. math::

       \theta = \operatorname{atan2}(\sin\Delta\lambda\cos\phi_2,\,
               \cos\phi_1\sin\phi_2 -
               \sin\phi_1\cos\phi_2\cos\Delta\lambda)

    Angles are converted to degrees and wrapped to ``[0, 360)``.

    Examples
    --------
    >>> from pycsamt.site.location import bearing, Coord
    >>> bearing(Coord(0, 0, 0), Coord(1, 0, 0))
    0.0
    >>> bearing(Coord(0, 0, 0), Coord(0, 1, 0))
    90.0

    See Also
    --------
    pycsamt.site.location.distance
    pycsamt.site.location.chainage_along
    """

    la1, lo1 = (a.lat, a.lon) if isinstance(a, Coord) else a
    la2, lo2 = (b.lat, b.lon) if isinstance(b, Coord) else b
    la1 = float(assert_lat_value(la1))
    lo1 = float(assert_lon_value(lo1))
    la2 = float(assert_lat_value(la2))
    lo2 = float(assert_lon_value(lo2))

    m = mode.lower()
    if m == "geodetic":
        phi1 = _rad(la1)
        phi2 = _rad(la2)
        dlambda = _rad(lo2 - lo1)
        y = math.sin(dlambda) * math.cos(phi2)
        x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(
            phi2
        ) * math.cos(dlambda)
        theta = math.degrees(math.atan2(y, x))
        return (theta + 360.0) % 360.0
    if m == "flat":
        dx, dy = _flat_offsets_m(la1, lo1, la2, lo2)
        θ = math.degrees(math.atan2(dx, dy))
        return (θ + 360.0) % 360.0
    if m == "utm":
        epsg = (
            _infer_utm_epsg((la1 + la2) * 0.5, (lo1 + lo2) * 0.5)
            if crs_to is None
            else crs_to
        )
        X, Y = project(
            [(lo1, la1), (lo2, la2)], crs_from="EPSG:4326", crs_to=epsg
        )
        dx = float(X[1] - X[0])
        dy = float(Y[1] - Y[0])
        theta = math.degrees(math.atan2(dx, dy))
        return (theta + 360.0) % 360.0
    raise ValueError("mode must be 'geodetic','flat','utm'")


def chainage_along(
    origin: Coord | tuple[float, float],
    azimuth: float,
    pts: Sequence[Coord | tuple[float, float]] | Coord | tuple[float, float],
) -> np.ndarray | float:
    r"""
    Project points onto a profile axis and return chainages.

    Chainage is the signed distance along the axis defined by an
    origin and an azimuth (0 deg north, 90 deg east). A local
    flat-Earth approximation is used, with 1 deg roughly equal
    to 111 km scaled by cosine of latitude for the east axis.

    Parameters
    ----------
    origin : Coord or (float, float)
        Profile origin as ``Coord`` or ``(lat, lon)``.
    azimuth : float
        Axis azimuth in degrees, 0 deg north, 90 deg east.
    pts : sequence of Coord or (float, float) or single
        Single point or a sequence of points to be projected.

    Returns
    -------
    numpy.ndarray or float
        Chainage(s) in meters. Returns a scalar for a single
        point, or a 1-D array for multiple points.

    Notes
    -----
    Let local offsets be ``dx`` east and ``dy`` north from the
    origin. The chainage uses

    .. math::

       s = dx \sin A + dy \cos A

    where :math:`A` is the azimuth in radians and
    :math:`dx, dy` are derived from degree differences using a
    local metric scale.

    Examples
    --------
    >>> from pycsamt.site.location import chainage_along
    >>> s = chainage_along((0.0, 0.0), 90.0, (0.0, 1.0))
    >>> s > 100000.0
    True

    See Also
    --------
    pycsamt.site.location.bearing
    pycsamt.site.location.distance
    """

    if isinstance(origin, Coord):
        la0, lo0 = origin.lat, origin.lon
    else:
        la0, lo0 = origin
    la0 = float(assert_lat_value(la0))
    lo0 = float(assert_lon_value(lo0))
    az = math.radians(float(azimuth))
    mperdeg = 111_000.0

    def _one(p: Coord | tuple[float, float]) -> float:
        if isinstance(p, Coord):
            la, lo = p.lat, p.lon
        else:
            la, lo = p
        la = float(assert_lat_value(la))
        lo = float(assert_lon_value(lo))
        dy = (la - la0) * mperdeg
        dx = (lo - lo0) * mperdeg * math.cos(_rad(la0))
        # project onto axis (north-clockwise azimuth)
        return dx * math.sin(az) + dy * math.cos(az)

    if isinstance(pts, (tuple, Coord)):
        return float(_one(pts))
    out = [_one(p) for p in pts]  # type: ignore[arg-type]
    return np.asarray(out, float)


# --------------helpers -----------------------
def _lower_cols(df: Any) -> Iterable[str]:
    try:
        cols = getattr(df, "columns", None)
        if cols is not None:
            return [str(c).strip() for c in cols]
        # accept Index/iterables (Series.index, lists, tuples)
        return [str(c).strip() for c in list(df)]
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
    except:
        return None

    idc = _pick_col(df, ("station", "site", "dataid", "id", "name"))
    if not idc:
        return None

    try:
        sref = str(sid).strip().upper()
        col = df[idc]
        # force to string, strip, and upper for robust match
        coln = col.astype(str).str.strip().str.upper()
        hits = df[coln == sref]
        if getattr(hits, "empty", True):
            return None
        return hits.iloc[0]
    except:
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


def _rad(x: float) -> float:
    return math.radians(float(x))


def _infer_utm_epsg(lat: float, lon: float) -> int:
    z = int(math.floor((float(lon) + 180.0) / 6.0) + 1)
    return (32600 + z) if float(lat) >= 0.0 else (32700 + z)


def _flat_offsets_m(
    la1: float,
    lo1: float,
    la2: float,
    lo2: float,
) -> tuple[float, float]:
    la_mid = math.radians((la1 + la2) * 0.5)
    m_lat = 111_000.0
    m_lon = 111_000.0 * math.cos(la_mid)
    dy = (la2 - la1) * m_lat
    dx = (lo2 - lo1) * m_lon
    return dx, dy


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
