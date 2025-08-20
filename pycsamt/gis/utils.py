# -*- coding: utf-8 -*-
# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

r"""
GIS utilities for *pycsamt* (v2). Provides coordinate
transformations between geographic (lat/lon) and UTM using
GDAL (preferred) or a PROJ fallback (``pyproj``).

Core features
-------------
- Geographic ⇄ UTM conversions (scalar and vectorized).
- Flexible latitude/longitude parsing (decimal and DMS).
- Small helpers for UTM zones and EPSG lookup.

Examples
--------
>>> from pycsamt.gis.utils import project_point_ll2utm
>>> project_point_ll2utm(-118.34, 34.05)  # doctest: +SKIP

.. note::

   Parts of this module are adapted from the **MTPy** project
   (https://github.com/MTgeophysics/mtpy-v2). Original work:
   Krieger, L., & Peacock, J. (2014), *MTpy: A Python toolbox
   for magnetotellurics*, Computers & Geosciences, 72, 167–175,
   https://doi.org/10.1016/j.cageo.2014.07.013.

   The implementation here has been revised, fixed, and
   extended for *pycsamt* v2, with additional validation,
   deprecations, and compatibility layers.
"""

import numpy as np
from typing import Optional, Tuple

from ..decorators import Deprecated
from ..log.logger import get_logger
from .config import (
    HAS_GDAL,
    EPSG_DICT,
    GDALMissingError,
    GisError,
)
from .constants import (
    DEG2RAD,
    RAD2DEG,
    ELLIPSOIDS,
    _EQUATORIAL_RADIUS_IDX,
    _ECC_SQUARED_IDX,
    utm_letter_designator,
)

if HAS_GDAL:
    from osgeo import osr
    from osgeo.ogr import OGRERR_NONE
else:
    import pyproj  # type: ignore

logger = get_logger(__name__)


__all__ = [
    'project_point_ll2utm',
    'project_point_utm2ll',
    'll_to_utm',
    'utm_to_ll',
    'get_utm_string_from_sr',
    'get_utm_zone',
    'utm_zone_to_epsg',
    'get_epsg',
    'assert_xy_coordinate_system', 
    'decimal_to_dms', 
    'dms_to_decimal',
    'GisError'
]

def assert_xy_coordinate_system(x, y) -> str:
    r"""
    Infer the coordinate system of paired ``x`` and ``y`` arrays.

    Heuristics detect one of three systems:

    * ``'ll'``  — longitude/latitude in decimal degrees.
    * ``'dms'`` — degree-minute-second strings ``DD:MM:SS``.
    * ``'utm'`` — any other numeric system (fallback).

    Parameters
    ----------
    x, y : array-like (1D)
        Arrays of horizontal (``x``) and vertical (``y``)
        positions. Elements may be numeric or DMS strings.

    Returns
    -------
    cs : {'utm', 'dms', 'll'}
        Inferred coordinate system token.

    Notes
    -----
    - DMS is detected if any array consists entirely of
      strings containing the ``':'`` separator.
    - Decimal degrees are detected when either orientation
      matches bounds:

      * ``|x| <= 180`` and ``|y| <= 90`` (lon, lat), or
      * ``|x| <= 90`` and ``|y| <= 180`` (lat, lon).

    - If neither DMS nor decimal-degree bounds match,
      the function returns ``'utm'`` by default.

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> x, y = np.random.rand(7) * 1.0, np.arange(7) * 0.1
    >>> assert_xy_coordinate_system(x, y)
    'll'
    >>> assert_xy_coordinate_system(x * 1000, y * 1000)
    'utm'
    >>> x = ['28:24:43.08', '28:24:42.69', '28:24:42.31']
    >>> y = ['109:19:58.34', '109:19:58.93', '109:19:59.51']
    >>> assert_xy_coordinate_system(x, y)
    'dms'
    """
    def _is_dms(a) -> bool:
        a = np.asarray(a, dtype=object).ravel()
        if a.size == 0:
            return False
        return all(":" in str(v) for v in a)

    # Normalize to numpy arrays
    xa = np.asarray(x, dtype=object).ravel()
    ya = np.asarray(y, dtype=object).ravel()

    # DMS detection
    if _is_dms(xa) or _is_dms(ya):
        return "dms"

    # Try numeric interpretation
    try:
        xf = np.asarray(xa, dtype=float)
        yf = np.asarray(ya, dtype=float)
    except (ValueError, TypeError):
        return "utm"

    # Decimal-degree bounds (either orientation)
    x180_y90 = (
        np.all(np.abs(xf) <= 180.0)
        and np.all(np.abs(yf) <= 90.0)
    )
    x90_y180 = (
        np.all(np.abs(xf) <= 90.0)
        and np.all(np.abs(yf) <= 180.0)
    )
    if x180_y90 or x90_y180:
        return "ll"

    return "utm"

def _assert_minutes(minutes: float) -> float:
    r"""
    Validate a minutes component in DMS notation.

    Ensures that ``minutes`` lies in the half-open interval
    ``[0, 60)``. Returns the input when valid, otherwise raises
    ``ValueError``.

    Parameters
    ----------
    minutes : float
        Minutes component to validate.

    Returns
    -------
    float
        The validated minutes value.

    Raises
    ------
    ValueError
        If ``minutes`` is not in ``[0, 60)``.

    Examples
    --------
    >>> _assert_minutes(59.9)
    59.9
    """
    if not (0 <= minutes < 60):
        msg = f"Minutes must be in [0, 60): {minutes}"
        raise ValueError(msg)
    return minutes


def _assert_seconds(seconds: float) -> float:
    r"""
    Validate a seconds component in DMS notation.

    Ensures that ``seconds`` lies in the half-open interval
    ``[0, 60)``. Returns the input when valid, otherwise raises
    ``ValueError``.

    Parameters
    ----------
    seconds : float
        Seconds component to validate.

    Returns
    -------
    float
        The validated seconds value.

    Raises
    ------
    ValueError
        If ``seconds`` is not in ``[0, 60)``.

    Examples
    --------
    >>> _assert_seconds(12.345)
    12.345
    """
    if not (0 <= seconds < 60):
        msg = f"Seconds must be in [0, 60): {seconds}"
        raise ValueError(msg)
    return seconds


def _rollover_dms(
    unit: float,
    value: float,
) -> tuple[float, float]:
    r"""
    Carry over excess seconds or minutes into the higher unit.

    If ``value >= 60``, this function computes the carry and
    remainder, returning ``(unit + carry, remainder)``. If no
    carry is needed, returns ``(unit, value)``.

    Parameters
    ----------
    unit : float
        The higher unit (minutes for seconds rollover, or
        degrees for minutes rollover).
    value : float
        The lower unit to check and normalize.

    Returns
    -------
    tuple of float
        A pair ``(new_unit, remainder)`` after rollover.

    Examples
    --------
    >>> _rollover_dms(10, 120)   # 120 sec -> 2 min
    (12, 0)
    """
    carry = int(value // 60)
    remainder = value % 60
    if carry:
        return (unit + carry, remainder)
    return (unit, remainder)


def dms_to_decimal(position_str):
    r"""
    Convert a DMS string or decimal degrees to a float.

    This is a convenience wrapper around
    :func:`convert_position_str2float`. It accepts either a
    decimal degrees string/number or a sexagesimal string in
    ``'DD:MM:SS.sss'`` format.

    Parameters
    ----------
    position_str : str or float
        Coordinate in decimal degrees or DMS string.

    Returns
    -------
    float or None
        Decimal degrees value, or ``None`` for invalid input.

    Examples
    --------
    >>> dms_to_decimal("34:03:00")
    34.05
    >>> dms_to_decimal("-118.34")
    -118.34
    """
    return convert_position_str2float(position_str)


def convert_position_str2float(
    position_str: str,
) -> float | None:
    r"""
    Convert DMS string ``'DD:MM:SS.sss'`` or decimal string to
    a float in decimal degrees.

    The function first attempts a direct ``float()`` cast. If
    that fails, it parses the DMS components, validates them,
    and applies rollover so that values like ``'10:60:00'``
    become ``'11:00:00'``.

    Parameters
    ----------
    position_str : str
        Input coordinate as decimal degrees or DMS string.

    Returns
    -------
    float or None
        Decimal degrees value. Returns ``None`` for inputs
        that are ``None`` or the string ``'None'``.

    Raises
    ------
    ValueError
        If the string is not a valid decimal or DMS in the
        form ``DD:MM:SS`` (with optional decimals).

    Notes
    -----
    - Sign is taken from the degrees component only. Minutes
      and seconds are treated as non-negative magnitudes.
    - Rollover is applied to keep minutes and seconds within
      ``[0, 60)``.

    Examples
    --------
    >>> convert_position_str2float("34:03:00")
    34.05
    >>> convert_position_str2float("-118:20:24")
    -118.34
    >>> convert_position_str2float("-118.34")
    -118.34
    """
    if position_str in (None, "None"):
        return None

    try:
        return float(position_str)
    except (TypeError, ValueError):
        pass

    parts = str(position_str).split(":")
    if len(parts) != 3:
        msg = (
            "Invalid DMS format, expected 'DD:MM:SS': "
            f"{position_str}"
        )
        raise ValueError(msg)

    deg = float(parts[0])
    minutes = _assert_minutes(float(parts[1]))
    seconds = _assert_seconds(float(parts[2]))

    minutes, seconds = _rollover_dms(minutes, seconds)
    deg, minutes = _rollover_dms(deg, minutes)

    sign = -1 if deg < 0 else 1
    return sign * (abs(deg) + minutes / 60.0 + seconds / 3600.0)



def assert_lat_value(latitude) -> Optional[float]:
    r"""
    Validate and coerce a latitude to decimal degrees.

    Accepts a float-like value or a sexagesimal string and
    returns the value in decimal degrees. Returns ``None`` for
    ``None``-like inputs. Raises ``ValueError`` if the coerced
    value is out of valid latitude bounds.

    Parameters
    ----------
    latitude : float or str or None
        Latitude in decimal degrees or as a string in
        ``"DD:MM:SS.sss"`` format. ``None`` and ``"None"``
        are treated as missing.

    Returns
    -------
    float or None
        Decimal degrees in ``[-90, 90)`` or ``None`` when the
        input is missing.

    Notes
    -----
    Sexagesimal parsing is delegated to
    :func:`convert_position_str2float`.

    Examples
    --------
    >>> assert_lat_value("34:03:00")
    34.05
    """
    if latitude in [None, "None"]:
        return None
    try:
        lat_value = float(latitude)
    except TypeError:
        return None
    except ValueError:
        lat_value = convert_position_str2float(latitude)

    if abs(lat_value) >= 90:
        print("==> The lat_value =", lat_value)
        raise ValueError(
            f"|Latitude| > 90, unacceptable!: {lat_value!r}"
        )
    return lat_value


def assert_lon_value(longitude) -> Optional[float]:
    r"""
    Validate and coerce a longitude to decimal degrees.

    Accepts a float-like value or a sexagesimal string and
    returns the value in decimal degrees. Returns ``None`` for
    ``None``-like inputs. Raises ``ValueError`` if the coerced
    value is out of valid longitude bounds.

    Parameters
    ----------
    longitude : float or str or None
        Longitude in decimal degrees or as a string in
        ``"DD:MM:SS.sss"`` format. ``None`` and ``"None"``
        are treated as missing.

    Returns
    -------
    float or None
        Decimal degrees in ``[-180, 180)`` or ``None`` when
        the input is missing.

    Notes
    -----
    Sexagesimal parsing is delegated to
    :func:`convert_position_str2float`.

    Examples
    --------
    >>> assert_lon_value("-118:20:24")
    -118.34
    """
    if longitude in [None, "None"]:
        return None
    try:
        lon_value = float(longitude)
    except TypeError:
        return None
    except ValueError:
        lon_value = convert_position_str2float(longitude)

    if abs(lon_value) >= 180:
        raise ValueError(
            f"|Longitude| > 180, unacceptable!: {lon_value!r}"
        )
    return lon_value


def assert_elevation_value(elevation) -> float:
    r"""
    Validate and coerce an elevation to a floating number.

    Attempts to cast ``elevation`` to ``float``. If coercion
    fails, returns ``0.0`` and logs a warning.

    Parameters
    ----------
    elevation : Any
        Elevation value expected to be numeric.

    Returns
    -------
    float
        Elevation as a floating number. ``0.0`` is returned on
        invalid input.

    Notes
    -----
    A warning is logged via the package logger when the input
    is not numeric.

    Examples
    --------
    >>> assert_elevation_value("12.5")
    12.5
    >>> assert_elevation_value("oops")
    0.0
    """
    try:
        elev_value = float(elevation)
    except (ValueError, TypeError):
        elev_value = 0.0
        logger.warning(
            "{} is not a number, setting elevation to 0".format(
                elevation
            )
        )
    return elev_value


def decimal_to_dms(position: float) -> str:
    r"""
    Convert decimal degrees to a ``DD:MM:SS.ss`` string.

    Parameters
    ----------
    position : float
        Decimal degrees for latitude or longitude.

    Returns
    -------
    str
        Sexagesimal string in ``"DD:MM:SS.ss"`` format.

    See Also
    --------
    convert_position_float2str
        Underlying converter used by this wrapper.

    Examples
    --------
    >>> decimal_to_dms(-118.34563)
    '-118:20:44.27'
    """
    return convert_position_float2str(position)


def convert_position_float2str(position: float) -> str:
    r"""
    Convert a decimal-degree value to ``DD:MM:SS.ss`` string.

    Parameters
    ----------
    position : float
        Decimal degrees of latitude or longitude.

    Returns
    -------
    str
        Sexagesimal string in ``"DD:MM:SS.ss"`` format.

    Notes
    -----
    The seconds field is rounded to 4 decimals to avoid
    carry-over artifacts from floating point precision. If
    rounding pushes seconds to ``60``, minutes are incremented
    and seconds reset to ``0``.

    Examples
    --------
    >>> convert_position_float2str(-118.34563)
    '-118:20:44.27'
    """
    assert type(position) is float, (
        "Given value is not a float"
    )

    deg = int(position)
    sign = -1 if deg < 0 else 1
    deg = abs(deg)

    minutes = (abs(position) - deg) * 60.0
    sec = np.round((minutes - int(minutes)) * 60.0, 4)
    if sec >= 60.0:
        minutes += 1.0
        sec = 0.0
    if int(minutes) == 60:
        deg += 1
        minutes = 0.0

    position_str = "{}:{:02.0f}:{:05.2f}".format(
        sign * int(deg),
        int(minutes),
        sec,
    )
    return position_str

@Deprecated(
    "GDAL SpatialReference → UTM string is deprecated; "
    "use 'get_utm_zone' for standard UTM formatting."
)
def get_utm_string_from_sr(
        spatial_ref: "osr.SpatialReference"
    ) -> str:
    r"""
    Return a UTM zone string (e.g., ``'11N'``) from a GDAL
    ``SpatialReference``.

    Parameters
    ----------
    spatial_ref : osr.SpatialReference
        GDAL spatial reference object with a UTM projection.

    Returns
    -------
    str
        UTM zone string such as ``'11N'`` or ``'55S'``. If a
        valid UTM zone is not encoded, returns ``'0'``.

    Notes
    -----
    This helper is deprecated. Prefer :func:`get_utm_zone`,
    which derives a zone from latitude/longitude using the
    standard UTM rules.

    Examples
    --------
    >>> # sr is an osr.SpatialReference set to a UTM CRS
    >>> # get_utm_string_from_sr(sr)  # doctest: +SKIP
    """
    zone = spatial_ref.GetUTMZone()
    if zone > 0:
        return f"{zone}N"
    if zone < 0:
        return f"{abs(zone)}S"
    return str(zone)


def get_utm_zone(
    latitude: float,
    longitude: float,
) -> Tuple[int, bool, str]:
    r"""
    Compute the UTM zone, hemisphere flag, and zone string for
    a geographic coordinate.

    Parameters
    ----------
    latitude : float
        Latitude in decimal degrees.
    longitude : float
        Longitude in decimal degrees.

    Returns
    -------
    zone_number : int
        UTM zone number in ``[1, 60]``.
    is_northern : bool
        ``True`` for the northern hemisphere, else ``False``.
    zone_string : str
        Concatenation of zone and latitude band letter, e.g.,
        ``'11N'``.

    Notes
    -----
    - Zone width is 6 degrees, starting at ``-180°``. Values
      are wrapped into ``[1, 60]``.
    - The latitude band letter is obtained via
      :func:`utm_letter_designator`.

    Examples
    --------
    >>> get_utm_zone(34.05, -118.34)[2]
    '11S'
    """
    zone_number = int((longitude + 180.0) / 6.0) + 1
    zone_number = ((zone_number - 1) % 60) + 1
    is_northern = latitude >= 0.0
    letter = utm_letter_designator(latitude)
    zone_string = f"{zone_number}{letter}"
    return zone_number, is_northern, zone_string


def utm_zone_to_epsg(
    zone_number: int,
    is_northern: bool,
) -> Optional[int]:
    r"""
    Resolve the WGS84 UTM EPSG code for a given UTM zone.

    Parameters
    ----------
    zone_number : int
        UTM zone number in ``[1, 60]``.
    is_northern : bool
        ``True`` for northern hemisphere, ``False`` for
        southern hemisphere.

    Returns
    -------
    int or None
        Matching EPSG code (e.g., ``32611`` or ``32755``) when
        found, else ``None``.

    Notes
    -----
    The lookup scans :data:`EPSG_DICT` for a PROJ string that
    matches ``+zone``, ``+datum=WGS84``, and an optional
    ``+south`` flag for the southern hemisphere.

    Examples
    --------
    >>> utm_zone_to_epsg(11, True) in (32611, None)
    True
    """
    datum_flag = "" if is_northern else "+south"
    for epsg_code, proj4 in EPSG_DICT.items():
        if (
            f"+zone={zone_number}" in proj4
            and "+datum=WGS84" in proj4
            and datum_flag in proj4
        ):
            return epsg_code
    return None


def get_epsg(
    latitude: float,
    longitude: float,
) -> Optional[int]:
    r"""
    Get the WGS84 UTM EPSG code for a geographic coordinate.

    Parameters
    ----------
    latitude : float
        Latitude in decimal degrees.
    longitude : float
        Longitude in decimal degrees.

    Returns
    -------
    int or None
        EPSG code for the inferred UTM CRS, or ``None`` if no
        match is found.

    Notes
    -----
    This is a convenience wrapper that calls
    :func:`get_utm_zone` and :func:`utm_zone_to_epsg`.

    Examples
    --------
    >>> epsg = get_epsg(34.05, -118.34)
    >>> epsg in (32611, 32711, None)
    True
    """
    zone_number, is_northern, _ = get_utm_zone(
        latitude,
        longitude,
    )
    return utm_zone_to_epsg(zone_number, is_northern)


def project_point_ll2utm(
    lat,
    lon,
    datum: str = "WGS84",
    utm_zone: Optional[str] = None,
    epsg: Optional[int] = None,
):
    r"""
    Transform one geographic point to UTM coordinates.

    Converts a single latitude/longitude (in decimal degrees
    or sexagesimal string) to UTM easting, northing, and zone.
    Uses GDAL/OSR when available, else falls back to PROJ via
    ``pyproj``.

    Parameters
    ----------
    lat : float or str
        Latitude in decimal degrees or sexagesimal string
        (``"DD:MM:SS.sss"``).
    lon : float or str
        Longitude in decimal degrees or sexagesimal string
        (``"DD:MM:SS.sss"``).
    datum : str, default "WGS84"
        Geodetic datum name (e.g., ``"WGS84"``, ``"NAD83"``).
        May also be an EPSG integer if GDAL is used.
    utm_zone : str, optional
        UTM zone designator like ``"55S"``. If omitted or set
        to a non-string (or ``"none"``), the zone is inferred
        from the input point.
    epsg : int, optional
        EPSG code for the projected CRS. When provided, it
        takes precedence over ``utm_zone`` and ``datum``.

    Returns
    -------
    easting : float or None
        UTM easting in meters, or ``None`` when inputs are
        ``None``.
    northing : float or None
        UTM northing in meters, or ``None`` when inputs are
        ``None``.
    zone : str or None
        UTM zone string (e.g., ``"55S"``). ``None`` when not
        determinable or inputs are ``None``.

    Notes
    -----
    - This variant accepts a single point. For vectorized
      inputs, use :func:`project_points_ll2utm`.
    - Inputs are validated with ``assert_lat_value`` and
      ``assert_lon_value`` before projection.
    - With GDAL, transforms use
      ``osr.CoordinateTransformation``; with PROJ they use
      ``pyproj.Proj``.

    Examples
    --------
    >>> e, n, z = project_point_ll2utm(34.05, -118.34)
    >>> z
    '11S'
    """
    if lat is None or lon is None:
        return (None, None, None)

    if np.iterable(lat) and np.iterable(lon):
        lat = np.array(
            [assert_lat_value(v) for v in lat]
        )
        lon = np.array(
            [assert_lon_value(v) for v in lon]
        )
        assert lat.size == lon.size
    else:
        lat = np.array([assert_lat_value(lat)])
        lon = np.array([assert_lon_value(lon)])

    if HAS_GDAL:
        ll_cs = osr.SpatialReference()
        if isinstance(datum, int):
            ogrerr = ll_cs.ImportFromEPSG(datum)
            if ogrerr != OGRERR_NONE:
                raise GisError(
                    "GDAL/OSR error code: {}".format(ogrerr)
                )
        elif isinstance(datum, str):
            ogrerr = ll_cs.SetWellKnownGeogCS(datum)
            if ogrerr != OGRERR_NONE:
                raise GisError(
                    "GDAL/OSR error code: {}".format(ogrerr)
                )
        else:
            raise GisError(
                "datum {!r} not understood; use EPSG int or "
                "well-known datum string".format(datum)
            )
        utm_cs = osr.SpatialReference()
    else:
        pp = None  # type: ignore[assignment]

    if isinstance(epsg, int):
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise GisError(
                    "GDAL/OSR error code: {}".format(ogrerr)
                )
        else:
            pp = pyproj.Proj("+init=EPSG:{}".format(epsg))
    elif epsg is None:
        if HAS_GDAL:
            ogrerr = utm_cs.CopyGeogCSFrom(ll_cs)
            if ogrerr != OGRERR_NONE:
                raise GisError(
                    "GDAL/OSR error code: {}".format(ogrerr)
                )
        if (
            utm_zone is None
            or not isinstance(utm_zone, str)
            or utm_zone.lower() == "none"
        ):
            zone_num, is_north, utm_zone = get_utm_zone(
                float(lat.mean()),
                float(lon.mean()),
            )
        else:
            zone_num = int(utm_zone[:-1])
            is_north = utm_zone[-1].lower() > "n"

        if HAS_GDAL:
            utm_cs.SetUTM(zone_num, is_north)
        else:
            projstring = (
                "+proj=utm +zone={} +{} +datum={}".format(
                    zone_num,
                    "north" if is_north else "south",
                    datum,
                )
            )
            pp = pyproj.Proj(projstring)

    proj = np.zeros_like(
        lat,
        dtype=[
            ("easting", np.float64),
            ("northing", np.float64),
            ("elev", np.float64),
            ("utm_zone", "U4"),
        ],
    )

    if HAS_GDAL:
        ll2utm = osr.CoordinateTransformation(
            ll_cs, utm_cs
        ).TransformPoint
    else:
        ll2utm = pp  # type: ignore[assignment]

    for ii in range(lat.size):
        if HAS_GDAL:
            x, y, z = ll2utm(
                float(lon[ii]),
                float(lat[ii]),
            )
            proj["easting"][ii] = x
            proj["northing"][ii] = y
            proj["elev"][ii] = z
        else:
            assert ll2utm is not None
            x, y = ll2utm(
                float(lon[ii]),
                float(lat[ii]),
            )
            proj["easting"][ii] = x
            proj["northing"][ii] = y
        if utm_zone is None:
            _, _, zstr = get_utm_zone(
                float(lat[ii]),
                float(lon[ii]),
            )
            proj["utm_zone"][ii] = zstr
        else:
            proj["utm_zone"][ii] = utm_zone

    if len(proj) == 1:
        return (
            proj["easting"][0],
            proj["northing"][0],
            proj["utm_zone"][0],
        )
    return np.rec.array(proj)


def project_point_utm2ll(
    easting: float,
    northing: float,
    utm_zone,
    datum: str = "WGS84",
    epsg: Optional[int] = 3149,
) -> tuple[float, float]:
    r"""
    Transform a UTM point to latitude/longitude.

    Converts UTM easting/northing to geographic latitude and
    longitude (decimal degrees). Uses GDAL/OSR when available,
    else falls back to PROJ via ``pyproj``.

    Parameters
    ----------
    easting : float
        Easting in meters.
    northing : float
        Northing in meters.
    utm_zone : str or int
        UTM zone designator. Either a string like ``"10S"`` or
        an integer UTM code (negative for south). Ignored when
        ``epsg`` is provided.
    datum : str, default "WGS84"
        Geodetic datum name (e.g., ``"WGS84"``, ``"NAD27"``).
    epsg : int, optional
        EPSG code defining the projected CRS. When provided, it
        takes precedence over ``utm_zone`` and ``datum``. The
        default is ``3149`` for historical compatibility.

    Returns
    -------
    lat : float
        Latitude in decimal degrees, rounded to 6 decimals.
    lon : float
        Longitude in decimal degrees, rounded to 6 decimals.

    Notes
    -----
    - With GDAL, the conversion uses
      ``osr.CoordinateTransformation`` to a geographic CRS
      cloned from the projected CRS.
    - With PROJ, the conversion uses ``pyproj.Proj`` with
      ``inverse=True``.
    - When ``epsg`` refers to a non-UTM CRS, ``utm_zone`` is
      not required.
    - The function validates ``easting`` and ``northing`` can
      be cast to ``float`` and raises ``GisError`` otherwise.

    Examples
    --------
    Using an explicit zone::

        >>> project_point_utm2ll(377274.0, 3762150.0, "11S")

    Using an EPSG code (preferred when known)::

        >>> project_point_utm2ll(
        ...     377274.0, 3762150.0, utm_zone="11S", epsg=32611
        ... )

    See Also
    --------
    project_points_ll2utm
        Forward transform from geographic to UTM.

    References
    ----------
    .. [1] PROJ documentation, https://proj.org/
    .. [2] GDAL/OSR, https://gdal.org/
    """
    try:
        easting = float(easting)
    except ValueError as exc:
        raise GisError("easting is not a float") from exc
    try:
        northing = float(northing)
    except ValueError as exc:
        raise GisError("northing is not a float") from exc

    if HAS_GDAL:
        utm_cs = osr.SpatialReference()
        utm_cs.SetWellKnownGeogCS(datum)
    else:
        pp = None  # type: ignore[assignment]

    if epsg is not None:
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise RuntimeError(
                    "GDAL/OSR error code: {}".format(ogrerr)
                )
        else:
            pp = pyproj.Proj("+init=EPSG:{}".format(epsg))
    elif isinstance(utm_zone, (str, np.bytes_)):
        if isinstance(utm_zone, np.bytes_):
            utm_zone = utm_zone.decode("UTF-8")
        try:
            zone_number = int(utm_zone[:-1])
            zone_letter = utm_zone[-1]
        except ValueError as exc:
            raise ValueError(
                "Zone number '{}' is not a number".format(
                    utm_zone[:-1]
                )
            ) from exc
        is_northern = zone_letter.lower() >= "n"
    elif isinstance(utm_zone, int):
        is_northern = utm_zone >= 0
        zone_number = abs(utm_zone)
    else:
        raise NotImplementedError(
            "utm_zone type '{}' not supported".format(
                type(utm_zone).__name__
            )
        )

    if epsg is None:
        if HAS_GDAL:
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = (
                "+proj=utm +zone={} +{} +datum={}".format(
                    zone_number,
                    "north" if is_northern else "south",
                    datum,
                )
            )
            pp = pyproj.Proj(projstring)

    if HAS_GDAL:
        ll_cs = utm_cs.CloneGeogCS()
        utm2ll = osr.CoordinateTransformation(
            utm_cs, ll_cs
        ).TransformPoint
        ll_point = list(utm2ll(easting, northing, 0.0))
        lon, lat = ll_point[0], ll_point[1]
    else:
        assert pp is not None
        lon, lat = pp(easting, northing, inverse=True)

    return (round(lat, 6), round(lon, 6))


def project_points_ll2utm(
    lat,
    lon,
    datum: str = "WGS84",
    utm_zone: Optional[str] = None,
    epsg: Optional[int] = None,
):
    r"""
    Transform latitude/longitude to UTM coordinates.

    Converts geographic coordinates to UTM eastings and
    northings. Uses GDAL/OSR when available, else falls back
    to PROJ via ``pyproj``. Accepts scalars or arrays; array
    inputs must share the same shape.

    Parameters
    ----------
    lat : float, str, or array-like
        Latitude(s) in decimal degrees or sexagesimal string.
        Sexagesimal uses ``"DD:MM:SS.sss"`` format. When an
        array is given, its shape must match ``lon``.
    lon : float, str, or array-like
        Longitude(s) in decimal degrees or sexagesimal string.
        Same formatting rules as for ``lat``.
    datum : str, default "WGS84"
        Geodetic datum name. Examples: ``"WGS84"``,
        ``"NAD83"``, ``"NAD27"``. Passed to the backend CRS.
    utm_zone : str, optional
        UTM zone designator (e.g., ``"55S"``). If omitted, a
        zone is inferred from the input centroid. Ignored when
        ``epsg`` is provided.
    epsg : int, optional
        EPSG code that fully defines the projected CRS. When
        set, it takes precedence over ``utm_zone``.

    Returns
    -------
    easting : float or ndarray
        UTM easting(s) in meters. Matches the input shape for
        array inputs.
    northing : float or ndarray
        UTM northing(s) in meters. Matches the input shape for
        array inputs.
    zone : str or None
        UTM zone string (e.g., ``"55S"``) when determinable.
        May be ``None`` for CRS that are not UTM.

    Notes
    -----
    - If either ``lat`` or ``lon`` is ``None``, the function
      returns ``(None, None, None)``.
    - Array inputs are flattened for computation and then
      reshaped back to the original layout.
    - With GDAL, the transform uses
      ``osr.CoordinateTransformation``; with PROJ it uses
      ``pyproj.Proj``.
    - When ``epsg`` refers to a non-UTM projected CRS, a UTM
      zone may not be derivable; ``zone`` can be ``None``.

    Examples
    --------
    Decimal degrees::

        >>> e, n, z = project_points_ll2utm(34.05, -118.34)
        >>> round(e, 1), round(n, 1), z
        (..., ..., ...)

    Sexagesimal strings::

        >>> project_points_ll2utm("34:03:00", "-118:20:24")

    Vectorized input::

        >>> lats = [34.00, 34.05]
        >>> lons = [-118.40, -118.34]
        >>> e, n, z = project_points_ll2utm(lats, lons)

    See Also
    --------
    project_point_utm2ll
        Inverse transform from UTM to geographic coords.

    References
    ----------
    .. [1] EPSG Registry, https://epsg.org/
    .. [2] PROJ documentation, https://proj.org/
    .. [3] GDAL/OSR, https://gdal.org/
    """

    lat = np.array(lat)
    lon = np.array(lon)

    if np.shape(lat) != np.shape(lon):
        raise ValueError(
            "latitude and longitude arrays are of "
            "different lengths"
        )

    flattened = False
    llshape = np.shape(lat)
    if llshape and llshape[0] > 1:
        flattened = True
        lat = lat.flatten()
        lon = lon.flatten()

    if lat is None or lon is None:
        return (None, None, None)

    if HAS_GDAL:
        utm_cs = osr.SpatialReference()
        utm_cs.SetWellKnownGeogCS(datum)
        ll_cs = utm_cs.CloneGeogCS()
        ll_cs.ExportToPrettyWkt()
    else:
        # pyproj will be used below via `pp`
        pp = None  # type: ignore[assignment]

    if epsg is not None:
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise RuntimeError(
                    "GDAL/osgeo ogr error code: {}".format(
                        ogrerr
                    )
                )
            zone_id = utm_cs.GetUTMZone()
            if zone_id and zone_id > 0:
                utm_cs.SetUTM(abs(zone_id), zone_id > 0)
        else:
            pp = pyproj.Proj("+init=EPSG:{}".format(epsg))
    else:
        if utm_zone is not None:
            zone_number = int(utm_zone[:-1])
            hemi_letter = utm_zone[-1].upper()
            is_northern = hemi_letter >= "N"
        else:
            latc = (np.nanmax(lat) + np.nanmin(lat)) / 2.0
            lonc = (np.nanmax(lon) + np.nanmin(lon)) / 2.0
            zone_number = int(
                np.floor((lonc + 180.0) / 6.0) + 1
            )
            letter = utm_letter_designator(latc)
            is_northern = letter >= "N"
            utm_zone = "{}{}".format(zone_number, letter)

        if HAS_GDAL:
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = (
                "+proj=utm +zone={:d} +{} +datum={}".format(
                    zone_number,
                    "north" if is_northern else "south",
                    datum,
                )
            )
            pp = pyproj.Proj(projstring)

    if HAS_GDAL:
        ll2utm = osr.CoordinateTransformation(
            ll_cs, utm_cs
        ).TransformPoints
        pts = np.array(ll2utm(np.array([lon, lat]).T))
        easting = pts[:, 0]
        northing = pts[:, 1]
    else:
        assert pp is not None
        easting, northing = pp(lon, lat)

    if flattened:
        easting = np.reshape(easting, llshape)
        northing = np.reshape(northing, llshape)

    return (easting, northing, utm_zone)

def ll_to_utm(
    reference_ellipsoid: int,
    lat: float,
    lon: float
) -> Tuple[str, float, float]:
    """
    Convert latitude/longitude to UTM coordinates.
    Equations from USGS Bulletin 1532.

    Parameters
    ----------
    reference_ellipsoid : int
        ID of ellipsoid from ELLIPSOIDS (e.g. 23 for WGS-84).
    lat, lon : float
        Geographic coordinates in decimal degrees.

    Returns
    -------
    utm_zone : str
        UTM zone string (e.g. '11N').
    easting : float
    northing : float
    """
    # Get ellipsoid parameters
    try:
        ell = next(e for e in ELLIPSOIDS if e[0] == reference_ellipsoid)
    except StopIteration:
        raise GisError(f"Unknown ellipsoid ID: {reference_ellipsoid}")
        
    a = ell[_EQUATORIAL_RADIUS_IDX]
    ecc_sq = ell[_ECC_SQUARED_IDX]
    k0 = 0.9996

    # normalize longitude into [-180, 180)
    long_temp = ((lon + 180) % 360) - 180

    lat_rad = lat * DEG2RAD
    long_rad = long_temp * DEG2RAD

    # determine UTM zone number
    zone_number = int((long_temp + 180) / 6) + 1
    # Norway exception
    if 56.0 <= lat < 64.0 and 3.0 <= long_temp < 12.0:
        zone_number = 32
    # Svalbard exceptions
    if 72.0 <= lat < 84.0:
        if   0.0  <= long_temp <  9.0: zone_number = 31
        elif  9.0 <= long_temp < 21.0: zone_number = 33
        elif 21.0 <= long_temp < 33.0: zone_number = 35
        elif 33.0 <= long_temp < 42.0: zone_number = 37

    long_origin = (zone_number - 1) * 6 - 180 + 3
    long_origin_rad = long_origin * DEG2RAD

    utm_zone = f"{zone_number}{utm_letter_designator(lat)}"

    ecc_prime_sq = ecc_sq / (1 - ecc_sq)
    N = a / np.sqrt(1 - ecc_sq * np.sin(lat_rad) ** 2)
    T = np.tan(lat_rad) ** 2
    C = ecc_prime_sq * np.cos(lat_rad) ** 2
    A = np.cos(lat_rad) * (long_rad - long_origin_rad)

    # meridional arc
    M = a * (
        (1 - ecc_sq/4 - 3*ecc_sq**2/64 - 5*ecc_sq**3/256) * lat_rad
        - (3*ecc_sq/8 + 3*ecc_sq**2/32 + 45*ecc_sq**3/1024) * np.sin(2*lat_rad)
        + (15*ecc_sq**2/256 + 45*ecc_sq**3/1024) * np.sin(4*lat_rad)
        - (35*ecc_sq**3/3072) * np.sin(6*lat_rad)
    )

    # easting
    easting = (
        k0 * N * (
            A
            + (1 - T + C) * A**3 / 6
            + (5 - 18*T + T**2 + 72*C - 58*ecc_prime_sq) * A**5 / 120
        )
        + 500_000.0
    )

    # northing
    northing = k0 * (
        M + N * np.tan(lat_rad) * (
            A**2 / 2
            + (5 - T + 9*C + 4*C**2) * A**4 / 24
            + (61 - 58*T + T**2 + 600*C - 330*ecc_prime_sq) * A**6 / 720
        )
    )
    # southern hemisphere offset
    if lat < 0:
        northing += 10_000_000.0

    return utm_zone, easting, northing

def utm_to_ll(
    reference_ellipsoid: int,
    northing: float,
    easting: float,
    zone: str
) -> Tuple[float, float]:
    """
    Convert UTM coordinates to latitude/longitude.
    Inverse of ll_to_utm, equations from USGS Bulletin 1532.

    Parameters
    ----------
    reference_ellipsoid : int
        ID of ellipsoid from ELLIPSOIDS.
    northing, easting : float
        UTM coordinates in meters.
    zone : str
        UTM zone string (e.g. '11N').

    Returns
    -------
    lat, lon : float
        Geographic coordinates in decimal degrees.
    """
    # Get ellipsoid parameters
    try:
        ell = next(e for e in ELLIPSOIDS if e[0] == reference_ellipsoid)
    except StopIteration:
        raise GisError(f"Unknown ellipsoid ID: {reference_ellipsoid}")
    a = ell[_EQUATORIAL_RADIUS_IDX]
    ecc_sq = ell[_ECC_SQUARED_IDX]
    k0 = 0.9996

    # Compute e1 parameter
    e1 = (1 - np.sqrt(1 - ecc_sq)) / (1 + np.sqrt(1 - ecc_sq))

    # Remove offsets
    x = easting - 500_000.0
    y = northing

    # Determine hemisphere
    zone_letter = zone[-1]
    if zone_letter.upper() < 'N':
        # southern hemisphere
        y -= 10_000_000.0

    # Central meridian of zone
    zone_number = int(zone[:-1])
    long_origin = (zone_number - 1) * 6 - 180 + 3

    # Calculate footpoint latitude
    M = y / k0
    mu = M / (a * (1 - ecc_sq/4 - 3*ecc_sq**2/64 - 5*ecc_sq**3/256))
    phi1_rad = (
        mu
        + (3*e1/2 - 27*e1**3/32) * np.sin(2*mu)
        + (21*e1**2/16 - 55*e1**4/32) * np.sin(4*mu)
        + (151*e1**3/96) * np.sin(6*mu)
    )

    # Precompute
    ecc_prime_sq = ecc_sq / (1 - ecc_sq)
    n1 = a / np.sqrt(1 - ecc_sq * np.sin(phi1_rad)**2)
    t1 = np.tan(phi1_rad)**2
    c1 = ecc_prime_sq * np.cos(phi1_rad)**2
    r1 = a * (1 - ecc_sq) / (1 - ecc_sq * np.sin(phi1_rad)**2)**1.5
    d  = x / (n1 * k0)

    # Latitude
    lat_rad = (
        phi1_rad
        - (n1 * np.tan(phi1_rad) / r1) * (
            d**2/2
            - (5 + 3*t1 + 10*c1 - 4*c1**2 - 9*ecc_prime_sq) * d**4/24
            + (61 + 90*t1 + 298*c1 + 45*t1**2 - 252*ecc_prime_sq - 3*c1**2) * d**6/720
        )
    )
    lat = lat_rad * RAD2DEG

    # Longitude
    lon_rad = (
        d
        - (1 + 2*t1 + c1) * d**3/6
        + (5 - 2*c1 + 28*t1 - 3*c1**2 + 8*ecc_prime_sq + 24*t1**2) * d**5/120
    ) / np.cos(phi1_rad)
    lon = long_origin + lon_rad * RAD2DEG

    return lat, lon


def epsg_project(
    x,
    y,
    epsg_from,
    epsg_to,
):
    r"""
    Project coordinates between two EPSG-defined CRSs.

    Leverages :mod:`pyproj` to transform ``(x, y)`` from the
    CRS identified by ``epsg_from`` to the CRS identified by
    ``epsg_to``. Inputs may be scalars or array-like; outputs
    match the broadcasted input shape.

    Parameters
    ----------
    x, y : float or array-like
        Coordinates to transform. When arrays are provided,
        they must be broadcastable to a common shape.
    epsg_from : int
        EPSG code of the source CRS.
    epsg_to : int
        EPSG code of the destination CRS.

    Returns
    -------
    x2, y2 : float or ndarray
        Transformed coordinates in the destination CRS. Returns
        ``None`` when :mod:`pyproj` is unavailable or EPSG
        codes are missing from :data:`EPSG_DICT`.

    Notes
    -----
    - Uses ``pyproj.Transformer`` under the hood. Older
      ``pyproj.transform`` calls are avoided.
    - This helper consults :data:`EPSG_DICT` for PROJ strings.
      If an EPSG code is absent from that dictionary, the
      function logs a warning and returns ``None``.

    Examples
    --------
    >>> # WGS84 lon/lat (EPSG:4326) to Web Mercator (EPSG:3857)
    >>> epsg_project(-118.34, 34.05, 4326, 3857)  # doctest: +SKIP

    References
    ----------
    .. [1] PROJ / pyproj documentation, https://pyproj.org/
    """
    try:
        import pyproj  # type: ignore
    except Exception:  # pragma: no cover
        logger.warning("Please install 'pyproj' to use "
                       "epsg_project.")
        return None

    if epsg_from is None or epsg_to is None:
        logger.warning("Both 'epsg_from' and 'epsg_to' must be "
                       "provided.")
        return None

    try:
        src = pyproj.CRS.from_user_input(EPSG_DICT[epsg_from])
        dst = pyproj.CRS.from_user_input(EPSG_DICT[epsg_to])
    except KeyError:
        logger.warning("EPSG code not in EPSG_DICT: from=%r to=%r",
                       epsg_from, epsg_to)
        return None

    transformer = pyproj.Transformer.from_crs(
        src,
        dst,
        always_xy=True,
    )
    x2, y2 = transformer.transform(x, y)
    return x2, y2


def utm_wgs84_conv(
    lat: float,
    lon: float,
) -> tuple[float, float, int, str]:
    r"""
    Convert WGS84 (lat, lon) to UTM and verify round-trip.

    Uses the :mod:`utm` package to convert a WGS84 geographic
    coordinate to UTM ``(easting, northing, zone, letter)``.
    Then converts back to lat/lon and checks numerical
    consistency.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees (WGS84).
    lon : float
        Longitude in decimal degrees (WGS84).

    Returns
    -------
    e : float
        UTM easting (meters).
    n : float
        UTM northing (meters).
    zone : int
        UTM zone number in ``[1, 60]``.
    letter : str
        UTM latitude band letter.

    Notes
    -----
    - Round-trip consistency is checked with a tolerance of
      ``1e-10`` degrees. If exceeded, a warning is logged.
    - Requires the external package ``utm`` (``pip install
      utm``).

    Examples
    --------
    >>> e, n, z, L = utm_wgs84_conv(34.05, -118.34)  # doctest: +SKIP
    >>> z, L  # doctest: +SKIP
    (11, 'S')

    References
    ----------
    .. [1] Turbo87/utm, https://github.com/Turbo87/utm
    """
    try:
        import utm  # type: ignore
    except Exception:  # pragma: no cover
        logger.warning("Please install 'utm' (pip install utm).")
        raise

    e, n, z, L = utm.from_latlon(lat, lon)
    new_lat, new_lon = utm.to_latlon(e, n, z, L)

    tol = 1e-10
    if abs(lat - new_lat) > tol:
        logger.warning("lat round-trip mismatch: in=%r out=%r",
                       lat, new_lat)
    if abs(lon - new_lon) > tol:
        logger.warning("lon round-trip mismatch: in=%r out=%r",
                       lon, new_lon)

    return e, n, z, L

@Deprecated(
    "Deprecated; will be removed in a future release. "
    "Use 'project_point_utm2ll' instead (GDAL/pyproj-backed)."
)
def transform_utm_to_ll(
    easting: float,
    northing: float,
    zone,
    reference_ellipsoid: str = "WGS84",
) -> Tuple[float, float, float]:
    r"""
    Transform a UTM point to (lon, lat, alt) using GDAL/OSR.

    Parameters
    ----------
    easting : float
        UTM easting in meters.
    northing : float
        UTM northing in meters.
    zone : str or int
        Zone string like ``'11N'`` / ``'55S'`` or a signed
        integer UTM code (negative for south).
    reference_ellipsoid : str, default "WGS84"
        Well-known geographic CRS name for GDAL.

    Returns
    -------
    lon : float
        Longitude in decimal degrees.
    lat : float
        Latitude in decimal degrees.
    alt : float
        Altitude (meters), usually ``0.0`` for 2D inputs.

    Raises
    ------
    GDALMissingError
        If GDAL/OSR is not available.
    ValueError
        If ``zone`` cannot be parsed.

    Notes
    -----
    - This function is **deprecated**. Prefer
      :func:`project_point_utm2ll`.
    """
    if not HAS_GDAL:
        raise GDALMissingError(
            "GDAL/OSR is required for 'transform_utm_to_ll'. "
            "Use 'project_point_utm2ll' for a GDAL/pyproj path."
        )

    utm_sr = osr.SpatialReference()
    utm_sr.SetWellKnownGeogCS(reference_ellipsoid)

    # parse zone
    if isinstance(zone, int):
        zone_num = abs(zone)
        is_north = zone >= 0
    elif isinstance(zone, str):
        try:
            zone_num = int(zone[:-1])
            zone_letter = zone[-1]
        except Exception as exc:  # noqa: BLE001
            raise ValueError(
                f"Invalid UTM zone string: {zone!r}"
            ) from exc
        is_north = zone_letter.lower() >= "n"
    else:
        raise ValueError(
            f"Unsupported zone type: {type(zone).__name__}"
        )

    utm_sr.SetUTM(zone_num, is_north)

    ll_sr = utm_sr.CloneGeogCS()
    to_ll = osr.CoordinateTransformation(utm_sr, ll_sr)

    lon, lat, alt = to_ll.TransformPoint(easting, northing, 0.0)
    return (lon, lat, alt)


@Deprecated(
    "Deprecated; will be removed in a future release. "
    "Use 'project_point_ll2utm' instead (GDAL/pyproj-backed)."
)
def transform_ll_to_utm(
    lon: float,
    lat: float,
    reference_ellipsoid: str = "WGS84",
):
    r"""
    Transform a geographic point (lon, lat) to UTM using GDAL.

    Parameters
    ----------
    lon : float
        Longitude in decimal degrees.
    lat : float
        Latitude in decimal degrees.
    reference_ellipsoid : str, default "WGS84"
        Well-known geographic CRS name for GDAL.

    Returns
    -------
    utm_sr : osr.SpatialReference
        UTM spatial reference describing the destination CRS.
    utm_point : tuple[float, float, float]
        ``(easting, northing, altitude)`` in meters.

    Raises
    ------
    GDALMissingError
        If GDAL/OSR is not available.

    Notes
    -----
    - This function is **deprecated**. Prefer
      :func:`project_point_ll2utm`, which can use PROJ when
      GDAL is unavailable.
    """
    if not HAS_GDAL:
        raise GDALMissingError(
            "GDAL/OSR is required for 'transform_ll_to_utm'. "
            "Use 'project_point_ll2utm' for a GDAL/pyproj path."
        )

    def _utm_zone_number(lon_deg: float) -> int:
        return int(1 + (lon_deg + 180.0) / 6.0)

    def _is_northern(lat_deg: float) -> bool:
        return lat_deg >= 0.0

    utm_sr = osr.SpatialReference()
    utm_sr.SetWellKnownGeogCS(reference_ellipsoid)
    utm_sr.SetUTM(_utm_zone_number(lon), _is_northern(lat))

    ll_sr = utm_sr.CloneGeogCS()
    to_utm = osr.CoordinateTransformation(ll_sr, utm_sr)

    easting, northing, alt = to_utm.TransformPoint(lon, lat, 0.0)
    return utm_sr, (easting, northing, alt)
