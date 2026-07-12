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
from __future__ import annotations

import re
from collections.abc import Sequence
from typing import (
    Any,
    Literal,
    overload,
)

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from ..decorators import Deprecated
from ..log.logger import get_logger
from .config import (
    EPSG_DICT,
    HAS_GDAL,
    ElevationAPIConfig,
    GDALMissingError,
    GisError,
)
from .constants import (
    _ECC_SQUARED_IDX,
    _EQUATORIAL_RADIUS_IDX,
    DEG2RAD,
    ELLIPSOIDS,
    RAD2DEG,
    utm_letter_designator,
)

if HAS_GDAL:
    from osgeo import osr
    from osgeo.ogr import OGRERR_NONE
else:
    import pyproj  # type: ignore

logger = get_logger(__name__)


Arr1D = NDArray[np.float64]
Obj1D = NDArray[object]

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
    "to_utm",
    "to_ll",
    "normalize_lat_lon",
    'GisError'
]

_HEMI_RE = re.compile(r"\s*([NSEW])\s*$", re.IGNORECASE)
_HEMI_SIGN = {"N": 1, "E": 1, "S": -1, "W": -1}


@overload
def normalize_lat_lon(
    a: float | str,
    b: float | str,
    *,
    assume: Literal["lonlat", "latlon", "auto"] = ...,
    error: Literal["ignore", "raise"] = ...,
    clip: bool = ...,
) -> tuple[float, float]:
    ...


@overload
def normalize_lat_lon(
    a: Sequence[float | str],
    b: Sequence[float | str],
    *,
    assume: Literal["lonlat", "latlon", "auto"] = ...,
    error: Literal["ignore", "raise"] = ...,
    clip: bool = ...,
) -> tuple[np.ndarray, np.ndarray]:
    ...

def normalize_lat_lon(
    a: Any,
    b: Any,
    *,
    assume: Literal["lonlat", "latlon", "auto"] = "lonlat",
    error: Literal["ignore", "raise"] = "ignore",
    clip: bool = False,
):
    r"""
    Resolve input pair to ``(lat, lon)`` regardless of order.

    Accepts values in either legacy order (lat, lon) or the
    new expected order (lon, lat). Returns a tuple in the
    canonical order ``(lat, lon)`` and avoids common
    ``|Latitude| > 90`` errors by swapping when clear.

    Parameters
    ----------
    a, b : float or str, or 1-D sequences
        Coordinate components. Strings may be decimal or DMS
        (``"DD:MM:SS"``). Sequences must be same length.
    assume : {"lonlat","latlon","auto"}, default "lonlat"
        Tie-break for ambiguous pairs (both ``|val| <= 90``).

        - ``"lonlat"``: treat input as (lon, lat).
        - ``"latlon"``: treat input as (lat, lon).
        - ``"auto"``: prefer (lon, lat) when ambiguous.

    error : {"ignore","raise"}, default "ignore"
        On impossible pairs (both ``|val| > 90`` or both
        ``|val| > 180``), either raise or return ``nan``s.
    clip : bool, default False
        If ``True``, clip lat to [-90, 90] and lon to
        [-180, 180] when slightly out-of-range.

    Returns
    -------
    lat, lon : float or ndarray
        Coordinates in canonical order.

    Notes
    -----
    Heuristics:

    - If one value is in (90, 180] and the other in [-90, 90],
      the former is lon and the latter is lat.
    - If both are ``|val| <= 90``, use ``assume``.
    - Values with ``|val| > 180`` are invalid unless ``clip``.
    """
    def _to_float(v):
        if v is None or v == "None":
            return np.nan
        if isinstance(v, str):
            try:
                return float(v)
            except Exception:
                try:
                    return convert_position_str2float(v)  # type: ignore  # noqa: E501
                except Exception:
                    return np.nan
        try:
            return float(v)
        except Exception:
            return np.nan

    def _coerce_pair(x, y):
        xv = _to_float(x)
        yv = _to_float(y)

        ax = abs(xv)
        ay = abs(yv)

        # invalid domain checks
        if ax > 1e9 or ay > 1e9:
            xv = np.nan
            yv = np.nan

        # quick invalid > 180
        if ax > 180 or ay > 180:
            if clip:
                xv = np.clip(xv, -180.0, 180.0)
                yv = np.clip(yv, -180.0, 180.0)
                ax = abs(xv)
                ay = abs(yv)
            else:
                if error == "raise":
                    raise ValueError("Values exceed 180.")
                return (np.nan, np.nan)

        # decisive cases
        if (90 < ax <= 180) and (ay <= 90):
            lon = xv
            lat = yv
            return (lat, lon)
        if (90 < ay <= 180) and (ax <= 90):
            lon = yv
            lat = xv
            return (lat, lon)

        # both within 90 → ambiguous
        if ax <= 90 and ay <= 90:
            if assume == "latlon":
                lat, lon = xv, yv
            else:
                # "lonlat" and "auto" prefer lon,lat input
                lat, lon = yv, xv
            return (lat, lon)

        # one slightly out of 90 but not decisive
        if clip:
            xv = np.clip(xv, -180.0, 180.0)
            yv = np.clip(yv, -180.0, 180.0)
            ax = abs(xv)
            ay = abs(yv)

        # fallback: try assume
        if assume == "latlon":
            lat, lon = xv, yv
        else:
            lat, lon = yv, xv

        # final range checks
        if abs(lat) > 90 or abs(lon) > 180:
            if error == "raise":
                raise ValueError("Unresolvable pair.")
            return (np.nan, np.nan)

        return (lat, lon)

    # scalar vs vector handling
    if np.isscalar(a) and np.isscalar(b):
        return _coerce_pair(a, b)

    a_arr = np.asarray(a, dtype=object)
    b_arr = np.asarray(b, dtype=object)

    if a_arr.shape != b_arr.shape:
        raise ValueError("Shapes of a and b must match.")

    lat_out = np.empty_like(a_arr, dtype=float)
    lon_out = np.empty_like(b_arr, dtype=float)

    it = np.nditer(
        [a_arr, b_arr, lat_out, lon_out],
        flags=["multi_index", "refs_ok"],
        op_flags=[["readonly"], ["readonly"],
                  ["writeonly"], ["writeonly"]],
    )
    for xa, xb, yl, yo in it:
        lat, lon = _coerce_pair(xa.item(), xb.item())
        yl[...] = lat
        yo[...] = lon

    return lat_out, lon_out

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
    Accepts optional trailing hemisphere letter:
    'N','S','E','W' (e.g., '26:00:00N', '10:00:00 W').

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

    # Fast path: plain decimal string/number
    try:
        return float(position_str)
    except (TypeError, ValueError):
        pass

    # Normalize and peel optional hemisphere at end
    s = str(position_str).strip().strip('"').strip("'")
    hemi_sign: int | None = None
    m = _HEMI_RE.search(s)
    if m:
        hemi = m.group(1).upper()
        hemi_sign = _HEMI_SIGN[hemi]
        s = s[: m.start()].rstrip()

    parts = s.split(":")
    if len(parts) != 3:
        msg = f"Invalid DMS, expected 'DD:MM:SS': {position_str}"
        raise ValueError(msg)

    # Degrees may carry sign; minutes/seconds must be non-negative
    deg = float(parts[0])
    try:
        minutes = abs(float(parts[1]))
    except ValueError as exc:
        raise ValueError(f"Invalid minutes: {parts[1]!r}") from exc
    try:
        seconds = abs(float(parts[2]))
    except ValueError as exc: # noqa
        # Handle cases like '00N' already stripped above, but be safe
        # by removing any trailing non-numeric tokens.
        t = re.sub(r"[^\d.+\-eE]", "", parts[2])
        try:
            seconds = abs(float(t))
        except Exception as exc2:
            raise ValueError(
                f"Invalid seconds: {parts[2]!r}") from exc2

    # Sign: hemisphere (if any) overrides degree sign. Capture it
    # before rollover so carries apply to the magnitude, not the
    # signed value ('-118:60:00' must become -119, not -117).
    sign = -1 if deg < 0 else 1
    if hemi_sign is not None:
        sign = hemi_sign
    deg = abs(deg)

    # Rollover SS->MM then MM->DD, then validate the remainders
    # (validating first would reject '10:60:00', which the rollover
    # exists to normalize to 11:00:00).
    minutes, seconds = _rollover_dms(minutes, seconds)
    deg, minutes = _rollover_dms(deg, minutes)
    minutes = _assert_minutes(minutes)
    seconds = _assert_seconds(seconds)

    return sign * (deg + minutes / 60.0 + seconds / 3600.0)


def assert_lat_value(latitude) -> float | None:
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
        raise ValueError(
            f"|Latitude| > 90, unacceptable!: {lat_value!r}"
        )
    return lat_value


def assert_lon_value(longitude) -> float | None:
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
            f"{elevation} is not a number, setting elevation to 0"
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

    position_str = f"{sign * int(deg)}:{int(minutes):02.0f}:{sec:05.2f}"
    return position_str


@Deprecated(
    "GDAL SpatialReference → UTM string is deprecated; "
    "use 'get_utm_zone' for standard UTM formatting."
)
def get_utm_string_from_sr(
        spatial_ref: osr.SpatialReference
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
) -> tuple[int, bool, str]:
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
) -> int | None:
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
) -> int | None:
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

@overload
def to_ll(
    easting: Any,
    northing: Any,
    zone: Any | None = ...,
    data: Any | None = ...,
    *,
    datum: str = ...,
    epsg: int | None = ...,
    as_frame: Literal[True],
) -> pd.DataFrame: ...
@overload
def to_ll(
    easting: float,
    northing: float,
    zone: Any | None = ...,
    data: Any | None = ...,
    *,
    datum: str = ...,
    epsg: int | None = ...,
    as_frame: Literal[False] = ...,
) -> tuple[float, float]: ...
@overload
def to_ll(
    easting: Sequence[float] | Arr1D,
    northing: Sequence[float] | Arr1D,
    zone: Sequence[Any] | Obj1D | None = ...,
    data: Any | None = ...,
    *,
    datum: str = ...,
    epsg: int | None = ...,
    as_frame: Literal[False] = ...,
) -> tuple[Arr1D, Arr1D]: ...

def to_ll(
    easting,
    northing,
    zone: Any | None = None,
    data: Any | None = None,
    *,
    datum: str = "WGS84",
    epsg: int | None = None,
    as_frame: bool = False,
) -> Any:
    r"""
    Convert UTM coordinates to geographic (lat, lon) using
    :func:`project_point_utm2ll`.

    Accepts flexible inputs (scalars, arrays, Series, or
    column names referencing a provided DataFrame). Returns a
    tuple or a DataFrame with original and converted fields.

    Parameters
    ----------
    easting, northing : float, array-like, or pd.Series
        UTM coordinates in meters. When strings are given,
        they are treated as column names in ``data``.
    zone : str, int, array-like, or pd.Series, optional
        UTM zone designator (e.g., ``"11S"``) or signed UTM
        code (negative for south). If strings are given, they
        are treated as column names in ``data``. May be
        omitted when ``epsg`` is provided.
    data : pd.DataFrame, optional
        Source of columns when any of the inputs are strings.
    datum : str, default "WGS84"
        Datum name for GDAL-backed paths.
    epsg : int, optional
        EPSG code that defines the projected CRS. When set,
        it takes precedence and ``zone`` is not required.
    as_frame : bool, default False
        If ``True``, return a DataFrame with columns
        ``['easting','northing','zone','lat','lon']``.
        Otherwise, return ``(lat, lon)`` as scalars or arrays.

    Returns
    -------
    DataFrame or tuple
        A DataFrame when ``as_frame=True``; otherwise a
        tuple ``(lat, lon)`` matching the input shape.

    Notes
    -----
    - Index is preserved when inputs originate from a Series
      or from ``data[<column>]`` with ``as_frame=True``.
    - When arrays are given and ``epsg`` is used, the same
      EPSG applies to all rows.
    - ``zone`` may be a scalar or an array/Series matching
      the shape of ``easting`` and ``northing``.

    Examples
    --------
    Scalars::

        >>> to_ll(377274.0, 3762150.0, "11S")  # doctest: +SKIP
        (34.05, -118.34)

    From a DataFrame::

        >>> # df has 'E', 'N', 'zone' columns
        >>> # to_ll('E','N','zone', data=df, as_frame=True)
        ...  # doctest: +SKIP
    """
    # --------------------------------------------------------------
    # Resolve inputs to 1-D numpy arrays; keep optional index
    # --------------------------------------------------------------

    idx =None
    def _extract(v, name: str):
        nonlocal idx
        if isinstance(v, str):
            if data is None:
                raise ValueError(
                    f"'{name}' is a column name but "
                    "'data' is not provided."
                )
            col = data[v]

            idx = getattr(col, "index", None)
            return np.asarray(col.to_numpy())
        if hasattr(v, "to_numpy"):  # pandas Series
            idx = getattr(v, "index", None)
            return np.asarray(v.to_numpy())
        if np.isscalar(v):
            return np.asarray([v])
        return np.asarray(v)

    e_arr = _extract(easting, "easting")
    n_arr = _extract(northing, "northing")

    if e_arr.shape != n_arr.shape:
        raise ValueError(
            "easting and northing shapes differ: "
            f"{e_arr.shape} != {n_arr.shape}"
        )

    if zone is not None:
        if isinstance (zone, str):
            z_arr = np.asarray([zone], dtype=object)
        else:
            z_arr = _extract(zone, "zone")
        if z_arr.size == 1 and e_arr.size > 1:
            z_arr = np.full(e_arr.shape, z_arr.item(), dtype=object)
        if z_arr.shape != e_arr.shape:
            raise ValueError(
                "zone shape must match easting/northing or be "
                "scalar."
            )
    else:
        # allow None when epsg is provided
        z_arr = np.full(e_arr.shape, None, dtype=object)

    # --------------------------------------------------------------
    # Compute lat/lon element-wise (project_point_utm2ll works
    # on scalars). Keep dtype as float64 for arrays.
    # --------------------------------------------------------------
    lat_list = []
    lon_list = []

    for i in range(e_arr.size):
        zval = z_arr.ravel()[i]
        lat_i, lon_i = project_point_utm2ll(  # type: ignore[name-defined]  # noqa: E501
            float(e_arr.ravel()[i]),
            float(n_arr.ravel()[i]),
            zval,
            datum=datum,
            epsg=epsg,
        )
        lat_list.append(lat_i)
        lon_list.append(lon_i)

    lat_arr = np.asarray(lat_list, dtype=float).reshape(e_arr.shape)
    lon_arr = np.asarray(lon_list, dtype=float).reshape(n_arr.shape)

    # --------------------------------------------------------------
    # Return per user request
    # --------------------------------------------------------------
    if as_frame:
        def _orig(v, arr):
            if np.isscalar(v):
                return np.asarray([v])
            if isinstance(v, str) and data is not None:
                return np.asarray(data[v].to_numpy())
            if hasattr(v, "to_numpy"):
                return np.asarray(v.to_numpy())
            return np.asarray(arr)

        e_col = _orig(easting, e_arr)
        n_col = _orig(northing, n_arr)
        z_col = (
            _orig(zone, z_arr)
            if zone is not None
            else np.asarray([None] * e_arr.size, dtype=object)
        )

        df = pd.DataFrame(
            {
                "easting": e_col,
                "northing": n_col,
                "zone": z_col,
                "lat": lat_arr,
                "lon": lon_arr,
            }
        )
        if idx is not None and len(df) == len(idx):
            df.index = idx
        return df

    # scalar squeeze if the inputs were scalars
    if np.isscalar(easting) and np.isscalar(northing) and (
        np.isscalar(zone) or zone is None
    ):
        return (lat_arr.item(), lon_arr.item())
    return (lat_arr, lon_arr)

# ---------- overloads for to_utm ----------

@overload
def to_utm(
    lat: Any,
    lon: Any,
    data: Any | None = ...,
    *,
    datum: str = ...,
    utm_zone: str | None = ...,
    epsg: int | None = ...,
    as_frame: Literal[True],
) -> pd.DataFrame: ...
@overload
def to_utm(
    lat: float,
    lon: float,
    data: Any | None = ...,
    *,
    datum: str = ...,
    utm_zone: str | None = ...,
    epsg: int | None = ...,
    as_frame: Literal[False] = ...,
) -> tuple[float, float, str | None]: ...
@overload
def to_utm(
    lat: Sequence[Any] | Arr1D,
    lon: Sequence[Any] | Arr1D,
    data: Any | None = ...,
    *,
    datum: str = ...,
    utm_zone: str | None = ...,
    epsg: int | None = ...,
    as_frame: Literal[False] = ...,
) -> tuple[Arr1D, Arr1D, Obj1D]: ...

def to_utm(
    lat: str | float | Sequence,
    lon: str | float | Sequence,
    data: Any | None = None,
    *,
    datum: str = "WGS84",
    utm_zone: str | None = None,
    epsg: int | None = None,
    as_frame: bool = False,
) -> Any:
    r"""
    Convert geographic coordinates to UTM using
    :func:`project_point_ll2utm`.

    Accepts flexible latitude/longitude inputs (scalars,
    arrays, pandas Series, or column names referencing a
    passed ``data`` frame). Returns either arrays/tuples or a
    ``DataFrame`` with original and UTM fields.

    Parameters
    ----------
    lat, lon : float, str, array-like, or pd.Series
        Geographic coordinates in decimal degrees or DMS
        strings. If a string is given, it is treated as a
        column name and ``data`` must be provided.
    data : pd.DataFrame, optional
        Source of columns when ``lat`` or ``lon`` are strings.
    datum : str, default "WGS84"
        Datum name or EPSG (int) accepted by the underlying
        transformer.
    utm_zone : str, optional
        UTM zone (e.g., ``"11S"``). If absent, it is inferred
        from the input centroid (or per-point for scalar).
    epsg : int, optional
        EPSG code for the projected CRS. Overrides
        ``utm_zone`` when provided.
    as_frame : bool, default False
        If ``True``, return a pandas ``DataFrame`` with
        columns ``['lat','lon','easting','northing','zone']``.
        If ``False``, return ``(easting, northing, zone)`` as
        arrays (or scalars for scalar inputs).

    Returns
    -------
    DataFrame or tuple
        When ``as_frame=True``, a frame with original and UTM
        columns. Otherwise a tuple ``(e, n, zone)`` where each
        element matches the input shape (scalar or 1-D array).

    Notes
    -----
    - Parsing of DMS strings and range checks are handled by
      the underlying projection helper.
    - For array inputs, this function preserves order and
      length; indices are preserved when built from
      ``data[<col>]`` and ``as_frame=True``.

    Examples
    --------
    Scalar input::

        >>> to_utm(34.05, -118.34)  # doctest: +SKIP
        (..., ..., '11S')

    From a DataFrame::

        >>> # df has columns 'lat_deg' and 'lon_deg'
        >>> # to_utm('lat_deg','lon_deg', data=df, as_frame=True)
        ...  # doctest: +SKIP
    """
    # ------------------------------------------------------------------
    # Resolve input into 1-D arrays and an optional index for frames
    # ------------------------------------------------------------------
    idx = None
    def _extract(v, name: str):
        nonlocal idx
        if isinstance(v, str):
            if data is None:
                raise ValueError(
                    f"'{name}' is a column name but 'data' "
                    "is not provided."
                )
            col = data[v]

            idx = getattr(col, "index", None)
            return np.asarray(col.to_numpy())
        if hasattr(v, "to_numpy"):  # pandas Series
            arr = v.to_numpy()  # type: ignore[attr-defined]
            idx = getattr(v, "index", None)
            return np.asarray(arr)
        if np.isscalar(v):
            return np.asarray([v])
        return np.asarray(v)

    lat_arr = _extract(lat, "lat")
    lon_arr = _extract(lon, "lon")

    if lat_arr.shape != lon_arr.shape:
        raise ValueError(
            "lat and lon shapes differ: "
            f"{lat_arr.shape} != {lon_arr.shape}"
        )
    # ------------------------------------------------------------------
    # Delegate to projection helper
    # - project_point_ll2utm supports scalar or array inputs and
    #   will return either a tuple (scalar) or a record array.
    # ------------------------------------------------------------------
    res = project_point_ll2utm(  # type: ignore[name-defined]
        lat_arr,
        lon_arr,
        datum=datum,
        utm_zone=utm_zone,
        epsg=epsg,
    )

    # Normalize output to arrays
    if isinstance(res, tuple):
        e, n, z = res
        e_arr = np.asarray([e])
        n_arr = np.asarray([n])
        z_arr = np.asarray([z], dtype=object)
    else:
        # expected a recarray with named fields
        e_arr = np.asarray(res["easting"])
        n_arr = np.asarray(res["northing"])
        z_arr = np.asarray(res["utm_zone"], dtype=object)

    # Return per user request
    if as_frame:
        # Build lat/lon columns reflecting original inputs
        if np.isscalar(lat):
            lat_col = np.asarray([lat])
        elif isinstance(lat, str) and data is not None:
            lat_col = np.asarray(data[lat].to_numpy())
        elif hasattr(lat, "to_numpy"):
            lat_col = np.asarray(lat.to_numpy())
        else:
            lat_col = np.asarray(lat_arr)

        if np.isscalar(lon):
            lon_col = np.asarray([lon])
        elif isinstance(lon, str) and data is not None:
            lon_col = np.asarray(data[lon].to_numpy())
        elif hasattr(lon, "to_numpy"):
            lon_col = np.asarray(lon.to_numpy())
        else:
            lon_col = np.asarray(lon_arr)

        frame = pd.DataFrame(
            {
                "lat": lat_col,
                "lon": lon_col,
                "easting": e_arr,
                "northing": n_arr,
                "zone": z_arr,
            }
        )
        if idx is not None and len(frame) == len(idx):
            frame.index = idx
        return frame

    # arrays/scalars: squeeze to scalar if input was scalar
    if np.isscalar(lat) and np.isscalar(lon):
        return (e_arr[0], n_arr[0], z_arr[0])
    return (e_arr, n_arr, z_arr)

def project_point_ll2utm(
    lat,
    lon,
    datum: str = "WGS84",
    utm_zone: str | None = None,
    epsg: int | None = None,
):
    r"""
    Transform one geographic point to UTM coordinates.

    Accepts a single point or 1-D arrays/Series of points in
    decimal degrees (or DMS strings) and returns UTM easting,
    northing, and zone. Uses GDAL/OSR when available; falls
    back to PROJ via ``pyproj`` with explicit XY order.

    Parameters
    ----------
    lat : float or str, array-like, or pd.Series
        Latitude in decimal degrees or sexagesimal string
        (``"DD:MM:SS.sss"``).
    lon : float, str, array-like, or pd.Series
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
    tuple or np.recarray
        For scalar input: ``(easting, northing, zone)``.
        For array input: record array with fields
        ``('easting','northing','elev','utm_zone')``.

        - easting : float or None
             UTM easting in meters, or ``None`` when inputs are
             ``None``.
        - northing : float or None
             UTM northing in meters, or ``None`` when inputs are
             ``None``.
        - zone : str or None
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

    Notes
    -----
    - Enforces traditional GIS order (lon, lat) in the GDAL
      path to avoid invalid-latitude errors with GDAL ≥ 3.
    - The ``pyproj`` path uses a Transformer with
      ``always_xy=True`` for the same reason.

    Examples
    --------
    >>> e, n, z = project_point_ll2utm(34.05, -118.34)
    >>> z
    '11S'
    """
    if lat is None or lon is None:
        return (None, None, None)

    # normalize to 1-D arrays; accept scalars/iterables
    if np.iterable(lat) and np.iterable(lon):
        lat_arr = np.asarray([assert_lat_value(v) for v in lat])
        lon_arr = np.asarray([assert_lon_value(v) for v in lon])
        if lat_arr.size != lon_arr.size:
            raise ValueError("lat and lon sizes differ.")
    else:
        lat_arr = np.asarray([assert_lat_value(lat)])
        lon_arr = np.asarray([assert_lon_value(lon)])

    n = lat_arr.size
    out = np.zeros(
        n,
        dtype=[
            ("easting", np.float64),
            ("northing", np.float64),
            ("elev", np.float64),
            ("utm_zone", "U4"),
        ],
    )

    if HAS_GDAL:
        # build source geographic SRS
        ll_cs = osr.SpatialReference()
        if isinstance(datum, int):
            ogrerr = ll_cs.ImportFromEPSG(datum)
        else:
            ogrerr = ll_cs.SetWellKnownGeogCS(str(datum))
        if ogrerr != OGRERR_NONE:
            raise GisError(
                f"GDAL/OSR error code: {ogrerr}"
            )

        for i in range(n):
            la = float(lat_arr[i])
            lo = float(lon_arr[i])

            utm_cs = osr.SpatialReference()
            if isinstance(epsg, int):
                ogrerr = utm_cs.ImportFromEPSG(epsg)
                if ogrerr != OGRERR_NONE:
                    raise GisError(
                        f"GDAL/OSR error code: {ogrerr}"
                    )
                # zone string for info (derive if missing)
                _, _, zone_str = get_utm_zone(la, lo)
            else:
                # tie projected CRS to same geographic part
                ogrerr = utm_cs.CopyGeogCSFrom(ll_cs)
                if ogrerr != OGRERR_NONE:
                    raise GisError(
                        f"GDAL/OSR error code: {ogrerr}"
                    )
                if (
                    utm_zone is None
                    or not isinstance(utm_zone, str)
                    or utm_zone.lower() == "none"
                ):
                    znum, znorth, zone_str = get_utm_zone(la, lo)
                else:
                    znum = int(utm_zone[:-1])
                    znorth = utm_zone[-1].lower() > "n"
                    zone_str = utm_zone
                utm_cs.SetUTM(znum, znorth)

            # enforce (lon, lat) order in GDAL ≥ 3
            if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):
                ll_cs.SetAxisMappingStrategy(
                    osr.OAMS_TRADITIONAL_GIS_ORDER
                )
                utm_cs.SetAxisMappingStrategy(
                    osr.OAMS_TRADITIONAL_GIS_ORDER
                )

            xform = osr.CoordinateTransformation(ll_cs, utm_cs)
            x, y, z = xform.TransformPoint(lo, la, 0.0)

            out["easting"][i] = x
            out["northing"][i] = y
            out["elev"][i] = z
            out["utm_zone"][i] = zone_str
    else:
        # pyproj fallback with explicit (lon, lat)
        import pyproj  # type: ignore

        # source CRS from datum (string or epsg-int both ok)
        src = pyproj.CRS.from_user_input(datum)

        for i in range(n):
            la = float(lat_arr[i])
            lo = float(lon_arr[i])

            if isinstance(epsg, int):
                dst = pyproj.CRS.from_epsg(epsg)
                _, _, zone_str = get_utm_zone(la, lo)
            else:
                if (
                    utm_zone is None
                    or not isinstance(utm_zone, str)
                    or utm_zone.lower() == "none"
                ):
                    znum, znorth, zone_str = get_utm_zone(la, lo)
                else:
                    znum = int(utm_zone[:-1])
                    znorth = utm_zone[-1].lower() > "n"
                    zone_str = utm_zone
                proj4 = (
                    "+proj=utm +zone={} +{} +datum={}".format(
                        znum,
                        "north" if znorth else "south",
                        datum,
                    )
                )
                dst = pyproj.CRS.from_string(proj4)

            tr = pyproj.Transformer.from_crs(
                src,
                dst,
                always_xy=True,
            )
            x, y = tr.transform(lo, la)

            out["easting"][i] = x
            out["northing"][i] = y
            out["elev"][i] = 0.0
            out["utm_zone"][i] = zone_str

    if n == 1:
        return (
            out["easting"][0],
            out["northing"][0],
            out["utm_zone"][0],
        )
    return np.rec.array(out)

def project_point_utm2ll(
    easting: float,
    northing: float,
    utm_zone,
    datum: str = "WGS84",
    epsg: int | None = 3149,
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
                    f"GDAL/OSR error code: {ogrerr}"
                )
        else:
            pp = pyproj.Proj(f"+init=EPSG:{epsg}")
    elif isinstance(utm_zone, (str, np.bytes_)):
        if isinstance(utm_zone, np.bytes_):
            utm_zone = utm_zone.decode("UTF-8")
        try:
            zone_number = int(utm_zone[:-1])
            zone_letter = utm_zone[-1]
        except ValueError as exc:
            raise ValueError(
                f"Zone number '{utm_zone[:-1]}' is not a number"
            ) from exc
        is_northern = zone_letter.lower() >= "n"
    elif isinstance(utm_zone, int):
        is_northern = utm_zone >= 0
        zone_number = abs(utm_zone)
    else:
        raise NotImplementedError(
            f"utm_zone type '{type(utm_zone).__name__}' not supported"
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
    # normalize_lat_lon resolve the classic axis order issue often
    # seen when interacting with different GIS libraries.
    lat, lon = normalize_lat_lon(lon, lat, assume="latlon")
    return (round(lat, 6), round(lon, 6))


def project_points_ll2utm(
    lat,
    lon,
    datum: str = "WGS84",
    utm_zone: str | None = None,
    epsg: int | None = None,
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
                    f"GDAL/osgeo ogr error code: {ogrerr}"
                )
            zone_id = utm_cs.GetUTMZone()
            if zone_id and zone_id > 0:
                utm_cs.SetUTM(abs(zone_id), zone_id > 0)
        else:
            pp = pyproj.Proj(f"+init=EPSG:{epsg}")
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
            utm_zone = f"{zone_number}{letter}"

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
) -> tuple[str, float, float]:
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
) -> tuple[float, float]:
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
) -> tuple[float, float, float]:
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


def _resolve_api_name(name: str | None) -> str:
    """
    Resolves the API name from a full name or a shorthand symbol.
    """
    # If no name is given, use the default from the config.
    if name is None:
        return ElevationAPIConfig.DEFAULT_API

    # Map shorthand symbols to their corresponding API names.
    symbol_map = {
        ',': 'open_meteo',       # Comma for comma_separated
        '|': 'open_topo_data',   # Pipe for pipe_separated
        'indiv': 'usgs_ned',
    }

    # First, check if the provided name is a direct valid key.
    if name in ElevationAPIConfig.APIS:
        return name

    # Next, check if it's a recognized shorthand symbol.
    if name in symbol_map:
        return symbol_map[name]

    # If the name is unrecognized, log a warning and fallback.
    logger.warning(
        f"Unknown API name or symbol: '{name}'. "
        f"Falling back to default API: "
        f"'{ElevationAPIConfig.DEFAULT_API}'."
    )
    return ElevationAPIConfig.DEFAULT_API

def _extract_value_from_nested_dict(
    data: dict, key_path: str
) -> Any:
    r"""Extract a value from a nested dict using a dot-path.

    This helper traverses a dictionary structure using a
    dot-separated string to locate and return a nested value.

    Parameters
    ----------
    data : dict
        The dictionary to search within.
    key_path : str
        A dot-separated path representing the keys to follow.
        For example, 'results.elevation' will access
        ``data['results']['elevation']``.

    Returns
    -------
    Any
        The value found at the specified path, or ``None`` if any
        key along the path does not exist.
    """
    # Split the path into individual keys.
    keys = key_path.split('.')
    value = data
    # Sequentially access each key in the path.
    for key in keys:
        if isinstance(value, list) and key.isdigit():
            try:
                # Handle list indexing if a key is a digit.
                value = value[int(key)]
            except IndexError:
                # Index is out of bounds.
                return None
        elif isinstance(value, dict):
            # Access the next level of the dictionary.
            value = value.get(key)
        else:
            # The path is invalid for the current data structure.
            return None

        if value is None:
            # A key was not found.
            return None

    return value

@overload
def get_elevation_from_api(
    latitude: float,
    longitude: float,
    api_name: str | None = None
) -> float:
    ...

@overload
def get_elevation_from_api(
    latitude: Sequence[float] | np.ndarray,
    longitude: Sequence[float] | np.ndarray,
    api_name: str | None = None
) -> np.ndarray:
    ...

def get_elevation_from_api(
    latitude: float | Sequence[float] | np.ndarray,
    longitude: float | Sequence[float] | np.ndarray,
    api_name: str | None = None
) -> float | np.ndarray:
    r"""Fetches elevation from sea level for geographic coords.

    This function queries an online service to retrieve elevation
    data (in meters) for one or more points on Earth's surface.
    It supports multiple APIs and includes a fallback mechanism.

    Parameters
    ----------
    latitude : float or array-like
        The latitude(s) of the point(s) in decimal degrees.
    longitude : float or array-like
        The longitude(s) of the point(s) in decimal degrees.
    api_name : str, optional
        The name of the elevation API to use (e.g., 'open_meteo',
        'open_topo_data'). Can also be a shorthand symbol like
        ',' or '|'. If ``None``, the default API is used.

    Returns
    -------
    float or numpy.ndarray
        - If inputs are scalars, returns a single float for the
          elevation in meters.
        - If inputs are arrays, returns a NumPy array of
          elevations.

    Raises
    ------
    GisError
        If the API request fails after all fallbacks or if the
        response is invalid.
    ImportError
        If the `requests` library is not installed.

    Notes
    -----
    This utility requires an active internet connection and the
    `requests` library. If a specified API fails, the function
    will automatically attempt to use the default API as a
    fallback. The underlying data quality depends on the chosen
    API provider [1]_.

    See Also
    --------
    get_elevation_from_utm : A similar function for UTM coords.
    pycsamt.gis.config.ElevationAPIConfig : API configurations.

    References
    ----------
    .. [1] Open-Meteo Elevation API Documentation.
           https://open-meteo.com/en/docs/elevation-api

    Examples
    --------
    >>> # Fetch elevation for a single point (Mt. Everest)
    >>> lat, lon = 27.9881, 86.9250
    >>> elev = get_elevation_from_api(lat, lon)
    >>> print(f"Elevation of Mount Everest: {elev:.2f} m")
    Elevation of Mount Everest: 8815.00 m

    >>> # Fetch for multiple points using an array
    >>> lats = [34.05, 40.71, 35.68]
    >>> lons = [-118.24, -74.00, 139.69]
    >>> elevs = get_elevation_from_api(lats, lons)
    >>> print(elevs)
    [  86. 3968. 3939.]
    """
    from ..utils._dependency import import_optional_dependency

    # Lazily import requests to avoid hard dependency.
    requests = import_optional_dependency(
        "requests",
        extra="'requests' is required to fetch elevation data."
    )
    # Resolve the API name from full name or symbol.
    resolved_api_name = _resolve_api_name(api_name)
    api_config = ElevationAPIConfig.get_api_config(
        resolved_api_name
    )
    api_url = api_config['url']
    params_format = api_config['params_format']
    response_key = api_config['response_key']

    is_scalar = np.isscalar(latitude)

    # Elevation APIs cap how many locations one request can carry
    # (Open-Meteo rejects overly long point lists with a 400). Batch
    # array input into chunks and concatenate, so a normal-sized
    # survey (dozens to a couple hundred stations) doesn't fail
    # outright just because it exceeds one request's limit.
    _MAX_BATCH = 90
    if not is_scalar and len(latitude) > _MAX_BATCH:
        lat_arr = np.asarray(latitude, dtype=float)
        lon_arr = np.asarray(longitude, dtype=float)
        chunks = [
            get_elevation_from_api(
                lat_arr[i:i + _MAX_BATCH],
                lon_arr[i:i + _MAX_BATCH],
                api_name,
            )
            for i in range(0, len(lat_arr), _MAX_BATCH)
        ]
        return np.concatenate(chunks)

    params = {}

    # Prepare parameters based on the API's required format.
    if params_format == 'comma_separated':
        params = {
            "latitude": (
                latitude if is_scalar
                else ",".join(map(str, latitude))
            ),
            "longitude": (
                longitude if is_scalar
                else ",".join(map(str, longitude))
            ),
        }
    elif params_format == 'pipe_separated':
        if is_scalar:
            locations = f"{latitude},{longitude}"
        else:
            locations = "|".join(
                [f"{lat},{lon}" for lat, lon in zip(
                    latitude, longitude)]
            )
        params = {"locations": locations}

    # For other formats like 'individual', handle them separately
    # if necessary, or add them to the ElevationAPIConfig.

    try:
        # Make the API call.
        response = requests.get(api_url, params=params)
        response.raise_for_status()  # Check for HTTP errors.
        data = response.json()

        # Extract elevation using the nested key path.
        elevations = _extract_value_from_nested_dict(
            data, response_key
        )

        if elevations is None:
            raise GisError(
                "Elevation data not found in API response."
            )

        if is_scalar:
            # Return a single float for scalar input.
            return float(
                elevations[0] if isinstance(elevations, list)
                else elevations
            )
        else:
            # Return a NumPy array for array input.
            return np.array(elevations, dtype=float)

    except requests.exceptions.RequestException as e:
        logger.error(
            f"API request to {resolved_api_name} failed: {e}"
        )
        # If the failed API was not the default, try fallback.
        if resolved_api_name != ElevationAPIConfig.DEFAULT_API:
            logger.info(
                "Trying fallback API: "
                f"{ElevationAPIConfig.DEFAULT_API}"
            )
            return get_elevation_from_api(
                latitude, longitude,
                ElevationAPIConfig.DEFAULT_API
            )

        # If fallback also fails or was the default, raise error.
        raise GisError(
            f"Failed to fetch elevation data from API: {e}"
        ) from e

@overload
def get_elevation_from_utm(
    easting: float,
    northing: float,
    zone: str,
    datum: str = "WGS84",
    **kws
) -> float:
    ...

@overload
def get_elevation_from_utm(
    easting: Sequence[float] | np.ndarray,
    northing: Sequence[float] | np.ndarray,
    zone: str,
    datum: str = "WGS84",
    **kws
) -> np.ndarray:
    ...

def get_elevation_from_utm(
    easting: float | Sequence[float] | np.ndarray,
    northing: float | Sequence[float] | np.ndarray,
    zone: str,
    datum: str = "WGS84",
    **kws
) -> float | np.ndarray:
    r"""Fetches elevation from sea level for UTM coordinates.

    This utility first converts UTM coordinates (easting,
    northing) to geographic latitude and longitude, then queries
    an online API to retrieve the elevation in meters.

    Parameters
    ----------
    easting : float or array-like
        The UTM easting coordinate(s) in meters.
    northing : float or array-like
        The UTM northing coordinate(s) in meters.
    zone : str
        The UTM zone designator for the coordinates, including
        the hemisphere letter (e.g., "11S", "32N").
    datum : str, default="WGS84"
        The geodetic datum of the UTM coordinates.

    Returns
    -------
    float or numpy.ndarray
        - If the input is a single point, returns a float
          representing the elevation in meters.
        - If the input is an array of points, returns a NumPy
          array of elevations.

    Raises
    ------
    GisError
        If the coordinate conversion or the subsequent API call
        fails.

    See Also
    --------
    get_elevation_from_api : The underlying function that
                             fetches elevation data.
    to_ll : The coordinate conversion function used internally.

    Examples
    --------
    >>> from pycsamt.gis.utils import get_elevation_from_utm
    >>> # Fetch elevation for a single UTM point in Los Angeles
    >>> easting, northing, zone = 377274.0, 3762150.0, "11S"
    >>> elevation = get_elevation_from_utm(
    ...     easting, northing, zone
    ... )
    >>> print(f"Elevation: {elevation:.2f} m")
    Elevation: 86.00 m
    """
    try:
        # Step 1: Convert UTM coordinates to Latitude/Longitude.
        # The `to_ll` function handles both scalar and array
        # inputs automatically.
        lat, lon = to_ll(
            easting=easting,
            northing=northing,
            zone=zone,
            datum=datum,
            as_frame=False
        )

        # Step 2: Fetch elevation using the geographic coords.
        # This function also handles scalar and array inputs.
        return get_elevation_from_api(lat, lon, **kws)

    except Exception as e:
        # Catch any error from conversion or the API call.
        logger.error(
            f"Failed to get elevation from UTM coordinates: {e}"
        )
        # Raise a specific error for better error handling.
        raise GisError(
            "Could not process UTM coordinates to get elevation"
        ) from e


def calculate_azimuth(
    easting: Sequence[float] | NDArray,
    northing: Sequence[float] | NDArray
) -> NDArray:
    r"""Calculates the forward azimuth between consecutive points.

    This function takes a sequence of UTM coordinates and computes
    the azimuth (compass direction) for the line segment starting
    at each point and ending at the next.

    Parameters
    ----------
    easting : array-like
        A sequence (list, NumPy array, etc.) of the easting
        coordinates in meters.
    northing : array-like
        A sequence of the northing coordinates in meters. Must be
        the same length as ``easting``.

    Returns
    -------
    numpy.ndarray
        An array of azimuths in degrees (0-360), measured
        clockwise from North. The last value is always ``NaN`` as
        there is no subsequent point to define a direction.

    Notes
    -----
    The azimuth at index ``i`` represents the direction of the
    vector from point ``i`` to point ``i+1``. The calculation is
    performed on the horizontal (2D) projection of the points;
    elevation is not used.

    Examples
    --------
    >>> import numpy as np
    >>> # A simple square path moving East, North, West
    >>> easting = [0, 100, 100, 0]
    >>> northing = [0, 0, 100, 100]
    >>> calculate_azimuth(easting, northing)
    array([ 90.,   0., 270.,  nan])
    """
    # Ensure inputs are NumPy arrays for vectorization
    easting = np.asarray(easting, dtype=float)
    northing = np.asarray(northing, dtype=float)

    # Azimuth requires at least two points to define a direction
    if easting.size < 2:
        return np.full_like(easting, np.nan)

    # Calculate the difference in coordinates between points
    delta_easting = np.diff(easting)
    delta_northing = np.diff(northing)

    # Calculate the angle in radians using arctan2 for quadrant
    # correctness. The arguments are (dy, dx) in traditional
    # math, so we use (delta_easting, delta_northing).
    angle_rad = np.arctan2(delta_easting, delta_northing)

    # Convert radians to degrees
    angle_deg = np.rad2deg(angle_rad)

    # Normalize the angle to a 0-360 degree azimuth
    azimuths = (angle_deg + 360) % 360

    # The result should have the same shape as the input. The
    # azimuth is assigned to the start point of each segment.
    # The last point has no forward segment, so it is NaN.
    result = np.full_like(easting, np.nan)
    result[:-1] = azimuths

    return result
