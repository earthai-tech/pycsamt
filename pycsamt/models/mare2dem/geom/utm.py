# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""UTM ↔ Longitude/Latitude coordinate conversion.

Port of:
  * ``LonLat2UTM.m``
  * ``UTM2LonLat.m``

When ``pyproj`` is installed the conversion delegates to it for full
accuracy and datum support.  Otherwise, the pure-Python WGS-84 Fortran-
origin implementation from the MATLAB code is used.

Reference
---------
Snyder, J. P. (1983). Map Projections Used by the US Geological Survey.
USGS Bulletin 1532, 2nd Edition, USGPO, Washington D.C.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "lonlat_to_utm",
    "utm_to_lonlat",
    "ELLIPSOIDS",
]

# ---------------------------------------------------------------------------
# Ellipsoid table  (semi-major axis in metres, inverse flattening)
# ---------------------------------------------------------------------------

ELLIPSOIDS: dict[str, tuple[float, float]] = {
    "wgs84":  (6_378_137.0,   1 / 298.257),
    "grs80":  (6_378_137.0,   1 / 298.257),
    "intl24": (6_378_388.0,   1 / 297.000),
    "ed50":   (6_378_388.0,   1 / 297.000),
    "sphere": (6_370_997.0,   0.0),
    "normal": (1.0,           0.0),
    "grs67":  (6_378_160.0,   1 / 247.247),
    "wgs72":  (6_378_135.0,   1 / 298.260),
    "wgs66":  (6_378_145.0,   1 / 298.250),
    "wgs60":  (6_378_165.0,   1 / 298.300),
    "clrk66": (6_378_206.4,   1 / 294.980),
    "clrk80": (6_378_249.1,   1 / 293.466),
    "intl67": (6_378_157.5,   1 / 298.250),
}

_K0 = 0.9996
_FALSE_EAST = 500_000.0
_FALSE_NORTH_SOUTH = 10_000_000.0
_DEG2RAD = math.pi / 180.0


# ---------------------------------------------------------------------------
# Pure-Python WGS-84 implementation (no external dependency)
# ---------------------------------------------------------------------------

def _ll_to_utm_pure(
    lat: np.ndarray,
    lon: np.ndarray,
    zone: int,
    south_hemi: bool,
    ellipsoid: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert Lat/Lon arrays to UTM E, N using Snyder's formulas."""
    a, f = ELLIPSOIDS[ellipsoid.lower()]
    e2 = 2 * f - f * f
    e4 = e2 * e2
    e6 = e4 * e2
    ep2 = e2 / (1 - e2)

    lat_r = lat * _DEG2RAD
    lon_r = lon * _DEG2RAD
    lambda0 = ((-180 + zone * 6) - 3) * _DEG2RAD

    sinL = np.sin(lat_r)
    tanL = np.tan(lat_r)
    cosL = np.cos(lat_r)
    T = tanL ** 2
    C = ep2 * cosL ** 2
    A = (lon_r - lambda0) * cosL
    A2 = A ** 2
    A4 = A2 ** 2
    S = sinL ** 2
    N = a / np.sqrt(1 - e2 * S)

    M0 = 1 - e2 * 0.25 - e4 * 0.046875 - e6 * 0.01953125
    M1 = e2 * 0.375 + e4 * 0.09375 + e6 * 0.043945313
    M2 = e4 * 0.05859375 + e6 * 0.043945313
    M3 = e6 * 0.011393229
    M = a * (
        M0 * lat_r
        - M1 * np.sin(2 * lat_r)
        + M2 * np.sin(4 * lat_r)
        - M3 * np.sin(6 * lat_r)
    )

    X0 = A4 * A / 120
    X1 = 5 - 18 * T + T ** 2 + 72 * C - 58 * ep2
    X2 = A2 * A / 6
    X3 = 1 - T + C
    easting = N * (A + X3 * X2 + X1 * X0)

    Y0 = 61 - 58 * T + T ** 2 + 600 * C - 330 * ep2
    Y1 = 5 - T + 9 * C + 4 * C ** 2
    northing = M + N * tanL * (A2 / 2 + Y1 * A4 / 24 + Y0 * A4 * A2 / 720)

    easting = easting * _K0 + _FALSE_EAST
    northing = northing * _K0
    if south_hemi:
        northing = northing + _FALSE_NORTH_SOUTH

    return easting, northing


def _utm_to_ll_pure(
    easting: np.ndarray,
    northing: np.ndarray,
    zone: int,
    south_hemi: bool,
    ellipsoid: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert UTM E, N arrays to Lat/Lon using Snyder's formulas."""
    a, f = ELLIPSOIDS[ellipsoid.lower()]
    e2 = 2 * f - f * f
    e4 = e2 * e2
    e6 = e4 * e2
    ep2 = e2 / (1 - e2)

    x = (easting - _FALSE_EAST) / _K0
    y = northing.copy()
    if south_hemi:
        y = y - _FALSE_NORTH_SOUTH
    y = y / _K0

    lambda0 = ((-180 + zone * 6) - 3) * _DEG2RAD

    M0 = 1 - e2 * 0.25 - e4 * 0.046875 - e6 * 0.01953125
    e2 * 0.375 + e4 * 0.09375 + e6 * 0.043945313
    e4 * 0.05859375 + e6 * 0.043945313
    e6 * 0.011393229

    mu = y / (a * M0)
    e1 = (1 - np.sqrt(1 - e2)) / (1 + np.sqrt(1 - e2))
    phi1 = (
        mu
        + (1.5 * e1 - 27 / 32 * e1 ** 3) * np.sin(2 * mu)
        + (21 / 16 * e1 ** 2 - 55 / 32 * e1 ** 4) * np.sin(4 * mu)
        + 151 / 96 * e1 ** 3 * np.sin(6 * mu)
        + 1097 / 512 * e1 ** 4 * np.sin(8 * mu)
    )

    sin_phi1 = np.sin(phi1)
    cos_phi1 = np.cos(phi1)
    tan_phi1 = np.tan(phi1)

    N1 = a / np.sqrt(1 - e2 * sin_phi1 ** 2)
    T1 = tan_phi1 ** 2
    C1 = ep2 * cos_phi1 ** 2
    R1 = a * (1 - e2) / (1 - e2 * sin_phi1 ** 2) ** 1.5
    D = x / (N1 * _K0)
    D2 = D ** 2
    D4 = D2 ** 2

    lat = phi1 - (N1 * tan_phi1 / R1) * (
        D2 / 2
        - (5 + 3 * T1 + 10 * C1 - 4 * C1 ** 2 - 9 * ep2) * D4 / 24
        + (61 + 90 * T1 + 298 * C1 + 45 * T1 ** 2 - 252 * ep2 - 3 * C1 ** 2)
        * D4 * D2 / 720
    )

    lon = lambda0 + (
        D
        - (1 + 2 * T1 + C1) * D2 * D / 6
        + (5 - 2 * C1 + 28 * T1 - 3 * C1 ** 2 + 8 * ep2 + 24 * T1 ** 2)
        * D4 * D / 120
    ) / cos_phi1

    return lat / _DEG2RAD, lon / _DEG2RAD


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def lonlat_to_utm(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    *,
    zone: int | None = None,
    south_hemi: bool | None = None,
    ellipsoid: str = "wgs84",
) -> tuple[np.ndarray, np.ndarray, int, bool]:
    """Convert geographic longitude / latitude to UTM.

    Port of ``LonLat2UTM.m``.  When ``pyproj`` is installed the
    conversion is delegated to it; otherwise the pure-Python Snyder
    formulas are used.

    Parameters
    ----------
    lon : float or array-like
        Longitude in decimal degrees (−180 to +180).
    lat : float or array-like
        Latitude in decimal degrees.
    zone : int, optional
        Force a specific UTM zone. Auto-computed from median longitude
        when omitted.
    south_hemi : bool, optional
        Force southern-hemisphere false northing. Auto-detected from
        median latitude when omitted.
    ellipsoid : str, default "wgs84"
        Ellipsoid name (see :data:`ELLIPSOIDS`).

    Returns
    -------
    easting : numpy.ndarray
        UTM easting in metres.
    northing : numpy.ndarray
        UTM northing in metres.
    zone : int
        UTM zone used.
    south_hemi : bool
        Hemisphere flag used.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.geom.utm import lonlat_to_utm
    >>> e, n, zone, sh = lonlat_to_utm(-70.0, 42.0)
    >>> zone
    19
    """
    lon_arr = np.atleast_1d(np.asarray(lon, dtype=float))
    lat_arr = np.atleast_1d(np.asarray(lat, dtype=float))

    # normalise lon to [-180, 180]
    lon_arr = lon_arr % 360
    lon_arr[lon_arr > 180] -= 360

    if zone is None:
        zone = int(1 + (np.median(lon_arr) + 180) / 6)
    if south_hemi is None:
        south_hemi = bool(np.median(lat_arr) < 0)

    try:
        from pyproj import CRS, Transformer
        crs_ll = CRS.from_epsg(4326)
        crs_utm = CRS.from_dict({
            "proj": "utm", "zone": zone, "south": south_hemi,
            "datum": "WGS84", "units": "m",
        })
        tf = Transformer.from_crs(crs_ll, crs_utm, always_xy=True)
        east, north = tf.transform(lon_arr, lat_arr)
        return np.asarray(east), np.asarray(north), zone, south_hemi
    except ImportError:
        pass

    east, north = _ll_to_utm_pure(lat_arr, lon_arr, zone, south_hemi, ellipsoid)
    return east, north, zone, south_hemi


def utm_to_lonlat(
    easting: float | np.ndarray,
    northing: float | np.ndarray,
    zone: int,
    south_hemi: bool | str = False,
    *,
    ellipsoid: str = "wgs84",
) -> tuple[np.ndarray, np.ndarray]:
    """Convert UTM coordinates to geographic longitude / latitude.

    Port of ``UTM2LonLat.m``.

    Parameters
    ----------
    easting : float or array-like
        UTM easting in metres.
    northing : float or array-like
        UTM northing in metres.
    zone : int
        UTM zone number.
    south_hemi : bool or str, default False
        ``True`` / ``'S'`` for southern hemisphere.
    ellipsoid : str, default "wgs84"
        Ellipsoid name (see :data:`ELLIPSOIDS`).

    Returns
    -------
    lon : numpy.ndarray
        Longitude in decimal degrees.
    lat : numpy.ndarray
        Latitude in decimal degrees.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.geom.utm import utm_to_lonlat
    >>> lon, lat = utm_to_lonlat(330000.0, 4650000.0, 19, False)
    """
    if isinstance(south_hemi, str):
        south_hemi = south_hemi.upper().startswith("S")
    south_hemi = bool(south_hemi)

    e_arr = np.atleast_1d(np.asarray(easting, dtype=float))
    n_arr = np.atleast_1d(np.asarray(northing, dtype=float))

    try:
        from pyproj import CRS, Transformer
        crs_ll = CRS.from_epsg(4326)
        crs_utm = CRS.from_dict({
            "proj": "utm", "zone": zone, "south": south_hemi,
            "datum": "WGS84", "units": "m",
        })
        tf = Transformer.from_crs(crs_utm, crs_ll, always_xy=True)
        lon, lat = tf.transform(e_arr, n_arr)
        return np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    except ImportError:
        pass

    lat, lon = _utm_to_ll_pure(e_arr, n_arr, zone, south_hemi, ellipsoid)
    return lon, lat
