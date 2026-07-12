# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.gis.utils` coordinate helpers.

Covers DMS <-> decimal conversion, coordinate validation, UTM zone
resolution, and the pure USGS Bulletin 1532 UTM <-> lat/lon
projections (no pyproj/GDAL required).
"""

from __future__ import annotations

import pytest

from pycsamt.gis.utils import (
    GisError,
    assert_elevation_value,
    assert_lat_value,
    assert_lon_value,
    convert_position_float2str,
    convert_position_str2float,
    decimal_to_dms,
    dms_to_decimal,
    get_utm_zone,
    ll_to_utm,
    utm_to_ll,
    utm_zone_to_epsg,
)

WGS84 = 23  # ELLIPSOIDS id used across the code base


# ----------------------- DMS <-> decimal conversion -------------------

@pytest.mark.parametrize(
    "text, expected",
    [
        ("34:03:00", 34.05),
        ("-118:20:24", -118.34),
        ("-118.34", -118.34),
        ("26:00:00N", 26.0),
        ("010:00:00 W", -10.0),
        ("10:60:00", 11.0),        # minute rollover
        ("-118:60:00", -119.0),    # rollover applies to magnitude
        ("10:00:90", 10.025),      # second rollover
    ],
)
def test_convert_position_str2float(text, expected):
    assert convert_position_str2float(text) == pytest.approx(expected)


def test_convert_position_str2float_none_and_invalid():
    assert convert_position_str2float(None) is None
    assert convert_position_str2float("None") is None
    with pytest.raises(ValueError):
        convert_position_str2float("12:34")  # not DD:MM:SS


def test_dms_to_decimal_wrapper():
    assert dms_to_decimal("34:03:00") == pytest.approx(34.05)
    assert dms_to_decimal("-118.34") == pytest.approx(-118.34)


def test_float_to_dms_roundtrip():
    text = convert_position_float2str(-118.34563)
    assert text.startswith("-118:20:")
    back = convert_position_str2float(text)
    assert back == pytest.approx(-118.34563, abs=1e-6)
    assert decimal_to_dms(34.05) == convert_position_float2str(34.05)
    with pytest.raises(AssertionError):
        convert_position_float2str("not-a-float")


# --------------------------- validators --------------------------------

def test_assert_lat_value_paths():
    assert assert_lat_value("34:03:00") == pytest.approx(34.05)
    assert assert_lat_value(45.0) == 45.0
    assert assert_lat_value(None) is None
    with pytest.raises(ValueError):
        assert_lat_value(95.0)


def test_assert_lon_value_paths():
    assert assert_lon_value("-118:20:24") == pytest.approx(-118.34)
    assert assert_lon_value(-118.34) == -118.34
    assert assert_lon_value(None) is None
    with pytest.raises(ValueError):
        assert_lon_value(190.0)


def test_assert_elevation_value_fallback():
    assert assert_elevation_value("12.5") == 12.5
    assert assert_elevation_value("oops") == 0.0


# ------------------------------ UTM zones ------------------------------

def test_get_utm_zone_los_angeles_and_sydney():
    zn, north, zstr = get_utm_zone(34.05, -118.34)
    assert (zn, north, zstr) == (11, True, "11S")
    zn, north, zstr = get_utm_zone(-33.86, 151.21)
    assert (zn, north) == (56, False)
    assert zstr.startswith("56")


def test_utm_zone_to_epsg_wgs84_codes():
    north = utm_zone_to_epsg(11, True)
    south = utm_zone_to_epsg(56, False)
    if north is not None:
        assert north == 32611
    if south is not None:
        assert south == 32756


# ----------------------- UTM projection roundtrip ----------------------

def test_ll_to_utm_reference_point():
    zone, easting, northing = ll_to_utm(WGS84, 34.05, -118.34)
    assert zone == "11S"
    # reference values from standard UTM converters
    assert easting == pytest.approx(376322.2, abs=1.0)
    assert northing == pytest.approx(3768509.8, abs=1.0)


def test_utm_ll_roundtrip_north_and_south():
    for lat, lon in [
        (34.05, -118.34),      # northern hemisphere
        (-33.86, 151.21),      # southern hemisphere
        (0.5, 6.6),            # near equator
        (63.0, 10.4),          # high latitude
    ]:
        zone, easting, northing = ll_to_utm(WGS84, lat, lon)
        lat2, lon2 = utm_to_ll(WGS84, northing, easting, zone)
        assert lat2 == pytest.approx(lat, abs=1e-5)
        assert lon2 == pytest.approx(lon, abs=1e-5)


def test_ll_to_utm_norway_exception():
    zone, _, _ = ll_to_utm(WGS84, 60.0, 5.0)
    assert zone.startswith("32")


def test_ll_to_utm_unknown_ellipsoid_raises():
    with pytest.raises(GisError):
        ll_to_utm(999, 10.0, 10.0)
    with pytest.raises(GisError):
        utm_to_ll(999, 1000.0, 500000.0, "31N")
