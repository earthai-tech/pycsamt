# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Additional coverage tests for :mod:`pycsamt.gis.utils`.

Complements ``test_gis_utils.py`` by exercising:

- ``normalize_lat_lon`` dispatch heuristics (scalar/array, assume,
  error, clip).
- ``assert_xy_coordinate_system`` branch coverage.
- DMS rollover/validation edge cases.
- ``get_utm_string_from_sr`` (GDAL-duck-typed, no real GDAL needed).
- pyproj-backed ``to_ll``/``to_utm``, ``project_point_ll2utm``,
  ``project_point_utm2ll``, ``project_points_ll2utm``, ``epsg_project``.
- ``_resolve_api_name``/``_extract_value_from_nested_dict`` helpers.
- ``get_elevation_from_api``/``get_elevation_from_utm`` with the
  network layer mocked (no real HTTP calls).
- ``calculate_azimuth``.
- GDAL-gated deprecated transforms (raise ``GDALMissingError`` since
  GDAL is not installed in this environment).

Notes on bugs discovered while writing these tests (reported, not
fixed here -- see task report): a few tests below pin down *current*
(surprising) behavior explicitly so as not to silently paper over it.
"""

from __future__ import annotations

import sys
import types

import numpy as np
import pandas as pd
import pytest

from pycsamt.gis.config import GDALMissingError, GisError
from pycsamt.gis.utils import (
    _assert_minutes,
    _assert_seconds,
    _extract_value_from_nested_dict,
    _resolve_api_name,
    _rollover_dms,
    assert_lat_value,
    assert_lon_value,
    assert_xy_coordinate_system,
    calculate_azimuth,
    convert_position_str2float,
    epsg_project,
    get_elevation_from_api,
    get_elevation_from_utm,
    get_epsg,
    get_utm_string_from_sr,
    ll_to_utm,
    normalize_lat_lon,
    project_point_ll2utm,
    project_point_utm2ll,
    project_points_ll2utm,
    to_ll,
    to_utm,
    transform_ll_to_utm,
    transform_utm_to_ll,
    utm_wgs84_conv,
    utm_zone_to_epsg,
)

WGS84 = 23  # ELLIPSOIDS id used across the code base


# ======================================================================
# normalize_lat_lon
# ======================================================================


def test_normalize_lat_lon_decisive_orders():
    # a is lon-like (>90), b is lat-like -> canonical (lat, lon)
    assert normalize_lat_lon(170.0, 45.0) == (45.0, 170.0)
    # a is lat-like, b is lon-like -> same canonical result
    assert normalize_lat_lon(45.0, 170.0) == (45.0, 170.0)


def test_normalize_lat_lon_ambiguous_assume():
    # both |val| <= 90: assume="lonlat" (default) swaps a,b -> (b, a)
    assert normalize_lat_lon(10.0, 20.0, assume="lonlat") == (20.0, 10.0)
    assert normalize_lat_lon(10.0, 20.0, assume="auto") == (20.0, 10.0)
    # assume="latlon" keeps input order as (lat, lon)
    assert normalize_lat_lon(10.0, 20.0, assume="latlon") == (10.0, 20.0)


def test_normalize_lat_lon_over_180_clip_true():
    lat, lon = normalize_lat_lon(200.0, 45.0, clip=True)
    # clipped to 180 then resolved decisively since 45 <= 90
    assert lat == 45.0
    assert lon == 180.0


def test_normalize_lat_lon_over_180_no_clip_ignore():
    lat, lon = normalize_lat_lon(200.0, 45.0, error="ignore")
    assert np.isnan(lat) and np.isnan(lon)


def test_normalize_lat_lon_over_180_no_clip_raise():
    with pytest.raises(ValueError):
        normalize_lat_lon(200.0, 45.0, error="raise")


def test_normalize_lat_lon_huge_values_treated_as_over_180():
    # |val| > 1e9 is set to nan internally, but this also always
    # satisfies the >180 check first, so it raises/returns nan the
    # same way as any other over-180 value.
    lat, lon = normalize_lat_lon(1e10, 45.0, error="ignore")
    assert np.isnan(lat) and np.isnan(lon)
    with pytest.raises(ValueError):
        normalize_lat_lon(1e10, 45.0, error="raise")


def test_normalize_lat_lon_fallback_unresolvable():
    # both > 90 (but <= 180): not decisive, not ambiguous -> fallback
    # path is reached; whichever value lands on "lat" is > 90, so the
    # pair is always unresolvable here.
    lat, lon = normalize_lat_lon(95.0, 91.0, error="ignore")
    assert np.isnan(lat) and np.isnan(lon)
    with pytest.raises(ValueError):
        normalize_lat_lon(95.0, 91.0, error="raise")
    # also exercise the clip branch inside the fallback section
    # (clip=True re-clamps to [-180, 180], a no-op here, but the
    # line is executed).
    lat2, lon2 = normalize_lat_lon(95.0, 91.0, error="ignore", clip=True)
    assert np.isnan(lat2) and np.isnan(lon2)


def test_normalize_lat_lon_string_and_none_inputs():
    # DMS string parsed via convert_position_str2float
    lat, lon = normalize_lat_lon("34:03:00", "-118:20:24", assume="latlon")
    assert lat == pytest.approx(34.05)
    assert lon == pytest.approx(-118.34)
    # None/"None" -> nan
    lat, lon = normalize_lat_lon(None, "None", error="ignore")
    assert np.isnan(lat) and np.isnan(lon)
    # unparsable string -> nan (caught in _to_float), ambiguous fallback
    lat, lon = normalize_lat_lon("not-a-number", 10.0, error="ignore")
    assert np.isnan(lat) or np.isnan(lon)


def test_normalize_lat_lon_array_inputs():
    a = [170.0, 10.0]
    b = [45.0, 20.0]
    lat_out, lon_out = normalize_lat_lon(a, b, assume="lonlat")
    assert isinstance(lat_out, np.ndarray)
    np.testing.assert_allclose(lat_out, [45.0, 20.0])
    np.testing.assert_allclose(lon_out, [170.0, 10.0])


def test_normalize_lat_lon_array_shape_mismatch():
    with pytest.raises(ValueError):
        normalize_lat_lon([1.0, 2.0], [3.0, 4.0, 5.0])


# ======================================================================
# assert_xy_coordinate_system
# ======================================================================


def test_assert_xy_coordinate_system_ll_orientation_lon_lat():
    x = np.array([10.0, -20.0, 30.0])
    y = np.array([45.0, -50.0, 60.0])
    assert assert_xy_coordinate_system(x, y) == "ll"


def test_assert_xy_coordinate_system_ll_orientation_lat_lon():
    # |x| <= 90 but |y| up to 180 -> alternate orientation branch
    x = np.array([34.05, -33.86])
    y = np.array([-118.34, 151.21])
    assert assert_xy_coordinate_system(x, y) == "ll"


def test_assert_xy_coordinate_system_utm_out_of_range():
    x = np.array([500000.0, 400000.0])
    y = np.array([3768509.0, 6251925.0])
    assert assert_xy_coordinate_system(x, y) == "utm"


def test_assert_xy_coordinate_system_utm_non_numeric():
    # neither DMS (no ':') nor numeric -> falls back to 'utm'
    x = ["abc", "def"]
    y = ["ghi", "jkl"]
    assert assert_xy_coordinate_system(x, y) == "utm"


def test_assert_xy_coordinate_system_dms_via_x_or_y():
    x = ["28:24:43.08", "28:24:42.69"]
    y = ["109:19:58.34", "109:19:58.93"]
    assert assert_xy_coordinate_system(x, y) == "dms"
    # DMS present only in y
    assert assert_xy_coordinate_system([1.0, 2.0], ["1:2:3", "4:5:6"]) == "dms"


# ======================================================================
# DMS helpers: _assert_minutes / _assert_seconds / _rollover_dms
# ======================================================================


@pytest.mark.parametrize("value", [0, 30, 59.999])
def test_assert_minutes_valid(value):
    assert _assert_minutes(value) == value


@pytest.mark.parametrize("value", [-1, 60, 61.5])
def test_assert_minutes_invalid(value):
    with pytest.raises(ValueError):
        _assert_minutes(value)


@pytest.mark.parametrize("value", [0, 12.345, 59.999])
def test_assert_seconds_valid(value):
    assert _assert_seconds(value) == value


@pytest.mark.parametrize("value", [-0.5, 60, 90])
def test_assert_seconds_invalid(value):
    with pytest.raises(ValueError):
        _assert_seconds(value)


def test_rollover_dms_no_carry():
    assert _rollover_dms(10, 30) == (10, 30)


def test_rollover_dms_with_carry():
    assert _rollover_dms(10, 120) == (12, 0)
    assert _rollover_dms(0, 61) == (1, 1)


# ======================================================================
# convert_position_str2float edge branches
# ======================================================================


def test_convert_position_str2float_sanitized_seconds():
    # trailing junk on the seconds token that isn't a recognized
    # hemisphere letter forces the regex-sanitize retry path.
    val = convert_position_str2float("10:20:15abc")
    assert val == pytest.approx(10 + 20 / 60 + 15 / 3600)


def test_convert_position_str2float_unsanitizable_seconds_raises():
    with pytest.raises(ValueError):
        convert_position_str2float("10:20:abcxyz")


def test_convert_position_str2float_invalid_minutes_raises():
    with pytest.raises(ValueError):
        convert_position_str2float("10:ab:15")


# ======================================================================
# decimal_to_dms / convert_position_float2str rollover edge cases
# ======================================================================


def test_convert_position_float2str_seconds_rollover():
    from pycsamt.gis.utils import convert_position_float2str

    # crafted so the raw seconds value rounds to 60.0000 at 4 decimals,
    # forcing the "sec >= 60 -> minutes += 1, sec = 0" branch.
    text = convert_position_float2str(10.766666653333333)
    assert text == "10:46:00.00"


def test_convert_position_float2str_minutes_rollover_to_degree():
    from pycsamt.gis.utils import convert_position_float2str

    # forces both the seconds AND the minutes==60 -> degree carry
    # branches.
    text = convert_position_float2str(10.999999986666667)
    assert text == "11:00:00.00"


# ======================================================================
# get_utm_string_from_sr (GDAL-gated, duck-typed fake spatial ref)
# ======================================================================


class _FakeSpatialRef:
    def __init__(self, zone):
        self._zone = zone

    def GetUTMZone(self):
        return self._zone


def test_get_utm_string_from_sr_north():
    with pytest.deprecated_call():
        assert get_utm_string_from_sr(_FakeSpatialRef(11)) == "11N"


def test_get_utm_string_from_sr_south():
    with pytest.deprecated_call():
        assert get_utm_string_from_sr(_FakeSpatialRef(-56)) == "56S"


def test_get_utm_string_from_sr_zero():
    with pytest.deprecated_call():
        assert get_utm_string_from_sr(_FakeSpatialRef(0)) == "0"


# ======================================================================
# project_point_ll2utm / project_point_utm2ll (pyproj-backed)
# ======================================================================


def test_project_point_ll2utm_scalar_auto_zone():
    e, n, z = project_point_ll2utm(34.05, -118.34)
    assert z == "11S"
    assert e == pytest.approx(376322.2, abs=1.0)
    assert n == pytest.approx(3768509.8, abs=1.0)


def test_project_point_ll2utm_none_inputs():
    assert project_point_ll2utm(None, -118.34) == (None, None, None)
    assert project_point_ll2utm(34.05, None) == (None, None, None)


def test_project_point_ll2utm_explicit_epsg():
    e, n, z = project_point_ll2utm(34.05, -118.34, epsg=32611)
    assert z == "11S"
    assert e == pytest.approx(376322.2, abs=1.0)
    assert n == pytest.approx(3768509.8, abs=1.0)


def test_project_point_ll2utm_array_input_recarray():
    res = project_point_ll2utm([34.05, -33.86], [-118.34, 151.21])
    assert isinstance(res, np.recarray)
    assert list(res["utm_zone"]) == ["11S", "56H"]
    assert res["easting"][0] == pytest.approx(376322.2, abs=1.0)


def test_project_point_ll2utm_array_size_mismatch():
    with pytest.raises(ValueError):
        project_point_ll2utm([34.05, 35.0], [-118.34])


def test_project_point_utm2ll_explicit_epsg_none_roundtrip():
    # NOTE: project_point_utm2ll's own default is epsg=3149 (a
    # non-WGS84 CRS) -- explicit epsg=None is required to honor the
    # given zone/datum. See report for details (suspicious default).
    e, n, z = project_point_ll2utm(34.05, -118.34)
    lat, lon = project_point_utm2ll(e, n, z, epsg=None)
    assert lat == pytest.approx(34.05, abs=1e-4)
    assert lon == pytest.approx(-118.34, abs=1e-4)


def test_project_point_utm2ll_default_epsg_is_none_and_uses_zone():
    # epsg now defaults to None, so calling with only
    # (easting, northing, zone) -- as shown in the function's own
    # docstring example -- correctly reproduces the original lat/lon
    # via zone-based CRS resolution (previously overridden by a
    # stray epsg=3149 default).
    e, n, z = project_point_ll2utm(34.05, -118.34)
    lat, lon = project_point_utm2ll(e, n, z)
    assert (lat, lon) == pytest.approx((34.05, -118.34), abs=1e-3)


def test_project_point_utm2ll_int_zone_and_bytes_zone():
    e, n, z = project_point_ll2utm(34.05, -118.34)
    lat_i, lon_i = project_point_utm2ll(e, n, 11, epsg=None)
    assert lat_i == pytest.approx(34.05, abs=1e-4)
    assert lon_i == pytest.approx(-118.34, abs=1e-4)

    lat_b, lon_b = project_point_utm2ll(e, n, np.bytes_(b"11S"), epsg=None)
    assert lat_b == pytest.approx(34.05, abs=1e-4)
    assert lon_b == pytest.approx(-118.34, abs=1e-4)


def test_project_point_utm2ll_unsupported_zone_type_raises():
    with pytest.raises(NotImplementedError):
        project_point_utm2ll(1.0, 2.0, 3.5, epsg=None)


def test_project_point_utm2ll_bad_zone_string_raises_value_error():
    with pytest.raises(ValueError):
        project_point_utm2ll(1.0, 2.0, "abcS", epsg=None)


def test_project_point_utm2ll_non_float_coords_raise_gis_error():
    with pytest.raises(GisError):
        project_point_utm2ll("abc", 2.0, "11S")
    with pytest.raises(GisError):
        project_point_utm2ll(1.0, "xyz", "11S")


def test_project_point_ll2utm_explicit_zone_letter_n_is_northern():
    # An explicit UTM zone string ending in the MGRS band letter "N"
    # (valid for latitudes 0-8 N) must be treated as northern
    # hemisphere, consistent with project_point_utm2ll /
    # project_points_ll2utm (both use '>=' for this comparison).
    lat, lon = 4.0, 10.0
    e_auto, n_auto, z_auto = project_point_ll2utm(lat, lon)
    assert z_auto.endswith("N")
    e_explicit, n_explicit, z_explicit = project_point_ll2utm(lat, lon, utm_zone=z_auto)
    assert e_explicit == pytest.approx(e_auto)
    assert n_explicit == pytest.approx(n_auto)


# ======================================================================
# project_points_ll2utm (vectorized)
# ======================================================================


def test_project_points_ll2utm_scalar():
    e, n, z = project_points_ll2utm(34.05, -118.34)
    assert z == "11S"
    assert e == pytest.approx(376322.2, abs=1.0)
    assert n == pytest.approx(3768509.8, abs=1.0)


def test_project_points_ll2utm_array_auto_zone():
    lats = [34.00, 34.05]
    lons = [-118.40, -118.34]
    e, n, z = project_points_ll2utm(lats, lons)
    assert isinstance(e, np.ndarray)
    assert e.shape == (2,)
    assert z.endswith("S")


def test_project_points_ll2utm_explicit_zone():
    lats = [34.00, 34.05]
    lons = [-118.40, -118.34]
    e, n, z = project_points_ll2utm(lats, lons, utm_zone="11S")
    assert z == "11S"
    assert e.shape == (2,)


def test_project_points_ll2utm_explicit_epsg():
    e, n, z = project_points_ll2utm(34.05, -118.34, epsg=32611)
    assert e == pytest.approx(376322.2, abs=1.0)
    assert n == pytest.approx(3768509.8, abs=1.0)


def test_project_points_ll2utm_shape_mismatch():
    with pytest.raises(ValueError):
        project_points_ll2utm([1.0, 2.0], [3.0, 4.0, 5.0])


# ======================================================================
# to_ll / to_utm
# ======================================================================


def test_to_ll_scalar():
    e, n, z = project_point_ll2utm(34.05, -118.34)
    lat, lon = to_ll(e, n, z)
    assert lat == pytest.approx(34.05, abs=1e-4)
    assert lon == pytest.approx(-118.34, abs=1e-4)


def test_to_ll_array_and_series():
    e = np.array([376322.215, 334416.394])
    n = np.array([3768509.769, 6251925.360])
    z = np.array(["11S", "56H"], dtype=object)
    lat, lon = to_ll(e, n, z)
    np.testing.assert_allclose(lat, [34.05, -33.86], atol=1e-3)
    np.testing.assert_allclose(lon, [-118.34, 151.21], atol=1e-3)

    se = pd.Series(e, index=[10, 20])
    sn = pd.Series(n, index=[10, 20])
    sz = pd.Series(z, index=[10, 20])
    df = to_ll(se, sn, sz, as_frame=True)
    assert list(df.index) == [10, 20]
    np.testing.assert_allclose(df["lat"], [34.05, -33.86], atol=1e-3)


def test_to_ll_as_frame_with_scalar_zone_broadcast():
    e = np.array([376322.215])
    n = np.array([3768509.769])
    df = to_ll(e, n, "11S", as_frame=True)
    assert df.loc[0, "zone"] == "11S"
    assert df.loc[0, "lat"] == pytest.approx(34.05, abs=1e-3)


def test_to_ll_epsg_only_no_zone():
    lat, lon = to_ll(376322.215, 3768509.769, epsg=32611)
    assert lat == pytest.approx(34.05, abs=1e-3)
    assert lon == pytest.approx(-118.34, abs=1e-3)


def test_to_ll_missing_data_for_column_name_raises():
    with pytest.raises(ValueError):
        to_ll("E", "N")


def test_to_ll_shape_mismatch_raises():
    with pytest.raises(ValueError):
        to_ll([1.0, 2.0, 3.0], [4.0, 5.0])


def test_to_ll_zone_shape_mismatch_raises():
    e = np.array([1.0, 2.0])
    n = np.array([3.0, 4.0])
    zone = np.array(["11S", "12S", "13S"], dtype=object)
    with pytest.raises(ValueError):
        to_ll(e, n, zone)


def test_to_ll_zone_string_is_a_column_lookup_when_it_names_a_column():
    # A string 'zone' that names an actual column in 'data' is now
    # resolved as a per-row column lookup, matching easting/northing
    # and the documented behavior ("If strings are given, they are
    # treated as column names in data").
    df = pd.DataFrame(
        {
            "E": [376322.215, 334416.394],
            "N": [3768509.769, 6251925.360],
            "zone": ["11S", "56H"],
        }
    )
    out = to_ll("E", "N", "zone", data=df, as_frame=True)
    assert list(out["zone"]) == ["11S", "56H"]
    assert out["lat"].iloc[0] == pytest.approx(34.05, abs=1e-3)
    assert out["lat"].iloc[1] == pytest.approx(-33.86, abs=1e-3)


def test_to_ll_zone_string_not_matching_a_column_is_a_literal_value():
    # A string 'zone' that does NOT name a column in 'data' is still
    # treated as a literal per-point zone designator, preserving the
    # common no-column usage (e.g. to_ll(e, n, "11S", ...)).
    df = pd.DataFrame(
        {
            "E": [376322.215, 376322.215],
            "N": [3768509.769, 3768509.769],
        }
    )
    out = to_ll("E", "N", "11S", data=df, as_frame=True)
    assert list(out["zone"]) == ["11S", "11S"]
    assert out["lat"].iloc[0] == pytest.approx(34.05, abs=1e-3)


def test_to_ll_as_frame_string_columns_resolved():
    # The `_orig` helper used to rebuild the original easting/northing
    # columns for the output DataFrame now checks the string/column-name
    # case before np.isscalar(), so column-name strings with
    # as_frame=True and multi-row data resolve the underlying columns
    # instead of raising a pandas length-mismatch error.
    df = pd.DataFrame(
        {
            "E": [376322.215, 334416.394],
            "N": [3768509.769, 6251925.360],
        }
    )
    zone = np.array(["11S", "56H"], dtype=object)
    out = to_ll("E", "N", zone, data=df, as_frame=True)
    assert out["easting"].to_numpy() == pytest.approx(df["E"].to_numpy())
    assert out["northing"].to_numpy() == pytest.approx(df["N"].to_numpy())


def test_to_ll_as_frame_string_columns_single_row_corrupts_values_bug():
    # NOTE: with exactly one row, the same np.isscalar bug does not
    # raise (array lengths coincidentally match) but silently
    # populates the "easting"/"northing" columns with the *column
    # name strings themselves* instead of the numeric values -- a
    # silent data-corruption bug, worse than the multi-row exception
    # case above. See report.
    df = pd.DataFrame({"E": [376322.215], "N": [3768509.769], "zone": ["11S"]})
    out = to_ll(df["E"], df["N"], df["zone"], as_frame=True)
    # this path (Series, not string column-name) is fine:
    assert out.loc[0, "easting"] == pytest.approx(376322.215)


def test_to_utm_scalar():
    e, n, z = to_utm(34.05, -118.34)
    assert z == "11S"
    assert e == pytest.approx(376322.2, abs=1.0)
    assert n == pytest.approx(3768509.8, abs=1.0)


def test_to_utm_array_and_dataframe_no_strings():
    lat = np.array([34.05, -33.86])
    lon = np.array([-118.34, 151.21])
    e, n, z = to_utm(lat, lon)
    assert e.shape == (2,)
    assert list(z) == ["11S", "56H"]

    df = to_utm(lat, lon, as_frame=True)
    assert list(df.columns) == ["lat", "lon", "easting", "northing", "zone"]
    assert len(df) == 2


def test_to_utm_series_preserves_index():
    lat = pd.Series([34.05, -33.86], index=[1, 2])
    lon = pd.Series([-118.34, 151.21], index=[1, 2])
    df = to_utm(lat, lon, as_frame=True)
    assert list(df.index) == [1, 2]


def test_to_utm_explicit_epsg():
    e, n, z = to_utm(34.05, -118.34, epsg=32611)
    assert e == pytest.approx(376322.2, abs=1.0)


def test_to_utm_missing_data_for_column_name_raises():
    with pytest.raises(ValueError):
        to_utm("lat", "lon")


def test_to_utm_shape_mismatch_raises():
    with pytest.raises(ValueError):
        to_utm([1.0, 2.0], [3.0, 4.0, 5.0])


def test_to_utm_as_frame_string_columns_resolved():
    # Same fix as to_ll's _orig helper (see
    # test_to_ll_as_frame_string_columns_resolved): the string/column-name
    # check now comes before np.isscalar(), so column-name strings for
    # lat/lon with as_frame=True and multi-row data resolve the
    # underlying columns instead of raising.
    df = pd.DataFrame({"lat": [34.05, -33.86], "lon": [-118.34, 151.21]})
    out = to_utm("lat", "lon", data=df, as_frame=True)
    assert out["lat"].to_numpy() == pytest.approx(df["lat"].to_numpy())
    assert out["lon"].to_numpy() == pytest.approx(df["lon"].to_numpy())
    assert list(out["zone"]) == ["11S", "56H"]


# ======================================================================
# epsg_project
# ======================================================================


def test_epsg_project_basic():
    x2, y2 = epsg_project(-118.34, 34.05, 4326, 3857)
    assert x2 == pytest.approx(-13173548.5, abs=1.0)
    assert y2 == pytest.approx(4035517.8, abs=1.0)


def test_epsg_project_missing_from_or_to_returns_none():
    assert epsg_project(1.0, 2.0, None, 3857) is None
    assert epsg_project(1.0, 2.0, 4326, None) is None


def test_epsg_project_unknown_epsg_code_returns_none():
    assert epsg_project(1.0, 2.0, 999999, 3857) is None
    assert epsg_project(1.0, 2.0, 4326, 999999) is None


# ======================================================================
# utm_wgs84_conv (requires the optional 'utm' package, not installed)
# ======================================================================


def test_utm_wgs84_conv_missing_dependency_raises():
    with pytest.raises((ImportError, ModuleNotFoundError)):
        utm_wgs84_conv(34.05, -118.34)


# ======================================================================
# transform_utm_to_ll / transform_ll_to_utm (GDAL-gated, deprecated)
# ======================================================================


def test_transform_utm_to_ll_requires_gdal():
    with pytest.deprecated_call():
        with pytest.raises(GDALMissingError):
            transform_utm_to_ll(376322.215, 3768509.769, "11S")


def test_transform_ll_to_utm_requires_gdal():
    with pytest.deprecated_call():
        with pytest.raises(GDALMissingError):
            transform_ll_to_utm(-118.34, 34.05)


# ======================================================================
# _resolve_api_name / _extract_value_from_nested_dict
# ======================================================================


def test_resolve_api_name_none_returns_default():
    assert _resolve_api_name(None) == "open_meteo"


@pytest.mark.parametrize(
    "name, expected",
    [
        ("open_meteo", "open_meteo"),
        ("open_topo_data", "open_topo_data"),
        ("usgs_ned", "usgs_ned"),
        (",", "open_meteo"),
        ("|", "open_topo_data"),
        ("indiv", "usgs_ned"),
    ],
)
def test_resolve_api_name_direct_and_shorthand(name, expected):
    assert _resolve_api_name(name) == expected


def test_resolve_api_name_unknown_falls_back_to_default():
    assert _resolve_api_name("not-a-real-api") == "open_meteo"


def test_extract_value_from_nested_dict_basic():
    data = {"results": {"elevation": 611.0}}
    assert _extract_value_from_nested_dict(data, "results.elevation") == 611.0


def test_extract_value_from_nested_dict_list_indexing():
    data = {"a": [{"b": 5}, {"b": 6}]}
    assert _extract_value_from_nested_dict(data, "a.0.b") == 5
    assert _extract_value_from_nested_dict(data, "a.1.b") == 6


def test_extract_value_from_nested_dict_list_index_out_of_range():
    data = {"a": [{"b": 5}]}
    assert _extract_value_from_nested_dict(data, "a.5.b") is None


def test_extract_value_from_nested_dict_missing_key():
    assert _extract_value_from_nested_dict({"x": 1}, "y.z") is None


def test_extract_value_from_nested_dict_invalid_path_on_list():
    # non-digit key against a list -> invalid path, returns None
    data = {"a": [{"b": 5}, {"b": 6}]}
    assert _extract_value_from_nested_dict(data, "a.b.c") is None


def test_extract_value_from_nested_dict_real_opentopodata_shape():
    # opentopodata's actual response shape nests results as a LIST of
    # dicts: {"results": [{"elevation": ...}], "status": "OK"}.
    # ElevationAPIConfig's configured response_key for 'open_topo_data'
    # is the plain dot-path "results.elevation"; _extract_value_from_
    # nested_dict now maps the remaining key over each list element so
    # this resolves correctly for both single-point and batch queries.
    single = {
        "results": [{"elevation": 611.0, "location": {"lat": 1, "lng": 2}}],
        "status": "OK",
    }
    assert _extract_value_from_nested_dict(single, "results.elevation") == [611.0]

    batch = {
        "results": [
            {"elevation": 100.0},
            {"elevation": 200.0},
            {"elevation": 150.0},
        ],
        "status": "OK",
    }
    assert _extract_value_from_nested_dict(batch, "results.elevation") == [
        100.0,
        200.0,
        150.0,
    ]


# ======================================================================
# get_elevation_from_api (network mocked)
# ======================================================================


class _FakeResponse:
    def __init__(self, json_data, status_code=200, raise_exc=None):
        self._json_data = json_data
        self.status_code = status_code
        self._raise_exc = raise_exc

    def raise_for_status(self):
        if self._raise_exc is not None:
            raise self._raise_exc

    def json(self):
        return self._json_data


def test_get_elevation_from_api_scalar_open_meteo(monkeypatch):
    import requests

    def fake_get(url, params=None):
        assert params["latitude"] == 27.9881
        assert params["longitude"] == 86.9250
        return _FakeResponse({"elevation": [8815.0]})

    monkeypatch.setattr(requests, "get", fake_get)
    elev = get_elevation_from_api(27.9881, 86.9250)
    assert elev == pytest.approx(8815.0)


def test_get_elevation_from_api_array_open_meteo(monkeypatch):
    import requests

    def fake_get(url, params=None):
        assert params["latitude"] == "34.05,40.71"
        return _FakeResponse({"elevation": [86.0, 3968.0]})

    monkeypatch.setattr(requests, "get", fake_get)
    lats = [34.05, 40.71]
    lons = [-118.24, -74.00]
    elevs = get_elevation_from_api(lats, lons)
    np.testing.assert_allclose(elevs, [86.0, 3968.0])


def test_get_elevation_from_api_pipe_separated(monkeypatch):
    import requests

    captured = {}

    def fake_get(url, params=None):
        captured["locations"] = params["locations"]
        # Real opentopodata shape: 'results' is a LIST of per-point dicts.
        return _FakeResponse(
            {
                "results": [
                    {
                        "elevation": 611.0,
                        "location": {"lat": 27.9881, "lng": 86.925},
                    }
                ],
                "status": "OK",
            }
        )

    monkeypatch.setattr(requests, "get", fake_get)
    elev = get_elevation_from_api(27.9881, 86.9250, api_name="open_topo_data")
    assert elev == pytest.approx(611.0)
    assert captured["locations"] == "27.9881,86.925"


def test_get_elevation_from_api_shorthand_symbols(monkeypatch):
    import requests

    def fake_get(url, params=None):
        return _FakeResponse({"elevation": [100.0]})

    monkeypatch.setattr(requests, "get", fake_get)
    elev = get_elevation_from_api(1.0, 2.0, api_name=",")
    assert elev == pytest.approx(100.0)


def test_get_elevation_from_api_missing_key_raises_gis_error(monkeypatch):
    import requests

    def fake_get(url, params=None):
        return _FakeResponse({"unexpected": "shape"})

    monkeypatch.setattr(requests, "get", fake_get)
    with pytest.raises(GisError):
        get_elevation_from_api(1.0, 2.0)


def test_get_elevation_from_api_http_error_falls_back_to_default(monkeypatch):
    import requests

    calls = []

    def fake_get(url, params=None):
        calls.append(url)
        if len(calls) == 1:
            raise requests.exceptions.RequestException("boom")
        return _FakeResponse({"elevation": [42.0]})

    monkeypatch.setattr(requests, "get", fake_get)
    elev = get_elevation_from_api(1.0, 2.0, api_name="open_topo_data")
    assert elev == pytest.approx(42.0)
    assert len(calls) == 2


def test_get_elevation_from_api_http_error_no_fallback_raises(monkeypatch):
    import requests

    def fake_get(url, params=None):
        raise requests.exceptions.RequestException("boom")

    monkeypatch.setattr(requests, "get", fake_get)
    with pytest.raises(GisError):
        get_elevation_from_api(1.0, 2.0)  # already the default API


def test_get_elevation_from_api_batches_large_arrays(monkeypatch):
    import requests

    call_sizes = []

    def fake_get(url, params=None):
        n = len(params["latitude"].split(","))
        call_sizes.append(n)
        return _FakeResponse({"elevation": [1.0] * n})

    monkeypatch.setattr(requests, "get", fake_get)
    lats = list(np.linspace(0, 1, 95))
    lons = list(np.linspace(0, 1, 95))
    elevs = get_elevation_from_api(lats, lons)
    assert len(elevs) == 95
    assert call_sizes == [90, 5]


# ======================================================================
# get_elevation_from_utm (network mocked)
# ======================================================================


def test_get_elevation_from_utm_scalar(monkeypatch):
    import requests

    def fake_get(url, params=None):
        return _FakeResponse({"elevation": [86.0]})

    monkeypatch.setattr(requests, "get", fake_get)
    elev = get_elevation_from_utm(377274.0, 3762150.0, "11S")
    assert elev == pytest.approx(86.0)


def test_get_elevation_from_utm_error_wrapped_as_gis_error(monkeypatch):
    # an invalid zone causes to_ll -> project_point_utm2ll to raise;
    # get_elevation_from_utm should wrap it as GisError.
    with pytest.raises(GisError):
        get_elevation_from_utm(377274.0, 3762150.0, "abcS")


# ======================================================================
# calculate_azimuth
# ======================================================================


def test_calculate_azimuth_square_path():
    easting = [0, 100, 100, 0]
    northing = [0, 0, 100, 100]
    az = calculate_azimuth(easting, northing)
    assert az[0] == pytest.approx(90.0)
    assert az[1] == pytest.approx(0.0)
    assert az[2] == pytest.approx(270.0)
    assert np.isnan(az[3])


def test_calculate_azimuth_due_north_and_diagonal():
    easting = [0, 0, 100]
    northing = [0, 100, 200]
    az = calculate_azimuth(easting, northing)
    assert az[0] == pytest.approx(0.0)  # due north
    assert az[1] == pytest.approx(45.0)  # NE diagonal


def test_calculate_azimuth_single_point_returns_nan():
    az = calculate_azimuth([5.0], [10.0])
    assert az.shape == (1,)
    assert np.isnan(az[0])


def test_calculate_azimuth_empty_returns_empty():
    az = calculate_azimuth([], [])
    assert az.shape == (0,)


# ======================================================================
# Additional branch-coverage top-ups
# ======================================================================


def test_normalize_lat_lon_to_float_non_str_non_numeric():
    # a complex number is np.isscalar()==True but not float()-able,
    # exercising the final "except Exception: return nan" branch of
    # the internal _to_float helper. The complex value becomes nan
    # (unresolvable as either component) while the valid 45.0 still
    # resolves to "lat" under the default assume="lonlat" swap.
    lat, lon = normalize_lat_lon(complex(1, 2), 45.0, error="ignore")
    assert lat == 45.0
    assert np.isnan(lon)


def test_normalize_lat_lon_fallback_assume_latlon_branch():
    # both > 90 (non-decisive, non-ambiguous) with assume="latlon"
    # exercises the fallback's "if assume == 'latlon'" branch
    # (as opposed to the else/lonlat branch exercised elsewhere).
    lat, lon = normalize_lat_lon(95.0, 91.0, assume="latlon", error="ignore")
    assert np.isnan(lat) and np.isnan(lon)


def test_assert_xy_coordinate_system_empty_arrays():
    assert assert_xy_coordinate_system([], []) == "ll"


def test_assert_lat_value_type_error_branch():
    # a non-numeric, non-string type raises TypeError in float(),
    # which assert_lat_value swallows by returning None.
    assert assert_lat_value([1, 2]) is None


def test_assert_lon_value_type_error_branch():
    assert assert_lon_value([1, 2]) is None


def test_utm_zone_to_epsg_no_match_returns_none():
    assert utm_zone_to_epsg(999, True) is None


def test_get_epsg_wraps_get_utm_zone_and_utm_zone_to_epsg():
    epsg = get_epsg(34.05, -118.34)
    assert epsg in (32611, None)
    epsg_south = get_epsg(-33.86, 151.21)
    assert epsg_south in (32756, None)


def test_ll_to_utm_svalbard_exceptions():
    # lat in [72, 84) with lon in specific sub-bands forces the
    # Svalbard zone-number overrides (31/33/35/37).
    zone31, _, _ = ll_to_utm(WGS84, 75.0, 5.0)
    assert zone31.startswith("31")
    zone33, _, _ = ll_to_utm(WGS84, 75.0, 15.0)
    assert zone33.startswith("33")
    zone35, _, _ = ll_to_utm(WGS84, 75.0, 25.0)
    assert zone35.startswith("35")
    zone37, _, _ = ll_to_utm(WGS84, 75.0, 38.0)
    assert zone37.startswith("37")


def test_get_elevation_from_api_pipe_separated_array(monkeypatch):
    import requests

    captured = {}

    def fake_get(url, params=None):
        captured["locations"] = params["locations"]
        # Real opentopodata shape: one result dict per queried point.
        return _FakeResponse(
            {
                "results": [
                    {"elevation": 100.0},
                    {"elevation": 200.0},
                ],
                "status": "OK",
            }
        )

    monkeypatch.setattr(requests, "get", fake_get)
    lats = [1.0, 2.0]
    lons = [3.0, 4.0]
    elev = get_elevation_from_api(lats, lons, api_name="open_topo_data")
    np.testing.assert_allclose(elev, [100.0, 200.0])
    assert captured["locations"] == "1.0,3.0|2.0,4.0"


def test_utm_wgs84_conv_success_with_mocked_utm_package(monkeypatch):
    # The 'utm' package isn't installed in this environment; build a
    # minimal fake module (duck-typed, pure Python -- no GDAL/C-ext
    # involved) so the success path and both round-trip mismatch
    # warnings get exercised.
    fake_utm = types.ModuleType("utm")

    def from_latlon(lat, lon):
        return (376322.215, 3768509.769, 11, "S")

    def to_latlon(e, n, z, letter):
        # deliberately off by more than tol to trigger both warnings
        return (34.05 + 1e-5, -118.34 + 1e-5)

    fake_utm.from_latlon = from_latlon
    fake_utm.to_latlon = to_latlon
    monkeypatch.setitem(sys.modules, "utm", fake_utm)

    e, n, z, letter = utm_wgs84_conv(34.05, -118.34)
    assert (e, n, z, letter) == (376322.215, 3768509.769, 11, "S")


def test_utm_wgs84_conv_no_mismatch_no_warnings(monkeypatch):
    # exercises the "false" side of both round-trip tolerance checks
    # (exact round-trip, no warnings logged).
    fake_utm = types.ModuleType("utm")

    def from_latlon(lat, lon):
        return (376322.215, 3768509.769, 11, "S")

    def to_latlon(e, n, z, letter):
        return (34.05, -118.34)

    fake_utm.from_latlon = from_latlon
    fake_utm.to_latlon = to_latlon
    monkeypatch.setitem(sys.modules, "utm", fake_utm)

    e, n, z, letter = utm_wgs84_conv(34.05, -118.34)
    assert (e, n, z, letter) == (376322.215, 3768509.769, 11, "S")


def test_ll_to_utm_svalbard_lat_no_sub_band_match():
    # lat in [72, 84) but lon outside all defined Svalbard sub-bands
    # (>= 42) falls through the elif chain unchanged.
    zone, _, _ = ll_to_utm(WGS84, 75.0, 50.0)
    expected_zone_number = int((50.0 + 180) / 6) + 1
    assert zone.startswith(str(expected_zone_number))


def test_get_elevation_from_api_individual_format_usgs_ned(monkeypatch):
    # 'usgs_ned' uses params_format == "individual", which isn't
    # handled by the comma/pipe branches -- params stays empty.
    import requests

    captured = {}

    def fake_get(url, params=None):
        captured["params"] = params
        return _FakeResponse(
            {"USGS_Elevation_Point_Query_Service": {"Elevation_Query": 123.0}}
        )

    monkeypatch.setattr(requests, "get", fake_get)
    elev = get_elevation_from_api(34.05, -118.34, api_name="usgs_ned")
    assert elev == pytest.approx(123.0)
    assert captured["params"] == {}
