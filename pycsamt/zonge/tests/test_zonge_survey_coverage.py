# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import ProcessingError, StationError
from pycsamt.zonge.survey import Station, Topography


def _utm_df():
    return pd.DataFrame(
        {
            "station": [0, 1, 2, 3],
            "easting": [500000.0, 500100.0, 500200.0, 500300.0],
            "northing": [4000000.0, 4000000.0, 4000000.0, 4000000.0],
            "elevation": [10.0, 11.0, 12.0, 13.0],
        }
    )


# ─────────────────────────────────────────────────────────────────────────
# Topography.read: type validation
# ─────────────────────────────────────────────────────────────────────────


def test_read_accepts_dataframe_directly():
    topo = Topography().read(_utm_df())
    assert not topo._frame.empty
    assert len(topo._frame) == 4


def test_read_raises_typeerror_for_unsupported_source():
    with pytest.raises(TypeError, match="file path or DataFrame"):
        Topography().read(12345)


def test_normalize_stn_columns_recognizes_heading_pitch_roll():
    df = pd.DataFrame(
        {
            "Station": [0, 1],
            "Easting": [1.0, 2.0],
            "Northing": [1.0, 2.0],
            "Elev": [1.0, 2.0],
            "Heading": [10.0, 20.0],
            "Pitch": [1.0, 2.0],
            "Roll": [0.5, 0.6],
        }
    )
    topo = Topography().read(df)
    for col in ("heading", "pitch", "roll"):
        assert col in topo._frame.columns


def test_normalize_stn_columns_raises_when_required_missing():
    df = pd.DataFrame({"station": [1, 2], "easting": [1.0, 2.0]})
    with pytest.raises(ProcessingError, match="missing required columns"):
        Topography().read(df)


# ─────────────────────────────────────────────────────────────────────────
# convert_coords
# ─────────────────────────────────────────────────────────────────────────


def test_convert_coords_raises_when_empty():
    with pytest.raises(ProcessingError, match="not been loaded"):
        Topography().convert_coords()


def test_convert_coords_noop_when_target_matches_current_system(caplog):
    topo = Topography(data=_utm_df(), verbose=True)
    out = topo.convert_coords(to="utm")
    assert out is None  # no-op; already utm


def test_median_inc_returns_none_for_fewer_than_two_values():
    from pycsamt.zonge.survey import _median_inc

    assert _median_inc(np.array([])) is None
    assert _median_inc(np.array([1.0])) is None


def test_convert_coords_to_ll_without_zone_or_epsg_logs_error_and_zeros():
    topo = Topography(data=_utm_df())
    topo.convert_coords(to="ll", inplace=True)
    assert np.allclose(topo.latitude, 0.0)
    assert np.allclose(topo.longitude, 0.0)


def test_convert_coords_to_ll_with_epsg_updates_inplace():
    topo = Topography(data=_utm_df(), epsg=32612)  # UTM zone 12N
    out = topo.convert_coords(to="ll", inplace=True)
    assert out is None
    assert "latitude" in topo._frame.columns
    assert "longitude" in topo._frame.columns
    assert not np.allclose(topo.latitude, 0.0)


def test_convert_coords_to_ll_not_inplace_returns_new_frame():
    topo = Topography(data=_utm_df(), epsg=32612)
    out = topo.convert_coords(to="ll", inplace=False)
    assert isinstance(out, pd.DataFrame)
    assert "latitude" in out.columns
    # original frame untouched
    assert "latitude" not in topo._frame.columns


def test_convert_coords_to_utm_from_ll():
    df = pd.DataFrame(
        {
            "station": [0, 1],
            "easting": [40.0, 40.001],  # actually lat when detected as ll
            "northing": [-110.0, -110.0],
            "elevation": [0.0, 0.0],
        }
    )
    topo = Topography(data=df, epsg=32612)
    out = topo.convert_coords(to="utm", inplace=False)
    assert "easting" in out.columns and "utm_zone" in out.columns


def test_convert_coords_raises_for_invalid_target():
    topo = Topography(data=_utm_df())
    with pytest.raises(ValueError, match="Invalid target system"):
        topo.convert_coords(to="bogus")


# ─────────────────────────────────────────────────────────────────────────
# to_grid / get_step
# ─────────────────────────────────────────────────────────────────────────


def test_to_grid_raises_when_empty():
    with pytest.raises(ProcessingError, match="not been loaded"):
        Topography().to_grid()


def test_get_step_empty_when_no_data():
    out = Topography().get_step()
    assert out.empty


def test_get_step_empty_when_single_row():
    df = _utm_df().iloc[:1]
    topo = Topography(data=df)
    assert topo.get_step().empty


# ─────────────────────────────────────────────────────────────────────────
# correct_coords
# ─────────────────────────────────────────────────────────────────────────


def test_correct_coords_warns_when_insufficient_data():
    with pytest.warns(UserWarning, match="Not enough data"):
        out = Topography().correct_coords()
    assert out is None


def test_correct_coords_auto_step_verbose_logs():
    topo = Topography(data=_utm_df(), verbose=True)
    out = topo.correct_coords()  # step=None -> auto-detected
    assert out is None
    assert len(topo._frame) == 4


def test_correct_coords_not_inplace_returns_new_frame():
    topo = Topography(data=_utm_df())
    original_easting = topo.easting.copy()
    out = topo.correct_coords(step=50.0, inplace=False)
    assert isinstance(out, pd.DataFrame)
    assert np.allclose(topo.easting, original_easting)  # untouched


# ─────────────────────────────────────────────────────────────────────────
# generate
# ─────────────────────────────────────────────────────────────────────────


def test_generate_raises_for_nonpositive_n_stations():
    with pytest.raises(ValueError, match="must be positive"):
        Topography.generate(
            start_coord=(0.0, 0.0), n_stations=0, step=10.0, azimuth=0.0
        )


def test_generate_ll_pathway_with_geopy_stubbed(monkeypatch):
    import sys
    import types

    fake_geopy = types.ModuleType("geopy")
    fake_distance = types.ModuleType("geopy.distance")

    class _FakePoint:
        def __init__(self, latitude, longitude):
            self.latitude = latitude
            self.longitude = longitude

    class _FakeGeodesic:
        def __init__(self, meters):
            self.meters = meters

        def destination(self, point, bearing):
            lat, lon = point
            return _FakePoint(lat + 0.001, lon)

    fake_distance.geodesic = _FakeGeodesic
    fake_geopy.distance = fake_distance
    monkeypatch.setitem(sys.modules, "geopy", fake_geopy)
    monkeypatch.setitem(sys.modules, "geopy.distance", fake_distance)

    from pycsamt.zonge import survey as survey_mod

    monkeypatch.setattr(
        survey_mod, "import_optional_dependency", lambda *a, **k: fake_geopy
    )

    topo = Topography.generate(
        start_coord=(40.0, -110.0),
        n_stations=3,
        step=100.0,
        azimuth=0.0,
        coord_type="ll",
    )
    assert isinstance(topo, Topography)
    assert len(topo._frame) == 3
    assert "easting" in topo._frame.columns


# ─────────────────────────────────────────────────────────────────────────
# write
# ─────────────────────────────────────────────────────────────────────────


def test_topography_write_empty_returns_empty_list():
    assert Topography().write() == []


def test_topography_write_nonempty_returns_header_and_rows():
    topo = Topography(data=_utm_df())
    lines = topo.write()
    assert lines[0] == "station,easting,northing,elevation"
    assert len(lines) == 1 + len(topo._frame)


# ─────────────────────────────────────────────────────────────────────────
# get_elevation_from
# ─────────────────────────────────────────────────────────────────────────


def test_get_elevation_from_raises_when_empty():
    with pytest.raises(ProcessingError, match="not been loaded"):
        Topography().get_elevation_from()


def test_get_elevation_from_utm_requires_zone():
    topo = Topography(data=_utm_df())
    with pytest.raises(ValueError, match="zone.*must be provided"):
        topo.get_elevation_from(from_="utm")


def test_get_elevation_from_utm_uses_zone_and_calls_helper(monkeypatch):
    from pycsamt.zonge import survey as survey_mod

    monkeypatch.setattr(
        survey_mod,
        "get_elevation_from_utm",
        lambda easting, northing, zone, datum: np.array([1.0, 2.0, 3.0, 4.0]),
    )
    topo = Topography(data=_utm_df(), utm_zone="12N")
    out = topo.get_elevation_from(from_="utm")
    assert np.allclose(out, [1.0, 2.0, 3.0, 4.0])


def test_get_elevation_from_api_converts_from_utm_when_no_latlon(monkeypatch):
    from pycsamt.zonge import survey as survey_mod

    monkeypatch.setattr(
        survey_mod,
        "get_elevation_from_api",
        lambda latitude, longitude: np.zeros(len(latitude)),
    )
    topo = Topography(data=_utm_df(), epsg=32612, verbose=True)
    out = topo.get_elevation_from(from_="api")
    assert len(out) == 4


def test_get_elevation_from_api_uses_existing_latlon(monkeypatch):
    from pycsamt.zonge import survey as survey_mod

    called = {}

    def _fake(latitude, longitude):
        called["lat"] = latitude
        return np.ones(len(latitude))

    monkeypatch.setattr(survey_mod, "get_elevation_from_api", _fake)
    df = _utm_df()
    df["latitude"] = [1.0, 2.0, 3.0, 4.0]
    df["longitude"] = [5.0, 6.0, 7.0, 8.0]
    topo = Topography(data=df)
    out = topo.get_elevation_from(from_="api")
    assert len(out) == 4
    assert "lat" in called


def test_get_elevation_from_raises_for_invalid_from():
    topo = Topography(data=_utm_df())
    with pytest.raises(ValueError, match="Invalid 'from_'"):
        topo.get_elevation_from(from_="bogus")


# ─────────────────────────────────────────────────────────────────────────
# get_azimuth / azimuth / bearing / properties
# ─────────────────────────────────────────────────────────────────────────


def test_get_azimuth_none_mode_returns_full_array():
    topo = Topography(data=_utm_df())
    out = topo.get_azimuth()
    assert isinstance(out, np.ndarray)


def test_get_azimuth_mean_and_median():
    topo = Topography(data=_utm_df())
    mean_az = topo.get_azimuth(mode="mean")
    median_az = topo.get_azimuth(mode="median")
    assert isinstance(mean_az, float)
    assert isinstance(median_az, float)


def test_get_azimuth_invalid_mode_raises():
    topo = Topography(data=_utm_df())
    with pytest.raises(ValueError, match="Invalid mode"):
        topo.get_azimuth(mode="bogus")


def test_get_azimuth_all_nan_returns_nan_for_mean():
    df = _utm_df()
    topo = Topography(data=df.iloc[:1])  # single row -> empty azimuth array
    out = topo.get_azimuth(mode="mean")
    assert np.isnan(out) or (isinstance(out, np.ndarray) and out.size == 0)


def test_azimuth_property_cached_and_single_row_case():
    topo = Topography(data=_utm_df().iloc[:1])
    assert topo.azimuth.size == 0
    # bearing is an alias
    assert topo.bearing.size == 0


def test_azimuth_property_caches_after_first_access():
    topo = Topography(data=_utm_df())
    first = topo.azimuth
    second = topo.azimuth
    assert first is second  # cached


def test_elevation_alias_property():
    topo = Topography(data=_utm_df())
    assert np.allclose(topo.elevation, topo.elevations)


def test_latitude_longitude_default_zeros_when_absent():
    topo = Topography(data=_utm_df())
    assert np.allclose(topo.latitude, 0.0)
    assert np.allclose(topo.longitude, 0.0)


# ─────────────────────────────────────────────────────────────────────────
# Station: remaining branches
# ─────────────────────────────────────────────────────────────────────────


def test_station_read_raises_when_station_column_missing_in_dataframe():
    df = pd.DataFrame({"freq": [1.0, 2.0]})
    with pytest.raises(StationError, match="column 'station' missing"):
        Station().read(df)


def test_station_read_raises_when_array_empty():
    with pytest.raises(StationError, match="empty station array"):
        Station().read(np.array([]))


def test_station_read_raises_when_no_valid_values():
    df = pd.DataFrame({"station": [np.nan, np.nan]})
    with pytest.raises(StationError, match="no valid station values"):
        Station().read(df)


def test_station_read_km_unit_conversion():
    st = Station()
    st.read([1.0, 2.0], unit="km")
    assert np.allclose(sorted(st._frame["station_m"]), [1000.0, 2000.0])


def test_station_read_raises_when_names_length_mismatch():
    df = pd.DataFrame({"station": [0.0, 1.0, 2.0]})
    with pytest.raises(StationError, match="length must match"):
        Station().read(df, names=["only_one"])


def test_station_write_empty_returns_empty_list():
    assert Station().write() == []


def test_station_span_none_when_empty():
    assert Station().span is None


def test_station_to_keywords_empty_when_no_values():
    assert Station().to_keywords() == {}
