# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.format.topo_source.

Coordinate conversion is checked against ground-truth pyproj output
directly (not just internal consistency), since this module's own
development caught a real lat/lon-swap bug in
pycsamt.gis.utils.project_point_utm2ll (see
fix_gis_to_ll_latlon_swap_bug in project memory / the git history) --
a real UTM-33N (most of Europe) test point is exactly what surfaced it.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from pycsamt.format.topo_source import (
    TopoAttribution,
    TopoTable,
    attribute_topo,
    read_topo_file,
    resolve_topo,
    topo_from_sites,
)

# UTM zone 33N (EPSG:32633), a point squarely inside "most of Europe" --
# the exact class of point the to_ll swap bug affected. Ground truth via
# a direct pyproj.Transformer check (see the bug-fix memory entry).
_UTM33N_EPSG = 32633
_UTM33N_POINTS = [
    (500000.0, 4500000.0, 15.0, 40.650857),
    (500100.0, 4500000.0, 15.001183, 40.650857),
    (500200.0, 4500000.0, 15.002366, 40.650856),
]


def _write_bln(path, points, *, flag=1):
    lines = [f"{len(points)},{flag}"]
    lines += [",".join(str(v) for v in row) for row in points]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class TestReadBln:
    def test_utm_points_convert_to_correct_lon_lat(self, tmp_path):
        rows = [
            (e, n, 100.0 + 10 * i) for i, (e, n, _, _) in enumerate(_UTM33N_POINTS)
        ]
        path = _write_bln(tmp_path / "line.bln", rows)
        table = read_topo_file(path, epsg=_UTM33N_EPSG)
        assert table.names is None
        assert table.n == 3
        expected_lon = [p[2] for p in _UTM33N_POINTS]
        expected_lat = [p[3] for p in _UTM33N_POINTS]
        np.testing.assert_allclose(table.lon, expected_lon, atol=1e-5)
        np.testing.assert_allclose(table.lat, expected_lat, atol=1e-5)
        np.testing.assert_allclose(table.elevation, [100.0, 110.0, 120.0])

    def test_latlon_true_skips_conversion(self, tmp_path):
        path = _write_bln(tmp_path / "line.bln", [(7.5, 45.1), (7.6, 45.2)])
        table = read_topo_file(path, latlon=True)
        np.testing.assert_allclose(table.lon, [7.5, 7.6])
        np.testing.assert_allclose(table.lat, [45.1, 45.2])
        assert table.elevation is None

    def test_missing_epsg_and_utm_zone_raises(self, tmp_path):
        path = _write_bln(tmp_path / "line.bln", [(500000.0, 4500000.0)])
        with pytest.raises(ValueError, match="epsg"):
            read_topo_file(path)

    def test_header_count_mismatch_raises(self, tmp_path):
        path = tmp_path / "bad.bln"
        path.write_text("3,1\n1,2\n3,4\n", encoding="utf-8")
        with pytest.raises(ValueError, match="declares 3"):
            read_topo_file(path, latlon=True)

    def test_malformed_header_raises(self, tmp_path):
        path = tmp_path / "bad.bln"
        path.write_text("not-a-number,1\n1,2\n", encoding="utf-8")
        with pytest.raises(ValueError, match="malformed"):
            read_topo_file(path, latlon=True)

    def test_unrecognised_extension_raises(self, tmp_path):
        path = tmp_path / "topo.xyz"
        path.write_text("1,2\n", encoding="utf-8")
        with pytest.raises(ValueError, match="unrecognised topo file extension"):
            read_topo_file(path)


class TestReadCsv:
    def test_leading_comment_lines_are_skipped(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text(
            "# a comment line\n# another one\nstation,lat,lon\nS0,45.1,7.5\n",
            encoding="utf-8",
        )
        table = read_topo_file(path)
        assert table.names == ["S0"]
        np.testing.assert_allclose(table.lon, [7.5])

    def test_named_latlon_columns(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text(
            "station,lat,lon,elev\nS0,45.10,7.50,1200\nS1,45.11,7.51,1195\n",
            encoding="utf-8",
        )
        table = read_topo_file(path)
        assert table.names == ["S0", "S1"]
        np.testing.assert_allclose(table.lon, [7.50, 7.51])
        np.testing.assert_allclose(table.lat, [45.10, 45.11])
        np.testing.assert_allclose(table.elevation, [1200, 1195])

    def test_named_easting_northing_columns_need_conversion(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text(
            "name,easting,northing\n"
            f"S0,{_UTM33N_POINTS[0][0]},{_UTM33N_POINTS[0][1]}\n"
            f"S1,{_UTM33N_POINTS[1][0]},{_UTM33N_POINTS[1][1]}\n",
            encoding="utf-8",
        )
        table = read_topo_file(path, epsg=_UTM33N_EPSG)
        np.testing.assert_allclose(
            table.lon, [_UTM33N_POINTS[0][2], _UTM33N_POINTS[1][2]], atol=1e-5
        )
        np.testing.assert_allclose(
            table.lat, [_UTM33N_POINTS[0][3], _UTM33N_POINTS[1][3]], atol=1e-5
        )

    def test_no_name_column_is_positional(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text("lat,lon\n45.1,7.5\n45.2,7.6\n", encoding="utf-8")
        table = read_topo_file(path)
        assert table.names is None
        assert table.n == 2

    def test_no_coordinate_columns_raises(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text("station,note\nS0,foo\n", encoding="utf-8")
        with pytest.raises(ValueError, match="could not find"):
            read_topo_file(path)


class TestReadStn:
    def test_real_zonge_stn_fixture(self):
        from pathlib import Path

        stn_path = Path(__file__).parents[3] / "data" / "avg" / "K1.stn"
        if not stn_path.exists():
            pytest.skip(f"bundled .stn fixture not found: {stn_path}")
        table = read_topo_file(stn_path, epsg=_UTM33N_EPSG)
        assert table.names is not None
        assert table.n > 0
        assert table.elevation is not None
        # every real K1.stn row has a non-degenerate elevation
        assert np.all(np.isfinite(table.elevation))


class TestAttributeTopo:
    def test_named_source_matches_by_id_normalized_fallback(self):
        table = TopoTable(
            lon=np.array([7.5, 7.6]),
            lat=np.array([45.1, 45.2]),
            elevation=np.array([1200.0, 1100.0]),
            names=["23-18-001", "23_18_002"],
            source="test",
        )
        attr = attribute_topo(table, ["23-18-001", "23-18-002", "23-18-999"])
        assert attr.lon["23-18-001"] == pytest.approx(7.5)
        assert attr.lon["23-18-002"] == pytest.approx(7.6)  # normalized match
        assert "23-18-999" in attr.unmatched_stations

    def test_positional_source_exact_count(self):
        table = TopoTable(lon=np.array([1.0, 2.0]), lat=np.array([3.0, 4.0]), source="test")
        attr = attribute_topo(table, ["A", "B"])
        assert attr.lon == {"A": 1.0, "B": 2.0}
        assert attr.lat == {"A": 3.0, "B": 4.0}
        assert attr.unmatched_stations == []

    def test_positional_count_mismatch_raises_by_default(self):
        table = TopoTable(lon=np.array([1.0, 2.0, 3.0]), lat=np.array([1.0, 2.0, 3.0]), source="test.bln")
        with pytest.raises(ValueError, match="has 3 point"):
            attribute_topo(table, ["A", "B"])

    def test_positional_count_mismatch_warns_and_truncates(self):
        table = TopoTable(lon=np.array([1.0, 2.0, 3.0]), lat=np.array([1.0, 2.0, 3.0]), source="test.bln")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            attr = attribute_topo(table, ["A", "B"], on_mismatch="warn")
        assert any(issubclass(w.category, UserWarning) for w in caught)
        assert attr.lon == {"A": 1.0, "B": 2.0}

    def test_invalid_on_mismatch_raises(self):
        table = TopoTable(lon=np.array([1.0]), lat=np.array([1.0]), source="test")
        with pytest.raises(ValueError, match="on_mismatch"):
            attribute_topo(table, ["A"], on_mismatch="bogus")


class TestTopoFromSites:
    def test_extracts_from_station_record_like_objects(self):
        from types import SimpleNamespace

        stations = [
            SimpleNamespace(id="S0", longitude=7.5, latitude=45.1, elevation=1200.0),
            SimpleNamespace(id="S1", longitude=7.6, latitude=45.2, elevation=None),
            SimpleNamespace(id="S2", longitude=None, latitude=None, elevation=None),
        ]
        table = topo_from_sites(stations)
        assert table.names == ["S0", "S1"]
        np.testing.assert_allclose(table.lon, [7.5, 7.6])
        assert np.isnan(table.elevation[1])

    def test_mapdata_like_wrapper_via_stations_attribute(self):
        from types import SimpleNamespace

        mapdata = SimpleNamespace(
            stations=[SimpleNamespace(id="S0", longitude=1.0, latitude=2.0, elevation=3.0)]
        )
        table = topo_from_sites(mapdata)
        assert table.names == ["S0"]

    def test_no_geolocated_stations_raises(self):
        from types import SimpleNamespace

        with pytest.raises(ValueError, match="no geo-located"):
            topo_from_sites([SimpleNamespace(id="S0", longitude=None, latitude=None)])


class TestResolveTopoFlat:
    def test_none_returns_empty_attribution(self):
        attr = resolve_topo(None, ["A", "B"])
        assert isinstance(attr, TopoAttribution)
        assert attr.lon == {}
        assert attr.unmatched_stations == ["A", "B"]

    def test_plain_mapping_passthrough(self):
        attr = resolve_topo({"A": (7.5, 45.1, 1200.0)}, ["A", "B"])
        assert attr.lon == {"A": 7.5}
        assert attr.elevation == {"A": 1200.0}
        assert attr.unmatched_stations == ["B"]

    def test_file_path_source(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text("station,lat,lon\nA,45.1,7.5\n", encoding="utf-8")
        attr = resolve_topo(path, ["A", "B"])
        assert attr.lon == {"A": 7.5}
        assert attr.unmatched_stations == ["B"]


class TestResolveTopoMultiline:
    def test_line_keyed_mapping_of_files(self, tmp_path):
        p1 = tmp_path / "l1.csv"
        p1.write_text("station,lat,lon\nA,45.1,7.5\nB,45.2,7.6\n", encoding="utf-8")
        p2 = tmp_path / "l2.csv"
        p2.write_text("station,lat,lon\nC,46.1,8.5\n", encoding="utf-8")
        attr = resolve_topo(
            {"L1": p1, "L2": p2}, {"L1": ["A", "B"], "L2": ["C", "D"]}
        )
        assert attr.lon == {"A": 7.5, "B": 7.6, "C": 8.5}
        assert attr.unmatched_stations == ["D"]

    def test_line_keyed_mapping_missing_line_marks_unmatched(self, tmp_path):
        p1 = tmp_path / "l1.csv"
        p1.write_text("station,lat,lon\nA,45.1,7.5\n", encoding="utf-8")
        attr = resolve_topo({"L1": p1}, {"L1": ["A"], "L2": ["C"]})
        assert attr.lon == {"A": 7.5}
        assert "C" in attr.unmatched_stations

    def test_positional_list_one_bln_per_line(self, tmp_path):
        p1 = _write_bln(tmp_path / "l1.bln", [(7.5, 45.1), (7.6, 45.2)])
        p2 = _write_bln(tmp_path / "l2.bln", [(8.5, 46.1)])
        attr = resolve_topo(
            [p1, p2], {"L1": ["A", "B"], "L2": ["C"]}, latlon=True
        )
        assert attr.lon == {"A": 7.5, "B": 7.6, "C": 8.5}

    def test_positional_list_length_mismatch_raises(self, tmp_path):
        p1 = _write_bln(tmp_path / "l1.bln", [(7.5, 45.1)], flag=1)
        with pytest.raises(ValueError, match="got 1 topo source"):
            resolve_topo([p1], {"L1": ["A"], "L2": ["B"]}, latlon=True)

    def test_single_named_source_spans_all_lines(self, tmp_path):
        path = tmp_path / "all.csv"
        path.write_text(
            "station,lat,lon\nA,45.1,7.5\nB,45.2,7.6\nC,46.1,8.5\n", encoding="utf-8"
        )
        attr = resolve_topo(path, {"L1": ["A", "B"], "L2": ["C"]})
        assert attr.lon == {"A": 7.5, "B": 7.6, "C": 8.5}

    def test_none_topo_multiline_all_unmatched(self):
        attr = resolve_topo(None, {"L1": ["A", "B"], "L2": ["C"]})
        assert attr.lon == {}
        assert set(attr.unmatched_stations) == {"A", "B", "C"}
