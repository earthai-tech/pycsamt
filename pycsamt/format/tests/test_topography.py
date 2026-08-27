# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 6 tests — per_station topography wired through pycsamt.map.topo.

Includes a real-data check against the bundled ``data/AMT/WILLY_DATA``
EDI survey (128 real, EDI-derived station elevations).
"""

from __future__ import annotations

import base64
from pathlib import Path

import numpy as np
import pytest

from pycsamt.format import read_pcsf, write_pcsf
from pycsamt.format.schema import (
    Grid2DGeometry,
    PCSFModel,
    TopographyPerStation,
    TopographyRaster,
)
from pycsamt.format.topography import (
    topography_from_elevation_file,
    topography_from_grid,
    topography_from_map_data,
    topography_raster_to_grid,
    topography_to_elev_map,
)
from pycsamt.map._core import MapData, StationRecord

_WILLY_DIR = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA"
_SKIP_WILLY = pytest.mark.skipif(
    not _WILLY_DIR.exists(), reason=f"bundled WILLY_DATA not found: {_WILLY_DIR}"
)


def _map_data() -> MapData:
    return MapData(
        sites=[],
        stations=(
            StationRecord("A0", 1.0, 2.0, 500.0, "L1", 0),
            StationRecord("A1", 1.1, 2.1, 650.0, "L1", 1),
            StationRecord("A2", 1.2, 2.2, None, "L1", 2),
        ),
    )


class TestTopographyFromMapData:
    def test_only_stations_with_elevation_are_kept(self):
        topo = topography_from_map_data(_map_data())
        assert topo.station_id == ["A0", "A1"]
        np.testing.assert_allclose(topo.elevation, [500.0, 650.0])

    def test_returns_none_when_no_elevation_anywhere(self):
        data = MapData(
            sites=[],
            stations=(StationRecord("A0", 1.0, 2.0, None, "L1", 0),),
        )
        assert topography_from_map_data(data) is None

    @_SKIP_WILLY
    def test_real_willy_data_elevations(self):
        from pycsamt.map import load_lines

        data = load_lines(_WILLY_DIR, detect="folder")
        topo = topography_from_map_data(data)
        assert topo is not None
        assert len(topo.station_id) == len(data.stations)
        assert np.all(np.isfinite(topo.elevation))
        # A real, specific station from this bundled survey.
        assert "18-001A" in topo.station_id
        idx = topo.station_id.index("18-001A")
        assert topo.elevation[idx] == pytest.approx(99.0)


class TestTopographyFromElevationFile:
    def test_parses_csv_upload(self):
        csv = "station,elevation\nA0,500\nA1,650\n"
        b64 = "data:text/csv;base64," + base64.b64encode(csv.encode()).decode()
        topo = topography_from_elevation_file(b64, "topo.csv")
        assert topo.station_id == ["A0", "A1"]
        np.testing.assert_allclose(topo.elevation, [500.0, 650.0])

    def test_recognises_station_names_column(self):
        # Regression guard: map3d.py's own former parser did not
        # recognise this column name; pycsamt.map.topo's always did.
        csv = "station_names,elevation\nA0,500\n"
        b64 = "data:text/csv;base64," + base64.b64encode(csv.encode()).decode()
        topo = topography_from_elevation_file(b64, "topo.csv")
        assert topo.station_id == ["A0"]

    def test_returns_none_for_unparseable_file(self):
        csv = "foo,bar\n1,2\n"
        b64 = "data:text/csv;base64," + base64.b64encode(csv.encode()).decode()
        assert topography_from_elevation_file(b64, "topo.csv") is None


class TestTopographyToElevMap:
    def test_round_trips_to_a_plain_dict(self):
        topo = TopographyPerStation(
            station_id=["A0", "A1"], elevation=np.array([500.0, 650.0])
        )
        assert topography_to_elev_map(topo) == {"A0": 500.0, "A1": 650.0}


class TestPcsfTopographyRoundTrip:
    def test_write_read_preserves_topography(self, tmp_path):
        topo = topography_from_map_data(_map_data())
        model = PCSFModel(
            geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([1.0, 2.0])),
            resistivity=np.ones((2, 2)) * 100.0,
            topography=topo,
        )
        path = write_pcsf(model, tmp_path / "topo.pcsf")
        restored = read_pcsf(path)
        assert restored.topography.station_id == topo.station_id
        np.testing.assert_array_equal(restored.topography.elevation, topo.elevation)


class TestTopographyFromGrid:
    def test_builds_expected_shape(self):
        x = np.linspace(0.0, 500.0, 6)
        y = np.linspace(0.0, 300.0, 4)
        elevation = 100.0 + 0.01 * np.add.outer(y, x)
        topo = topography_from_grid(x, y, elevation)
        assert topo.kind == "raster"
        assert topo.elevation.shape == (4, 6)
        np.testing.assert_allclose(topo.x, x)
        np.testing.assert_allclose(topo.y, y)

    def test_round_trips_through_grid_helpers(self):
        x = np.array([0.0, 10.0, 20.0])
        y = np.array([0.0, 5.0])
        elevation = np.array([[100.0, 101.0, 102.0], [103.0, 104.0, 105.0]])
        topo = topography_from_grid(x, y, elevation)
        rx, ry, relev = topography_raster_to_grid(topo)
        np.testing.assert_allclose(rx, x)
        np.testing.assert_allclose(ry, y)
        np.testing.assert_allclose(relev, elevation)


class TestPcsfTopographyRasterRoundTrip:
    def test_write_read_preserves_raster_topography(self, tmp_path):
        x = np.array([0.0, 250.0, 500.0])
        y = np.array([0.0, 300.0])
        elevation = np.array([[100.0, 110.0, 120.0], [130.0, 140.0, 150.0]])
        topo = topography_from_grid(x, y, elevation)
        model = PCSFModel(
            geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([1.0, 2.0])),
            resistivity=np.ones((2, 2)) * 100.0,
            topography=topo,
        )
        path = write_pcsf(model, tmp_path / "topo_raster.pcsf")
        restored = read_pcsf(path)

        assert isinstance(restored.topography, TopographyRaster)
        np.testing.assert_allclose(restored.topography.x, x)
        np.testing.assert_allclose(restored.topography.y, y)
        np.testing.assert_allclose(restored.topography.elevation, elevation)


class TestMap3DUploadParserConsolidation:
    """map3d.py's own topo-upload parser must delegate to
    pycsamt.map.topo, not keep an independent implementation — the
    actual Phase 6 DoD."""

    def test_delegates_and_recognises_station_names(self):
        dash = pytest.importorskip("dash")
        from pycsamt.app.web.callbacks.map3d import _parse_topo_upload

        csv = "station_names,elevation\nA0,500\nA1,650\n"
        b64 = "data:text/csv;base64," + base64.b64encode(csv.encode()).decode()
        records = _parse_topo_upload(b64, "topo.csv")
        assert records == [
            {"station": "A0", "elev": 500.0},
            {"station": "A1", "elev": 650.0},
        ]

    def test_matches_pycsamt_map_topo_exactly(self):
        pytest.importorskip("dash")
        from pycsamt.app.web.callbacks.map3d import _parse_topo_upload
        from pycsamt.map.topo import parse_elevation_file

        csv = "ID,Z\nA0,12.5\nA1,13.0\n"
        b64 = base64.b64encode(csv.encode()).decode()
        expected = parse_elevation_file(b64, "x.csv")
        actual = {r["station"]: r["elev"] for r in _parse_topo_upload(b64, "x.csv")}
        assert actual == expected
