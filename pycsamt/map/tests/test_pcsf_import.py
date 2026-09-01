# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 7 tests — pycsamt.map.inversion.load_pcsf_lines / MapView.from_pcsf.

Validated against the real bundled Occam2D dataset (data/occam2D, 47
stations) so the ``rho`` column extracted per station is genuinely the
nearest-column value from a real inverted grid, not synthetic data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.map._core import StationRecord
from pycsamt.map.inversion import load_pcsf_lines
from pycsamt.map.view import MapView

_OCCAM_DIR = Path(__file__).parents[3] / "data" / "occam2D"
_SKIP_OCCAM = pytest.mark.skipif(
    not _OCCAM_DIR.exists(), reason=f"bundled occam2D data not found: {_OCCAM_DIR}"
)


def _write_grid2d_pcsf(path, *, station_elevations=None, station_lonlat=None):
    from pycsamt.format import write_pcsf
    from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
    from pycsamt.models.occam2d.results import InversionResult

    result = InversionResult(workdir=_OCCAM_DIR)
    model = occam2d_to_pcsf(
        result,
        station_elevations=station_elevations,
        station_lonlat=station_lonlat,
    )
    return write_pcsf(model, path), result


def _write_multiline_pcsf(path):
    from pycsamt.format import write_pcsf
    from pycsamt.format.multiline import build_multiline_pcsf

    profiles = {
        "L1": {
            "x": np.array([0.0, 100.0, 200.0]),
            "z": np.array([10.0, 50.0]),
            "rho": np.array([[100.0, 110.0, 120.0], [50.0, 55.0, 60.0]]),
            "sta_x": [0.0, 100.0, 200.0],
            "sta_names": ["A0", "A1", "A2"],
            "sta_elev": [10.0, 12.0, 9.0],
        },
        "L2": {
            "x": np.array([0.0, 100.0, 200.0]),
            "z": np.array([10.0, 50.0]),
            "rho": np.array([[200.0, 210.0, 220.0], [90.0, 95.0, 99.0]]),
            "sta_x": [0.0, 100.0, 200.0],
            "sta_names": ["B0", "B1", "B2"],
            "sta_elev": [11.0, 13.0, 8.0],
        },
    }
    model = build_multiline_pcsf(profiles)
    return write_pcsf(model, path)


@_SKIP_OCCAM
class TestLoadPcsfLinesGrid2D:
    def test_one_line_all_real_stations(self, tmp_path):
        path, result = _write_grid2d_pcsf(tmp_path / "grid2d.pcsf")
        data = load_pcsf_lines(path, fetch_elevation=False)
        assert data.lines == ("line1",)
        assert len(data.stations) == result.data.n_sites

    def test_rho_columns_are_real_nearest_mesh_columns(self, tmp_path):
        from pycsamt.format import read_pcsf

        path, result = _write_grid2d_pcsf(tmp_path / "grid2d.pcsf")
        model = read_pcsf(path)
        data = load_pcsf_lines(path, fetch_elevation=False)
        section = data.metadata["sections"]["line1"]

        assert section["rho"].shape == (
            model.geometry.z.shape[0],
            len(model.stations.name),
        )
        # Spot-check the first station's column really is its nearest
        # column in the real inverted grid.
        first_name = str(model.stations.name[0])
        col_idx = int(
            np.argmin(np.abs(model.geometry.x - model.stations.x[0]))
        )
        expected = model.resistivity[:, col_idx]
        actual_idx = list(section["stations"]).index(first_name)
        np.testing.assert_array_equal(
            section["rho"][:, actual_idx], expected
        )

    def test_elevations_flow_through(self, tmp_path):
        elev = {"S00": 1200.0, "S01": 1195.0}
        path, _ = _write_grid2d_pcsf(
            tmp_path / "grid2d.pcsf", station_elevations=elev
        )
        data = load_pcsf_lines(path, fetch_elevation=False)
        by_id = {s.id: s.elevation for s in data.stations}
        assert by_id["S00"] == pytest.approx(1200.0)
        assert by_id["S01"] == pytest.approx(1195.0)
        assert by_id["S02"] is None  # not fabricated

    def test_known_stations_supply_geo_coordinates(self, tmp_path):
        path, _ = _write_grid2d_pcsf(tmp_path / "grid2d.pcsf")
        known = (
            StationRecord("S00", latitude=10.0, longitude=20.0, elevation=500.0),
        )
        data = load_pcsf_lines(
            path, known_stations=known, fetch_elevation=False
        )
        s00 = next(s for s in data.stations if s.id == "S00")
        assert s00.latitude == pytest.approx(10.0)
        assert s00.longitude == pytest.approx(20.0)

    def test_map3d_renders_from_a_real_pcsf_file(self, tmp_path):
        path, _ = _write_grid2d_pcsf(tmp_path / "grid2d.pcsf")
        mv = MapView.from_pcsf(path, fetch_elevation=False)
        fig = mv.map3d(mode="fence")
        assert fig is not None

    def test_file_own_lon_lat_used_without_known_stations(self, tmp_path):
        # A single Occam2D line, self-georeferenced via
        # station_lonlat -- no EDI known_stations match needed at all
        # to place it on a real basemap.
        path, result = _write_grid2d_pcsf(
            tmp_path / "grid2d.pcsf",
            station_lonlat={"S00": (7.5, 45.1), "S01": (7.6, 45.2)},
        )
        data = load_pcsf_lines(path, fetch_elevation=False)
        by_id = {s.id: s for s in data.stations}
        assert by_id["S00"].longitude == pytest.approx(7.5)
        assert by_id["S00"].latitude == pytest.approx(45.1)
        # Unmatched stations stay unknown, not fabricated.
        assert by_id["S02"].longitude is None

    def test_known_stations_take_priority_over_file_own_lon_lat(self, tmp_path):
        path, _ = _write_grid2d_pcsf(
            tmp_path / "grid2d.pcsf",
            station_lonlat={"S00": (7.5, 45.1)},
        )
        known = (
            StationRecord("S00", latitude=99.0, longitude=88.0, elevation=None),
        )
        data = load_pcsf_lines(
            path, known_stations=known, fetch_elevation=False
        )
        s00 = next(s for s in data.stations if s.id == "S00")
        assert s00.longitude == pytest.approx(88.0)
        assert s00.latitude == pytest.approx(99.0)


class TestLoadPcsfLinesMultiline:
    def test_two_lines_grouped_by_line_id(self, tmp_path):
        path = _write_multiline_pcsf(tmp_path / "multiline.pcsf")
        data = load_pcsf_lines(path, fetch_elevation=False)
        assert set(data.lines) == {"L1", "L2"}
        assert len(data.stations) == 6
        sections = data.metadata["sections"]
        assert sections["L1"]["rho"].shape == (2, 3)
        assert sections["L2"]["rho"].shape == (2, 3)

    def test_mapview_from_pcsf_end_to_end(self, tmp_path):
        path = _write_multiline_pcsf(tmp_path / "multiline.pcsf")
        mv = MapView.from_pcsf(path, fetch_elevation=False)
        assert mv.n_stations == 6
        assert set(mv.lines) == {"L1", "L2"}


def _write_grid3d_pcsf(path, *, n_air=2, station_line_id=None, station_lonlat=None):
    """A synthetic grid3d PCSF with a per-cell-unique resistivity value
    (``1000*iz + 10*iy + ix``) so a curtain's sampled values can be
    checked against hand-computed expectations, not just shapes --
    and station coordinates in a *different*, zero-centred local frame
    (mirroring ModEM's own ``.dat`` convention) from the grid's own
    node-based frame, so the registration offset is genuinely
    exercised rather than accidentally trivial (offset == 0)."""
    from pycsamt.format import PCSFModel, StationTable, write_pcsf
    from pycsamt.format.schema import Grid3DGeometry

    nz_earth, ny, nx = 5, 6, 7
    x_widths = np.full(nx, 100.0)
    y_widths = np.full(ny, 100.0)
    z_widths = np.concatenate([np.full(n_air, 50.0), np.full(nz_earth, 80.0)])

    x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])
    y_nodes = np.concatenate([[0.0], np.cumsum(y_widths)])
    z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])
    x_c = (x_nodes[:-1] + x_nodes[1:]) / 2
    y_c = (y_nodes[:-1] + y_nodes[1:]) / 2
    z_c = (z_nodes[:-1] + z_nodes[1:]) / 2

    nz = n_air + nz_earth
    rho = np.fromfunction(
        lambda iz, iy, ix: 1000 * iz + 10 * iy + ix, (nz, ny, nx)
    ).astype(float)

    geometry = Grid3DGeometry(
        x=x_c,
        y=y_c,
        z=z_c,
        x_nodes=x_nodes,
        y_nodes=y_nodes,
        z_nodes=z_nodes,
        n_air=n_air,
    )
    names = ["23-01-001", "23-01-002", "23-01-003", "23-01-004"]
    # Grid centre sits at (350, 300); station-local frame is centred at
    # (0, 0) -- exercises the symmetric-padding registration offset.
    lon = lat = None
    if station_lonlat is not None:
        lon = np.array([station_lonlat.get(n, (np.nan, np.nan))[0] for n in names])
        lat = np.array([station_lonlat.get(n, (np.nan, np.nan))[1] for n in names])
    stations = StationTable(
        name=names,
        x=np.array([-300.0, -100.0, 100.0, 300.0]),
        y=np.zeros(4),
        z=np.array([10.0, 11.0, 12.0, 13.0]),
        line_id=station_line_id,
        lon=lon,
        lat=lat,
    )
    model = PCSFModel(
        geometry=geometry,
        resistivity=rho,
        source_backend="modem3d",
        stations=stations,
    )
    return write_pcsf(model, path), (x_c, y_c, z_c, n_air, rho)


class TestLoadPcsfLinesGrid3D:
    def test_curtain_samples_real_cells_with_registration_offset(self, tmp_path):
        path, (x_c, y_c, z_c, n_air, rho) = _write_grid3d_pcsf(
            tmp_path / "grid3d.pcsf"
        )
        data = load_pcsf_lines(path, fetch_elevation=False)
        assert data.lines == ("01",)
        assert len(data.stations) == 4

        section = data.metadata["sections"]["01"]
        # The earth depth axis is re-zeroed to the top of the earth
        # domain (0 at the first earth cell's top edge), matching
        # pycsamt.models.modem.section.station_curtain -- an n_air > 0
        # model like this one carries topography in its own air layers,
        # so no per-station datum shift is applied. 80 m earth cells ->
        # centres 40, 120, 200, 280, 360.
        np.testing.assert_allclose(
            section["z"], [40.0, 120.0, 200.0, 280.0, 360.0]
        )

        # Hand-computed expectation: real x = local_x + 350 (grid centre),
        # real y = local_y + 300 -> nearest (iy, ix) per station, then the
        # earth-only column at that cell (unchanged: n_air > 0 -> no
        # datum re-reference).
        expected_ix = [0, 2, 4, 6]  # real x = 50, 250, 450, 650
        expected_iy = 2  # real y = 300 -> nearest y_c is 250 (tie -> lower)
        for k, ix in enumerate(expected_ix):
            expected_col = rho[n_air:, expected_iy, ix]
            np.testing.assert_array_equal(section["rho"][:, k], expected_col)

    def test_residual_air_fill_is_masked_for_pre_mask_files(self, tmp_path):
        """A grid3d PCSF written before the adapter learned to mask air
        fill (n_air=0, top layers still ~1e12 ohm.m) must still render
        with the air dropped -- load_pcsf_lines nan-masks cells above
        _AIR_FILL_THRESHOLD defensively."""
        from pycsamt.format import PCSFModel, StationTable, write_pcsf
        from pycsamt.format.schema import Grid3DGeometry

        nz, ny, nx = 6, 4, 5
        rho = np.full((nz, ny, nx), 200.0)
        rho[:2, :, :] = 1e12  # two un-flagged air layers
        z_c = np.arange(nz) * 50.0 + 25.0
        geometry = Grid3DGeometry(
            x=np.arange(nx) * 100.0,
            y=np.arange(ny) * 100.0,
            z=z_c,
            n_air=0,
        )
        stations = StationTable(
            name=["23-07-001", "23-07-002"],
            x=np.array([100.0, 300.0]),
            y=np.array([150.0, 150.0]),
            z=np.array([0.0, 0.0]),
        )
        path = write_pcsf(
            PCSFModel(
                geometry=geometry,
                resistivity=rho,
                source_backend="modem3d",
                stations=stations,
            ),
            tmp_path / "grid3d_premask.pcsf",
        )
        data = load_pcsf_lines(path, fetch_elevation=False)
        section = data.metadata["sections"]["07"]
        col = section["rho"][:, 0]
        assert np.isnan(col[:2]).all()  # air rows dropped
        np.testing.assert_allclose(col[2:], 200.0)  # real earth kept

    def test_depth_axis_is_referenced_to_each_station_surface(self, tmp_path):
        """n_air == 0 + a recorded a.s.l. datum -> every column is
        re-referenced to its own ground surface, so ``section['z']`` is a
        true depth-below-surface (not depth-below-model-top)."""
        from pycsamt.format import PCSFModel, StationTable, write_pcsf
        from pycsamt.format.schema import Grid3DGeometry

        # 10 uniform 40 m earth cells; a linear rho ramp with depth so the
        # re-reference is easy to read back.
        nz = 10
        z_nodes = np.arange(nz + 1) * 40.0
        z_c = (z_nodes[:-1] + z_nodes[1:]) / 2.0
        rho = np.tile(z_c[:, None, None], (1, 3, 3))  # rho == depth-below-top
        geometry = Grid3DGeometry(
            x=np.arange(3) * 100.0,
            y=np.arange(3) * 100.0,
            z=z_c,
            z_nodes=z_nodes,
            n_air=0,
        )
        stations = StationTable(
            name=["23-05-001", "23-05-002"],
            x=np.array([100.0, 100.0]),
            y=np.array([100.0, 100.0]),
            z=np.array([90.0, 60.0]),  # a.s.l. -> 110 / 140 m below top
        )
        model = PCSFModel(
            geometry=geometry,
            resistivity=rho,
            source_backend="modem3d",
            stations=stations,
            metadata={"station_z": {"datum_masl": 200.0}},
        )
        path = write_pcsf(model, tmp_path / "grid3d_topo.pcsf")
        section = load_pcsf_lines(path, fetch_elevation=False).metadata[
            "sections"
        ]["05"]

        # rho == depth-below-model-top, so a column read at
        # "depth-below-surface d" must come back as d + (200 - elev).
        z = section["z"]
        s1 = np.interp(120.0, z, section["rho"][:, 0])  # elev 90 -> +110
        s2 = np.interp(120.0, z, section["rho"][:, 1])  # elev 60 -> +140
        assert s1 == pytest.approx(230.0, abs=5.0)
        assert s2 == pytest.approx(260.0, abs=5.0)

    def test_deep_bc_padding_is_trimmed(self, tmp_path):
        """The deep, geometrically-growing boundary-condition cells are
        dropped from a sliced curtain so they don't bloat 3-D views."""
        from pycsamt.format import PCSFModel, StationTable, write_pcsf
        from pycsamt.format.schema import Grid3DGeometry

        # 12 x 40 m core cells, then 8 cells doubling each time (deep BC).
        core = np.full(12, 40.0)
        pad = 80.0 * 2.0 ** np.arange(8)
        widths = np.concatenate([core, pad])
        z_nodes = np.concatenate([[0.0], np.cumsum(widths)])
        z_c = (z_nodes[:-1] + z_nodes[1:]) / 2.0
        rho = np.full((z_c.size, 3, 3), 100.0)
        geometry = Grid3DGeometry(
            x=np.arange(3) * 100.0, y=np.arange(3) * 100.0,
            z=z_c, z_nodes=z_nodes, n_air=0,
        )
        stations = StationTable(
            name=["23-06-001", "23-06-002"],
            x=np.array([100.0, 100.0]), y=np.array([100.0, 100.0]),
            z=np.array([0.0, 0.0]),
        )
        path = write_pcsf(
            PCSFModel(geometry=geometry, resistivity=rho,
                      source_backend="modem3d", stations=stations),
            tmp_path / "grid3d_deep.pcsf",
        )
        section = load_pcsf_lines(path, fetch_elevation=False).metadata[
            "sections"
        ]["06"]
        # core zone ends at 480 m; the trim keeps the core plus a cell or
        # two, never the multi-km padding.
        assert section["z"].max() < 2000.0
        assert section["z"].size < z_c.size

    def test_multiline_grouping_falls_back_to_name_heuristic(self, tmp_path):
        path, _ = _write_grid3d_pcsf(tmp_path / "grid3d_2lines.pcsf")
        # Overwrite with two distinct line tokens via the station-name
        # convention (no explicit line_id in the file).
        from pycsamt.format import PCSFModel, StationTable, read_pcsf, write_pcsf

        model = read_pcsf(path)
        model.stations = StationTable(
            name=["23-01-001", "23-01-002", "23-02-001", "23-02-002"],
            x=model.stations.x,
            y=model.stations.y,
            z=model.stations.z,
        )
        write_pcsf(model, path)

        data = load_pcsf_lines(path, fetch_elevation=False)
        assert set(data.lines) == {"01", "02"}
        assert len(data.stations) == 4

    def test_elevations_and_known_stations_flow_through(self, tmp_path):
        path, _ = _write_grid3d_pcsf(tmp_path / "grid3d.pcsf")
        known = (
            StationRecord(
                "23-01-001", latitude=10.0, longitude=20.0, elevation=500.0
            ),
        )
        data = load_pcsf_lines(
            path, known_stations=known, fetch_elevation=False
        )
        by_id = {s.id: s for s in data.stations}
        assert by_id["23-01-001"].latitude == pytest.approx(10.0)
        assert by_id["23-01-001"].longitude == pytest.approx(20.0)
        # Station table's own z (a real ModEM-recorded value here, 11.0)
        # is used as-is for a station with no known-station match.
        assert by_id["23-01-002"].elevation == pytest.approx(11.0)

    def test_file_own_lon_lat_used_without_known_stations(self, tmp_path):
        # Same real-basemap capability as grid2d, exercised for a
        # native ModEM grid3d volume (e.g. modem3d_to_pcsf's own
        # GG_Lat/GG_Lon passthrough) -- no known_stations needed.
        path, _ = _write_grid3d_pcsf(
            tmp_path / "grid3d.pcsf",
            station_lonlat={"23-01-001": (119.1, 32.1), "23-01-002": (119.2, 32.2)},
        )
        data = load_pcsf_lines(path, fetch_elevation=False)
        by_id = {s.id: s for s in data.stations}
        assert by_id["23-01-001"].longitude == pytest.approx(119.1)
        assert by_id["23-01-001"].latitude == pytest.approx(32.1)
        assert by_id["23-01-003"].longitude is None

    def test_mapview_from_pcsf_end_to_end(self, tmp_path):
        path, _ = _write_grid3d_pcsf(tmp_path / "grid3d.pcsf")
        mv = MapView.from_pcsf(path, fetch_elevation=False)
        assert mv.n_stations == 4
        fig = mv.map3d(mode="fence")
        assert fig is not None


def _write_split_region_mesh_pcsf(path):
    """A hand-built rectangular mesh split into a left half
    (rho=100) and a right half (rho=500) by two triangulated quads
    each, so point-location sampling can be checked against a known
    answer, not just a plausible-looking shape."""
    from pycsamt.format import PCSFModel, StationTable, write_pcsf
    from pycsamt.format.schema import UnstructuredMeshGeometry

    nodes = np.array(
        [
            [0, 0], [100, 0], [200, 0],
            [0, 50], [100, 50], [200, 50],
            [0, 100], [100, 100], [200, 100],
        ],
        dtype=float,
    )
    tris = np.array(
        [
            [0, 1, 4], [0, 4, 3],
            [1, 2, 5], [1, 5, 4],
            [3, 4, 7], [3, 7, 6],
            [4, 5, 8], [4, 8, 7],
        ],
        dtype=np.int64,
    )
    region_ids = np.array([1, 1, 2, 2, 1, 1, 2, 2], dtype=np.int32)
    resistivity = np.where(region_ids == 1, 100.0, 500.0)

    geometry = UnstructuredMeshGeometry(
        nodes=nodes, connectivity=tris, region_ids=region_ids, plane="xz"
    )
    stations = StationTable(
        name=["S0", "S1"], x=np.array([50.0, 150.0]), y=np.zeros(2), z=np.zeros(2)
    )
    model = PCSFModel(
        geometry=geometry,
        resistivity=resistivity,
        source_backend="mare2dem",
        stations=stations,
    )
    return write_pcsf(model, path)


class TestLoadPcsfLinesMeshUnstructured:
    def test_curtain_samples_the_correct_region_per_station(self, tmp_path):
        path = _write_split_region_mesh_pcsf(tmp_path / "mesh.pcsf")
        data = load_pcsf_lines(path, fetch_elevation=False, mesh_z_samples=5)
        assert data.lines == ("line1",)
        assert len(data.stations) == 2

        section = data.metadata["sections"]["line1"]
        assert section["z"].shape == (5,)
        assert section["rho"].shape == (5, 2)
        stations = list(section["stations"])
        s0_col = section["rho"][:, stations.index("S0")]
        s1_col = section["rho"][:, stations.index("S1")]
        np.testing.assert_allclose(s0_col, 100.0)
        np.testing.assert_allclose(s1_col, 500.0)

    def test_query_outside_the_mesh_is_nan_not_fabricated(self, tmp_path):
        from pycsamt.format import PCSFModel, StationTable, read_pcsf, write_pcsf

        path = _write_split_region_mesh_pcsf(tmp_path / "mesh.pcsf")
        model = read_pcsf(path)
        # A station far outside the mesh's x range (mesh spans [0, 200]).
        model.stations = StationTable(
            name=["Sfar"], x=np.array([5000.0]), y=np.array([0.0]), z=np.array([0.0])
        )
        write_pcsf(model, path)
        data = load_pcsf_lines(path, fetch_elevation=False, mesh_z_samples=5)
        section = data.metadata["sections"]["line1"]
        assert np.all(np.isnan(section["rho"]))

    def test_mapview_from_pcsf_end_to_end(self, tmp_path):
        path = _write_split_region_mesh_pcsf(tmp_path / "mesh.pcsf")
        mv = MapView.from_pcsf(path, fetch_elevation=False, mesh_z_samples=5)
        assert mv.n_stations == 2
        fig = mv.map3d(mode="fence")
        assert fig is not None

    def test_compact_per_region_resistivity_is_rejected(self, tmp_path):
        # A model.validate()-legal file (resistivity in the compact,
        # region-collapsed form PCSFModel itself accepts for
        # mesh_unstructured) but not yet expanded onto each of this
        # mesh's 2 triangles -- load_pcsf_lines needs the per-triangle
        # form (every adapter in pycsamt.format.adapters already
        # expands it; this simulates a hand-built file that skipped
        # that step).
        from pycsamt.format import PCSFModel, StationTable, write_pcsf
        from pycsamt.format.schema import UnstructuredMeshGeometry

        geometry = UnstructuredMeshGeometry(
            nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64),
            region_ids=np.array([1, 1], dtype=np.int32),  # both triangles -> 1 region
        )
        stations = StationTable(
            name=["S0"], x=np.array([0.3]), y=np.array([0.0]), z=np.array([0.0])
        )
        model = PCSFModel(
            geometry=geometry,
            resistivity=np.array([100.0]),  # compact (n_regions=1,) form
            stations=stations,
        )
        path = write_pcsf(model, tmp_path / "mesh_bad.pcsf")
        with pytest.raises(ValueError, match="per-triangle"):
            load_pcsf_lines(path)


@pytest.mark.skipif(
    pytest.importorskip("triangle", reason="triangle package not installed") is None,
    reason="triangle package not installed",
)
class TestLoadPcsfLinesMeshUnstructuredRealData:
    """Same point-location sampling, exercised against the real bundled
    MARE2DEM ``demo_mt_inversion`` mesh (6540 real triangulated
    regions) — not just a hand-built synthetic split."""

    _DATA_DIR = Path(__file__).parents[3] / "data" / "mare2dem" / "demo_mt_inversion"

    def _real_model(self, tmp_path):
        import triangle

        from pycsamt.format import StationTable, write_pcsf
        from pycsamt.format.adapters.mare2dem import mare2dem_to_pcsf
        from pycsamt.forward.maxwell.contracts_tri import TriMesh
        from pycsamt.models.mare2dem.iotools.poly import read_poly
        from pycsamt.models.mare2dem.results import InversionResult

        poly_path = self._DATA_DIR / "demo.poly"
        if not poly_path.exists():
            pytest.skip(f"bundled MARE2DEM data not found: {poly_path}")

        pf = read_poly(poly_path)
        pslg = {"vertices": pf.nodes, "segments": pf.segments - 1, "regions": pf.regions}
        out = triangle.triangulate(pslg, "pA")
        region_ids = np.round(out["triangle_attributes"].ravel()).astype(np.int64)
        mesh = TriMesh(
            nodes_m=out["vertices"], triangles=out["triangles"], region_ids=region_ids
        )
        result = InversionResult(workdir=self._DATA_DIR)
        model = mare2dem_to_pcsf(result, mesh, created_by="pytest")

        x = model.geometry.nodes[:, 0]
        xs = np.linspace(
            x.min() + (x.max() - x.min()) * 0.3,
            x.min() + (x.max() - x.min()) * 0.7,
            4,
        )
        model.stations = StationTable(
            name=[f"S{i}" for i in range(4)], x=xs, y=np.zeros(4), z=np.zeros(4)
        )
        return write_pcsf(model, tmp_path / "mare2dem_real.pcsf")

    def test_curtain_stays_inside_the_real_mesh(self, tmp_path):
        path = self._real_model(tmp_path)
        data = load_pcsf_lines(path, fetch_elevation=False, mesh_z_samples=40)
        section = data.metadata["sections"]["line1"]
        assert section["rho"].shape == (40, 4)
        # Stations were placed well inside the mesh's real x extent, so
        # every query should land inside some real triangle.
        assert np.all(np.isfinite(section["rho"]))

    def test_mapview_end_to_end_on_real_mesh(self, tmp_path):
        path = self._real_model(tmp_path)
        mv = MapView.from_pcsf(path, fetch_elevation=False, mesh_z_samples=40)
        assert mv.n_stations == 4
        fig = mv.map3d(mode="fence")
        assert fig is not None


class TestLoadPcsfLinesErrors:
    def test_no_station_table_raises(self, tmp_path):
        from pycsamt.format import Grid2DGeometry, PCSFModel, write_pcsf

        geometry = Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0]))
        model = PCSFModel(geometry=geometry, resistivity=np.ones((2, 2)))
        path = write_pcsf(model, tmp_path / "nostations.pcsf")
        with pytest.raises(ValueError, match="no station table"):
            load_pcsf_lines(path)
