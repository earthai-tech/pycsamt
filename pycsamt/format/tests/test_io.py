# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Round-trip tests for pycsamt.format.io — write_pcsf / read_pcsf."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from pycsamt.format.io import read_pcsf, write_pcsf
from pycsamt.format.schema import (
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSF_VERSION,
    PCSFModel,
    RESISTIVITY_UNIT,
    StationTable,
    TopographyPerStation,
    UnstructuredMeshGeometry,
)


def _grid2d_model() -> PCSFModel:
    x = np.array([0.0, 100.0, 200.0, 300.0])
    z = np.array([10.0, 50.0, 150.0])
    rho_log10 = np.array(
        [
            [2.0, 2.1, 2.2, 2.3],
            [2.5, 2.6, 2.7, 2.8],
            [3.0, 3.1, 3.2, 3.3],
        ]
    )
    geometry = Grid2DGeometry(
        x=x, z=z, x_nodes=np.array([-50.0, 50.0, 150.0, 250.0, 350.0]),
        z_nodes=np.array([0.0, 30.0, 100.0, 200.0]),
        origin=np.array([500000.0, 4500000.0]), azimuth_deg=45.0,
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=10.0**rho_log10,
        resistivity_native=rho_log10,
        resistivity_native_encoding="log10",
        uncertainty=np.full_like(rho_log10, 0.1),
        stations=StationTable(
            name=["S0", "S1", "S2", "S3"], x=x, y=np.zeros(4), z=np.array([10.0, 12.0, 9.0, 11.0]),
            lon=np.array([7.50, 7.51, 7.52, 7.53]),
            lat=np.array([45.10, 45.11, 45.12, 45.13]),
        ),
        topography=TopographyPerStation(
            station_id=["S0", "S1", "S2", "S3"],
            elevation=np.array([10.0, 12.0, 9.0, 11.0]),
        ),
        survey={"name": "demo survey", "n_stations": 4, "bbox": [0.0, 0.0, 300.0, 0.0]},
        history={"rms": np.array([3.2, 2.1, 1.4, 1.05]), "lambda": np.array([10.0, 5.0, 2.0, 1.0])},
        source_backend="occam2d",
        created_by="pytest",
        crs="EPSG:32633",
        description="synthetic round-trip fixture",
        metadata={"note": "unit test"},
    )


def _grid3d_model() -> PCSFModel:
    geometry = Grid3DGeometry(
        x=np.array([0.0, 100.0, 200.0]),
        y=np.array([0.0, 100.0]),
        z=np.array([10.0, 50.0, 150.0, 400.0]),
        origin=np.array([500000.0, 4500000.0, 0.0]),
        rotation_deg=12.5,
        n_air=3,
    )
    rho = np.full(geometry.resistivity_shape, 100.0)
    return PCSFModel(geometry=geometry, resistivity=rho, source_backend="modem3d")


def _mesh_model() -> PCSFModel:
    geometry = UnstructuredMeshGeometry(
        nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
        region_ids=np.array([0, 1]),
        plane="xz",
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=np.array([120.0, 340.0]),
        resistivity_by_region=np.array([120.0, 340.0]),
        resistivity_by_node=np.array([100.0, 150.0, 200.0, 250.0]),
        source_backend="mare2dem",
    )


def _multiline_model() -> PCSFModel:
    line_geo = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
    lines = [
        LineEntry(
            line_id="L1", geometry=line_geo, resistivity=np.array([[100.0, 110.0], [50.0, 55.0]]),
            offset_y=0.0, offset_kind="real", azimuth_deg=0.0,
        ),
        LineEntry(
            line_id="L2", geometry=line_geo, resistivity=np.array([[200.0, 210.0], [90.0, 95.0]]),
            offset_y=250.0, offset_kind="real", azimuth_deg=0.0,
        ),
    ]
    derived = DerivedVolume(
        grid=Grid3DGeometry(x=np.array([0.0, 100.0]), y=np.array([0.0, 125.0, 250.0]), z=np.array([10.0, 50.0])),
        resistivity=np.full((2, 3, 2), 150.0),
        derivation_method="linear_interp",
        derived_from=["L1", "L2"],
        synthesized=True,
    )
    geometry = MultilineGeometry(lines=lines, derived_volume=derived)
    return PCSFModel(geometry=geometry, source_backend="generic")


class TestGrid2DRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsf(model, tmp_path / "grid2d.pcsf")
        restored = read_pcsf(path)

        assert restored.kind == "grid2d"
        assert restored.source_backend == "occam2d"
        np.testing.assert_allclose(restored.geometry.x, model.geometry.x)
        np.testing.assert_allclose(restored.geometry.z, model.geometry.z)
        np.testing.assert_allclose(restored.geometry.x_nodes, model.geometry.x_nodes)
        np.testing.assert_allclose(restored.geometry.origin, model.geometry.origin)
        assert restored.geometry.azimuth_deg == pytest.approx(45.0)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        np.testing.assert_allclose(restored.resistivity_native, model.resistivity_native)
        assert restored.resistivity_native_encoding == "log10"
        np.testing.assert_allclose(restored.uncertainty, model.uncertainty)
        assert restored.stations.name == model.stations.name
        np.testing.assert_allclose(restored.stations.z, model.stations.z)
        np.testing.assert_allclose(restored.stations.lon, model.stations.lon)
        np.testing.assert_allclose(restored.stations.lat, model.stations.lat)
        assert restored.topography.station_id == model.topography.station_id
        np.testing.assert_allclose(restored.topography.elevation, model.topography.elevation)
        assert restored.survey == model.survey
        np.testing.assert_allclose(restored.history["rms"], model.history["rms"])
        assert restored.crs == "EPSG:32633"
        assert restored.description == "synthetic round-trip fixture"
        assert restored.metadata == {"note": "unit test"}

    def test_root_attrs(self, tmp_path):
        path = write_pcsf(_grid2d_model(), tmp_path / "grid2d.pcsf")
        with h5py.File(path, "r") as fh:
            assert fh.attrs["pcsf_version"] == PCSF_VERSION
            assert fh.attrs["resistivity_unit"] == RESISTIVITY_UNIT
            assert fh.attrs["source_backend"] == "occam2d"

    def test_lon_lat_datasets_present_when_set(self, tmp_path):
        path = write_pcsf(_grid2d_model(), tmp_path / "grid2d.pcsf")
        with h5py.File(path, "r") as fh:
            assert "lon" in fh["stations"]
            assert "lat" in fh["stations"]

    def test_stations_without_lon_lat_round_trip_as_none(self, tmp_path):
        model = _grid2d_model()
        model.stations = StationTable(
            name=model.stations.name,
            x=model.stations.x,
            y=model.stations.y,
            z=model.stations.z,
        )
        path = write_pcsf(model, tmp_path / "grid2d_nolonlat.pcsf")
        with h5py.File(path, "r") as fh:
            assert "lon" not in fh["stations"]
            assert "lat" not in fh["stations"]
        restored = read_pcsf(path)
        assert restored.stations.lon is None
        assert restored.stations.lat is None


class TestGrid3DRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _grid3d_model()
        path = write_pcsf(model, tmp_path / "grid3d.pcsf")
        restored = read_pcsf(path)

        assert restored.kind == "grid3d"
        assert restored.geometry.n_air == 3
        assert restored.geometry.rotation_deg == pytest.approx(12.5)
        np.testing.assert_allclose(restored.geometry.origin, model.geometry.origin)
        # (n_z, n_y, n_x) axis order preserved end to end.
        assert restored.resistivity.shape == (4, 2, 3)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)


class TestMeshUnstructuredRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _mesh_model()
        path = write_pcsf(model, tmp_path / "mesh.pcsf")
        restored = read_pcsf(path)

        assert restored.kind == "mesh_unstructured"
        np.testing.assert_allclose(restored.geometry.nodes, model.geometry.nodes)
        np.testing.assert_array_equal(restored.geometry.connectivity, model.geometry.connectivity)
        np.testing.assert_array_equal(restored.geometry.region_ids, model.geometry.region_ids)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        np.testing.assert_allclose(restored.resistivity_by_region, model.resistivity_by_region)
        np.testing.assert_allclose(restored.resistivity_by_node, model.resistivity_by_node)


class TestMultilineRoundTrip:
    def test_round_trip_preserves_line_order_and_derived_volume(self, tmp_path):
        model = _multiline_model()
        path = write_pcsf(model, tmp_path / "multiline.pcsf")
        restored = read_pcsf(path)

        assert restored.kind == "multiline"
        assert [line.line_id for line in restored.geometry.lines] == ["L1", "L2"]
        np.testing.assert_allclose(
            restored.geometry.lines[1].resistivity, model.geometry.lines[1].resistivity
        )
        assert restored.geometry.lines[1].offset_y == pytest.approx(250.0)
        assert restored.geometry.lines[1].offset_kind == "real"
        assert restored.geometry.derived_volume is not None
        assert restored.geometry.derived_volume.synthesized is True
        assert restored.geometry.derived_volume.derived_from == ["L1", "L2"]
        np.testing.assert_allclose(
            restored.geometry.derived_volume.resistivity, model.geometry.derived_volume.resistivity
        )
        assert restored.resistivity is None


class TestErrors:
    def test_write_rejects_non_pcsf_model(self, tmp_path):
        with pytest.raises(TypeError):
            write_pcsf(object(), tmp_path / "bad.pcsf")

    def test_write_rejects_invalid_model(self, tmp_path):
        model = PCSFModel(geometry=Grid2DGeometry(x=np.array([0.0]), z=np.array([0.0])))
        with pytest.raises(ValueError, match="resistivity is required"):
            write_pcsf(model, tmp_path / "invalid.pcsf")

    def test_read_rejects_missing_geometry_group(self, tmp_path):
        path = tmp_path / "no_geometry.pcsf"
        with h5py.File(path, "w") as fh:
            fh.attrs["pcsf_version"] = PCSF_VERSION
        with pytest.raises(ValueError, match="geometry"):
            read_pcsf(path)

    def test_read_rejects_unknown_geometry_kind(self, tmp_path):
        path = tmp_path / "bad_kind.pcsf"
        with h5py.File(path, "w") as fh:
            fh.attrs["pcsf_version"] = PCSF_VERSION
            fh.create_group("geometry").attrs["kind"] = "not_a_real_kind"
        with pytest.raises(ValueError, match="unknown geometry kind"):
            read_pcsf(path)

    def test_write_creates_missing_parent_directories(self, tmp_path):
        model = _grid2d_model()
        nested = tmp_path / "a" / "b" / "c.pcsf"
        path = write_pcsf(model, nested)
        assert path.exists()

    def test_read_rejects_missing_pcsf_version(self, tmp_path):
        path = tmp_path / "no_version.pcsf"
        with h5py.File(path, "w") as fh:
            fh.create_group("geometry").attrs["kind"] = "grid2d"
        with pytest.raises(ValueError, match="pcsf_version"):
            read_pcsf(path)

    def test_read_rejects_malformed_pcsf_version(self, tmp_path):
        path = tmp_path / "bad_version.pcsf"
        with h5py.File(path, "w") as fh:
            fh.attrs["pcsf_version"] = "not-a-version"
            fh.create_group("geometry").attrs["kind"] = "grid2d"
        with pytest.raises(ValueError, match="malformed pcsf_version"):
            read_pcsf(path)

    def test_read_rejects_unrecognised_major_version(self, tmp_path):
        path = write_pcsf(_grid2d_model(), tmp_path / "future_major.pcsf")
        with h5py.File(path, "a") as fh:
            fh.attrs["pcsf_version"] = "99.0.0"
        with pytest.raises(ValueError, match="unsupported pcsf_version"):
            read_pcsf(path)

    def test_read_warns_on_newer_minor_version(self, tmp_path):
        current_major, current_minor, _ = (
            int(part) for part in PCSF_VERSION.split(".")
        )
        path = write_pcsf(_grid2d_model(), tmp_path / "future_minor.pcsf")
        with h5py.File(path, "a") as fh:
            fh.attrs["pcsf_version"] = f"{current_major}.{current_minor + 1}.0"
        with pytest.warns(UserWarning, match="newer than this reader's"):
            restored = read_pcsf(path)
        assert restored.kind == "grid2d"

    def test_read_accepts_older_minor_version(self, tmp_path):
        path = write_pcsf(_grid2d_model(), tmp_path / "older_minor.pcsf")
        with h5py.File(path, "a") as fh:
            fh.attrs["pcsf_version"] = "0.0.0"
        # No warning, no error: an older MINOR is fully understood by
        # a newer, backward-compatible reader.
        restored = read_pcsf(path)
        assert restored.kind == "grid2d"
