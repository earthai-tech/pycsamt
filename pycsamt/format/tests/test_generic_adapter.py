# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.format.adapters.generic` — the array-only
entry point any AI/DL inversion result (UNet, GCN, ResNet, or a third
party's own model) uses to become a PCSF/PCSM file, with no dependency
on any pycsamt solver result class."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.format import (
    read_pcsf,
    read_pcsm,
    write_pcsf,
    write_pcsm,
)
from pycsamt.format.adapters.generic import (
    grid2d_to_pcsf,
    grid3d_to_pcsf,
    mesh_to_pcsf,
)
from pycsamt.format.provenance import ModelProvenance
from pycsamt.format.schema import StationTable


class TestGrid2DToPcsf:
    def test_linear_encoding_passthrough(self):
        rho = np.array([[100.0, 110.0], [50.0, 55.0]])
        model = grid2d_to_pcsf(rho, x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
        model.validate()
        assert model.kind == "grid2d"
        np.testing.assert_allclose(model.resistivity, rho)
        assert model.resistivity_native is None
        assert model.resistivity_native_encoding is None
        assert model.source_backend == "ai"

    def test_log10_encoding_is_converted_and_native_is_preserved(self):
        log10_rho = np.array([[2.0, 2.1], [1.7, 1.8]])
        model = grid2d_to_pcsf(
            log10_rho, x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]),
            encoding="log10", source_backend="unet",
        )
        model.validate()
        np.testing.assert_allclose(model.resistivity, 10.0**log10_rho)
        np.testing.assert_allclose(model.resistivity_native, log10_rho)
        assert model.resistivity_native_encoding == "log10"
        assert model.source_backend == "unet"

    def test_unknown_encoding_raises(self):
        with pytest.raises(ValueError):
            grid2d_to_pcsf(
                np.zeros((2, 2)), x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0]),
                encoding="bogus",
            )

    def test_origin_and_azimuth_are_the_offset_mechanism(self):
        model = grid2d_to_pcsf(
            np.zeros((2, 2)), x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]),
            origin=(500.0, 1000.0), azimuth_deg=30.0,
        )
        np.testing.assert_allclose(model.geometry.origin, [500.0, 1000.0])
        assert model.geometry.azimuth_deg == 30.0

    def test_station_elevations_and_lonlat_populate_stations_and_topography(self):
        model = grid2d_to_pcsf(
            np.zeros((1, 2)), x=np.array([0.0, 100.0]), z=np.array([25.0]),
            station_names=["s1", "s2"], station_x=[0.0, 100.0],
            station_elevations={"s1": 10.0, "s2": 12.0},
            station_lonlat={"s1": (10.0, 45.0), "s2": (10.1, 45.1)},
        )
        model.validate()
        assert model.stations.name == ["s1", "s2"]
        np.testing.assert_allclose(model.stations.lon, [10.0, 10.1])
        np.testing.assert_allclose(model.stations.lat, [45.0, 45.1])
        assert model.topography is not None
        np.testing.assert_allclose(model.topography.elevation, [10.0, 12.0])

    def test_prebuilt_stations_bypass_the_builder(self):
        stations = StationTable(
            name=["a"], x=np.array([0.0]), y=np.array([0.0]), z=np.array([5.0])
        )
        model = grid2d_to_pcsf(
            np.zeros((1, 1)), x=np.array([0.0]), z=np.array([5.0]),
            stations=stations, station_names=["ignored"],
        )
        assert model.stations is stations

    def test_uncertainty_sensitivity_history_and_provenance_round_trip(self, tmp_path):
        rho = np.full((2, 2), 100.0)
        model = grid2d_to_pcsf(
            rho, x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]),
            uncertainty=np.full((2, 2), 0.1),
            sensitivity=np.full((2, 2), 0.5),
            history={"epoch": [0, 1, 2], "val_loss": [0.5, 0.3, 0.2]},
            provenance=ModelProvenance(
                architecture="UNet", framework="pytorch", random_seed=42,
            ),
        )
        model.validate()
        path = write_pcsf(model, tmp_path / "unet_result.pcsf")
        round_tripped = read_pcsf(path)
        np.testing.assert_allclose(round_tripped.resistivity, rho)
        np.testing.assert_allclose(round_tripped.uncertainty, 0.1)
        np.testing.assert_allclose(round_tripped.history["val_loss"], [0.5, 0.3, 0.2])
        assert round_tripped.metadata["model_provenance"]["architecture"] == "UNet"
        assert round_tripped.metadata["model_provenance"]["random_seed"] == 42

    def test_pcsm_round_trip(self, tmp_path):
        model = grid2d_to_pcsf(
            np.array([[100.0, 110.0], [50.0, 55.0]]),
            x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]),
            source_backend="unet",
        )
        pcsm_path = write_pcsm(model, tmp_path / "unet_result.pcsm")
        round_tripped = read_pcsm(pcsm_path)
        np.testing.assert_allclose(round_tripped.resistivity, model.resistivity)
        assert round_tripped.source_backend == "unet"


class TestGrid3DToPcsf:
    def test_shape_and_kind(self):
        rho = np.full((2, 2, 2), 100.0)
        model = grid3d_to_pcsf(
            rho, x=np.array([0.0, 100.0]), y=np.array([0.0, 100.0]),
            z=np.array([10.0, 50.0]), source_backend="resnet",
        )
        model.validate()
        assert model.kind == "grid3d"
        assert model.resistivity.shape == (2, 2, 2)
        assert model.source_backend == "resnet"

    def test_origin_rotation_and_n_air_passthrough(self):
        model = grid3d_to_pcsf(
            np.zeros((1, 1, 1)), x=np.array([0.0]), y=np.array([0.0]),
            z=np.array([10.0]), origin=(100.0, 200.0, 0.0),
            rotation_deg=15.0, n_air=3,
        )
        np.testing.assert_allclose(model.geometry.origin, [100.0, 200.0, 0.0])
        assert model.geometry.rotation_deg == 15.0
        assert model.geometry.n_air == 3

    def test_station_xyz_and_topo_wiring(self):
        model = grid3d_to_pcsf(
            np.zeros((1, 1, 1)), x=np.array([0.0]), y=np.array([0.0]),
            z=np.array([10.0]),
            station_names=["s1"], station_x=[10.0], station_y=[20.0],
            station_z=[30.0], station_lonlat={"s1": (11.0, 46.0)},
        )
        assert model.stations.z[0] == 30.0
        np.testing.assert_allclose(model.stations.lon, [11.0])


class TestMeshToPcsf:
    _NODES = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    _CONNECTIVITY = np.array([[0, 1, 2], [1, 3, 2]])

    def test_requires_at_least_one_resistivity_source(self):
        with pytest.raises(ValueError):
            mesh_to_pcsf(self._NODES, self._CONNECTIVITY)

    def test_per_cell_resistivity_used_as_is(self):
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY, resistivity=np.array([100.0, 200.0]),
        )
        model.validate()
        assert model.kind == "mesh_unstructured"
        np.testing.assert_allclose(model.resistivity, [100.0, 200.0])
        assert model.resistivity_by_node is None

    def test_gcn_node_values_are_averaged_per_triangle(self):
        # Triangle 0 = nodes (0,1,2) -> mean(1,2,3); triangle 1 = nodes
        # (1,3,2) -> mean(2,4,3) -- a hand-computed reference for the
        # documented arithmetic-mean projection.
        node_values = np.array([1.0, 2.0, 3.0, 4.0])
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY, resistivity_by_node=node_values,
            source_backend="gcn",
        )
        model.validate()
        np.testing.assert_allclose(model.resistivity, [2.0, 3.0])
        np.testing.assert_allclose(model.resistivity_by_node, node_values)
        assert model.source_backend == "gcn"

    def test_log10_encoding_applies_to_the_primary_source(self):
        node_values_log10 = np.array([1.0, 2.0, 3.0, 4.0])
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY,
            resistivity_by_node=node_values_log10, encoding="log10",
        )
        expected_node_linear = 10.0**node_values_log10
        np.testing.assert_allclose(model.resistivity_by_node, expected_node_linear)
        np.testing.assert_allclose(
            model.resistivity,
            expected_node_linear[self._CONNECTIVITY].mean(axis=1),
        )

    def test_region_table_expansion_is_zero_based(self):
        region_ids = np.array([0, 1])
        by_region = np.array([100.0, 200.0])
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY, region_ids=region_ids,
            resistivity_by_region=by_region,
        )
        model.validate()
        np.testing.assert_allclose(model.resistivity, [100.0, 200.0])
        assert model.geometry.region_ids.tolist() == [0, 1]

    def test_out_of_range_region_id_raises(self):
        with pytest.raises(ValueError):
            mesh_to_pcsf(
                self._NODES, self._CONNECTIVITY, region_ids=np.array([0, 5]),
                resistivity_by_region=np.array([100.0, 200.0]),
            )

    def test_missing_region_ids_default_to_single_region(self):
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY, resistivity=np.array([100.0, 200.0]),
        )
        assert model.geometry.region_ids.tolist() == [0, 0]

    def test_resistivity_by_node_round_trips_through_pcsf(self, tmp_path):
        node_values = np.array([1.0, 2.0, 3.0, 4.0])
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY, resistivity_by_node=node_values,
        )
        path = write_pcsf(model, tmp_path / "gcn_result.pcsf")
        round_tripped = read_pcsf(path)
        np.testing.assert_allclose(round_tripped.resistivity_by_node, node_values)
        np.testing.assert_allclose(round_tripped.resistivity, model.resistivity)

    def test_resistivity_by_node_round_trips_through_pcsm(self, tmp_path):
        node_values = np.array([1.0, 2.0, 3.0, 4.0])
        model = mesh_to_pcsf(
            self._NODES, self._CONNECTIVITY, resistivity_by_node=node_values,
        )
        path = write_pcsm(model, tmp_path / "gcn_result.pcsm")
        round_tripped = read_pcsm(path)
        np.testing.assert_allclose(round_tripped.resistivity_by_node, node_values)
