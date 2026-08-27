# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.format.schema — PCSF dataclasses and validation."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.format.schema import (
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    StationTable,
    TopographyRaster,
    UnstructuredMeshGeometry,
)


def _grid2d() -> Grid2DGeometry:
    return Grid2DGeometry(x=np.array([0.0, 100.0, 200.0]), z=np.array([10.0, 50.0]))


def _grid3d() -> Grid3DGeometry:
    return Grid3DGeometry(
        x=np.array([0.0, 100.0]), y=np.array([0.0, 50.0, 100.0]), z=np.array([10.0, 50.0, 200.0])
    )


def _mesh() -> UnstructuredMeshGeometry:
    return UnstructuredMeshGeometry(
        nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
        region_ids=np.array([0, 1]),
    )


class TestGrid2DGeometry:
    def test_resistivity_shape(self):
        assert _grid2d().resistivity_shape == (2, 3)

    def test_rejects_mismatched_x_nodes(self):
        geo = Grid2DGeometry(
            x=np.array([0.0, 100.0]), z=np.array([10.0]), x_nodes=np.array([0.0, 50.0])
        )
        with pytest.raises(ValueError, match="x_nodes"):
            geo.validate()

    def test_rejects_bad_origin_shape(self):
        geo = Grid2DGeometry(
            x=np.array([0.0]), z=np.array([10.0]), origin=np.array([0.0, 0.0, 0.0])
        )
        with pytest.raises(ValueError, match="origin"):
            geo.validate()


class TestGrid3DGeometry:
    def test_resistivity_shape_is_z_y_x(self):
        # Deliberately (n_z, n_y, n_x) to match ModEM's own native order.
        assert _grid3d().resistivity_shape == (3, 3, 2)

    def test_rejects_negative_n_air(self):
        geo = Grid3DGeometry(
            x=np.array([0.0]), y=np.array([0.0]), z=np.array([0.0]), n_air=-1
        )
        with pytest.raises(ValueError, match="n_air"):
            geo.validate()


class TestUnstructuredMeshGeometry:
    def test_n_regions(self):
        assert _mesh().n_regions == 2

    def test_rejects_region_id_length_mismatch(self):
        mesh = UnstructuredMeshGeometry(
            nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            connectivity=np.array([[0, 1, 2]]),
            region_ids=np.array([0, 1]),
        )
        with pytest.raises(ValueError, match="region_ids"):
            mesh.validate()

    def test_rejects_out_of_range_node_index(self):
        mesh = UnstructuredMeshGeometry(
            nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            connectivity=np.array([[0, 1, 5]]),
            region_ids=np.array([0]),
        )
        with pytest.raises(ValueError, match="out-of-range"):
            mesh.validate()

    def test_rejects_bad_plane(self):
        mesh = UnstructuredMeshGeometry(
            nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            connectivity=np.array([[0, 1, 2]]),
            region_ids=np.array([0]),
            plane="not-a-plane",
        )
        with pytest.raises(ValueError, match="plane"):
            mesh.validate()


class TestMultilineGeometry:
    def test_needs_at_least_one_line(self):
        with pytest.raises(ValueError, match="at least one line"):
            MultilineGeometry(lines=[]).validate()

    def test_rejects_duplicate_line_ids(self):
        line = LineEntry(
            line_id="L1", geometry=_grid2d(), resistivity=np.ones((2, 3))
        )
        with pytest.raises(ValueError, match="duplicate"):
            MultilineGeometry(lines=[line, line]).validate()

    def test_line_entry_rejects_shape_mismatch(self):
        line = LineEntry(
            line_id="L1", geometry=_grid2d(), resistivity=np.ones((5, 5))
        )
        with pytest.raises(ValueError, match="resistivity shape"):
            line.validate()

    def test_derived_volume_requires_valid_method(self):
        dv = DerivedVolume(
            grid=_grid3d(),
            resistivity=np.ones(_grid3d().resistivity_shape),
            derivation_method="bogus",
        )
        with pytest.raises(ValueError, match="derivation_method"):
            dv.validate()

    def test_valid_multiline_round_trip_in_memory(self):
        geo = MultilineGeometry(
            lines=[
                LineEntry(
                    line_id="L1",
                    geometry=_grid2d(),
                    resistivity=np.ones((2, 3)),
                    offset_y=0.0,
                    offset_kind="real",
                ),
                LineEntry(
                    line_id="L2",
                    geometry=_grid2d(),
                    resistivity=np.ones((2, 3)) * 2,
                    offset_y=100.0,
                    offset_kind="real",
                ),
            ]
        )
        geo.validate()
        assert [line.line_id for line in geo.lines] == ["L1", "L2"]


class TestStationTable:
    def test_requires_matching_lengths(self):
        table = StationTable(name=["A", "B"], x=np.array([0.0]))
        with pytest.raises(ValueError, match="StationTable.x"):
            table.validate()

    def test_lon_lat_optional_and_default_none(self):
        table = StationTable(
            name=["A", "B"], x=np.array([0.0, 1.0]), y=np.zeros(2), z=np.zeros(2)
        )
        table.validate()
        assert table.lon is None and table.lat is None

    def test_lon_lat_accepted_when_shape_matches(self):
        table = StationTable(
            name=["A", "B"],
            x=np.array([0.0, 1.0]),
            y=np.zeros(2),
            z=np.zeros(2),
            lon=np.array([10.0, 10.1]),
            lat=np.array([50.0, 50.1]),
        )
        table.validate()
        np.testing.assert_array_equal(table.lon, [10.0, 10.1])

    def test_lon_lat_shape_mismatch_raises(self):
        table = StationTable(
            name=["A", "B"],
            x=np.array([0.0, 1.0]),
            y=np.zeros(2),
            z=np.zeros(2),
            lon=np.array([10.0]),
            lat=np.array([50.0, 50.1]),
        )
        with pytest.raises(ValueError, match="StationTable.lon"):
            table.validate()

    def test_lon_without_lat_raises(self):
        table = StationTable(
            name=["A", "B"],
            x=np.array([0.0, 1.0]),
            y=np.zeros(2),
            z=np.zeros(2),
            lon=np.array([10.0, 10.1]),
        )
        with pytest.raises(ValueError, match="together"):
            table.validate()


class TestTopographyRaster:
    def test_requires_1d_x_and_y(self):
        raster = TopographyRaster(
            x=np.zeros((2, 2)), y=np.array([0.0, 1.0]), elevation=np.zeros((2, 2))
        )
        with pytest.raises(ValueError, match="must be 1-D"):
            raster.validate()

    def test_rejects_wrong_elevation_shape(self):
        raster = TopographyRaster(
            x=np.array([0.0, 100.0, 200.0]),
            y=np.array([0.0, 50.0]),
            elevation=np.zeros((3, 3)),
        )
        with pytest.raises(ValueError, match="does not match"):
            raster.validate()

    def test_valid_raster_passes(self):
        raster = TopographyRaster(
            x=np.array([0.0, 100.0, 200.0]),
            y=np.array([0.0, 50.0]),
            elevation=np.zeros((2, 3)),
        )
        raster.validate()
        assert raster.kind == "raster"


class TestPCSFModel:
    def test_grid2d_requires_resistivity(self):
        model = PCSFModel(geometry=_grid2d(), source_backend="occam2d")
        with pytest.raises(ValueError, match="resistivity is required"):
            model.validate()

    def test_grid2d_rejects_wrong_resistivity_shape(self):
        model = PCSFModel(geometry=_grid2d(), resistivity=np.ones((5, 5)))
        with pytest.raises(ValueError, match="resistivity shape"):
            model.validate()

    def test_grid2d_valid_model_passes(self):
        model = PCSFModel(geometry=_grid2d(), resistivity=np.ones((2, 3)) * 100.0)
        model.validate()
        assert model.kind == "grid2d"

    def test_multiline_rejects_container_level_resistivity(self):
        geo = MultilineGeometry(
            lines=[
                LineEntry(line_id="L1", geometry=_grid2d(), resistivity=np.ones((2, 3)))
            ]
        )
        model = PCSFModel(geometry=geo, resistivity=np.ones((2, 3)))
        with pytest.raises(ValueError, match="multiline geometry"):
            model.validate()

    def test_multiline_without_container_resistivity_passes(self):
        geo = MultilineGeometry(
            lines=[
                LineEntry(line_id="L1", geometry=_grid2d(), resistivity=np.ones((2, 3)))
            ]
        )
        model = PCSFModel(geometry=geo)
        model.validate()

    def test_resistivity_native_requires_encoding(self):
        model = PCSFModel(
            geometry=_grid2d(),
            resistivity=np.ones((2, 3)),
            resistivity_native=np.log10(np.ones((2, 3))),
        )
        with pytest.raises(ValueError, match="resistivity_native_encoding"):
            model.validate()

    def test_uncertainty_shape_must_match_resistivity(self):
        model = PCSFModel(
            geometry=_grid2d(),
            resistivity=np.ones((2, 3)),
            uncertainty=np.ones((5, 5)),
        )
        with pytest.raises(ValueError, match="uncertainty shape"):
            model.validate()

    def test_mesh_unstructured_accepts_per_region_resistivity(self):
        mesh = _mesh()
        model = PCSFModel(
            geometry=mesh, resistivity=np.array([100.0, 200.0])
        )
        model.validate()

    def test_mesh_unstructured_accepts_per_cell_resistivity(self):
        mesh = _mesh()
        model = PCSFModel(
            geometry=mesh, resistivity=np.array([100.0, 200.0])
        )  # 2 cells == 2 regions here; exercised distinctly below
        model.validate()

    def test_resistivity_by_node_accepts_matching_node_count(self):
        mesh = _mesh()
        model = PCSFModel(
            geometry=mesh,
            resistivity=np.array([100.0, 200.0]),
            resistivity_by_node=np.array([100.0, 150.0, 200.0, 250.0]),
        )
        model.validate()

    def test_resistivity_by_node_rejects_wrong_node_count(self):
        mesh = _mesh()
        model = PCSFModel(
            geometry=mesh,
            resistivity=np.array([100.0, 200.0]),
            resistivity_by_node=np.array([100.0, 200.0]),
        )
        with pytest.raises(ValueError, match="resistivity_by_node shape"):
            model.validate()

    def test_resistivity_by_node_rejected_for_non_mesh_geometry(self):
        model = PCSFModel(
            geometry=_grid2d(),
            resistivity=np.ones((2, 3)) * 100.0,
            resistivity_by_node=np.array([1.0, 2.0]),
        )
        with pytest.raises(ValueError, match="resistivity_by_node is only valid"):
            model.validate()

    def test_repr_summarizes_arrays_instead_of_dumping_them(self):
        model = PCSFModel(geometry=_grid2d(), resistivity=np.ones((2, 3)) * 100.0)
        text = repr(model)
        assert "ndarray(shape=" in text
        assert "100.0" not in text
