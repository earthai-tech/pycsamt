"""Tests for geology-to-Maxwell mesh construction."""

import numpy as np
import pytest

from pycsamt.ai.geology import GeologyGrid, TopographicSurface
from pycsamt.forward.maxwell import (
    MeshDesign,
    ReceiverSet,
    SolverMeshModel,
    build_solver_mesh,
    skin_depth_m,
)


def _design(**changes):
    values = dict(
        horizontal_padding_cells=1,
        bottom_padding_cells=1,
        air_layers=2,
        padding_expansion=1.2,
        air_expansion=1.1,
    )
    values.update(changes)
    return MeshDesign(**values)


def test_skin_depth_broadcasting_and_validation():
    assert float(skin_depth_m(100, 1)) == pytest.approx(5032.921, rel=1e-5)
    assert skin_depth_m([100, 400], 1).shape == (2,)
    with pytest.raises(ValueError, match="resistivity"):
        skin_depth_m(0, 1)
    with pytest.raises(ValueError, match="frequency"):
        skin_depth_m(100, 0)


def test_design_normalization_roundtrip_and_validation():
    design = _design(horizontal_padding_cells=(2, 3))
    assert design.horizontal_padding == (2, 3)
    assert MeshDesign.from_dict(design.to_dict()) == design
    with pytest.raises(ValueError, match="two values"):
        MeshDesign(horizontal_padding_cells=(1, 2, 3))
    with pytest.raises(ValueError, match="greater than"):
        MeshDesign(padding_expansion=0.9)


def test_build_2d_padding_air_core_and_resistivity_mapping():
    grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=100, dz_m=50)
    resistivity = np.array([[100, 200, 300], [400, 500, 600]])
    model = build_solver_mesh(
        grid,
        resistivity_ohm_m=resistivity,
        frequencies_hz=[10, 1],
        design=_design(),
    )
    assert model.mesh.shape == (5, 5)
    np.testing.assert_allclose(
        model.conductivity_s_m[model.core_slices], 1 / resistivity
    )
    assert np.all(model.air_mask[:2])
    assert model.quality.cell_count == 25
    with pytest.raises(ValueError):
        model.conductivity_s_m[0, 0] = 1


def test_2d_topography_changes_air_earth_boundary():
    grid = GeologyGrid.regular_2d(nx=3, nz=3, dx_m=100, dz_m=50)
    topography = TopographicSurface(grid, [100, 50, 100], 100)
    model = build_solver_mesh(
        grid,
        conductivity_s_m=np.ones(grid.shape) * 0.01,
        frequencies_hz=[1],
        topography=topography,
        design=_design(air_layers=1),
    )
    core_earth = model.earth_mask[model.core_slices]
    assert not core_earth[0, 1]
    assert core_earth[1, 1]
    assert (
        model.conductivity_s_m[model.core_slices][0, 1]
        == model.design.air_conductivity_s_m
    )


def test_build_3d_shape_padding_and_topography():
    grid = GeologyGrid.regular_3d(nx=2, ny=3, nz=2, dx_m=100, dy_m=80, dz_m=50)
    elevation = np.array([[100, 100], [90, 90], [80, 80]])
    topography = TopographicSurface(grid, elevation, 100)
    model = build_solver_mesh(
        grid,
        resistivity_ohm_m=np.full(grid.shape, 100),
        frequencies_hz=[1],
        topography=topography,
        design=_design(air_layers=1),
    )
    assert model.mesh.shape == (4, 5, 4)
    assert model.conductivity_s_m[model.core_slices].shape == grid.shape
    assert model.mesh.dimension == 3


def test_mesh_model_archive_hash_and_problem_conversion(tmp_path):
    grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=100, dz_m=50)
    model = build_solver_mesh(
        grid,
        resistivity_ohm_m=np.full(grid.shape, 100),
        frequencies_hz=[1],
        design=_design(),
    )
    restored = SolverMeshModel.from_npz(model.to_npz(tmp_path / "mesh.npz"))
    assert restored.model_hash == model.model_hash
    receivers = ReceiverSet([[150, 0]], ["S00"])
    problem = restored.to_problem([1], receivers, mark_air_inactive=True)
    assert np.array_equal(problem.active_cells, restored.earth_mask)
    assert problem.metadata["mesh_model_hash"] == restored.model_hash


def test_receiver_and_source_validation():
    grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=10, dz_m=10)
    model = build_solver_mesh(
        grid,
        conductivity_s_m=np.ones(grid.shape),
        frequencies_hz=[1],
        design=_design(),
    )
    outside = ReceiverSet([[1e6, 0]], ["far"])
    assert "outside" in model.assess_receivers(outside)[0]
    with pytest.raises(ValueError, match="outside"):
        model.to_problem([1], outside)
    with pytest.raises(ValueError, match="exactly one"):
        build_solver_mesh(grid, frequencies_hz=[1])
    with pytest.raises(ValueError, match="exactly one"):
        build_solver_mesh(
            grid,
            conductivity_s_m=np.ones(grid.shape),
            resistivity_ohm_m=np.ones(grid.shape),
            frequencies_hz=[1],
        )
    below = ReceiverSet([[5, 5]], ["buried"])
    assert "below local terrain" in model.assess_receivers(below)[0]


def test_shifted_depth_origin_is_rejected():
    grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=10, dz_m=10, z_origin_m=100)
    with pytest.raises(ValueError, match="depth zero"):
        build_solver_mesh(
            grid,
            conductivity_s_m=np.ones(grid.shape),
            frequencies_hz=[1],
            design=_design(),
        )
