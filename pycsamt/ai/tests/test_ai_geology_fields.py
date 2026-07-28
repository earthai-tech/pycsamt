"""Contracts for :mod:`pycsamt.ai.geology.fields`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.geology import (
    CorrelatedField,
    DirectionalVariogram,
    GaussianCorrelation,
    GeologyGrid,
    directional_variogram,
    generate_gaussian_field,
)


def test_regular_grid_2d_axes_extent_serialization_and_readonly():
    grid = GeologyGrid.regular_2d(
        nx=4,
        nz=3,
        dx_m=100,
        dz_m=50,
        x_origin_m=1000,
        z_origin_m=0,
        crs="EPSG:32630",
    )
    assert grid.dimension == 2
    assert grid.shape == (3, 4)
    assert grid.spacing_m == (50.0, 100.0)
    assert grid.extent_m == {"x": (1000.0, 1400.0), "z": (0.0, 150.0)}
    assert not grid.x_m.flags.writeable
    restored = GeologyGrid.from_dict(grid.to_dict())
    np.testing.assert_array_equal(restored.x_m, grid.x_m)
    np.testing.assert_array_equal(restored.z_m, grid.z_m)


def test_regular_grid_3d_uses_canonical_zyx_shape():
    grid = GeologyGrid.regular_3d(nx=5, ny=4, nz=3, dx_m=100, dy_m=200, dz_m=50)
    assert grid.dimension == 3
    assert grid.shape == (3, 4, 5)
    assert grid.spacing_m == (50.0, 200.0, 100.0)


def test_grid_rejects_irregular_axes_for_spectral_spacing():
    with pytest.raises(ValueError, match="regularly spaced"):
        GeologyGrid([0, 1, 3], [0, 1, 2])
    with pytest.raises(ValueError, match="at least two"):
        GeologyGrid.regular_2d(nx=1, nz=2, dx_m=1, dz_m=1)


def test_correlation_validates_dimensions_and_normalizes_azimuth():
    model = GaussianCorrelation(1000, 100, length_y_m=500, azimuth_deg=210)
    assert model.azimuth_deg == pytest.approx(30)
    assert model.anisotropy_xz == pytest.approx(10)
    assert GaussianCorrelation.from_dict(model.to_dict()) == model
    grid3d = GeologyGrid.regular_3d(nx=2, ny=2, nz=2, dx_m=1, dy_m=1, dz_m=1)
    with pytest.raises(ValueError, match="length_y_m"):
        GaussianCorrelation(2, 1).validate_grid(grid3d)


@pytest.mark.parametrize("boundary", ["periodic", "reflect"])
def test_2d_field_is_seed_deterministic_standardized_and_hashed(boundary):
    grid = GeologyGrid.regular_2d(nx=32, nz=16, dx_m=100, dz_m=50)
    model = GaussianCorrelation(400, 100)
    first = generate_gaussian_field(grid, model, seed=7, boundary=boundary)
    second = generate_gaussian_field(grid, model, seed=7, boundary=boundary)
    np.testing.assert_array_equal(first.values, second.values)
    assert np.mean(first.values) == pytest.approx(0, abs=1e-12)
    assert np.std(first.values) == pytest.approx(1, abs=1e-12)
    assert first.field_hash == second.field_hash
    assert not first.values.flags.writeable


def test_3d_field_azimuth_generation_and_rescale():
    grid = GeologyGrid.regular_3d(nx=12, ny=10, nz=8, dx_m=100, dy_m=100, dz_m=50)
    model = GaussianCorrelation(500, 100, length_y_m=200, azimuth_deg=30)
    field = generate_gaussian_field(grid, model, seed=2)
    assert field.values.shape == (8, 10, 12)
    scaled = field.rescale(mean=2.0, standard_deviation=0.5)
    assert np.mean(scaled.values) == pytest.approx(2.0)
    assert np.std(scaled.values) == pytest.approx(0.5)
    assert not scaled.standardized


def test_field_npz_roundtrip_preserves_identity(tmp_path):
    field = generate_gaussian_field(
        GeologyGrid.regular_2d(nx=12, nz=8, dx_m=100, dz_m=50),
        GaussianCorrelation(300, 100),
        seed=9,
    )
    restored = CorrelatedField.from_npz(field.to_npz(tmp_path / "field.npz"))
    np.testing.assert_array_equal(restored.values, field.values)
    assert restored.provenance() == field.provenance()
    assert restored.field_hash == field.field_hash


def test_directional_variogram_shapes_counts_and_axis_validation():
    field = generate_gaussian_field(
        GeologyGrid.regular_2d(nx=16, nz=8, dx_m=100, dz_m=50),
        GaussianCorrelation(300, 100),
        seed=3,
    )
    result = directional_variogram(field, "x", max_lag_cells=4)
    assert isinstance(result, DirectionalVariogram)
    assert result.n_lags == 4
    np.testing.assert_array_equal(result.lag_m, [100, 200, 300, 400])
    assert np.all(np.diff(result.pair_count) < 0)
    assert np.all(result.semivariance >= 0)
    with pytest.raises(ValueError, match="axis must"):
        directional_variogram(field, "y")


def test_generation_contract_rejects_invalid_seed_boundary_and_constant_case():
    grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=1, dz_m=1)
    model = GaussianCorrelation(1, 1)
    with pytest.raises(TypeError, match="seed"):
        generate_gaussian_field(grid, model, seed=1.5)
    with pytest.raises(ValueError, match="boundary"):
        generate_gaussian_field(grid, model, seed=1, boundary="wrap")
    with pytest.raises(ValueError, match="max_lag_cells"):
        directional_variogram(
            generate_gaussian_field(grid, model, seed=1),
            "x",
            max_lag_cells=4,
        )
