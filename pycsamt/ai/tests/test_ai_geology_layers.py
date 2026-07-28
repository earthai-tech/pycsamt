"""Contracts for :mod:`pycsamt.ai.geology.layers`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.geology import (
    ElectricalLayer,
    GaussianCorrelation,
    GeologyGrid,
    LayeredGeology,
    generate_layered_geology,
)


def _units_2d():
    return (
        ElectricalLayer("cover", 20),
        ElectricalLayer("sediments", 100),
        ElectricalLayer("basement", 1000),
    )


def test_electrical_layer_validation_serialization_and_log_property():
    layer = ElectricalLayer(
        "basement",
        1000,
        log10_std=0.2,
        heterogeneity=GaussianCorrelation(500, 100),
        resistivity_bounds_ohm_m=(100, 5000),
    )
    assert layer.log10_resistivity == pytest.approx(3)
    assert ElectricalLayer.from_dict(layer.to_dict()) == layer
    with pytest.raises(ValueError, match="requires"):
        ElectricalLayer("invalid", 100, log10_std=0.1)
    with pytest.raises(ValueError, match="increasing"):
        ElectricalLayer("invalid", 100, resistivity_bounds_ohm_m=(100, 10))


def test_flat_2d_layers_have_expected_indices_and_resistivity():
    grid = GeologyGrid.regular_2d(nx=8, nz=8, dx_m=100, dz_m=50)
    model = generate_layered_geology(
        grid, _units_2d(), [100, 250], seed=1, interface_policy="raise"
    )
    assert model.interface_depth_m.shape == (2, 8)
    np.testing.assert_array_equal(model.interface(0), 100)
    assert set(np.unique(model.layer_index)) == {0, 1, 2}
    for index, expected in enumerate((20, 100, 1000)):
        np.testing.assert_array_equal(
            model.resistivity_ohm_m[model.layer_mask(index)], expected
        )
    assert model.adjusted_interface_fraction == 0
    assert model.summary()["n_layers"] == 3
    assert model.generation_config["mean_interface_depth_m"] == (100.0, 250.0)
    with pytest.raises(TypeError):
        model.generation_config["new"] = "value"


def test_correlated_interfaces_are_seed_deterministic_and_vary_laterally():
    grid = GeologyGrid.regular_2d(nx=32, nz=20, dx_m=100, dz_m=50)
    kwargs = dict(
        grid=grid,
        layers=_units_2d(),
        mean_interface_depth_m=[250, 650],
        seed=9,
        interface_relief_std_m=[30, 60],
        interface_correlation=GaussianCorrelation(500, 100),
    )
    first = generate_layered_geology(**kwargs)
    second = generate_layered_geology(**kwargs)
    np.testing.assert_array_equal(first.interface_depth_m, second.interface_depth_m)
    np.testing.assert_array_equal(first.resistivity_ohm_m, second.resistivity_ohm_m)
    assert np.std(first.interface(0)) > 0
    assert first.model_hash == second.model_hash


def test_projection_records_adjustments_and_strict_policy_rejects():
    grid = GeologyGrid.regular_2d(nx=24, nz=10, dx_m=100, dz_m=50)
    kwargs = dict(
        grid=grid,
        layers=_units_2d(),
        mean_interface_depth_m=[100, 350],
        seed=2,
        interface_relief_std_m=[150, 150],
        interface_correlation=GaussianCorrelation(300, 100),
        minimum_thickness_m=50,
    )
    projected = generate_layered_geology(**kwargs, interface_policy="project")
    assert projected.adjusted_interface_fraction > 0
    boundaries = np.concatenate(
        [
            np.zeros((1, grid.shape[1])),
            projected.interface_depth_m,
            np.full((1, grid.shape[1]), 500),
        ]
    )
    assert np.all(np.diff(boundaries, axis=0) >= 50 - 1e-9)
    with pytest.raises(ValueError, match="violate"):
        generate_layered_geology(**kwargs, interface_policy="raise")


def test_heterogeneous_layer_is_bounded_and_other_units_stay_constant():
    grid = GeologyGrid.regular_2d(nx=24, nz=16, dx_m=100, dz_m=50)
    units = (
        ElectricalLayer("cover", 10),
        ElectricalLayer(
            "variable",
            100,
            log10_std=0.4,
            heterogeneity=GaussianCorrelation(400, 100),
            resistivity_bounds_ohm_m=(30, 300),
        ),
        ElectricalLayer("basement", 1000),
    )
    model = generate_layered_geology(grid, units, [150, 500], seed=4)
    variable = model.resistivity_ohm_m[model.layer_mask("variable")]
    assert np.std(variable) > 0
    assert np.min(variable) >= 30 and np.max(variable) <= 300
    np.testing.assert_array_equal(
        model.resistivity_ohm_m[model.layer_mask("cover")], 10
    )


def test_3d_layered_volume_and_rotated_interface_correlation():
    grid = GeologyGrid.regular_3d(nx=12, ny=10, nz=8, dx_m=100, dy_m=100, dz_m=50)
    units = (ElectricalLayer("cover", 20), ElectricalLayer("basement", 800))
    model = generate_layered_geology(
        grid,
        units,
        [150],
        seed=5,
        interface_relief_std_m=25,
        interface_correlation=GaussianCorrelation(
            500, 100, length_y_m=200, azimuth_deg=30
        ),
    )
    assert model.interface_depth_m.shape == (1, 10, 12)
    assert model.resistivity_ohm_m.shape == (8, 10, 12)


def test_single_layer_model_has_empty_interfaces():
    grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
    model = generate_layered_geology(
        grid, [ElectricalLayer("halfspace", 100)], [], seed=0
    )
    assert model.interface_depth_m.shape == (0, 4)
    assert model.adjusted_interface_fraction == 0
    np.testing.assert_array_equal(model.resistivity_ohm_m, 100)


def test_layered_npz_roundtrip_preserves_model_hash(tmp_path):
    model = generate_layered_geology(
        GeologyGrid.regular_2d(nx=8, nz=8, dx_m=100, dz_m=50),
        _units_2d(),
        [100, 250],
        seed=7,
    )
    restored = LayeredGeology.from_npz(model.to_npz(tmp_path / "layers.npz"))
    assert restored.model_hash == model.model_hash
    np.testing.assert_array_equal(restored.interface_depth_m, model.interface_depth_m)


def test_generation_rejects_count_correlation_and_feasibility_errors():
    grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
    with pytest.raises(ValueError, match="contain 2"):
        generate_layered_geology(grid, _units_2d(), [100], seed=0)
    with pytest.raises(ValueError, match="requires interface_correlation"):
        generate_layered_geology(
            grid,
            [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
            [100],
            seed=0,
            interface_relief_std_m=10,
        )
    with pytest.raises(ValueError, match="infeasible"):
        generate_layered_geology(
            grid,
            _units_2d(),
            [50, 100],
            seed=0,
            minimum_thickness_m=100,
        )
