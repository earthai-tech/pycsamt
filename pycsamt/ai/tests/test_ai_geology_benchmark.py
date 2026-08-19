"""Contracts for seeded publication-benchmark geology families."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.geology import (
    ID_BENCHMARK_FAMILIES,
    OOD_BENCHMARK_FAMILIES,
    BenchmarkGeology,
    GeologyGrid,
    generate_benchmark_geology,
)

# A self-contained stand-in for the frozen DUHI-paper geology_v1.yaml
# config (that file lives in the private, gitignored DUHI-paper/ tree and
# is not distributed with the package). Same schema and parameter ranges
# generate_benchmark_geology() expects, so the public contract is tested
# without depending on paper-reproduction material outside the repo.
CONFIG = {
    "base_layered": {
        "cover_resistivity_ohm_m": [20.0, 80.0],
        "host_resistivity_ohm_m": [90.0, 350.0],
        "basement_resistivity_ohm_m": [600.0, 1500.0],
        "first_interface_depth_m": [300.0, 510.0],
        "second_interface_depth_m": [900.0, 1200.0],
        "interface_length_x_m": [600.0, 1800.0],
        "interface_length_z_m": 300.0,
        "minimum_thickness_m": 150.0,
    },
    "id_families": {
        "layered": {
            "interface_relief_std_m": [15.0, 90.0],
        },
        "intrusion": {
            "center_x_fraction": [0.30, 0.70],
            "center_depth_m": [450.0, 1050.0],
            "radius_x_m": [450.0, 1050.0],
            "radius_z_m": [150.0, 360.0],
            "resistivity_ohm_m": [3.0, 18.0],
            "dip_deg": [-35.0, 35.0],
            "transition_fraction": [0.05, 0.25],
        },
        "intrusion_halo": {
            "center_x_fraction": [0.30, 0.70],
            "center_depth_m": [450.0, 1050.0],
            "halo_radius_x_m": [750.0, 1350.0],
            "halo_radius_z_m": [300.0, 540.0],
            "halo_resistivity_ohm_m": [18.0, 60.0],
            "core_radius_fraction": [0.40, 0.65],
            "core_resistivity_ohm_m": [2.0, 10.0],
            "dip_deg": [-35.0, 35.0],
        },
        "dipping_fault": {
            "dip_deg": [45.0, 75.0],
            "trace_x_fraction": [0.30, 0.60],
            "throw_cells": [1, 3],
        },
        "multiple_body": {
            "body_count": [2, 4],
            "radius_x_m": [240.0, 480.0],
            "radius_z_m": [120.0, 270.0],
            "conductive_resistivity_ohm_m": [3.0, 20.0],
            "resistive_resistivity_ohm_m": [1600.0, 3000.0],
        },
        "correlated_heterogeneous": {
            "correlation_length_x_m": [600.0, 1800.0],
            "correlation_length_z_m": [150.0, 450.0],
            "log10_resistivity_mean": [1.7, 2.5],
            "log10_resistivity_std": [0.25, 0.55],
        },
    },
    "ood_families": {
        "extreme_dip": {
            "dip_deg": [20.0, 35.0],
            "trace_x_fraction": [0.25, 0.65],
            "throw_cells": [2, 4],
        },
        "deep_small_conductor": {
            "center_x_fraction": [0.25, 0.75],
            "center_depth_m": [1350.0, 1650.0],
            "radius_x_m": [150.0, 300.0],
            "radius_z_m": [90.0, 180.0],
            "resistivity_ohm_m": [1.0, 6.0],
        },
        "overlapping_multibody": {
            "body_count": [3, 5],
            "centre_spread_x_m": 360.0,
            "centre_spread_z_m": 180.0,
            "radius_x_m": [450.0, 750.0],
            "radius_z_m": [210.0, 390.0],
            "resistivity_ohm_m": [2.0, 25.0],
        },
        "high_contrast_halo": {
            "center_x_fraction": [0.30, 0.70],
            "center_depth_m": [450.0, 1050.0],
            "halo_radius_x_m": [900.0, 1500.0],
            "halo_radius_z_m": [360.0, 600.0],
            "halo_resistivity_ohm_m": [5.0, 15.0],
            "core_radius_fraction": [0.35, 0.55],
            "core_resistivity_ohm_m": [0.3, 1.5],
            "dip_deg": [-45.0, 45.0],
        },
        "anisotropic_correlation": {
            "correlation_length_z_m": [120.0, 240.0],
            "horizontal_vertical_ratio": [8.0, 14.0],
            "log10_resistivity_mean": [1.7, 2.5],
            "log10_resistivity_std": [0.35, 0.65],
        },
        "rugged_interfaces": {
            "interface_relief_std_m": [150.0, 260.0],
            "interface_length_x_m": [450.0, 900.0],
        },
    },
}


def _grid() -> GeologyGrid:
    return GeologyGrid.regular_2d(nx=20, nz=12, dx_m=300, dz_m=150)


@pytest.mark.parametrize(
    "regime,families",
    [
        ("id", ID_BENCHMARK_FAMILIES),
        ("ood", OOD_BENCHMARK_FAMILIES),
    ],
)
def test_all_frozen_families_generate_valid_models(regime, families):
    for index, family in enumerate(families):
        model = generate_benchmark_geology(
            _grid(),
            family,
            seed=100 + index,
            configuration=CONFIG,
            regime=regime,
        )
        assert isinstance(model, BenchmarkGeology)
        assert model.resistivity_ohm_m.shape == (12, 20)
        assert np.all(np.isfinite(model.resistivity_ohm_m))
        assert np.all(model.resistivity_ohm_m > 0)
        assert len(model.model_hash) == 64


def test_benchmark_geology_is_seed_deterministic():
    first = generate_benchmark_geology(
        _grid(),
        "intrusion_halo",
        seed=42,
        configuration=CONFIG,
    )
    second = generate_benchmark_geology(
        _grid(),
        "intrusion_halo",
        seed=42,
        configuration=CONFIG,
    )
    np.testing.assert_array_equal(
        first.resistivity_ohm_m, second.resistivity_ohm_m
    )
    assert first.parameters == second.parameters
    assert first.model_hash == second.model_hash


def test_deep_small_conductor_snaps_only_when_raster_would_be_empty():
    model = generate_benchmark_geology(
        _grid(),
        "deep_small_conductor",
        seed=3237934865,
        configuration=CONFIG,
        regime="ood",
    )
    lens = model.parameters["lenses"][0]
    assert model.parameters["center_snapped_to_grid"] is True
    assert lens["center_x_m"] in model.grid.x_m
    assert lens["center_z_m"] in model.grid.z_m
    assert np.any(
        np.isclose(
            model.resistivity_ohm_m,
            lens["resistivity_ohm_m"],
        )
    )


def test_benchmark_geology_rejects_cross_regime_family():
    with pytest.raises(ValueError, match="incompatible"):
        generate_benchmark_geology(
            _grid(),
            "extreme_dip",
            seed=0,
            configuration=CONFIG,
            regime="id",
        )
