"""Contracts for seeded publication-benchmark geology families."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from pycsamt.ai.geology import (
    ID_BENCHMARK_FAMILIES,
    OOD_BENCHMARK_FAMILIES,
    BenchmarkGeology,
    GeologyGrid,
    generate_benchmark_geology,
)

ROOT = Path(__file__).resolve().parents[3]
CONFIG = yaml.safe_load(
    (ROOT / "DUHI-paper" / "configs" / "geology_v1.yaml").read_text(
        encoding="utf-8"
    )
)


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
