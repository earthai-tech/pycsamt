"""Tests for empirical field-calibrated survey corruption."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.data.contracts import SurveyData
from pycsamt.ai.domain_gap import apply_empirical_corruption


def _survey() -> SurveyData:
    return SurveyData(
        np.full((3, 2, 2), 10.0 + 5.0j),
        [10.0, 1.0],
        ["A", "B", "C"],
        ["zxy", "zyx"],
        [[0, 0], [1, 0], [2, 0]],
    )


def _kwargs() -> dict:
    return {
        "field_frequencies_hz": [100.0, 10.0, 1.0],
        "static_log10_resistivity_profiles": [
            [0.2, 0.2, 0.2],
            [-0.2, -0.2, -0.2],
        ],
        "measurement_reliability_profiles": [
            [0.8, 0.7, 0.6],
            [0.6, 0.5, 0.4],
        ],
        "dimensionality_reliability_profiles": [
            [1.0, 0.5, 0.2],
            [0.8, 0.4, 0.1],
        ],
        "observation_reliability_profiles": [
            [0.8, 0.35, 0.12],
            [0.48, 0.2, 0.04],
        ],
        "relative_error_quantiles": {
            "zxy": {"levels": [0.0, 1.0], "values": [0.01, 0.05]},
            "zyx": {"levels": [0.0, 1.0], "values": [0.02, 0.08]},
        },
    }


def test_empirical_corruption_is_seed_deterministic_and_non_mutating():
    clean = _survey()
    original = np.array(clean.impedance)
    first = apply_empirical_corruption(clean, seed=4, **_kwargs())
    second = apply_empirical_corruption(clean, seed=4, **_kwargs())
    np.testing.assert_array_equal(
        first.survey.impedance, second.survey.impedance
    )
    np.testing.assert_array_equal(clean.impedance, original)
    assert first.record_hash == second.record_hash


def test_empirical_corruption_returns_all_latent_arrays():
    result = apply_empirical_corruption(_survey(), seed=2, **_kwargs())
    assert result.survey.shape == (3, 2, 2)
    assert result.static_log10_resistivity_factor.shape == (3, 2)
    assert result.noise_realization.shape == result.survey.shape
    assert result.observation_reliability.shape == result.survey.shape
    assert result.survey.impedance_error is not None
    assert np.all(result.relative_error_fraction[:, :, 1] >= 0.02)


def test_static_resistivity_factor_uses_impedance_square_root():
    kwargs = _kwargs()
    kwargs["static_log10_resistivity_profiles"] = [[2.0, 2.0, 2.0]]
    kwargs["measurement_reliability_profiles"] = [[1.0, 1.0, 1.0]]
    kwargs["dimensionality_reliability_profiles"] = [[1.0, 1.0, 1.0]]
    kwargs["observation_reliability_profiles"] = [[1.0, 1.0, 1.0]]
    kwargs["relative_error_quantiles"] = {
        component: {"levels": [0.0, 1.0], "values": [0.0, 0.0]}
        for component in ("zxy", "zyx")
    }
    result = apply_empirical_corruption(_survey(), seed=0, **kwargs)
    np.testing.assert_allclose(result.survey.impedance, 10.0 * _survey().impedance)


def test_empirical_corruption_applies_component_missing_rates():
    result = apply_empirical_corruption(
        _survey(),
        seed=0,
        missing_rate_by_component={"zxy": 1.0, "zyx": 0.0},
        **_kwargs(),
    )
    assert np.all(result.missing_mask[:, :, 0])
    assert not np.any(result.missing_mask[:, :, 1])
    assert result.survey.n_valid == 6


def test_empirical_corruption_rejects_unaligned_profiles():
    kwargs = _kwargs()
    kwargs["measurement_reliability_profiles"] = [[0.5, 0.5]]
    with pytest.raises(ValueError, match="measurement_reliability_profiles"):
        apply_empirical_corruption(_survey(), seed=0, **kwargs)
