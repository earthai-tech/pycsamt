"""Contracts for :mod:`pycsamt.ai.domain_gap.report`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.data.contracts import SurveyData
from pycsamt.ai.domain_gap.report import (
    compare_feature_distributions,
    compare_survey_distributions,
)
from pycsamt.ai.domain_gap.simulator import apply_corruption_suite


def _survey(n_station=6, n_frequency=10, error=None, seed=0):
    rng = np.random.default_rng(seed)
    magnitude = rng.uniform(50, 200, size=(n_station, n_frequency, 1))
    phase = np.deg2rad(rng.uniform(20, 70, size=(n_station, n_frequency, 1)))
    z = magnitude * np.exp(1j * phase)
    return SurveyData(
        z,
        np.linspace(1000.0, 1.0, n_frequency),
        [f"S{i}" for i in range(n_station)],
        ["xy"],
        np.column_stack([np.arange(n_station) * 100.0, np.zeros(n_station)]),
        impedance_error=error,
    )


def test_identical_samples_have_zero_ks_statistic_and_mean_difference():
    rng = np.random.default_rng(0)
    values = rng.normal(0, 1, 300)
    comparison = compare_feature_distributions(values, values, feature="demo")
    assert comparison.ks_statistic == pytest.approx(0.0)
    assert comparison.mean_difference == pytest.approx(0.0)
    assert comparison.std_ratio == pytest.approx(1.0)


def test_clearly_different_samples_have_large_ks_statistic_and_low_pvalue():
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, 300)
    b = rng.normal(10, 1, 300)
    comparison = compare_feature_distributions(a, b, feature="demo")
    assert comparison.ks_statistic > 0.9
    assert comparison.ks_pvalue < 0.01
    assert comparison.mean_difference == pytest.approx(-10.0, abs=0.5)


def test_empty_samples_yield_nan_statistics_without_raising():
    comparison = compare_feature_distributions(np.array([]), np.array([1.0, 2.0]), feature="demo")
    assert np.isnan(comparison.ks_statistic)
    assert np.isnan(comparison.ks_pvalue)
    assert comparison.simulated_stats["count"] == 0


def test_compare_survey_distributions_covers_default_features():
    survey = _survey()
    report = compare_survey_distributions(survey, survey)
    assert set(report.comparisons) == {
        "log_impedance_magnitude",
        "phase_deg",
        "error_to_magnitude_ratio",
    }
    # identical survey against itself: every finite comparison is ~0
    for feature in ("log_impedance_magnitude", "phase_deg"):
        assert report.comparisons[feature].ks_statistic == pytest.approx(0.0)


def test_error_ratio_feature_is_empty_without_declared_errors():
    survey = _survey()
    report = compare_survey_distributions(survey, survey)
    assert report.comparisons["error_to_magnitude_ratio"].simulated_stats["count"] == 0


def test_worst_feature_identifies_the_most_corrupted_feature():
    survey = _survey(n_station=15, n_frequency=15)
    corrupted, _ = apply_corruption_suite(survey, severity="held_out_corruption", seed=0)
    report = compare_survey_distributions(corrupted, survey)
    assert report.worst_feature() in report.comparisons


def test_worst_feature_raises_when_every_statistic_is_nan():
    empty_survey = _survey(n_station=1, n_frequency=1)
    # error feature will be empty (no declared error); force the others empty
    # by comparing a survey with no valid data against itself.
    from dataclasses import replace

    invalid = replace(
        empty_survey, impedance=np.full_like(empty_survey.impedance, np.nan)
    )
    report = compare_survey_distributions(invalid, invalid)
    with pytest.raises(ValueError, match="finite KS statistic"):
        report.worst_feature()


def test_compare_feature_distributions_rejects_unknown_feature_via_survey_path():
    survey = _survey()
    with pytest.raises(ValueError, match="unknown feature"):
        compare_survey_distributions(survey, survey, features=("not_a_feature",))


def test_to_dict_is_json_shaped():
    survey = _survey()
    report = compare_survey_distributions(survey, survey)
    payload = report.to_dict()
    assert isinstance(payload, dict)
    assert payload["phase_deg"]["feature"] == "phase_deg"
