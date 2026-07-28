"""Contracts for :mod:`pycsamt.ai.domain_gap.simulator`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.data.contracts import SurveyData
from pycsamt.ai.domain_gap.simulator import (
    SEVERITY_PRESETS,
    CorruptionConfig,
    CorruptionRecord,
    add_heteroscedastic_noise,
    apply_corruption_suite,
    apply_dropout,
    apply_error_floor,
    apply_galvanic_distortion,
    apply_static_shift,
    inject_outliers,
    perturb_coordinates,
)


def _survey(n_station=4, n_frequency=6, components=("xy", "yx"), error=None):
    z = np.full((n_station, n_frequency, len(components)), 100 + 50j, dtype=complex)
    return SurveyData(
        impedance=z,
        frequencies_hz=np.linspace(1000.0, 1.0, n_frequency),
        station_names=[f"S{i}" for i in range(n_station)],
        components=list(components),
        coordinates_m=np.column_stack(
            [np.arange(n_station) * 100.0, np.zeros(n_station), np.full(n_station, 50.0)]
        ),
        impedance_error=error,
    )


# --------------------------------------------------------------------------- #
# CorruptionConfig / CorruptionRecord
# --------------------------------------------------------------------------- #


def test_clean_preset_is_the_all_zero_default():
    assert SEVERITY_PRESETS["clean"] == CorruptionConfig()
    config = CorruptionConfig()
    assert config.noise_level_range == (0.0, 0.0)
    assert config.station_dropout_rate == 0.0


def test_config_round_trips_through_dict_and_hashes_deterministically():
    config = CorruptionConfig(noise_level_range=(0.01, 0.05), station_dropout_rate=0.1)
    restored = CorruptionConfig.from_dict(config.to_dict())
    assert restored == config
    assert restored.config_hash() == config.config_hash()
    assert config.config_hash() != CorruptionConfig().config_hash()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"noise_level_range": (0.5, 0.1)},
        {"noise_level_range": (-0.1, 0.1)},
        {"station_dropout_rate": 1.5},
        {"station_dropout_rate": -0.1},
        {"static_shift_log10_sigma": -1.0},
    ],
)
def test_config_rejects_invalid_ranges_and_rates(kwargs):
    with pytest.raises(ValueError):
        CorruptionConfig(**kwargs)


def test_record_to_dict_embeds_config_hash():
    record = CorruptionRecord(CorruptionConfig(), seed=7, sampled={"n_outliers": 0})
    payload = record.to_dict()
    assert payload["config_hash"] == CorruptionConfig().config_hash()
    assert payload["seed"] == 7


# --------------------------------------------------------------------------- #
# Individual corruption steps
# --------------------------------------------------------------------------- #


def test_heteroscedastic_noise_is_a_noop_for_zero_range():
    survey = _survey()
    out = add_heteroscedastic_noise(survey, level_range=(0.0, 0.0), rng=np.random.default_rng(0))
    np.testing.assert_array_equal(out.impedance, survey.impedance)
    assert out.impedance_error is None


def test_heteroscedastic_noise_perturbs_and_sets_error():
    survey = _survey()
    out = add_heteroscedastic_noise(
        survey, level_range=(0.05, 0.05), rng=np.random.default_rng(0)
    )
    assert not np.array_equal(out.impedance, survey.impedance)
    assert out.impedance_error is not None
    assert np.all(out.impedance_error[out.valid] > 0)
    assert out.n_valid == survey.n_valid


def test_heteroscedastic_noise_combines_existing_error_in_quadrature():
    error = np.full((4, 6, 2), 3.0)
    survey = _survey(error=error)
    out = add_heteroscedastic_noise(
        survey, level_range=(0.05, 0.05), rng=np.random.default_rng(0)
    )
    assert np.all(out.impedance_error[out.valid] >= 3.0)


def test_error_floor_is_noop_at_zero_and_raises_existing_error():
    survey = _survey(error=np.ones((4, 6, 2)))
    unchanged = apply_error_floor(survey, floor_fraction=0.0)
    np.testing.assert_array_equal(unchanged.impedance_error, survey.impedance_error)

    floored = apply_error_floor(survey, floor_fraction=0.5)
    magnitude = np.abs(survey.impedance)
    np.testing.assert_allclose(floored.impedance_error, 0.5 * magnitude)


def test_error_floor_initializes_error_when_absent_without_invalidating():
    survey = _survey()
    floored = apply_error_floor(survey, floor_fraction=0.02)
    assert floored.impedance_error is not None
    assert floored.n_valid == survey.n_valid


def test_static_shift_scales_impedance_and_error_by_same_factor():
    survey = _survey(error=np.ones((4, 6, 2)))
    shifted, info = apply_static_shift(
        survey, log10_sigma=0.2, rng=np.random.default_rng(0), return_info=True
    )
    factor = info["static_shift_factor"]
    np.testing.assert_allclose(
        shifted.impedance, survey.impedance * factor[:, None, None]
    )
    np.testing.assert_allclose(
        shifted.impedance_error, survey.impedance_error * factor[:, None, None]
    )


def test_static_shift_zero_sigma_is_noop():
    survey = _survey()
    out = apply_static_shift(survey, log10_sigma=0.0, rng=np.random.default_rng(0))
    np.testing.assert_array_equal(out.impedance, survey.impedance)


def test_galvanic_distortion_noop_when_all_params_zero():
    survey = _survey(components=("xx", "xy", "yx", "yy"))
    out = apply_galvanic_distortion(survey, rng=np.random.default_rng(0))
    np.testing.assert_array_equal(out.impedance, survey.impedance)


def test_galvanic_distortion_perturbs_full_tensor_and_preserves_shape():
    survey = _survey(components=("xx", "xy", "yx", "yy"), error=np.ones((4, 6, 4)))
    out, sampled = apply_galvanic_distortion(
        survey,
        gain_log10_sigma=0.05,
        twist_deg_sigma=15.0,
        shear_sigma=0.2,
        anisotropy_sigma=0.1,
        rng=np.random.default_rng(1),
        return_info=True,
    )
    assert out.shape == survey.shape
    assert not np.array_equal(out.impedance, survey.impedance)
    assert sampled["twist_deg"].shape == (4,)
    # error scales with gain only (first-order approximation), stays positive
    assert np.all(out.impedance_error > 0)


def test_galvanic_distortion_handles_missing_diagonal_components():
    survey = _survey(components=("xy", "yx"))
    out = apply_galvanic_distortion(
        survey, twist_deg_sigma=10.0, rng=np.random.default_rng(2)
    )
    assert out.components == ("xy", "yx")
    assert out.shape == survey.shape


def test_dropout_station_rate_one_invalidates_everything():
    survey = _survey()
    out = apply_dropout(survey, station_rate=1.0, rng=np.random.default_rng(0))
    assert out.n_valid == 0


def test_dropout_zero_rates_are_noop():
    survey = _survey()
    out = apply_dropout(survey, rng=np.random.default_rng(0))
    assert out.n_valid == survey.n_valid


def test_dropout_frequency_rate_drops_whole_frequency_columns():
    survey = _survey(n_station=6, n_frequency=8)
    out, info = apply_dropout(
        survey, frequency_rate=1.0, rng=np.random.default_rng(0), return_info=True
    )
    assert out.n_valid == 0
    assert info["dropped_frequency_count"] == 8


def test_inject_outliers_zero_rate_is_noop():
    survey = _survey()
    out = inject_outliers(survey, rate=0.0, rng=np.random.default_rng(0))
    np.testing.assert_array_equal(out.impedance, survey.impedance)


def test_inject_outliers_perturbs_a_subset_and_keeps_valid_mask():
    survey = _survey(n_station=2, n_frequency=20, components=("xy",))
    out, info = inject_outliers(
        survey, rate=0.5, rng=np.random.default_rng(0), return_info=True
    )
    assert info["n_outliers"] == int(round(0.5 * survey.n_valid))
    assert out.n_valid == survey.n_valid  # outliers stay "valid"
    n_changed = np.count_nonzero(~np.isclose(out.impedance, survey.impedance))
    assert n_changed == info["n_outliers"]


def test_perturb_coordinates_leaves_nan_elevation_untouched():
    survey = _survey()
    coordinates = np.array(survey.coordinates_m)
    coordinates[0, 2] = np.nan
    survey = SurveyData(
        survey.impedance,
        survey.frequencies_hz,
        survey.station_names,
        survey.components,
        coordinates,
    )
    moved = perturb_coordinates(
        survey, coordinate_sigma_m=10.0, elevation_sigma_m=5.0, rng=np.random.default_rng(0)
    )
    assert np.isnan(moved.coordinates_m[0, 2])
    assert not np.allclose(moved.coordinates_m[1:, 2], survey.coordinates_m[1:, 2])
    assert not np.allclose(moved.coordinates_m[:, :2], survey.coordinates_m[:, :2])


# --------------------------------------------------------------------------- #
# apply_corruption_suite
# --------------------------------------------------------------------------- #


def test_clean_severity_is_exact_passthrough():
    survey = _survey()
    out, record = apply_corruption_suite(survey, severity="clean", seed=0)
    np.testing.assert_array_equal(out.impedance, survey.impedance)
    assert out.n_valid == survey.n_valid
    assert record.severity == "clean"


def test_same_seed_is_exactly_reproducible():
    survey = _survey()
    out1, record1 = apply_corruption_suite(survey, severity="severe", seed=123)
    out2, record2 = apply_corruption_suite(survey, severity="severe", seed=123)
    np.testing.assert_array_equal(out1.impedance, out2.impedance)
    np.testing.assert_array_equal(out1.valid, out2.valid)
    assert record1.to_dict() == record2.to_dict()


def test_different_seeds_diverge():
    survey = _survey()
    out1, _ = apply_corruption_suite(survey, severity="severe", seed=1)
    out2, _ = apply_corruption_suite(survey, severity="severe", seed=2)
    assert not np.array_equal(out1.impedance, out2.impedance)


def test_severity_presets_increase_corruption_monotonically():
    survey = _survey(n_station=20, n_frequency=20)
    coverages = {}
    for severity in ("clean", "in_distribution", "severe", "held_out_corruption"):
        out, _ = apply_corruption_suite(survey, severity=severity, seed=0)
        coverages[severity] = out.coverage().overall
    order = ["clean", "in_distribution", "severe", "held_out_corruption"]
    values = [coverages[name] for name in order]
    assert values == sorted(values, reverse=True)


def test_config_and_severity_are_mutually_exclusive():
    survey = _survey()
    with pytest.raises(ValueError, match="exactly one"):
        apply_corruption_suite(survey, seed=0)
    with pytest.raises(ValueError, match="exactly one"):
        apply_corruption_suite(survey, CorruptionConfig(), severity="clean", seed=0)


def test_unknown_severity_raises():
    survey = _survey()
    with pytest.raises(ValueError, match="unknown severity"):
        apply_corruption_suite(survey, severity="nonexistent", seed=0)


def test_masks_and_errors_stay_consistent_after_full_suite():
    survey = _survey(n_station=10, n_frequency=12, error=np.ones((10, 12, 2)))
    out, record = apply_corruption_suite(survey, severity="severe", seed=5)
    # SurveyData's own invariants guarantee finiteness/positivity; re-validating
    # via replace() would raise if they were violated, so simply checking the
    # object exists and reporting on it is a meaningful regression guard.
    assert out.valid.shape == out.impedance.shape
    if out.impedance_error is not None:
        assert np.all(out.impedance_error[out.valid] > 0)
    assert isinstance(record.to_dict()["config_hash"], str)
