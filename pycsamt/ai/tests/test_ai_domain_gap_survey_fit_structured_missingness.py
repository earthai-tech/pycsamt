# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage for the StructuredMissingness half of survey_fit.py.

test_ai_domain_gap_survey_fit.py (the real WILLY-line contract tests)
covers ``survey_data_from_sites``/``fit_corruption_config``/
``fit_distortion_priors_from_sites``. This file fills the remaining gap:
``_sum_zero_contrast``, ``StructuredMissingness`` (``dropout_probability``,
``to_dict``/``from_dict`` round trip, schema-version guard), and
``fit_structured_missingness`` itself (the L-BFGS-B fit, its sum-to-zero
identifiability constraint, and its ``l2`` validation).
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.data.contracts import SurveyData
from pycsamt.ai.domain_gap.survey_fit import (
    StructuredMissingness,
    _project_local_metres,
    _sum_zero_contrast,
    fit_structured_missingness,
)


# ---------------------------------------------------------------------------
# _project_local_metres
# ---------------------------------------------------------------------------


def test_project_local_metres_uses_equirectangular_projection_for_real_coords():
    lats = np.array([10.0, 10.001, 10.002])
    lons = np.array([20.0, 20.001, 20.002])
    xy = _project_local_metres(lats, lons, station_spacing=500.0)
    assert xy.shape == (3, 2)
    assert np.ptp(xy[:, 0]) > 10.0 or np.ptp(xy[:, 1]) > 10.0


def test_project_local_metres_falls_back_to_uniform_grid_when_coords_collapse():
    lats = np.full(4, np.nan)
    lons = np.full(4, np.nan)
    xy = _project_local_metres(lats, lons, station_spacing=250.0)
    assert xy.shape == (4, 2)
    # 2x2 uniform grid at 250 m spacing
    assert set(xy[:, 0]) == {0.0, 250.0}
    assert set(xy[:, 1]) == {0.0, 250.0}


def test_project_local_metres_fallback_when_all_below_threshold():
    lats = np.array([10.0, 10.0, 10.0])
    lons = np.array([20.0, 20.0, 20.0])
    xy = _project_local_metres(lats, lons, station_spacing=100.0)
    assert xy.shape == (3, 2)


# ---------------------------------------------------------------------------
# _sum_zero_contrast
# ---------------------------------------------------------------------------


def test_sum_zero_contrast_shape_and_sum_to_zero_property():
    C = _sum_zero_contrast(4)
    assert C.shape == (4, 3)
    theta = np.array([1.0, -2.0, 0.5])
    effect = C @ theta
    assert effect.sum() == pytest.approx(0.0, abs=1e-12)


def test_sum_zero_contrast_n_equal_1_has_no_free_params():
    C = _sum_zero_contrast(1)
    assert C.shape == (1, 0)


def test_sum_zero_contrast_rejects_non_positive_n():
    with pytest.raises(ValueError, match="n must be positive"):
        _sum_zero_contrast(0)


# ---------------------------------------------------------------------------
# StructuredMissingness
# ---------------------------------------------------------------------------


def _toy_model():
    return StructuredMissingness(
        alpha=0.5,
        station_effect=(0.2, -0.2),
        frequency_effect=(0.1, 0.0, -0.1),
        component_effect=(0.3, -0.3),
        station_names=("A", "B"),
        frequencies_hz=(10.0, 5.0, 1.0),
        components=("xy", "yx"),
        l2=2.0,
    )


def test_dropout_probability_shape_and_range():
    model = _toy_model()
    p = model.dropout_probability()
    assert p.shape == (2, 3, 2)
    assert np.all((p >= 0.0) & (p <= 1.0))


def test_dropout_probability_matches_manual_sigmoid():
    model = StructuredMissingness(
        alpha=0.0,
        station_effect=(0.0,),
        frequency_effect=(0.0,),
        component_effect=(0.0,),
        station_names=("A",),
        frequencies_hz=(1.0,),
        components=("xy",),
    )
    np.testing.assert_allclose(model.dropout_probability(), [[[0.5]]])


def test_to_dict_and_from_dict_round_trip():
    model = _toy_model()
    d = model.to_dict()
    assert d["schema_version"] == 1
    restored = StructuredMissingness.from_dict(d)
    assert restored == model


def test_from_dict_rejects_unsupported_schema_version():
    d = _toy_model().to_dict()
    d["schema_version"] = 2
    with pytest.raises(ValueError, match="unsupported"):
        StructuredMissingness.from_dict(d)


def test_to_dict_returns_plain_python_types():
    model = _toy_model()
    d = model.to_dict()
    assert isinstance(d["station_effect"], list)
    assert isinstance(d["alpha"], float)


# ---------------------------------------------------------------------------
# fit_structured_missingness
# ---------------------------------------------------------------------------


def _survey_with_structured_missingness(seed=0):
    rng = np.random.default_rng(seed)
    n_station, n_freq, n_comp = 6, 8, 1
    z = (
        rng.normal(size=(n_station, n_freq, n_comp))
        + 1j * rng.normal(size=(n_station, n_freq, n_comp))
    )
    # Station 0 is missing at every frequency; frequency 0 is missing at
    # every station -- a real structured (not random) missingness pattern.
    z[0, :, :] = np.nan
    z[:, 0, :] = np.nan
    return SurveyData(
        z,
        np.linspace(100.0, 1.0, n_freq),
        [f"S{i}" for i in range(n_station)],
        ["xy"],
        np.zeros((n_station, 2)),
    )


def test_fit_structured_missingness_returns_correct_shape_and_identifiability():
    survey = _survey_with_structured_missingness()
    model = fit_structured_missingness(survey, max_iter=200)

    assert isinstance(model, StructuredMissingness)
    assert model.dropout_probability().shape == survey.shape
    assert abs(sum(model.station_effect)) < 1e-6
    assert abs(sum(model.frequency_effect)) < 1e-6
    assert abs(sum(model.component_effect)) < 1e-6


def test_fit_structured_missingness_station_0_has_high_dropout():
    survey = _survey_with_structured_missingness()
    model = fit_structured_missingness(survey, max_iter=300, l2=0.1)
    p = model.dropout_probability()
    # station 0 is missing everywhere -- its fitted probability should be
    # markedly higher than a station with full coverage.
    assert p[0].mean() > p[2:].mean()


def test_fit_structured_missingness_rejects_negative_l2():
    survey = _survey_with_structured_missingness()
    with pytest.raises(ValueError, match="l2 must be finite"):
        fit_structured_missingness(survey, l2=-1.0)


def test_fit_structured_missingness_rejects_non_finite_l2():
    survey = _survey_with_structured_missingness()
    with pytest.raises(ValueError, match="l2 must be finite"):
        fit_structured_missingness(survey, l2=float("nan"))


def test_fit_structured_missingness_provenance_matches_survey():
    survey = _survey_with_structured_missingness()
    model = fit_structured_missingness(survey, max_iter=100)
    assert model.station_names == tuple(survey.station_names)
    assert model.components == tuple(survey.components)
    assert len(model.frequencies_hz) == survey.n_frequencies
