"""Contracts for :mod:`pycsamt.ai.domain_gap.survey_fit`, exercised against
the real, tracked WILLY L18PLT line (28 stations, 53 shared frequencies).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.ai.domain_gap.simulator import CorruptionConfig
from pycsamt.ai.domain_gap.survey_fit import (
    fit_corruption_config,
    fit_distortion_priors_from_sites,
    survey_data_from_sites,
)

_L18 = Path("data/AMT/WILLY_DATA/L18PLT")

pytestmark = pytest.mark.skipif(
    not _L18.exists(), reason="WILLY_DATA/L18PLT is not available in this checkout."
)


@pytest.fixture(scope="module")
def willy_survey():
    return survey_data_from_sites(str(_L18), recursive=False, verbose=0)


def test_bridge_reads_the_full_l18_line_with_full_tensor_and_errors(willy_survey):
    assert willy_survey.n_stations == 28
    assert willy_survey.components == ("xx", "xy", "yx", "yy")
    assert willy_survey.n_frequencies == 53
    assert willy_survey.impedance_error is not None
    assert willy_survey.coverage().overall == pytest.approx(1.0)


def test_bridge_projects_coordinates_to_finite_local_metres(willy_survey):
    coordinates = willy_survey.coordinates_m
    assert np.all(np.isfinite(coordinates[:, :2]))
    # elevation was reported for every WILLY station
    assert np.all(np.isfinite(coordinates[:, 2]))
    # stations must not all collapse onto a single point
    assert np.ptp(coordinates[:, 0]) > 10.0 or np.ptp(coordinates[:, 1]) > 10.0


class _FakeSite:
    def __init__(self, name, freq, z, z_err=None, coords=(0.0, 0.0, 0.0)):
        self.name = name
        self.freq = freq
        self.z = z
        self.z_err = z_err
        self.coords = coords


def test_bridge_rejects_mismatched_frequency_grids(monkeypatch):
    # survey_fit imports ensure_sites at module level (`from ...emtools
    # ._core import ensure_sites`), so the name to patch is the one
    # bound in survey_fit's own namespace, not emtools._core's.
    import pycsamt.ai.domain_gap.survey_fit as survey_fit

    a = _FakeSite("A", np.array([10.0, 1.0]), np.ones((2, 2, 2), dtype=complex))
    b = _FakeSite("B", np.array([10.0, 2.0]), np.ones((2, 2, 2), dtype=complex))
    monkeypatch.setattr(
        survey_fit, "ensure_sites", lambda *args, **kwargs: [a, b]
    )
    with pytest.raises(ValueError, match="does not share the frequency grid"):
        survey_data_from_sites([a, b])


def test_bridge_reads_same_line_twice_identically():
    survey_a = survey_data_from_sites(str(_L18), recursive=False, verbose=0)
    survey_b = survey_data_from_sites(str(_L18), recursive=False, verbose=0)
    np.testing.assert_array_equal(survey_a.frequencies_hz, survey_b.frequencies_hz)


def test_fit_corruption_config_from_real_errors_yields_plausible_ranges(willy_survey):
    config = fit_corruption_config(willy_survey)
    assert isinstance(config, CorruptionConfig)
    lo, hi = config.noise_level_range
    assert 0.0 < lo <= hi < 1.0
    assert 0.0 < config.error_floor_fraction < 1.0
    # L18PLT has complete coverage, so there is nothing to learn about dropout
    assert config.station_dropout_rate == 0.0
    assert config.frequency_dropout_rate == 0.0


def test_fit_corruption_config_severity_scale_scales_ranges(willy_survey):
    base = fit_corruption_config(willy_survey, severity_scale=1.0)
    scaled = fit_corruption_config(willy_survey, severity_scale=2.0)
    assert scaled.noise_level_range[1] == pytest.approx(base.noise_level_range[1] * 2.0)


def test_fit_corruption_config_rejects_survey_without_errors():
    from pycsamt.ai.data.contracts import SurveyData

    z = np.full((2, 3, 1), 1 + 1j)
    survey = SurveyData(z, [3.0, 2.0, 1.0], ["A", "B"], ["xy"], np.zeros((2, 2)))
    with pytest.raises(ValueError, match="impedance_error"):
        fit_corruption_config(survey)


def test_fit_distortion_priors_from_real_willy_line_are_finite_and_nonnegative():
    priors = fit_distortion_priors_from_sites(str(_L18), recursive=False, verbose=0)
    expected_keys = {
        "static_shift_log10_sigma",
        "distortion_gain_log10_sigma",
        "distortion_twist_deg_sigma",
        "distortion_shear_sigma",
        "distortion_anisotropy_sigma",
    }
    assert set(priors) == expected_keys
    for value in priors.values():
        assert np.isfinite(value)
        assert value >= 0.0
    # A real, non-degenerate survey line should show some twist spread.
    assert priors["distortion_twist_deg_sigma"] > 0.0
