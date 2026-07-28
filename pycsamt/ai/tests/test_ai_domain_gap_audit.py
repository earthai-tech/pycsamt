"""Contracts for :mod:`pycsamt.ai.domain_gap.audit`, exercised against the
real, tracked WILLY L18PLT line (28 stations, 53 shared frequencies) and
against small synthetic fixtures for the exclusion/mismatch paths that a
complete, non-degenerate field survey cannot exhibit.
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pytest

from pycsamt.ai.domain_gap import audit as audit_module
from pycsamt.ai.domain_gap.audit import (
    DimensionalitySummary,
    FrequencyGridReport,
    StationExclusion,
    SurveyAuditReport,
    audit_survey,
)

_L18 = Path("data/AMT/WILLY_data/L18PLT")

pytestmark = pytest.mark.skipif(
    not _L18.exists(), reason="WILLY_data/L18PLT is not available in this checkout."
)


@pytest.fixture(scope="module")
def willy_report() -> SurveyAuditReport:
    return audit_survey(str(_L18), recursive=False, verbose=0)


def test_real_l18_reports_full_accounting_with_no_exclusions(willy_report):
    assert willy_report.n_stations_input == 28
    assert willy_report.n_stations_included == 28
    assert willy_report.excluded_stations == ()
    assert willy_report.frequency_grid.matched
    assert willy_report.coverage is not None
    assert willy_report.coverage.overall == pytest.approx(1.0)
    assert willy_report.frequency_range_hz == pytest.approx((1.008, 10400.0))


def test_real_l18_error_ratio_percentiles_are_ordered_and_finite(willy_report):
    assert willy_report.error_ratio_p05 is not None
    assert (
        0.0
        < willy_report.error_ratio_p05
        <= willy_report.error_ratio_p50
        <= willy_report.error_ratio_p95
    )


def test_real_l18_station_spacing_is_positive_and_ordered(willy_report):
    spacing = willy_report.station_spacing_m
    assert spacing is not None
    assert 0.0 < spacing["min"] <= spacing["median"] <= spacing["max"]


def test_real_l18_elevation_is_fully_reported(willy_report):
    assert willy_report.elevation_coverage == pytest.approx(1.0)


def test_real_l18_dimensionality_fractions_are_a_partition(willy_report):
    dim = willy_report.dimensionality
    assert dim.n_samples == 28 * 53
    assert dim.frac_1d + dim.frac_2d + dim.frac_3d == pytest.approx(1.0, abs=1e-6)
    assert dim.strike_consensus_deg is not None
    assert dim.strike_consensus_iqr_deg is not None
    # this is a real, complex survey line: not every station is 1-D.
    assert dim.frac_3d > 0.0


def test_real_l18_distortion_indicators_show_real_spread(willy_report):
    assert willy_report.static_shift_log10_sigma > 0.0
    assert willy_report.distortion_twist_deg_sigma > 0.0


def test_real_l18_crs_is_reported_as_not_declared(willy_report):
    # coordinates are locally projected from lat/lon; no formal CRS is
    # ever inferred by this bridge.
    assert willy_report.crs_declared is False


def test_summary_is_nonempty_text_containing_station_counts(willy_report):
    text = willy_report.summary()
    assert "28 input" in text
    assert "28 included" in text


def test_json_round_trip_preserves_every_field(willy_report):
    with TemporaryDirectory() as directory:
        path = willy_report.write_json(Path(directory) / "audit.json")
        reloaded = SurveyAuditReport.read_json(path)
    assert reloaded.to_dict() == willy_report.to_dict()


class _FakeSite:
    def __init__(self, name, freq, z, z_err=None, coords=(0.0, 0.0, 0.0)):
        self.name = name
        self.freq = freq
        self.z = z
        self.z_err = z_err
        self.coords = coords


def _patch_ensure_sites(monkeypatch, stations):
    # audit.py imports ensure_sites at module level (`from ...emtools._core
    # import ensure_sites`), so the name to patch is the one bound in
    # audit's own namespace, matching the gotcha documented for
    # survey_fit.py's equivalent test.
    monkeypatch.setattr(
        audit_module, "ensure_sites", lambda *args, **kwargs: stations
    )


def test_stations_missing_freq_or_z_are_excluded_and_accounted_for(monkeypatch):
    good = _FakeSite(
        "GOOD", np.array([10.0, 1.0]), np.ones((2, 2, 2), dtype=complex)
    )
    missing_freq = _FakeSite("NO-FREQ", None, np.ones((2, 2, 2), dtype=complex))
    missing_z = _FakeSite("NO-Z", np.array([10.0, 1.0]), None)
    _patch_ensure_sites(monkeypatch, [good, missing_freq, missing_z])

    # Only "GOOD" survives exclusion, so audit_survey proceeds to the
    # canonical bridge; that bridge's own re-normalization through the
    # real (unpatched) ensure_sites cannot recognize a bare list of
    # these lightweight fakes, which is a pre-existing ensure_sites
    # limitation unrelated to what this test targets. Stub the bridge
    # itself so this test exercises audit_survey's own exclusion and
    # frequency-grid logic in isolation.
    from pycsamt.ai.data.contracts import SurveyData

    stub_survey = SurveyData(
        np.ones((1, good.freq.size, 4), dtype=complex),
        good.freq,
        [good.name],
        ["xx", "xy", "yx", "yy"],
        [[0.0, 0.0]],
    )
    monkeypatch.setattr(
        audit_module, "survey_data_from_sites", lambda *a, **k: stub_survey
    )

    report = audit_survey([good, missing_freq, missing_z])

    assert report.n_stations_input == 3
    assert report.n_stations_included == 1
    assert set(item.station for item in report.excluded_stations) == {
        "NO-FREQ",
        "NO-Z",
    }
    assert all(
        item.reason == "missing freq or z array" for item in report.excluded_stations
    )
    assert report.frequency_grid.matched
    assert report.frequency_grid.reference_station == "GOOD"


def test_all_stations_missing_data_raises(monkeypatch):
    missing_freq = _FakeSite("NO-FREQ", None, np.ones((2, 2, 2), dtype=complex))
    _patch_ensure_sites(monkeypatch, [missing_freq])

    with pytest.raises(ValueError, match="no station"):
        audit_survey([missing_freq])


def test_mismatched_frequency_grid_is_reported_not_raised(monkeypatch):
    a = _FakeSite("A", np.array([10.0, 1.0]), np.ones((2, 2, 2), dtype=complex))
    b = _FakeSite("B", np.array([10.0, 2.0]), np.ones((2, 2, 2), dtype=complex))
    _patch_ensure_sites(monkeypatch, [a, b])

    report = audit_survey([a, b])

    assert report.n_stations_input == 2
    assert report.n_stations_included == 2
    assert not report.frequency_grid.matched
    assert report.frequency_grid.mismatched_stations == ("B",)
    assert report.frequency_grid.n_frequencies_by_station == {"A": 2, "B": 2}
    # the canonical bridge cannot run without a shared grid.
    assert report.coverage is None
    assert report.frequency_range_hz is None
    assert report.error_ratio_p50 is None


def test_frequency_grid_report_rejects_inconsistent_matched_flag():
    with pytest.raises(ValueError, match="matched must be true exactly when"):
        FrequencyGridReport(True, "A", {"A": 1}, ("B",))


def test_survey_audit_report_rejects_excess_exclusions():
    with pytest.raises(ValueError, match="cannot exceed n_stations_input"):
        SurveyAuditReport(
            1,
            (
                StationExclusion("A", "missing freq or z array"),
                StationExclusion("B", "missing freq or z array"),
            ),
            FrequencyGridReport(True, "A", {"A": 1}, ()),
            None,
            None,
            None,
            None,
            None,
            None,
            0.0,
            False,
            DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            "2026-01-01T00:00:00Z",
        )


def test_dimensionality_summary_rejects_negative_n_samples():
    with pytest.raises(ValueError, match="n_samples cannot be negative"):
        DimensionalitySummary(-1, 0.0, 0.0, 0.0, None, None)
