"""Tests for method-aware AMT edge QC: the ``method`` argument to
:func:`~pycsamt.iot.edge_amt.amt_edge_report` wires the per-method
``powerline_sensitive`` flag and target bands into the diagnostics.
"""

from __future__ import annotations

import numpy as np

from pycsamt.iot import amt_edge_report, amt_edge_table

FS = 2048.0


def _contaminated_signal(mains_hz: float = 50.0, seed: int = 0) -> np.ndarray:
    t = np.arange(0, 4, 1.0 / FS)
    rng = np.random.default_rng(seed)
    return 0.5 * np.sin(2 * np.pi * mains_hz * t) + rng.standard_normal(
        t.size
    )


# ---------------------------------------------------------------------------
# backward compatibility
# ---------------------------------------------------------------------------
def test_no_method_preserves_default_behaviour():
    report = amt_edge_report(_contaminated_signal(), FS)
    assert report["method"] is None
    assert report["powerline_applicable"] is True
    assert report["powerline"] is not None
    # no target bands imposed -> coverage fraction is undefined
    assert np.isnan(report["frequency_coverage"]["coverage_fraction"])


# ---------------------------------------------------------------------------
# powerline_sensitive gating
# ---------------------------------------------------------------------------
def test_powerline_sensitive_method_runs_detection():
    report = amt_edge_report(_contaminated_signal(), FS, method="amt")
    assert report["powerline_applicable"] is True
    assert report["powerline"] is not None
    assert report["powerline"]["contaminated"] is True  # 50 Hz present


def test_non_powerline_sensitive_method_skips_detection():
    # TDEM is a time-domain method; powerline harmonics are not applicable
    report = amt_edge_report(_contaminated_signal(), FS, method="tdem")
    assert report["powerline_applicable"] is False
    assert report["powerline"] is None


def test_csem_and_mt_are_powerline_sensitive():
    for method in ("mt", "csamt", "csem"):
        report = amt_edge_report(_contaminated_signal(), FS, method=method)
        assert report["powerline_applicable"] is True


# ---------------------------------------------------------------------------
# unknown / unrecognised methods fall back to the default
# ---------------------------------------------------------------------------
def test_unknown_method_defaults_to_detection_on():
    report = amt_edge_report(_contaminated_signal(), FS, method="unknown")
    assert report["powerline_applicable"] is True


def test_unrecognised_method_does_not_raise():
    report = amt_edge_report(_contaminated_signal(), FS, method="seismic")
    assert report["powerline_applicable"] is True  # safe default
    assert report["powerline"] is not None


# ---------------------------------------------------------------------------
# target bands feed coverage scoring
# ---------------------------------------------------------------------------
def test_method_populates_coverage_fraction():
    # the tone resolves above the noise floor, so coverage is scored
    # against the AMT target bands (a real number, not the undefined nan)
    report = amt_edge_report(_contaminated_signal(), FS, method="amt")
    frac = report["frequency_coverage"]["coverage_fraction"]
    assert frac is not None
    assert not np.isnan(frac)  # target bands were applied
    assert 0.0 <= frac <= 1.0

    # without a method, no target bands -> coverage fraction stays undefined
    plain = amt_edge_report(_contaminated_signal(), FS)
    assert np.isnan(plain["frequency_coverage"]["coverage_fraction"])


# ---------------------------------------------------------------------------
# table tolerates mixed powerline-applicable rows
# ---------------------------------------------------------------------------
def test_table_handles_skipped_powerline():
    amt = amt_edge_report(_contaminated_signal(), FS, method="amt")
    tdem = amt_edge_report(_contaminated_signal(), FS, method="tdem")
    table = amt_edge_table({"ex": amt, "hz": tdem})
    assert "powerline_applicable" in table.columns
    assert "method" in table.columns
    rows = {r["channel"]: r for r in table.to_dict("records")}
    assert rows["ex"]["powerline_applicable"] is True
    assert rows["hz"]["powerline_applicable"] is False
    # a skipped-powerline row carries no contamination verdict
    assert rows["hz"]["powerline_contaminated"] is None
