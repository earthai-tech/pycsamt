"""Tests for :mod:`pycsamt.iot.edge_csamt` — CSAMT / controlled-source
edge quality control: skin depth, field-zone classification, transmitter
frequency-comb detection, and source-signal stability.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.iot import (
    DeviceRole,
    FieldZone,
    assess_source_stability,
    classify_field_zones,
    csamt_edge_report,
    csamt_edge_table,
    detect_transmitter_frequencies,
    skin_depth_m,
)


# ---------------------------------------------------------------------------
# skin depth
# ---------------------------------------------------------------------------
def test_skin_depth_reference_value():
    # 503.29 * sqrt(100 / 1) = 5032.9 m
    assert float(skin_depth_m(100.0, 1.0)) == pytest.approx(5032.9, rel=1e-3)


def test_skin_depth_scales_with_frequency():
    # doubling frequency reduces skin depth by sqrt(2)
    d1 = float(skin_depth_m(100.0, 1.0))
    d2 = float(skin_depth_m(100.0, 2.0))
    assert d1 / d2 == pytest.approx(np.sqrt(2.0), rel=1e-6)


def test_skin_depth_invalid_is_nan():
    assert np.isnan(float(skin_depth_m(100.0, 0.0)))
    assert np.isnan(float(skin_depth_m(-1.0, 10.0)))


def test_skin_depth_array_broadcast():
    out = skin_depth_m(100.0, np.array([1.0, 4.0, 16.0]))
    assert out.shape == (3,)
    assert out[0] > out[1] > out[2]


# ---------------------------------------------------------------------------
# field zones
# ---------------------------------------------------------------------------
def test_classify_field_zones_orders_high_to_low():
    freqs = np.array([4096.0, 256.0, 16.0, 1.0])
    cov = classify_field_zones(freqs, 100.0, offset_m=5000.0)
    # high frequencies (shallow skin depth) are far field; low ones are not
    zmap = dict(zip([int(f) for f in cov.freq_hz], cov.zones))
    assert zmap[4096] == FieldZone.FAR.value
    assert zmap[1] == FieldZone.NEAR.value
    assert cov.n_far + cov.n_transition + cov.n_near == 4


def test_classify_field_zones_summary_flags():
    freqs = np.array([4096.0, 256.0, 16.0, 4.0, 1.0])
    cov = classify_field_zones(freqs, 100.0, offset_m=5000.0)
    assert not cov.all_far_field
    assert cov.correction_recommended
    assert 0.0 <= cov.far_fraction <= 1.0
    assert cov.first_far_field_hz() is not None


def test_classify_field_zones_all_far_when_offset_large():
    freqs = np.array([4096.0, 1024.0, 256.0])
    cov = classify_field_zones(freqs, 100.0, offset_m=500000.0)
    assert cov.all_far_field
    assert not cov.correction_recommended
    assert cov.far_fraction == pytest.approx(1.0)


def test_classify_field_zones_per_frequency_resistivity():
    freqs = np.array([100.0, 10.0])
    rho = np.array([50.0, 500.0])
    cov = classify_field_zones(freqs, rho, offset_m=3000.0)
    assert len(cov.zones) == 2


def test_classify_field_zones_bad_resistivity_shape():
    with pytest.raises(ValueError):
        classify_field_zones(
            np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0]), offset_m=1000.0
        )


def test_classify_field_zones_bad_ratios():
    with pytest.raises(ValueError):
        classify_field_zones(
            np.array([1.0]),
            100.0,
            offset_m=1000.0,
            near_ratio=5.0,
            far_ratio=1.0,
        )


# ---------------------------------------------------------------------------
# transmitter comb
# ---------------------------------------------------------------------------
def _comb_signal(fs: float, tx_present, duration: float = 4.0, seed: int = 0):
    t = np.arange(0, duration, 1.0 / fs)
    sig = sum(np.sin(2 * np.pi * f * t) for f in tx_present)
    sig = sig + 0.01 * np.random.default_rng(seed).standard_normal(t.size)
    return sig


def test_detect_transmitter_frequencies_finds_present_lines():
    fs = 2048.0
    present = [8.0, 32.0, 128.0]
    expected = [8.0, 32.0, 128.0, 512.0]
    sig = _comb_signal(fs, present)
    comb = detect_transmitter_frequencies(sig, fs, expected)
    assert comb.n_expected == 4
    assert comb.n_detected == 3
    assert comb.missing() == [512.0]
    assert comb.detection_fraction == pytest.approx(0.75)


def test_detect_transmitter_frequencies_above_nyquist_undetected():
    fs = 256.0
    comb = detect_transmitter_frequencies(_comb_signal(fs, [16.0]), fs, [16.0, 200.0])
    line = next(ln for ln in comb.lines if ln.frequency_hz == 200.0)
    assert not line.detected  # above Nyquist (128 Hz)


def test_detect_transmitter_frequencies_empty_signal():
    comb = detect_transmitter_frequencies([], 1000.0, [10.0, 20.0])
    assert comb.n_detected == 0
    assert comb.n_expected == 2


# ---------------------------------------------------------------------------
# source stability
# ---------------------------------------------------------------------------
def test_assess_source_stability_steady_source():
    rng = np.random.default_rng(1)
    current = np.full(300, 10.0) + 0.02 * rng.standard_normal(300)
    st = assess_source_stability(current, tx_voltage=np.full(300, 250.0))
    assert st.stable
    assert st.flags == []
    assert st.current_mean_a == pytest.approx(10.0, abs=0.1)
    assert st.voltage_mean_v == pytest.approx(250.0, abs=0.1)


def test_assess_source_stability_flags_unstable():
    rng = np.random.default_rng(2)
    current = 10.0 + 5.0 * rng.standard_normal(300)
    st = assess_source_stability(current)
    assert not st.stable
    assert "current_unstable" in st.flags


def test_assess_source_stability_on_off_keying():
    # on/off-keyed source: on-state statistics ignore the off samples
    current = np.concatenate([np.full(100, 10.0), np.zeros(40), np.full(100, 10.0)])
    st = assess_source_stability(current)
    assert st.on_fraction == pytest.approx(200 / 240, abs=1e-6)
    assert st.stable  # on-state current is steady


def test_assess_source_stability_insufficient_samples():
    st = assess_source_stability([1.0])
    assert not st.stable
    assert "insufficient_samples" in st.flags


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------
def test_csamt_edge_report_and_table():
    fs = 2048.0
    tx = [8.0, 32.0, 128.0, 512.0]
    sig = _comb_signal(fs, [8.0, 32.0, 128.0])
    current = np.full(200, 10.0)
    report = csamt_edge_report(
        sig,
        fs,
        tx_frequencies=tx,
        offset_m=5000.0,
        resistivity=100.0,
        tx_current=current,
    )
    assert "transmitter" in report
    assert "field_zones" in report
    assert "source_stability" in report

    table = csamt_edge_table({"ex": report})
    assert table.shape[0] == 1
    for col in (
        "n_tx_detected",
        "far_fraction",
        "correction_recommended",
        "source_stable",
    ):
        assert col in table.columns


def test_csamt_edge_report_without_source():
    report = csamt_edge_report(
        _comb_signal(2048.0, [32.0]),
        2048.0,
        tx_frequencies=[32.0],
        offset_m=5000.0,
        resistivity=100.0,
    )
    assert "source_stability" not in report


# ---------------------------------------------------------------------------
# transmitter role
# ---------------------------------------------------------------------------
def test_transmitter_device_role():
    assert DeviceRole.TRANSMITTER.value == "transmitter"
    from pycsamt.iot import DeviceConfig

    tx = DeviceConfig(device_id="tx-1", role="transmitter")
    assert tx.role is DeviceRole.TRANSMITTER
