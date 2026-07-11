"""Tests for CSEM support: the ``EMMethod.CSEM`` profile and the
offset-domain (MVO/PVO) edge QC in :mod:`pycsamt.iot.edge_csem`.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.iot import (
    EMMethod,
    MonitoringConfig,
    OffsetResponse,
    csem_edge_report,
    csem_edge_table,
    field_vs_offset,
    method_profile,
    target_bands_for_method,
)

OFFSETS = np.array([1000.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0])


def _decaying_amplitude(offsets=OFFSETS, scale=3000.0, a0=1e-9):
    return a0 * np.exp(-offsets / scale)


# ---------------------------------------------------------------------------
# CSEM as a first-class method
# ---------------------------------------------------------------------------
def test_csem_is_registered_method():
    assert EMMethod.CSEM.value == "csem"
    profile = method_profile("csem")
    assert profile.controlled_source
    assert profile.frequency_band_hz == (0.01, 100.0)


def test_csem_target_bands_and_config():
    bands = target_bands_for_method("csem")
    assert bands[0][0] == pytest.approx(0.01)
    assert bands[-1][1] == pytest.approx(100.0)
    cfg = MonitoringConfig.for_method("csem")
    assert cfg.method is EMMethod.CSEM
    assert cfg.required_channels == ["ex", "ey", "hx", "hy"]


# ---------------------------------------------------------------------------
# field_vs_offset (MVO / PVO)
# ---------------------------------------------------------------------------
def test_field_vs_offset_monotonic_decay():
    resp = field_vs_offset(OFFSETS, _decaying_amplitude(), frequency_hz=1.0)
    assert isinstance(resp, OffsetResponse)
    assert resp.n_offsets == OFFSETS.size
    assert resp.monotonic_decay
    assert resp.dynamic_range_db > 0
    assert resp.frequency_hz == pytest.approx(1.0)


def test_field_vs_offset_detects_bump():
    amp = _decaying_amplitude().copy()
    amp[3] = amp[2] * 1.5              # a receiver reading that increases
    resp = field_vs_offset(OFFSETS, amp)
    assert not resp.monotonic_decay


def test_field_vs_offset_noise_floor_limits_detectability():
    # amplitude falls below the floor at the far offsets
    amp = _decaying_amplitude(scale=1500.0)
    floor = 1e-11
    resp = field_vs_offset(OFFSETS, amp, noise_floor=floor)
    assert resp.n_detectable < resp.n_offsets
    assert resp.max_detectable_offset_m is not None
    # every offset beyond the detectability limit is below the floor
    for off, ok in zip(resp.offsets_m, resp.above_noise):
        if off > resp.max_detectable_offset_m:
            assert not ok


def test_field_vs_offset_sorts_by_offset():
    resp = field_vs_offset([4000.0, 1000.0, 2000.0], [1.0, 3.0, 2.0])
    assert resp.offsets_m == [1000.0, 2000.0, 4000.0]
    assert resp.amplitudes == [3.0, 2.0, 1.0]


def test_field_vs_offset_carries_phase():
    phases = -OFFSETS / 1000.0 * 30.0
    resp = field_vs_offset(OFFSETS, _decaying_amplitude(), phases_deg=phases)
    assert resp.phases_deg is not None
    assert len(resp.phases_deg) == OFFSETS.size


def test_field_vs_offset_drops_nonfinite():
    offsets = [1000.0, 2000.0, np.nan, 4000.0]
    amps = [3.0, 2.0, 1.0, 0.5]
    resp = field_vs_offset(offsets, amps)
    assert resp.n_offsets == 3            # NaN offset dropped


def test_field_vs_offset_single_point():
    resp = field_vs_offset([1000.0], [1e-9])
    assert resp.monotonic_decay
    assert resp.dynamic_range_db == 0.0


def test_field_vs_offset_empty():
    resp = field_vs_offset([], [])
    assert resp.n_offsets == 0
    assert resp.max_detectable_offset_m is None


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------
def test_field_vs_offset_length_mismatch():
    with pytest.raises(ValueError):
        field_vs_offset([1.0, 2.0, 3.0], [1.0, 2.0])


def test_field_vs_offset_rejects_negative_amplitude():
    with pytest.raises(ValueError):
        field_vs_offset([1.0, 2.0], [1.0, -2.0])


def test_field_vs_offset_phase_length_mismatch():
    with pytest.raises(ValueError):
        field_vs_offset([1.0, 2.0], [1.0, 2.0], phases_deg=[0.0])


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------
def test_csem_edge_report_and_table():
    report = csem_edge_report(
        OFFSETS, _decaying_amplitude(), noise_floor=1e-13, frequency_hz=1.0
    )
    assert "offset_response" in report
    table = csem_edge_table({"1.0Hz": report})
    assert table.shape[0] == 1
    for col in ("frequency_hz", "n_detectable", "max_detectable_offset_m",
                "monotonic_decay", "dynamic_range_db"):
        assert col in table.columns
