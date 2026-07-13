"""Tests for the static-shift indicator (:func:`estimate_static_shift`) and
the controlled-source transmitter timing-lock fields on ``SyncPayload``.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.iot import (
    StaticShift,
    estimate_static_shift,
    parse_payload,
    validate_payload,
)

FREQ = np.logspace(3, 0, 12)
BASE = 100.0 * np.ones_like(FREQ)
PHASE = 45.0 * np.ones_like(FREQ)


# ---------------------------------------------------------------------------
# static shift
# ---------------------------------------------------------------------------
def test_no_static_shift_for_equal_modes():
    res = estimate_static_shift(BASE, BASE, phase_xy=PHASE, phase_yx=PHASE)
    assert isinstance(res, StaticShift)
    assert res.shift_factor == pytest.approx(1.0)
    assert not res.static_shift


def test_static_shift_detected():
    res = estimate_static_shift(
        BASE * 3.0, BASE, phase_xy=PHASE, phase_yx=PHASE
    )
    assert res.static_shift
    assert res.shift_factor == pytest.approx(3.0, rel=1e-6)
    assert res.split_decades == pytest.approx(np.log10(3.0), rel=1e-6)
    assert res.consistency_std == pytest.approx(0.0, abs=1e-9)


def test_frequency_dependent_split_is_not_static():
    # a split that grows strongly with frequency is structure, not shift
    varying = BASE * np.linspace(1.0, 10.0, FREQ.size)
    res = estimate_static_shift(varying, BASE)
    assert not res.static_shift
    assert res.consistency_std > 0.15


def test_phase_disagreement_vetoes_static_shift():
    # a persistent resistivity split but disagreeing phases -> anisotropy
    res = estimate_static_shift(
        BASE * 3.0, BASE, phase_xy=PHASE + 20.0, phase_yx=PHASE
    )
    assert res.phase_diff_deg == pytest.approx(20.0, rel=1e-6)
    assert not res.static_shift


def test_static_shift_without_phases():
    res = estimate_static_shift(BASE * 2.0, BASE)
    assert res.static_shift  # phases optional; split is enough
    assert np.isnan(res.phase_diff_deg)


def test_static_shift_all_invalid():
    res = estimate_static_shift([np.nan, -1.0], [0.0, np.nan])
    assert not res.static_shift
    assert np.isnan(res.shift_factor)


def test_static_shift_as_dict():
    d = estimate_static_shift(BASE * 3.0, BASE).as_dict()
    assert set(d) == {
        "shift_factor",
        "split_decades",
        "consistency_std",
        "phase_diff_deg",
        "static_shift",
    }


# ---------------------------------------------------------------------------
# transmitter timing lock on SyncPayload
# ---------------------------------------------------------------------------
def test_sync_payload_tx_lock_fields():
    p = parse_payload(
        "sync",
        {
            "offset_ms": 1.2,
            "transmitter_locked": "yes",
            "tx_offset_ms": 0.3,
            "source_id": "TX1",
        },
    )
    assert p.tx_locked is True
    assert p.tx_sync_offset_ms == pytest.approx(0.3)
    assert p.tx_id == "TX1"


def test_sync_payload_tx_lock_canonical():
    canon = validate_payload("sync", {"tx_locked": True, "tx_id": "TX1"})
    assert canon["tx_locked"] is True
    assert canon["tx_id"] == "TX1"


def test_sync_payload_without_tx_fields_unaffected():
    # legacy sync packets with no transmitter fields still parse
    p = parse_payload("sync", {"offset_ms": 0.5, "gps_lock": True})
    assert p.tx_locked is None
    assert p.tx_sync_offset_ms is None
    assert p.tx_id is None
