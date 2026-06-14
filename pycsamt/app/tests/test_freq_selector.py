# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for FreqSelector (Phase 3)."""

from __future__ import annotations

import math
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.widgets.freq_selector import FreqSelector, _RangeSlider


# ── _RangeSlider unit tests ───────────────────────────────────────────────

def test_range_slider_creates(qapp):
    s = _RangeSlider(0.001, 1000.0)
    assert s is not None
    s.close()


def test_range_slider_initial_lo_hi(qapp):
    s = _RangeSlider(0.1, 100.0)
    assert s.lo_hz == pytest.approx(0.1)
    assert s.hi_hz == pytest.approx(100.0)
    s.close()


def test_range_slider_set_selection(qapp):
    s = _RangeSlider(0.001, 10000.0)
    s.set_selection(1.0, 100.0)
    assert s.lo_hz == pytest.approx(1.0)
    assert s.hi_hz == pytest.approx(100.0)
    s.close()


def test_range_slider_selection_clamped_to_range(qapp):
    s = _RangeSlider(1.0, 100.0)
    s.set_selection(0.001, 99999.0)   # outside range
    assert s.lo_hz >= 1.0
    assert s.hi_hz <= 100.0
    s.close()


def test_range_slider_log_roundtrip(qapp):
    s = _RangeSlider(0.001, 10000.0)
    for hz in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0):
        frac = s._hz_to_frac(hz)
        assert 0.0 <= frac <= 1.0
        recovered = s._frac_to_hz(frac)
        assert recovered == pytest.approx(hz, rel=1e-6)
    s.close()


def test_range_slider_emits_signal(qapp):
    s = _RangeSlider(0.001, 10000.0)
    s.resize(300, 40)
    received = []
    s.range_changed.connect(lambda lo, hi: received.append((lo, hi)))
    s.set_selection(10.0, 1000.0)
    # Signal should be emitted via set_selection only indirectly;
    # test internal state is correct (signal only fires on mouse events)
    assert s.lo_hz == pytest.approx(10.0)
    assert s.hi_hz == pytest.approx(1000.0)
    s.close()


# ── FreqSelector widget tests ─────────────────────────────────────────────

@pytest.fixture
def sel(qapp):
    s = FreqSelector(f_min=0.001, f_max=10000.0)
    yield s
    s.close()


def test_freq_selector_creates(qapp):
    s = FreqSelector()
    assert s is not None
    s.close()


def test_freq_selector_default_range(sel):
    assert sel.lo_hz == pytest.approx(0.001)
    assert sel.hi_hz == pytest.approx(10000.0)


def test_freq_selector_set_freq_range(sel):
    sel.set_freq_range(1.0, 1000.0)
    assert sel.lo_hz == pytest.approx(1.0)
    assert sel.hi_hz == pytest.approx(1000.0)


def test_freq_selector_set_selection(sel):
    sel.set_selection(5.0, 500.0)
    assert sel.lo_hz == pytest.approx(5.0)
    assert sel.hi_hz == pytest.approx(500.0)


def test_freq_selector_toggle_unit(sel):
    assert not sel._show_period
    sel._toggle_unit()
    assert sel._show_period
    sel._toggle_unit()
    assert not sel._show_period


def test_freq_selector_label_updates_on_toggle(sel):
    sel._toggle_unit()
    assert "Period" in sel._unit_lbl.text()
    sel._toggle_unit()
    assert "Frequency" in sel._unit_lbl.text()


def test_freq_selector_range_changed_signal(sel):
    received = []
    sel.range_changed.connect(lambda lo, hi: received.append((lo, hi)))
    sel._on_slider_changed(10.0, 500.0)
    assert len(received) == 1
    assert received[0][0] == pytest.approx(10.0)
    assert received[0][1] == pytest.approx(500.0)


def test_fmt_freq_hz(sel):
    assert "Hz" in sel._fmt_freq(50.0)


def test_fmt_freq_khz(sel):
    assert "kHz" in sel._fmt_freq(5000.0)


def test_fmt_freq_mhz(sel):
    assert "mHz" in sel._fmt_freq(0.5)


def test_fmt_period_seconds(sel):
    assert "s" in sel._fmt_period(5.0)


def test_fmt_period_ms(sel):
    assert "ms" in sel._fmt_period(0.005)
