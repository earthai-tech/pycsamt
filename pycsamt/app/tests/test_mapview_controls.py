# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.controls."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestFmtFreq:
    def test_none_returns_dash(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        assert _fmt_freq(None) == "—"

    def test_khz_range(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        assert _fmt_freq(1500) == "1.5 kHz"

    def test_hz_range(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        assert _fmt_freq(50) == "50 Hz"

    def test_sub_hz_range_includes_period(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        result = _fmt_freq(0.1)
        assert "Hz" in result and "s)" in result


class TestPresets:
    def test_depth_presets_cover_full_and_bands(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.callbacks.controls import _DEPTH_PRESETS

        assert _DEPTH_PRESETS[IDs.BTN_DEPTH_FULL] == (None, None)
        assert _DEPTH_PRESETS[IDs.BTN_DEPTH_500] == (0, 500)
        assert _DEPTH_PRESETS[IDs.TB3D_DEPTH_2K] == (0, 2000)

    def test_rho_presets_are_ordered_bands(self):
        from pycsamt.app.mapview.callbacks.controls import _RHO_PRESETS

        for lo, hi in _RHO_PRESETS.values():
            assert lo < hi


class TestRegisterControls:
    def test_register_controls_is_callable(self):
        from pycsamt.app.mapview.callbacks.controls import register_controls

        assert callable(register_controls)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_CONTROLS in cb_outputs
        assert IDs.CTL_DEPTH_LO in cb_outputs
        assert IDs.CTL_RHO_LO in cb_outputs
