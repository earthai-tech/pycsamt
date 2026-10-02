# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.settings."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestStatPill:
    def test_contains_value_and_label(self):
        from pycsamt.app.mapview.callbacks.settings import _stat_pill

        pill = _stat_pill("42", "visible")
        text = str(pill)
        assert "42" in text
        assert "visible" in text
        assert "mv-set-stat" in text


class TestRegisterSettings:
    def test_register_settings_is_callable(self):
        from pycsamt.app.mapview.callbacks.settings import register_settings

        assert callable(register_settings)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.CANVAS_SETTINGS in cb_outputs
        assert IDs.STORE_MASKED in cb_outputs
        assert IDs.SET_SUMMARY in cb_outputs
