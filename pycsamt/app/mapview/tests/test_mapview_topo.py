# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.topo."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestMsg:
    def test_ok_message_uses_check_icon(self):
        from pycsamt.app.mapview.callbacks.topo import _msg

        span = _msg("All good", ok=True)
        text = str(span)
        assert "bi-check-circle" in text
        assert "All good" in text

    def test_error_message_uses_warning_icon(self):
        from pycsamt.app.mapview.callbacks.topo import _msg

        span = _msg("Something broke", ok=False)
        text = str(span)
        assert "bi-exclamation-triangle" in text
        assert "Something broke" in text


class TestRegisterTopo:
    def test_register_topo_is_callable(self):
        from pycsamt.app.mapview.callbacks.topo import register_topo

        assert callable(register_topo)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.TOPO_UPLOAD_STORE in cb_outputs
        assert IDs.STORE_DATA in cb_outputs
        assert IDs.TOPO_EXPORT_DL in cb_outputs
