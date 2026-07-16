# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.chrome — seed/theme/help/dock."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestRegisterChrome:
    def test_register_chrome_is_callable(self):
        from pycsamt.app.mapview.callbacks.chrome import register_chrome

        assert callable(register_chrome)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_DATA in cb_outputs
        assert IDs.STORE_THEME in cb_outputs
        assert IDs.MODAL_HELP in cb_outputs

    def test_resize_snippet_dispatches_resize_event(self):
        from pycsamt.app.mapview.callbacks.chrome import _RESIZE

        assert "dispatchEvent" in _RESIZE
        assert "resize" in _RESIZE
