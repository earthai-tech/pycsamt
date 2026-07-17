# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.export."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestRegisterExport:
    def test_register_export_is_callable(self):
        from pycsamt.app.mapview.callbacks.export import register_export

        assert callable(register_export)

    def test_expected_output_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.EXPORT_DL in cb_outputs
        assert IDs.BTN_EXPORT in cb_outputs
