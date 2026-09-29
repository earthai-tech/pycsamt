# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.view — rail switch + render."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestRailMapping:
    def test_rail_maps_to_view_names(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.callbacks.view import _RAIL

        assert _RAIL[IDs.RAIL_MAP] == "map"
        assert _RAIL[IDs.RAIL_3D] == "map3d"


class TestRegisterView:
    def test_register_view_is_callable(self):
        from pycsamt.app.mapview.callbacks.view import register_view

        assert callable(register_view)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_VIEW in cb_outputs
        assert IDs.CANVAS_GRAPH in cb_outputs
