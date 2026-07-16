# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.toolbar."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestBasemapMapping:
    def test_all_basemap_buttons_mapped(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.callbacks.toolbar import _BM_BTN_STYLE

        assert _BM_BTN_STYLE[IDs.TB_BM_DARK] == "carto-darkmatter"
        assert _BM_BTN_STYLE[IDs.TB_BM_SAT] == "esri-satellite"
        assert len(_BM_BTN_STYLE) == 5


class TestMode3dMapping:
    def test_all_mode3d_buttons_mapped(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.callbacks.toolbar import _MODE3D_BTN

        assert _MODE3D_BTN[IDs.TB3D_MODE_FENCE] == "fence"
        assert _MODE3D_BTN[IDs.TB3D_MODE_BLOCK] == "block"
        assert _MODE3D_BTN[IDs.TB3D_MODE_DEPTH] == "depth"
        assert _MODE3D_BTN[IDs.TB3D_MODE_SURFACE] == "surface"


class TestRegisterToolbar:
    def test_register_toolbar_is_callable(self):
        from pycsamt.app.mapview.callbacks.toolbar import register_toolbar

        assert callable(register_toolbar)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.TB_INFO in cb_outputs
        assert IDs.CTL_BASEMAP in cb_outputs
        assert IDs.CRS_INFO in cb_outputs
