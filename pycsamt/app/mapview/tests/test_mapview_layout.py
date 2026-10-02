# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Smoke tests for pycsamt.app.mapview.layout — the Dash component tree."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestCreateLayout:
    def test_returns_html_div(self):
        from dash import html

        from pycsamt.app.mapview.layout import create_layout

        result = create_layout()
        assert isinstance(result, html.Div)

    def test_root_id_and_class(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.layout import create_layout

        result = create_layout()
        assert result.id == IDs.ROOT
        assert "mv-root" in result.className

    def test_contains_stores(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.layout import create_layout

        children_str = str(create_layout())
        for store_id in (
            IDs.STORE_DATA,
            IDs.STORE_THEME,
            IDs.STORE_VIEW,
            IDs.STORE_SELECTION,
        ):
            assert store_id in children_str

    def test_contains_load_and_help_modals(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.layout import create_layout

        children_str = str(create_layout())
        assert IDs.MODAL_LOAD in children_str
        assert IDs.MODAL_HELP in children_str

    def test_contains_canvas_and_inspector(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.layout import create_layout

        children_str = str(create_layout())
        assert IDs.CANVAS_GRAPH in children_str
        assert IDs.INSPECTOR in children_str

    def test_contains_view_rail(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.layout import create_layout

        children_str = str(create_layout())
        assert IDs.RAIL_MAP in children_str
        assert IDs.RAIL_3D in children_str

    def test_settings_and_session_canvases_present(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.layout import create_layout

        children_str = str(create_layout())
        assert IDs.CANVAS_SETTINGS in children_str
        assert IDs.CANVAS_SESSION in children_str

    def test_call_returns_a_fresh_tree_each_time(self):
        # a new session id is minted per call, so the trees are distinct
        # objects but structurally equivalent components.
        from pycsamt.app.mapview.layout import create_layout

        first = create_layout()
        second = create_layout()
        assert first is not second
        assert type(first) is type(second)


class TestDashTablePlaceholder:
    def test_returns_a_component(self):
        from pycsamt.app.mapview.layout import dash_table_placeholder

        placeholder = dash_table_placeholder()
        assert placeholder is not None
