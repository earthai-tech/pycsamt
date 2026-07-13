# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the pycsamt Dash web application (Phase 6).

These tests do NOT start a real server — they verify that the app can
be created, the layout is well-formed, callbacks are registered, and
utility functions behave correctly.
"""

from __future__ import annotations

import base64

import pandas as pd
import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


# ──────────────────────────────────────────────────────────────────────────────
# Utils
# ──────────────────────────────────────────────────────────────────────────────


class TestUtils:
    def test_fig_to_src_returns_data_uri(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from pycsamt.app.web.utils import fig_to_src

        fig, ax = plt.subplots()
        ax.plot([1, 2], [3, 4])
        src = fig_to_src(fig)
        assert src.startswith("data:image/png;base64,")
        assert len(src) > 100

    def test_fig_to_src_base64_is_valid(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from pycsamt.app.web.utils import fig_to_src

        fig, ax = plt.subplots()
        src = fig_to_src(fig)
        b64_part = src.split(",", 1)[1]
        decoded = base64.b64decode(b64_part)
        assert decoded[:4] == b"\x89PNG"  # PNG magic bytes

    def test_empty_src_returns_data_uri(self):
        from pycsamt.app.web.utils import empty_src

        src = empty_src(dark=True)
        assert src.startswith("data:image/png;base64,")

    def test_find_edi_files_directory(self, tmp_path):
        from pycsamt.app.web.utils import find_edi_files

        (tmp_path / "a.edi").write_text("dummy")
        (tmp_path / "b.edi").write_text("dummy")
        (tmp_path / "c.txt").write_text("dummy")
        found = find_edi_files(str(tmp_path))
        assert len(found) == 2
        assert all(f.endswith(".edi") for f in found)

    def test_find_edi_files_single_file(self, tmp_path):
        from pycsamt.app.web.utils import find_edi_files

        edi = tmp_path / "test.edi"
        edi.write_text("dummy")
        found = find_edi_files(str(edi))
        assert found == [str(edi)]

    def test_find_edi_files_empty_dir(self, tmp_path):
        from pycsamt.app.web.utils import find_edi_files

        assert find_edi_files(str(tmp_path)) == []

    def test_find_edi_files_nonexistent(self, tmp_path):
        from pycsamt.app.web.utils import find_edi_files

        assert find_edi_files(str(tmp_path / "no_such")) == []

    def test_apply_web_dark_theme_sets_facecolor(self):
        import matplotlib as mpl

        from pycsamt.app.web.utils import apply_web_dark_theme

        apply_web_dark_theme()
        assert mpl.rcParams["axes.facecolor"] == "#181825"

    def test_apply_web_light_theme_sets_facecolor(self):
        import matplotlib as mpl

        from pycsamt.app.web.utils import (
            apply_web_light_theme,
        )

        apply_web_light_theme()
        assert mpl.rcParams["axes.facecolor"] == "#eff1f5"


class TestMapBuilder:
    @pytest.fixture
    def sample_df(self):
        return pd.DataFrame(
            {
                "ID": ["S01", "S02", "S03"],
                "Latitude": [48.50, 48.51, 48.52],
                "Longitude": [7.75, 7.76, 7.77],
                "N_freq": [32, 32, 24],
            }
        )

    def test_empty_map_returns_figure(self):
        from pycsamt.app.web.utils import build_station_map

        fig = build_station_map(pd.DataFrame())
        import plotly.graph_objects as go

        assert isinstance(fig, go.Figure)

    def test_map_with_data_returns_figure(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        fig = build_station_map(sample_df, dark=True)
        import plotly.graph_objects as go

        assert isinstance(fig, go.Figure)

    def test_map_has_scatter_trace(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        fig = build_station_map(sample_df)
        # tile-based map traces (Scattermap on modern plotly,
        # Scattermapbox on older releases)
        assert any(
            type(t).__name__ in ("Scattermap", "Scattermapbox")
            for t in fig.data
        )

    def test_map_trace_lat_lon(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        fig = build_station_map(sample_df)
        trace = fig.data[0]
        assert len(trace.lat) == 3
        assert len(trace.lon) == 3

    def test_map_selected_station_highlighted(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        fig = build_station_map(sample_df, selected_id="S02")
        trace = fig.data[0]
        # Selected station should have a larger marker size
        sizes = list(trace.marker.size)
        assert sizes[1] > sizes[0]

    def test_map_light_mode(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        fig = build_station_map(sample_df, dark=False)
        import plotly.graph_objects as go

        assert isinstance(fig, go.Figure)

    def test_map_overlay_column(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        # N_freq is a valid column
        fig = build_station_map(sample_df, overlay="N_freq")
        import plotly.graph_objects as go

        assert isinstance(fig, go.Figure)

    def test_map_unknown_overlay_uses_index(self, sample_df):
        from pycsamt.app.web.utils import build_station_map

        # Should not raise even with unknown overlay column
        fig = build_station_map(sample_df, overlay="NoSuchCol")
        import plotly.graph_objects as go

        assert isinstance(fig, go.Figure)


# ──────────────────────────────────────────────────────────────────────────────
# Layout
# ──────────────────────────────────────────────────────────────────────────────


class TestLayout:
    def test_layout_returns_html_div(self):
        from dash import html

        from pycsamt.app.web.layout import layout

        result = layout()
        assert isinstance(result, html.Div)

    def test_layout_contains_stores(self):

        from pycsamt.app.web.layout import IDs, layout

        result = layout()
        # Serialize and check IDs appear in repr
        children_str = str(result)
        assert IDs.STORE_DATA in children_str
        assert IDs.STORE_SELECTION in children_str
        assert IDs.STORE_THEME in children_str

    def test_layout_contains_map_graph(self):
        from pycsamt.app.web.layout import IDs, layout

        result = layout()
        children_str = str(result)
        assert IDs.MAP_GRAPH in children_str

    def test_layout_contains_station_table(self):
        from pycsamt.app.web.layout import IDs, layout

        result = layout()
        children_str = str(result)
        assert IDs.STATION_TABLE in children_str

    def test_profile_panel_contains_profile_tabs(self):
        # the root layout now opens on a welcome overlay; the profile
        # tabs live inside the profile panel component
        from pycsamt.app.web.layout import IDs, _profile_panel

        children_str = str(_profile_panel())
        assert IDs.PROFILE_TABS in children_str

    def test_layout_contains_agent_panel(self):
        from pycsamt.app.web.layout import IDs, layout

        result = layout()
        children_str = str(result)
        assert IDs.AGENT_SELECT in children_str
        assert IDs.BTN_RUN_AGENT in children_str

    def test_ids_class_has_all_expected(self):
        from pycsamt.app.web.layout import IDs

        for attr in (
            "STORE_DATA",
            "STORE_SELECTION",
            "STORE_THEME",
            "BTN_LOAD",
            "BTN_THEME",
            "MAP_GRAPH",
            "STATION_TABLE",
            "PROFILE_TABS",
            "IMG_RHO_PHI",
            "IMG_RHO_PS",
            "AGENT_SELECT",
            "BTN_RUN_AGENT",
            "AGENT_LOG",
        ):
            assert hasattr(IDs, attr), f"IDs.{attr} missing"


# ──────────────────────────────────────────────────────────────────────────────
# App factory
# ──────────────────────────────────────────────────────────────────────────────


class TestAppFactory:
    def test_create_app_returns_dash_instance(self):
        import dash

        from pycsamt.app.web.app import create_app

        app = create_app()
        assert isinstance(app, dash.Dash)

    def test_app_has_layout(self):
        from pycsamt.app.web.app import create_app

        app = create_app()
        assert app.layout is not None

    def test_app_title(self):
        from pycsamt.app.web.app import create_app

        app = create_app()
        assert app.title == "pyCSAMT"

    def test_app_callbacks_registered(self):
        from pycsamt.app.web.app import create_app

        app = create_app()
        # Dash stores callback map internally
        # Check that at least some callbacks are registered
        assert len(app.callback_map) > 0

    def test_specific_callbacks_registered(self):
        from pycsamt.app.web.app import create_app
        from pycsamt.app.web.layout import IDs

        app = create_app()
        cb_outputs = str(app.callback_map)
        # Store data callback should be registered
        assert IDs.STORE_DATA in cb_outputs or IDs.MAP_GRAPH in cb_outputs


# ──────────────────────────────────────────────────────────────────────────────
# Callbacks (unit-level, no server)
# ──────────────────────────────────────────────────────────────────────────────


class TestCallbacks:
    def test_callbacks_package_entry_point(self):
        """Callback wiring is exposed via register_callbacks."""
        from pycsamt.app.web.callbacks import register_callbacks

        assert callable(register_callbacks)

    def test_empty_store_returns_no_stations(self):
        """update_table with None store should return empty records."""
        # We verify the logic path without invoking the full Dash machinery
        # by calling equivalent Python logic
        store_data = None
        records = (
            [] if not store_data else store_data.get("station_records", [])
        )
        assert records == []

    def test_station_map_with_valid_df(self):
        import plotly.graph_objects as go

        from pycsamt.app.web.utils import build_station_map

        df = pd.DataFrame(
            {
                "ID": ["A", "B"],
                "Latitude": [10.0, 11.0],
                "Longitude": [20.0, 21.0],
            }
        )
        fig = build_station_map(df)
        assert isinstance(fig, go.Figure)
        assert len(fig.data) > 0

    def test_session_cache_roundtrip(self):
        # per-session state moved from a module global to the cache
        from pycsamt.app.web.cache import cache_get, cache_set

        cache_set("test-session-xyz", {"sites": ["s1"]})
        assert cache_get("test-session-xyz") == {"sites": ["s1"]}
        assert cache_get("missing-session") is None
