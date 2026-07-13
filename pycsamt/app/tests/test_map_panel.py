# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MapPanel (Phase 2)."""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.map_panel import MapPanel


@pytest.fixture
def coords_df():
    return pd.DataFrame(
        {
            "ID": ["A01", "A02", "A03"],
            "Latitude": [48.50, 48.51, 48.52],
            "Longitude": [7.75, 7.76, 7.77],
        }
    )


@pytest.fixture
def panel(qapp):
    p = MapPanel()
    yield p
    p.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_map_panel_creates(qapp):
    p = MapPanel()
    assert p is not None
    p.close()


def test_canvas_created(panel):
    assert panel._canvas is not None


def test_panel_public_surface(panel):
    # the overlay combo moved out of the panel; the panel now exposes
    # data-binding and redraw entry points instead
    assert callable(panel.set_dataframe)
    assert callable(panel.set_sites)
    assert callable(panel.redraw)
    assert callable(panel.highlight_station)


# ── Data binding ──────────────────────────────────────────────────────────


def test_set_dataframe_does_not_raise(panel, coords_df):
    panel.set_dataframe(coords_df)


def test_set_dataframe_stores_data(panel, coords_df):
    panel.set_dataframe(coords_df)
    assert len(panel._df) == 3


def test_clear_empties_dataframe(panel, coords_df):
    panel.set_dataframe(coords_df)
    panel.clear()
    assert panel._df.empty


def test_scatter_created_after_set_dataframe(panel, coords_df):
    panel.set_dataframe(coords_df)
    assert panel._scatter is not None


# ── Highlight ────────────────────────────────────────────────────────────


def test_highlight_station_sets_selected_id(panel, coords_df):
    panel.set_dataframe(coords_df)
    panel.highlight_station("A02")
    assert panel._selected_id == "A02"


def test_highlight_unknown_station_does_not_raise(panel, coords_df):
    panel.set_dataframe(coords_df)
    panel.highlight_station("UNKNOWN")  # must not raise


# ── Signal ───────────────────────────────────────────────────────────────


def test_station_selected_signal_exists(panel):
    assert hasattr(panel, "station_selected")
