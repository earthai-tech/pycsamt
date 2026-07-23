# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for StationPanel.

Note: pycsamt/app/tests/test_station_model.py tests a different module
(pycsamt.app.desktop.models.station_model) despite the similar name — this
file targets the panel widget itself (pycsamt.app.desktop.panels.station_panel).
"""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.station_panel import StationPanel


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "ID": ["S01", "S02", "S03"],
            "Line": ["L1", "L1", "L2"],
            "Latitude": [48.50, 48.51, 48.52],
            "Longitude": [7.75, 7.76, 7.77],
            "Elevation": [200.0, 201.0, 202.0],
            "N_freq": [10, 12, 8],
            "Tipper": [True, False, True],
        }
    )


@pytest.fixture
def panel(qapp):
    p = StationPanel()
    yield p
    p.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_creates(qapp):
    p = StationPanel()
    assert p is not None
    p.close()


def test_initial_summary_label(panel):
    assert panel._summary_lbl.text() == "No stations loaded"


def test_table_created(panel):
    assert panel._table is not None


# ── set_dataframe ─────────────────────────────────────────────────────────


def test_set_dataframe_updates_summary_plural(panel, df):
    panel.set_dataframe(df)
    assert panel._summary_lbl.text() == "3 stations loaded"


def test_set_dataframe_singular(panel, df):
    panel.set_dataframe(df.iloc[:1])
    assert panel._summary_lbl.text() == "1 station loaded"


def test_set_dataframe_empty(panel, df):
    panel.set_dataframe(df.iloc[:0])
    assert panel._summary_lbl.text() == "0 stations loaded"


# ── clear ─────────────────────────────────────────────────────────────────


def test_clear_resets_summary(panel, df):
    panel.set_dataframe(df)
    panel.clear()
    assert panel._summary_lbl.text() == "No stations loaded"


# ── highlight_station ───────────────────────────────────────────────────


def test_highlight_station_delegates_to_table(panel, df, monkeypatch):
    from unittest import mock

    panel.set_dataframe(df)
    with mock.patch.object(panel._table, "select_station_id") as m:
        panel.highlight_station("S02")
        m.assert_called_once_with("S02")


# ── station_selected signal ──────────────────────────────────────────────


def test_rows_selected_emits_first_id(panel, df):
    panel.set_dataframe(df)
    received = []
    panel.station_selected.connect(received.append)
    panel._on_rows_selected(["S02", "S03"])
    assert received == ["S02"]


def test_rows_selected_empty_list_no_emit(panel):
    received = []
    panel.station_selected.connect(received.append)
    panel._on_rows_selected([])
    assert received == []


def test_table_selection_signal_wired(panel, df):
    """rows_selected from the real StationTable propagates through to
    station_selected (end-to-end signal wiring, not just the handler)."""
    panel.set_dataframe(df)
    received = []
    panel.station_selected.connect(received.append)
    panel._table.rows_selected.emit(["S01"])
    assert received == ["S01"]
