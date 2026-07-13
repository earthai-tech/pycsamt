# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for StationModel and StationTable (Phase 2)."""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import Qt

from pycsamt.app.desktop.models.station_model import (
    _COLUMNS,
    StationModel,
)

# ── Fixtures ──────────────────────────────────────────────────────────────


@pytest.fixture
def sample_df():
    return pd.DataFrame(
        {
            "ID": ["S01", "S02", "S03"],
            "Latitude": [48.5001, 48.5102, 48.5203],
            "Longitude": [7.7501, 7.7602, 7.7703],
            "Elevation": [200.0, 210.0, 220.0],
            "N_freq": [32, 32, 24],
            "Tipper": [True, False, True],
        }
    )


@pytest.fixture
def model(qapp, sample_df):
    m = StationModel()
    m.set_dataframe(sample_df)
    return m


# ── StationModel ──────────────────────────────────────────────────────────


def test_empty_model_row_count(qapp):
    assert StationModel().rowCount() == 0


def test_column_count(qapp):
    assert StationModel().columnCount() == len(_COLUMNS)


def test_row_count_after_set_dataframe(model):
    assert model.rowCount() == 3


def test_header_data_horizontal(model):
    for i, col in enumerate(_COLUMNS):
        assert model.headerData(i, Qt.Orientation.Horizontal) == col


def test_header_data_vertical_is_row_number(model):
    assert model.headerData(0, Qt.Orientation.Vertical) == "1"
    assert model.headerData(2, Qt.Orientation.Vertical) == "3"


def test_display_role_id(model):
    idx = model.index(0, 0)  # column 0 = ID
    assert model.data(idx, Qt.ItemDataRole.DisplayRole) == "S01"


def test_display_role_latitude_precision(model):
    idx = model.index(0, 1)  # Latitude
    val = model.data(idx, Qt.ItemDataRole.DisplayRole)
    assert "." in val
    # 5 decimal places
    assert len(val.split(".")[1]) == 5


def test_display_role_tipper_yes_no(model):
    idx_yes = model.index(0, 5)  # Tipper = True
    idx_no = model.index(1, 5)  # Tipper = False
    assert model.data(idx_yes, Qt.ItemDataRole.DisplayRole) == "yes"
    assert model.data(idx_no, Qt.ItemDataRole.DisplayRole) == "no"


def test_user_role_returns_raw_value(model):
    idx = model.index(0, 4)  # N_freq
    assert model.data(idx, Qt.ItemDataRole.UserRole) == 32


def test_invalid_index_returns_none(model):
    from PySide6.QtCore import QModelIndex

    assert model.data(QModelIndex()) is None


def test_station_id_at_row(model):
    assert model.station_id_at_row(1) == "S02"


def test_station_id_at_row_out_of_range(model):
    assert model.station_id_at_row(99) is None


def test_row_for_station_id(model):
    assert model.row_for_station_id("S03") == 2


def test_row_for_station_id_missing(model):
    assert model.row_for_station_id("NOPE") == -1


def test_clear_resets_row_count(model):
    model.clear()
    assert model.rowCount() == 0


# ── StationTable ──────────────────────────────────────────────────────────


def test_station_table_creates(qapp):
    from pycsamt.app.desktop.widgets.station_table import (
        StationTable,
    )

    t = StationTable()
    assert t is not None
    t.close()


def test_station_table_set_dataframe(qapp, sample_df):
    from pycsamt.app.desktop.widgets.station_table import (
        StationTable,
    )

    t = StationTable()
    t.set_dataframe(sample_df)
    assert t._model.rowCount() == 3
    t.close()


def test_station_table_emits_rows_selected(qapp, sample_df):
    from pycsamt.app.desktop.widgets.station_table import (
        StationTable,
    )

    t = StationTable()
    t.set_dataframe(sample_df)

    received = []
    t.rows_selected.connect(received.append)
    t.select_station_id("S02")
    # Signal was emitted with a list containing "S02"
    assert any("S02" in ids for ids in received)
    t.close()
