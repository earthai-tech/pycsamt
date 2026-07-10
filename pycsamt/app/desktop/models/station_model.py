# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationModel — QAbstractTableModel wrapping a pandas DataFrame of sites.

Phase 2: Columns exposed to StationTable:
  ID · Latitude · Longitude · Elevation · N_freq · Tipper

Supports sorting via QSortFilterProxyModel.  Data is refreshed with
beginResetModel / endResetModel so connected views redraw correctly.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
from PySide6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    Qt,
)

_COLUMNS = ["ID", "Latitude", "Longitude", "Elevation", "N_freq", "Tipper"]
_PRECISION = {"Latitude": 5, "Longitude": 5, "Elevation": 1}

# Unicode badge appended to station IDs that have been recomputed
_RECOMPUTED_BADGE = " ◈"

from PySide6.QtGui import QBrush, QColor


class StationModel(QAbstractTableModel):
    """Qt table model backed by a pandas DataFrame."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._df = pd.DataFrame(columns=_COLUMNS)
        self._recomputed_ids: set[str] = set()

    # ── Public mutators ───────────────────────────────────────────────

    def set_dataframe(self, df: pd.DataFrame) -> None:
        """Replace the backing DataFrame and refresh all views."""
        self._recomputed_ids.clear()   # clear before reset so views never see stale badges
        self.beginResetModel()
        self._df = df.reset_index(drop=True)
        self.endResetModel()

    def clear(self) -> None:
        self.set_dataframe(pd.DataFrame(columns=_COLUMNS))

    def mark_recomputed(self, station_ids: set[str]) -> None:
        """Mark a set of station IDs as recomputed (adds badge in the ID column)."""
        self._recomputed_ids = set(station_ids)
        # Notify the view that column 0 (ID) needs repainting.
        # Pass [] for roles so Qt treats all roles as changed (avoids PySide6
        # strict-int issues with enum values in QList<int>).
        if len(self._df) > 0:
            top    = self.index(0, 0)
            bottom = self.index(len(self._df) - 1, 0)
            self.dataChanged.emit(top, bottom, [])

    def clear_recomputed(self) -> None:
        """Remove all recomputed badges."""
        self.mark_recomputed(set())

    # ── QAbstractTableModel interface ─────────────────────────────────

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return len(self._df)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return len(_COLUMNS)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        col_name = _COLUMNS[index.column()]
        value = self._df.iloc[index.row()][col_name]

        if role == Qt.ItemDataRole.DisplayRole:
            if pd.isna(value):
                return "—"
            prec = _PRECISION.get(col_name)
            if prec is not None:
                return f"{float(value):.{prec}f}"
            if col_name == "Tipper":
                return "yes" if value else "no"
            text = str(value)
            if col_name == "ID" and text in self._recomputed_ids:
                return text + _RECOMPUTED_BADGE
            return text

        if role == Qt.ItemDataRole.UserRole:
            return value

        if role == Qt.ItemDataRole.ForegroundRole:
            if col_name == "ID" and str(value) in self._recomputed_ids:
                # Teal-green accent — visible in both dark and light themes
                return QBrush(QColor("#89dceb"))
            return None

        if role == Qt.ItemDataRole.ToolTipRole:
            if col_name == "ID" and str(value) in self._recomputed_ids:
                return "This station has been recomputed with EDIRecomputer (◈)"
            return None

        if role == Qt.ItemDataRole.TextAlignmentRole:
            if col_name == "ID":
                return int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        return None

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return _COLUMNS[section]
        return str(section + 1)

    # ── Helpers ───────────────────────────────────────────────────────

    def station_id_at_row(self, row: int) -> str | None:
        if 0 <= row < len(self._df):
            return str(self._df.iloc[row]["ID"])
        return None

    def row_for_station_id(self, station_id: str) -> int:
        """Return the row index for *station_id*, or -1 if not found."""
        matches = self._df.index[self._df["ID"] == station_id].tolist()
        return matches[0] if matches else -1
