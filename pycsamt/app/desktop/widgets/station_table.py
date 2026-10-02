# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationTable — QTableView backed by StationModel with sort proxy.

Phase 2: Emits rows_selected(list[str]) with station IDs whenever
the selection changes.  Programmatic selection via select_station_id()
scrolls to and highlights the matching row.
"""

from __future__ import annotations

import pandas as pd
from PySide6.QtCore import QSortFilterProxyModel, Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QTableView,
    QWidget,
)

from pycsamt.app.desktop.models.station_model import (
    StationModel,
)


class _StationFilter(QSortFilterProxyModel):
    """Matches the filter text in any shown column *or* the station's line
    (kept in the model's DataFrame, not shown as a column).  Comma-separated
    terms are alternatives: ``S01, S05`` or ``L22, L26``."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._terms: list[str] = []

    def set_terms(self, text: str) -> None:
        self._terms = [t.strip().lower() for t in (text or "").split(",")
                       if t.strip()]
        self.invalidateFilter()

    def filterAcceptsRow(self, row, parent) -> bool:  # noqa: N802
        if not self._terms:
            return True
        model = self.sourceModel()
        cells = [str(model.data(model.index(row, c, parent)) or "").lower()
                 for c in range(model.columnCount())]
        df = getattr(model, "_df", None)
        if df is not None and "Line" in df.columns and row < len(df):
            cells.append(str(df["Line"].iloc[row]).lower())
        return any(term in cell for term in self._terms for cell in cells)


class StationTable(QTableView):
    """Sortable station list table."""

    rows_selected = Signal(list)  # list[str] — station IDs

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._model = StationModel(self)
        self._proxy = _StationFilter(self)
        self._proxy.setSourceModel(self._model)
        self._proxy.setSortCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.setModel(self._proxy)
        self._setup_view()
        self.selectionModel().selectionChanged.connect(
            self._on_selection_changed
        )

    # ── View configuration ────────────────────────────────────────────

    def _setup_view(self) -> None:
        self.setSortingEnabled(True)
        self.sortByColumn(0, Qt.SortOrder.AscendingOrder)
        self.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self.setAlternatingRowColors(True)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setWordWrap(False)

    # ── Data binding ──────────────────────────────────────────────────

    def set_dataframe(self, df: pd.DataFrame) -> None:
        """Load a new station DataFrame into the model."""
        self._model.set_dataframe(df)
        self.resizeColumnsToContents()
        self.horizontalHeader().setStretchLastSection(True)

    def clear(self) -> None:
        self._model.clear()

    # ── Filtering ─────────────────────────────────────────────────────

    def filter(self, text: str) -> int:
        """Show only rows containing *text* in any column or their line
        (case-insensitive; comma-separated terms are alternatives);
        returns the number of rows shown."""
        self._proxy.set_terms(text)
        return self._proxy.rowCount()

    def visible_count(self) -> int:
        return self._proxy.rowCount()

    def total_count(self) -> int:
        return self._model.rowCount()

    # ── Programmatic selection ────────────────────────────────────────

    def select_station_id(self, station_id: str) -> None:
        """Highlight the row for *station_id* and scroll to it."""
        src_row = self._model.row_for_station_id(station_id)
        if src_row < 0:
            return
        proxy_idx = self._proxy.mapFromSource(self._model.index(src_row, 0))
        if not proxy_idx.isValid():  # hidden by the current filter
            return
        self.selectionModel().clearSelection()
        self.selectionModel().select(
            proxy_idx,
            self.selectionModel().SelectionFlag.Select
            | self.selectionModel().SelectionFlag.Rows,
        )
        self.scrollTo(proxy_idx, QAbstractItemView.ScrollHint.PositionAtCenter)

    # ── Selection signal ──────────────────────────────────────────────

    def _on_selection_changed(self, selected, deselected) -> None:
        rows = self.selectionModel().selectedRows()
        ids: list[str] = []
        for proxy_idx in rows:
            src_idx = self._proxy.mapToSource(proxy_idx)
            sid = self._model.station_id_at_row(src_idx.row())
            if sid:
                ids.append(sid)
        self.rows_selected.emit(ids)
