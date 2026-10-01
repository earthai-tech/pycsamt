# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationPanel — left dock widget showing the station list table.

Phase 2: Contains a StationTable and a summary label showing
"N stations loaded".  Propagates rows_selected → station_selected.
"""

from __future__ import annotations

import pandas as pd
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QLabel,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.station_table import (
    StationTable,
)


class StationPanel(QWidget):
    """Station list panel — contains the sortable station table."""

    station_selected = Signal(str)  # single station ID from click/keyboard

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        self._summary_lbl = QLabel("No stations loaded")
        self._summary_lbl.setObjectName("StationSummaryLabel")
        layout.addWidget(self._summary_lbl)

        self._table = StationTable(self)
        layout.addWidget(self._table)

        self._table.rows_selected.connect(self._on_rows_selected)

    # ── Data binding ──────────────────────────────────────────────────

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self._table.set_dataframe(df)
        self._update_summary()

    def clear(self) -> None:
        self._table.clear()
        self._summary_lbl.setText("No stations loaded")

    # ── Filtering (the main window's "Filter stations" box) ──────────

    def filter(self, text: str) -> int:
        """Show only stations whose row contains *text*; returns how many
        are shown.  (The search box called this method, which did not
        exist -- the error was swallowed, so filtering did nothing.)"""
        self._filter_text = (text or "").strip()
        shown = self._table.filter(self._filter_text)
        self._update_summary()
        return shown

    def _update_summary(self) -> None:
        total = self._table.total_count()
        shown = self._table.visible_count()
        if getattr(self, "_filter_text", "") and shown != total:
            self._summary_lbl.setText(
                f"{shown} of {total} stations shown")
        else:
            self._summary_lbl.setText(
                f"{total} station{'s' if total != 1 else ''} loaded")

    # ── Programmatic highlight (called by AppController) ─────────────

    def highlight_station(self, station_id: str) -> None:
        self._table.select_station_id(station_id)

    # ── Internal slot ─────────────────────────────────────────────────

    def _on_rows_selected(self, ids: list) -> None:
        if ids:
            self.station_selected.emit(ids[0])
