# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SearchableComboBox — a station-picker with an integrated search row in its
popup.

Closed state
  [ — select station (N) — ][ ▼ ]
  or
  [ 18-003A                ][ ▼ ]  (after selection)

Open state (custom popup panel)
  ┌──────────────────────────────────┐
  │ 🔍  Type to search…          [×] │
  ├──────────────────────────────────┤
  │ 18-001A                          │
  │ 18-002U                          │
  │ 18-003A                     ↕    │
  │ …  (up to max_visible rows)      │
  └──────────────────────────────────┘

Typing in the search field filters the list in real-time (substring,
case-insensitive).  Arrow keys move the highlight, Enter or a mouse click
confirms the selection.  Escape or clicking outside dismisses the popup.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import (
    QComboBox,
    QFrame,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
    QWidget,
)


class _StationPopup(QFrame):
    """Internal floating popup: search field + scrollable station list."""

    station_chosen = Signal(str)  # emitted when user confirms a station

    def __init__(self, max_visible: int = 12, parent: QWidget | None = None) -> None:
        # Qt.WindowType.Popup: auto-closes on Escape / click outside
        super().__init__(parent, Qt.WindowType.Popup)
        self.setObjectName("StationPopupFrame")
        self._max_visible = max_visible
        self._all_names: list[str] = []
        self._setup_ui()

    # ── UI construction ───────────────────────────────────────────────

    def _setup_ui(self) -> None:
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(2)

        # Search row
        self._search = QLineEdit(self)
        self._search.setObjectName("StationSearch")
        self._search.setPlaceholderText("🔍  Type to search…")
        self._search.setClearButtonEnabled(True)
        self._search.textChanged.connect(self._apply_filter)
        lay.addWidget(self._search)

        # Thin separator
        sep = QFrame(self)
        sep.setObjectName("StationPopupSep")
        sep.setFrameShape(QFrame.Shape.HLine)
        lay.addWidget(sep)

        # Station list
        self._list = QListWidget(self)
        self._list.setObjectName("StationList")
        self._list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._list.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._list.setEditTriggers(QListWidget.EditTrigger.NoEditTriggers)
        self._list.itemClicked.connect(lambda it: self.station_chosen.emit(it.text()))
        lay.addWidget(self._list)

    # ── Public API ────────────────────────────────────────────────────

    def set_names(self, names: list[str]) -> None:
        self._all_names = list(names)
        self._refresh_list(self._all_names)

    def reset_search(self) -> None:
        """Clear the search field and restore the full list."""
        self._search.blockSignals(True)
        self._search.clear()
        self._search.blockSignals(False)
        self._refresh_list(self._all_names)

    def preferred_height(self) -> int:
        """Height that accommodates the search row + up to max_visible items."""
        row_h = max(self._list.sizeHintForRow(0), 24) if self._list.count() else 26
        n = min(self._max_visible, max(1, len(self._all_names)))
        return (
            self._search.sizeHint().height()
            + 8  # separator + spacing
            + n * row_h
            + 12  # frame padding
        )

    # ── Internals ─────────────────────────────────────────────────────

    def _refresh_list(self, names: list[str]) -> None:
        self._list.blockSignals(True)
        self._list.clear()
        for name in names:
            self._list.addItem(QListWidgetItem(name))
        self._list.blockSignals(False)

    def _apply_filter(self, text: str) -> None:
        lo = text.lower().strip()
        filtered = (
            [n for n in self._all_names if lo in n.lower()] if lo else self._all_names
        )
        self._refresh_list(filtered)
        # Auto-highlight first match so Enter immediately confirms it
        if self._list.count():
            self._list.setCurrentRow(0)

    # ── Key navigation ────────────────────────────────────────────────

    def keyPressEvent(self, event: QKeyEvent) -> None:
        key = event.key()
        if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            sel = self._list.selectedItems()
            if sel:
                self.station_chosen.emit(sel[0].text())
            elif self._list.count():
                self.station_chosen.emit(self._list.item(0).text())
        elif key == Qt.Key.Key_Escape:
            self.hide()
        elif key == Qt.Key.Key_Down:
            # Move list selection down; if nothing selected pick row 0
            row = self._list.currentRow()
            if row < self._list.count() - 1:
                self._list.setCurrentRow(row + 1)
            self._list.setFocus()
        elif key == Qt.Key.Key_Up:
            row = self._list.currentRow()
            if row > 0:
                self._list.setCurrentRow(row - 1)
            self._list.setFocus()
        else:
            # Route all other keystrokes to the search field
            if not self._search.hasFocus():
                self._search.setFocus()
            super().keyPressEvent(event)


class SearchableComboBox(QComboBox):
    """
    Non-editable station picker with an integrated search field in its popup.

    Closed state shows the selected station name or a placeholder.
    The popup contains a live-filter search row and a scrollable list.

    Signals
    -------
    station_selected(str)
        Emitted when the user confirms a station choice (click or Enter).
        Not emitted by programmatic :meth:`select_station` calls.
    """

    station_selected = Signal(str)

    _PLACEHOLDER = "— select station —"

    def __init__(
        self,
        max_visible: int = 12,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("StationCombo")
        self.setEditable(False)

        # Seed with placeholder so the combo is never visually empty
        self.addItem(self._PLACEHOLDER)

        # Custom popup (created once, reused)
        self._popup = _StationPopup(max_visible=max_visible)
        self._popup.station_chosen.connect(self._on_popup_chosen)

    # ── Public API ────────────────────────────────────────────────────

    def set_names(self, names: list[str]) -> None:
        """Replace the station list.  Resets any current selection."""
        n = len(names)
        placeholder = f"— select station ({n}) —" if n > 0 else self._PLACEHOLDER
        self.blockSignals(True)
        self.clear()
        self.addItem(placeholder)  # index 0 — the "nothing selected" state
        for name in names:
            self.addItem(name)  # index 1..N
        self.setCurrentIndex(0)
        self.blockSignals(False)

        self._popup.set_names(names)

    def select_station(self, name: str) -> None:
        """Programmatically set the active station (no signal emitted)."""
        idx = self.findText(name)
        if idx > 0:
            self.blockSignals(True)
            self.setCurrentIndex(idx)
            self.blockSignals(False)

    def current_station(self) -> str:
        """Return the currently selected station name, or '' if none."""
        idx = self.currentIndex()
        if idx <= 0:
            return ""
        return self.currentText()

    # ── Qt overrides ──────────────────────────────────────────────────

    def showPopup(self) -> None:
        """Override: open our custom popup instead of the built-in dropdown."""
        pos = self.mapToGlobal(self.rect().bottomLeft())
        self._popup.reset_search()
        w = max(self.width(), 220)
        h = self._popup.preferred_height()
        self._popup.move(pos)
        self._popup.resize(w, h)
        self._popup.show()
        self._popup.raise_()
        self._popup._search.setFocus()

    def hidePopup(self) -> None:
        self._popup.hide()
        super().hidePopup()

    # ── Internal slot ─────────────────────────────────────────────────

    def _on_popup_chosen(self, name: str) -> None:
        self._popup.hide()
        self.select_station(name)
        self.station_selected.emit(name)
