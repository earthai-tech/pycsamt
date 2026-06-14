# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
CategoryComboBar — reusable two-level dropdown selector bar.

Layout:
    [Category ▾] [Item ▾] [↻] [⬆]

Mirrors the correction-window category/correction pattern so every panel
that groups content by category uses the same interaction model.

The item combo can be hidden via ``set_item_combo_visible(False)`` when a
panel only needs one level of selection (e.g. ProfilePanel view-type).
"""
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QWidget,
)


class CategoryComboBar(QWidget):
    """
    Two-level selector bar: [Category ▾] [Item ▾] [↻] [⬆]

    Signals
    -------
    category_changed(int)
        Fired when the category combo index changes.
    item_changed(int)
        Fired when the item combo index changes.
    refresh_clicked
        The ↻ button was pressed.
    export_clicked
        The ⬆ button was pressed.
    """

    category_changed = Signal(int)
    item_changed     = Signal(int)
    refresh_clicked  = Signal()
    export_clicked   = Signal()

    def __init__(
        self,
        category_label: str = "",
        item_label: str = "",
        show_item: bool = True,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setFixedHeight(34)

        bar = QHBoxLayout(self)
        bar.setContentsMargins(4, 2, 4, 2)
        bar.setSpacing(4)

        if category_label:
            lbl = QLabel(category_label)
            lbl.setObjectName("BarLabel")
            bar.addWidget(lbl)

        self._combo_cat = QComboBox()
        self._combo_cat.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_cat.currentIndexChanged.connect(self._on_cat_changed)
        bar.addWidget(self._combo_cat, 2)

        self._lbl_item = None
        if item_label:
            self._lbl_item = QLabel(item_label)
            self._lbl_item.setObjectName("BarLabel")
            bar.addWidget(self._lbl_item)

        self._combo_item = QComboBox()
        self._combo_item.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_item.currentIndexChanged.connect(self._on_item_changed)
        bar.addWidget(self._combo_item, 3)

        self._btn_refresh = QPushButton("↻")
        self._btn_refresh.setFixedWidth(28)
        self._btn_refresh.setToolTip("Refresh plot")
        self._btn_refresh.clicked.connect(self.refresh_clicked)
        bar.addWidget(self._btn_refresh)

        self._btn_export = QPushButton("⬆")
        self._btn_export.setFixedWidth(28)
        self._btn_export.setToolTip("Export figure…")
        self._btn_export.clicked.connect(self.export_clicked)
        bar.addWidget(self._btn_export)

        self.set_item_combo_visible(show_item)

    # ── Public API ────────────────────────────────────────────────────

    def set_categories(self, cats: list[str]) -> None:
        """Populate the category combo, suppressing intermediate signals."""
        self._combo_cat.blockSignals(True)
        self._combo_cat.clear()
        self._combo_cat.addItems(cats)
        self._combo_cat.blockSignals(False)

    def set_items(self, items: list[str]) -> None:
        """Populate the item combo, suppressing intermediate signals."""
        self._combo_item.blockSignals(True)
        self._combo_item.clear()
        self._combo_item.addItems(items)
        self._combo_item.blockSignals(False)

    def set_item_combo_visible(self, visible: bool) -> None:
        self._combo_item.setVisible(visible)
        if self._lbl_item is not None:
            self._lbl_item.setVisible(visible)

    def set_export_visible(self, visible: bool) -> None:
        self._btn_export.setVisible(visible)

    def set_refresh_visible(self, visible: bool) -> None:
        self._btn_refresh.setVisible(visible)

    def current_category(self) -> int:
        return self._combo_cat.currentIndex()

    def current_item(self) -> int:
        return self._combo_item.currentIndex()

    def current_category_text(self) -> str:
        return self._combo_cat.currentText()

    def current_item_text(self) -> str:
        return self._combo_item.currentText()

    def set_current_category(self, idx: int) -> None:
        self._combo_cat.setCurrentIndex(idx)

    def set_current_item(self, idx: int) -> None:
        self._combo_item.setCurrentIndex(idx)

    # ── Internal slots ────────────────────────────────────────────────

    def _on_cat_changed(self, idx: int) -> None:
        self.category_changed.emit(idx)

    def _on_item_changed(self, idx: int) -> None:
        self.item_changed.emit(idx)
