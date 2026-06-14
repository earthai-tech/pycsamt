# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
OverlayCombo — QComboBox for choosing the station colour-coding attribute.

Phase 2/3: Emits overlay_changed(str) when the selection changes so
MapPanel can redraw station colours.
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QWidget

OVERLAY_OPTIONS = [
    "Index",
    "Apparent Resistivity",
    "Phase",
    "Quality Score",
    "Z-strike",
    "Skew",
    "Tipper magnitude",
    "Ellipticity",
]


class OverlayCombo(QComboBox):
    """Drop-down that selects what attribute colours station markers."""

    overlay_changed = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.addItems(OVERLAY_OPTIONS)
        self.currentTextChanged.connect(self.overlay_changed)

    @property
    def current_overlay(self) -> str:
        return self.currentText()
