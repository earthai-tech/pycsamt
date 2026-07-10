# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TopographyPage — settings for PYCSAMT_TOPO (pycsamt.topo.config.TopoConfig).

Controls
--------
  • Enabled (checkbox)
  • Vertical exaggeration (spin box)
  • Marker pad fraction  (fine-tunes the gap between topo surface and station
    markers — expressed as a fraction of the axes height)
"""
from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QVBoxLayout,
)

from .base_page import SettingsPage


class TopographyPage(SettingsPage):
    """Settings page for the terrain-following geometry integration."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()
        self.populate()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        grp = QGroupBox("Topography Integration  (PYCSAMT_TOPO)")
        form = QFormLayout(grp)
        form.setSpacing(8)

        self._enabled_cb = QCheckBox("Enable terrain-following geometry")
        form.addRow("", self._enabled_cb)

        self._exag_spin = QDoubleSpinBox()
        self._exag_spin.setRange(0.1, 20.0)
        self._exag_spin.setSingleStep(0.5)
        self._exag_spin.setDecimals(1)
        form.addRow("Vertical exaggeration:", self._exag_spin)

        self._pad_spin = QDoubleSpinBox()
        self._pad_spin.setRange(0.001, 0.20)
        self._pad_spin.setSingleStep(0.001)
        self._pad_spin.setDecimals(3)
        form.addRow("Marker pad fraction:", self._pad_spin)

        note = QLabel(
            "<i>Topography is only active when elevation data is attached to the<br>"
            "survey (via the Map Viewer or loaded alongside EDI files).</i>"
        )
        note.setWordWrap(True)
        form.addRow(note)

        root.addWidget(grp)
        root.addStretch()

    # ── SettingsPage interface ────────────────────────────────────────────────

    def populate(self) -> None:
        try:
            from pycsamt.topo import PYCSAMT_TOPO as T
            self._enabled_cb.setChecked(T.enabled)
            self._exag_spin.setValue(T.exaggeration)
            self._pad_spin.setValue(T.marker_pad_fraction)
        except Exception:
            self._enabled_cb.setChecked(False)
            self._exag_spin.setValue(1.0)
            self._pad_spin.setValue(0.015)

    def collect(self) -> dict:
        return {
            "topography": {
                "enabled":     self._enabled_cb.isChecked(),
                "exaggeration": self._exag_spin.value(),
                "marker_pad":  self._pad_spin.value(),
            }
        }

    def reset(self) -> None:
        try:
            from pycsamt.topo import PYCSAMT_TOPO
            PYCSAMT_TOPO.reset()
        except Exception:
            pass
        self.populate()
