# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PseudosectionsPage — settings for PYCSAMT_STATION_RENDERING and PYCSAMT_SECTION.

Controls
--------
Station group
  • Side (top / bottom)
  • Show markers (checkbox)
  • Marker symbol (▽ v  ▲ ^  ● o  ■ s  ◆ D)
  • Marker size
  • Marker y-offset  (fraction of axes height above top edge)
  • Max station labels

Section group
  • Y-direction  (down = shallow at top / up = shallow at bottom)
  • Station side inherited from station group (synced automatically)
"""
from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QSpinBox,
    QVBoxLayout,
)

from .base_page import SettingsPage

_MARKER_CHOICES = [
    ("▽  (down-triangle)", "v"),
    ("△  (up-triangle)",   "^"),
    ("●  (circle)",        "o"),
    ("■  (square)",        "s"),
    ("◆  (diamond)",       "D"),
]

_Y_DIR_CHOICES = [
    ("Down  (shallow → top, standard)",  "down"),
    ("Up    (shallow → bottom)",          "up"),
    ("None  (no inversion applied)",      "none"),
]

_SIDE_CHOICES = [("Top", "top"), ("Bottom", "bottom")]


class PseudosectionsPage(SettingsPage):
    """Settings page for station markers and pseudosection axis style."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()
        self.populate()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        # ── Station group ─────────────────────────────────────────────────────
        grp_st = QGroupBox("Station Rendering  (PYCSAMT_STATION_RENDERING)")
        form_st = QFormLayout(grp_st)
        form_st.setSpacing(8)

        self._side_combo = QComboBox()
        for label, _ in _SIDE_CHOICES:
            self._side_combo.addItem(label)
        form_st.addRow("Station labels side:", self._side_combo)

        self._show_cb = QCheckBox("Show station markers")
        form_st.addRow("", self._show_cb)

        self._marker_combo = QComboBox()
        for label, _ in _MARKER_CHOICES:
            self._marker_combo.addItem(label)
        form_st.addRow("Marker symbol:", self._marker_combo)

        self._size_spin = QDoubleSpinBox()
        self._size_spin.setRange(1.0, 30.0)
        self._size_spin.setSingleStep(0.5)
        self._size_spin.setDecimals(1)
        form_st.addRow("Marker size (pt):", self._size_spin)

        self._offset_spin = QDoubleSpinBox()
        self._offset_spin.setRange(0.90, 1.15)
        self._offset_spin.setSingleStep(0.005)
        self._offset_spin.setDecimals(3)
        form_st.addRow("Marker y-offset:", self._offset_spin)

        self._maxlbl_spin = QSpinBox()
        self._maxlbl_spin.setRange(1, 100)
        form_st.addRow("Max station labels:", self._maxlbl_spin)

        root.addWidget(grp_st)

        # ── Section group ─────────────────────────────────────────────────────
        grp_sec = QGroupBox("Section Axis Style  (PYCSAMT_SECTION)")
        form_sec = QFormLayout(grp_sec)
        form_sec.setSpacing(8)

        self._ydir_combo = QComboBox()
        for label, _ in _Y_DIR_CHOICES:
            self._ydir_combo.addItem(label)
        form_sec.addRow("Y-direction:", self._ydir_combo)

        root.addWidget(grp_sec)
        root.addStretch()

    # ── SettingsPage interface ────────────────────────────────────────────────

    def populate(self) -> None:
        try:
            from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR
            ps = SR.pseudosection
            self._side_combo.setCurrentIndex(
                next((i for i, (_, v) in enumerate(_SIDE_CHOICES) if v == ps.side), 0)
            )
            self._show_cb.setChecked(ps.show_markers)
            self._marker_combo.setCurrentIndex(
                next((i for i, (_, v) in enumerate(_MARKER_CHOICES) if v == ps.marker.marker), 0)
            )
            self._size_spin.setValue(ps.marker.size)
            self._offset_spin.setValue(ps.marker.offset)
            self._maxlbl_spin.setValue(ps.max_labels)
        except Exception:
            pass

        try:
            from pycsamt.api.section import PYCSAMT_SECTION as SEC
            y_dir = SEC.pseudosection.axis.y_direction
            self._ydir_combo.setCurrentIndex(
                next((i for i, (_, v) in enumerate(_Y_DIR_CHOICES) if v == y_dir), 0)
            )
        except Exception:
            pass

    def collect(self) -> dict:
        station = {
            "side":          _SIDE_CHOICES[self._side_combo.currentIndex()][1],
            "show_markers":  self._show_cb.isChecked(),
            "marker_symbol": _MARKER_CHOICES[self._marker_combo.currentIndex()][1],
            "marker_size":   self._size_spin.value(),
            "marker_offset": self._offset_spin.value(),
            "max_labels":    self._maxlbl_spin.value(),
        }
        section = {
            "y_direction":  _Y_DIR_CHOICES[self._ydir_combo.currentIndex()][1],
            "station_side": _SIDE_CHOICES[self._side_combo.currentIndex()][1],
        }
        return {"station": station, "section": section}

    def reset(self) -> None:
        try:
            from pycsamt.api.station import PYCSAMT_STATION_RENDERING
            PYCSAMT_STATION_RENDERING.reset()
        except Exception:
            pass
        try:
            from pycsamt.api.section import PYCSAMT_SECTION
            PYCSAMT_SECTION.reset()
        except Exception:
            pass
        self.populate()
