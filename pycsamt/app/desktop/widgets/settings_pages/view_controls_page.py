# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ViewControlsPage — settings for PYCSAMT_CONTROL.

Controls
--------
Resistivity
  • Display scale  (log10 / linear)

Phase
  • Range min / max (degrees)
  • Unit  (degree / radian)
  • Wrap  (checkbox)

X-axis  (period / frequency axis on profile plots)
  • View  (log10_period / period / log10_frequency / frequency)
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QVBoxLayout,
)

from .base_page import SettingsPage

_RHO_VIEWS = [("log₁₀ ρ  (log10)", "log10"), ("ρ  (linear)", "linear")]
_PHASE_UNITS = [("Degree (°)", "degree"), ("Radian (rad)", "radian")]
_X_VIEWS = [
    ("log₁₀ Period  (s)  — recommended", "log10_period"),
    ("Period  (s)", "period"),
    ("log₁₀ Frequency  (Hz)", "log10_frequency"),
    ("Frequency  (Hz)", "frequency"),
]


class ViewControlsPage(SettingsPage):
    """Settings page for resistivity, phase, and x-axis view controls."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()
        self.populate()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        # ── Resistivity ───────────────────────────────────────────────────────
        grp_rho = QGroupBox("Apparent Resistivity  (ρ_a)")
        form_rho = QFormLayout(grp_rho)
        form_rho.setSpacing(8)

        self._rho_combo = QComboBox()
        for label, _ in _RHO_VIEWS:
            self._rho_combo.addItem(label)
        form_rho.addRow("Display scale:", self._rho_combo)
        root.addWidget(grp_rho)

        # ── Phase ─────────────────────────────────────────────────────────────
        grp_ph = QGroupBox("Phase  (φ)")
        form_ph = QFormLayout(grp_ph)
        form_ph.setSpacing(8)

        range_row = QHBoxLayout()
        self._ph_min = QDoubleSpinBox()
        self._ph_min.setRange(-360.0, 360.0)
        self._ph_min.setSingleStep(10.0)
        self._ph_min.setDecimals(1)
        self._ph_max = QDoubleSpinBox()
        self._ph_max.setRange(-360.0, 360.0)
        self._ph_max.setSingleStep(10.0)
        self._ph_max.setDecimals(1)
        range_row.addWidget(self._ph_min)
        range_row.addWidget(QLabel("→"))
        range_row.addWidget(self._ph_max)
        form_ph.addRow("Range (°):", range_row)

        self._ph_unit_combo = QComboBox()
        for label, _ in _PHASE_UNITS:
            self._ph_unit_combo.addItem(label)
        form_ph.addRow("Unit:", self._ph_unit_combo)

        self._ph_wrap_cb = QCheckBox("Wrap phase into range")
        form_ph.addRow("", self._ph_wrap_cb)
        root.addWidget(grp_ph)

        # ── X-axis ────────────────────────────────────────────────────────────
        grp_x = QGroupBox("X-axis (profile plots)")
        form_x = QFormLayout(grp_x)
        form_x.setSpacing(8)

        self._x_combo = QComboBox()
        for label, _ in _X_VIEWS:
            self._x_combo.addItem(label)
        form_x.addRow("View:", self._x_combo)
        root.addWidget(grp_x)

        root.addStretch()

    # ── SettingsPage interface ────────────────────────────────────────────────

    def populate(self) -> None:
        try:
            from pycsamt.api.control import (
                PYCSAMT_CONTROL as C,
            )

            self._rho_combo.setCurrentIndex(
                next(
                    (i for i, (_, v) in enumerate(_RHO_VIEWS) if v == C.rho.view),
                    0,
                )
            )
            lo, hi = C.phase.range
            self._ph_min.setValue(lo)
            self._ph_max.setValue(hi)
            self._ph_unit_combo.setCurrentIndex(
                next(
                    (i for i, (_, v) in enumerate(_PHASE_UNITS) if v == C.phase.unit),
                    0,
                )
            )
            self._ph_wrap_cb.setChecked(C.phase.wrap)
            self._x_combo.setCurrentIndex(
                next(
                    (i for i, (_, v) in enumerate(_X_VIEWS) if v == C.x.view),
                    0,
                )
            )
        except Exception:
            pass

    def collect(self) -> dict:
        return {
            "view_controls": {
                "rho_view": _RHO_VIEWS[self._rho_combo.currentIndex()][1],
                "phase_range": (self._ph_min.value(), self._ph_max.value()),
                "phase_unit": _PHASE_UNITS[self._ph_unit_combo.currentIndex()][1],
                "phase_wrap": self._ph_wrap_cb.isChecked(),
                "x_view": _X_VIEWS[self._x_combo.currentIndex()][1],
            }
        }

    def reset(self) -> None:
        try:
            from pycsamt.api.control import PYCSAMT_CONTROL

            PYCSAMT_CONTROL.reset()
        except Exception:
            pass
        self.populate()
