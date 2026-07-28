# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
DisplayPage — settings for PYCSAMT_STYLE (MT component colours and line widths).

Each MT component (XY, YX, XX, YY, TE, TM) exposes:
  • Color  — hex string, edited via a colour-picker button
  • Line width — spin box

Correction before/after and raw-data alpha are also exposed.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout,
)

from .base_page import SettingsPage


def _color_btn(color: str, parent) -> QPushButton:
    """Return a button that shows *color* and opens a QColorDialog on click."""
    btn = QPushButton()
    btn.setFixedSize(48, 22)
    btn.setObjectName("ColorPicker")

    def _apply(c: str) -> None:
        btn.setStyleSheet(f"background-color:{c}; border:1px solid #666;")
        btn._color = c

    _apply(color)

    def _pick() -> None:
        c = QColorDialog.getColor(QColor(btn._color), parent, "Pick colour")
        if c.isValid():
            _apply(c.name())

    btn.clicked.connect(_pick)
    return btn


def _lw_spin(value: float) -> QDoubleSpinBox:
    sp = QDoubleSpinBox()
    sp.setRange(0.2, 8.0)
    sp.setSingleStep(0.2)
    sp.setDecimals(1)
    sp.setValue(value)
    sp.setFixedWidth(60)
    return sp


_COMPONENTS = [
    ("XY  (TE mode / off-diag)", "xy"),
    ("YX  (TM mode / off-diag)", "yx"),
    ("XX  (diagonal)", "xx"),
    ("YY  (diagonal)", "yy"),
    ("TE  (apparent)", "te"),
    ("TM  (apparent)", "tm"),
]


class DisplayPage(SettingsPage):
    """Settings page for MT component colours and line widths."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()
        self.populate()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        # ── MT components ─────────────────────────────────────────────────────
        grp_mt = QGroupBox("MT Component Style  (PYCSAMT_STYLE.mt)")
        form_mt = QFormLayout(grp_mt)
        form_mt.setSpacing(6)

        self._comp_widgets: dict[str, tuple[QPushButton, QDoubleSpinBox]] = {}
        for label, key in _COMPONENTS:
            row = QHBoxLayout()
            row.setAlignment(Qt.AlignmentFlag.AlignLeft)
            cbtn = _color_btn("#1e66f5", self)
            lw = _lw_spin(1.5)
            row.addWidget(cbtn)
            row.addWidget(lw)
            form_mt.addRow(label + ":", row)
            self._comp_widgets[key] = (cbtn, lw)

        root.addWidget(grp_mt)

        # ── Correction colours ────────────────────────────────────────────────
        grp_corr = QGroupBox("Correction Style  (PYCSAMT_STYLE.correction)")
        form_corr = QFormLayout(grp_corr)
        form_corr.setSpacing(6)

        self._corr_before_btn = _color_btn("#005a9e", self)
        form_corr.addRow("Before correction:", self._corr_before_btn)
        self._corr_after_btn = _color_btn("#9e0a00", self)
        form_corr.addRow("After correction:", self._corr_after_btn)
        root.addWidget(grp_corr)

        root.addStretch()

    # ── SettingsPage interface ────────────────────────────────────────────────

    def populate(self) -> None:
        try:
            from pycsamt.api.style import PYCSAMT_STYLE as S

            for key, (cbtn, lw_spin) in self._comp_widgets.items():
                comp = getattr(S.mt, key, None)
                if comp is None:
                    continue
                color = getattr(comp, "color", "#1e66f5") or "#1e66f5"
                cbtn.setStyleSheet(f"background-color:{color}; border:1px solid #666;")
                cbtn._color = color
                lw_spin.setValue(getattr(comp, "lw", 1.5) or 1.5)
        except Exception:
            pass

        try:
            from pycsamt.api.style import PYCSAMT_STYLE as S

            bc = getattr(S.correction.before, "color", "#005a9e") or "#005a9e"
            ac = getattr(S.correction.after, "color", "#9e0a00") or "#9e0a00"
            self._corr_before_btn.setStyleSheet(
                f"background-color:{bc}; border:1px solid #666;"
            )
            self._corr_before_btn._color = bc
            self._corr_after_btn.setStyleSheet(
                f"background-color:{ac}; border:1px solid #666;"
            )
            self._corr_after_btn._color = ac
        except Exception:
            pass

    def collect(self) -> dict:
        changes: dict = {}
        try:
            from pycsamt.api.style import PYCSAMT_STYLE as S

            for key, (cbtn, lw_spin) in self._comp_widgets.items():
                comp = getattr(S.mt, key, None)
                if comp is None:
                    continue
                new_color = getattr(cbtn, "_color", None)
                if new_color and new_color != getattr(comp, "color", None):
                    comp.color = new_color
                new_lw = lw_spin.value()
                if abs(new_lw - (getattr(comp, "lw", 1.5) or 1.5)) > 0.05:
                    comp.lw = new_lw
            # correction
            bc = getattr(self._corr_before_btn, "_color", None)
            if bc:
                try:
                    S.correction.before.color = bc
                except Exception:
                    pass
            ac = getattr(self._corr_after_btn, "_color", None)
            if ac:
                try:
                    S.correction.after.color = ac
                except Exception:
                    pass
        except Exception:
            pass
        return changes

    def reset(self) -> None:
        try:
            from pycsamt.api.style import PYCSAMT_STYLE

            PYCSAMT_STYLE.reset()
        except Exception:
            pass
        self.populate()
