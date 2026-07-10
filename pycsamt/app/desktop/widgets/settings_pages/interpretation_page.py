# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InterpretationPage — settings for PYCSAMT_INTERP (hydrogeophysical style).

Controls
--------
Section style
  • Colormap (resistivity section)
  • Water-table marker style
  • Contour alpha

Profile style
  • Colormap
  • Reference line style
"""
from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QVBoxLayout,
)

from .base_page import SettingsPage

_CMAPS = [
    "viridis", "plasma", "magma", "inferno",
    "RdBu_r",  "seismic", "jet",  "turbo",
    "coolwarm", "bwr",    "Spectral_r",
]

_WT_MARKERS = [
    ("None (no marker)",  "none"),
    ("Dashed line  --",   "--"),
    ("Solid line   —",    "-"),
    ("Dotted line  ·",    ":"),
]


class InterpretationPage(SettingsPage):
    """Settings page for the hydrogeophysical interpretation style."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()
        self.populate()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        # ── Section style ─────────────────────────────────────────────────────
        grp_sec = QGroupBox("Section Style  (PYCSAMT_INTERP.pseudosection)")
        form_sec = QFormLayout(grp_sec)
        form_sec.setSpacing(8)

        self._sec_cmap = QComboBox()
        self._sec_cmap.addItems(_CMAPS)
        form_sec.addRow("Colormap:", self._sec_cmap)

        self._wt_combo = QComboBox()
        for label, _ in _WT_MARKERS:
            self._wt_combo.addItem(label)
        form_sec.addRow("Water-table marker:", self._wt_combo)

        self._alpha_spin = QDoubleSpinBox()
        self._alpha_spin.setRange(0.0, 1.0)
        self._alpha_spin.setSingleStep(0.05)
        self._alpha_spin.setDecimals(2)
        form_sec.addRow("Contour alpha:", self._alpha_spin)

        root.addWidget(grp_sec)

        # ── Profile style ─────────────────────────────────────────────────────
        grp_prof = QGroupBox("Profile Style  (PYCSAMT_INTERP.profile)")
        form_prof = QFormLayout(grp_prof)
        form_prof.setSpacing(8)

        self._prof_cmap = QComboBox()
        self._prof_cmap.addItems(_CMAPS)
        form_prof.addRow("Colormap:", self._prof_cmap)

        root.addWidget(grp_prof)
        root.addStretch()

    # ── SettingsPage interface ────────────────────────────────────────────────

    def populate(self) -> None:
        try:
            from pycsamt.api.interp import PYCSAMT_INTERP as I
            sec  = I.pseudosection
            prof = I.profile

            sec_cmap = getattr(sec, "cmap", None) or "viridis"
            if sec_cmap in _CMAPS:
                self._sec_cmap.setCurrentText(sec_cmap)

            wt_ls = getattr(sec, "wt_linestyle", "--") or "--"
            idx = next((i for i, (_, v) in enumerate(_WT_MARKERS) if v == wt_ls), 1)
            self._wt_combo.setCurrentIndex(idx)

            alpha = getattr(sec, "alpha", 0.85) or 0.85
            self._alpha_spin.setValue(alpha)

            prof_cmap = getattr(prof, "cmap", None) or "viridis"
            if prof_cmap in _CMAPS:
                self._prof_cmap.setCurrentText(prof_cmap)
        except Exception:
            self._alpha_spin.setValue(0.85)

    def collect(self) -> dict:
        try:
            from pycsamt.api.interp import PYCSAMT_INTERP as I
            sec_cmap  = self._sec_cmap.currentText()
            wt_ls     = _WT_MARKERS[self._wt_combo.currentIndex()][1]
            alpha     = self._alpha_spin.value()
            prof_cmap = self._prof_cmap.currentText()

            try:
                sec = I.pseudosection
                if hasattr(sec, "cmap"):
                    sec.cmap = sec_cmap
                if hasattr(sec, "wt_linestyle"):
                    sec.wt_linestyle = wt_ls
                if hasattr(sec, "alpha"):
                    sec.alpha = alpha
            except Exception:
                pass

            try:
                prof = I.profile
                if hasattr(prof, "cmap"):
                    prof.cmap = prof_cmap
            except Exception:
                pass
        except Exception:
            pass
        return {}

    def reset(self) -> None:
        try:
            from pycsamt.api.interp import PYCSAMT_INTERP
            PYCSAMT_INTERP.reset()
        except Exception:
            pass
        self.populate()
