"""Settings page for package-wide survey site ordering."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QVBoxLayout,
)

from .base_page import SettingsPage


class OrderingPage(SettingsPage):
    """Configure the v2.6 automatic site-ordering policy."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        root = QVBoxLayout(self)
        note = QLabel(
            "These defaults are used by loaders and processors when no explicit "
            "ordering is requested. They take effect on the next data load."
        )
        note.setWordWrap(True)
        root.addWidget(note)
        group = QGroupBox("Survey Site Ordering  (PYCSAMT_ORDERING)")
        form = QFormLayout(group)
        self._mode = QComboBox()
        self._mode.addItems(
            ["auto", "chainage", "input", "station", "latitude", "longitude"]
        )
        form.addRow("Ordering mode:", self._mode)
        self._linearity = self._ratio_spin()
        form.addRow("Minimum line linearity:", self._linearity)
        self._cross_track = self._ratio_spin()
        form.addRow("Maximum cross-track ratio:", self._cross_track)
        self._coord_fraction = self._ratio_spin()
        form.addRow("Minimum coordinate coverage:", self._coord_fraction)
        root.addWidget(group)
        root.addStretch()
        self.populate()

    @staticmethod
    def _ratio_spin() -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setRange(0.0, 1.0)
        spin.setDecimals(2)
        spin.setSingleStep(0.05)
        return spin

    def populate(self) -> None:
        try:
            from pycsamt.api.ordering import PYCSAMT_ORDERING as O

            self._mode.setCurrentText(O.mode)
            self._linearity.setValue(O.min_linearity)
            self._cross_track.setValue(O.max_cross_track_ratio)
            self._coord_fraction.setValue(O.min_coordinate_fraction)
        except Exception:
            pass

    def collect(self) -> dict:
        return {
            "ordering": {
                "mode": self._mode.currentText(),
                "min_linearity": self._linearity.value(),
                "max_cross_track_ratio": self._cross_track.value(),
                "min_coordinate_fraction": self._coord_fraction.value(),
            }
        }

    def reset(self) -> None:
        try:
            from pycsamt.api.ordering import PYCSAMT_ORDERING

            PYCSAMT_ORDERING.reset()
        except Exception:
            pass
        self.populate()
