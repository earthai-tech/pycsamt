"""Settings page for package-wide plot export defaults."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
)

from .base_page import SettingsPage


class OutputPage(SettingsPage):
    """Configure :data:`pycsamt.api.plot.PLOT_CONFIG`."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        root = QVBoxLayout(self)
        group = QGroupBox("Plot Export Defaults  (PLOT_CONFIG)")
        form = QFormLayout(group)

        self._fmt = QLineEdit()
        self._fmt.setPlaceholderText("png, +svg, +pdf")
        self._fmt.setToolTip(
            "A format such as png, or additive formats such as +svg,+pdf."
        )
        form.addRow("Output format(s):", self._fmt)

        self._dpi = QSpinBox()
        self._dpi.setRange(50, 2400)
        self._dpi.setSingleStep(25)
        form.addRow("Resolution (DPI):", self._dpi)

        self._bbox = QComboBox()
        self._bbox.addItems(["tight", "standard"])
        form.addRow("Bounding box:", self._bbox)

        self._savedir = QLineEdit()
        self._savedir.setPlaceholderText("Use the selected export folder")
        form.addRow("Default output folder:", self._savedir)

        self._transparent = QCheckBox("Use transparent figure backgrounds")
        form.addRow("", self._transparent)
        self._close = QCheckBox("Close figures after saving")
        form.addRow("", self._close)
        root.addWidget(group)
        root.addStretch()
        self.populate()

    def populate(self) -> None:
        try:
            from pycsamt.api.plot import PLOT_CONFIG as P

            fmt = P.fmt if isinstance(P.fmt, str) else ",".join(P.fmt)
            self._fmt.setText(fmt)
            self._dpi.setValue(int(P.dpi))
            self._bbox.setCurrentText(
                "tight" if P.bbox_inches == "tight" else "standard"
            )
            self._savedir.setText("" if P.savedir is None else str(P.savedir))
            self._transparent.setChecked(bool(P.transparent))
            self._close.setChecked(bool(P.close_after_save))
        except Exception:
            self._fmt.setText("png")
            self._dpi.setValue(150)

    def collect(self) -> dict:
        return {
            "plot": {
                "fmt": self._fmt.text().strip() or "png",
                "dpi": self._dpi.value(),
                "bbox_inches": (
                    "tight" if self._bbox.currentText() == "tight" else None
                ),
                "savedir": self._savedir.text().strip() or None,
                "transparent": self._transparent.isChecked(),
                "close_after_save": self._close.isChecked(),
            }
        }

    def reset(self) -> None:
        try:
            from pycsamt.api.plot import PLOT_CONFIG

            PLOT_CONFIG.reset()
        except Exception:
            pass
        self.populate()
