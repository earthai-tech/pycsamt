"""Settings page for v2.6 contour and mesh rendering defaults."""

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


class RenderingPage(SettingsPage):
    """Configure shared contour and mesh-review styles."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        root = QVBoxLayout(self)

        contours = QGroupBox("Contour Overlays  (PYCSAMT_CONTOUR.default)")
        contour_form = QFormLayout(contours)
        self._contour_enabled = QCheckBox("Draw contour overlays")
        contour_form.addRow("", self._contour_enabled)
        self._levels = QSpinBox()
        self._levels.setRange(2, 100)
        contour_form.addRow("Number of levels:", self._levels)
        self._contour_width = QDoubleSpinBox()
        self._contour_width.setRange(0.05, 10.0)
        self._contour_width.setSingleStep(0.1)
        contour_form.addRow("Line width:", self._contour_width)
        self._contour_alpha = QDoubleSpinBox()
        self._contour_alpha.setRange(0.0, 1.0)
        self._contour_alpha.setSingleStep(0.05)
        contour_form.addRow("Opacity:", self._contour_alpha)
        self._labels = QCheckBox("Label contour lines")
        contour_form.addRow("", self._labels)
        root.addWidget(contours)

        mesh = QGroupBox("Mesh Review Style  (PYCSAMT_MESH.review)")
        mesh_form = QFormLayout(mesh)
        self._mesh_edges = QCheckBox("Show cell boundaries")
        mesh_form.addRow("", self._mesh_edges)
        self._mesh_width = QDoubleSpinBox()
        self._mesh_width.setRange(0.05, 10.0)
        self._mesh_width.setSingleStep(0.05)
        mesh_form.addRow("Boundary width:", self._mesh_width)
        self._mesh_alpha = QDoubleSpinBox()
        self._mesh_alpha.setRange(0.0, 1.0)
        self._mesh_alpha.setSingleStep(0.05)
        mesh_form.addRow("Boundary opacity:", self._mesh_alpha)
        self._mesh_style = QComboBox()
        self._mesh_style.addItems(["-", "--", ":", "-."])
        mesh_form.addRow("Boundary style:", self._mesh_style)
        root.addWidget(mesh)
        root.addStretch()
        self.populate()

    def populate(self) -> None:
        try:
            from pycsamt.api.contour import PYCSAMT_CONTOUR as C

            style = C.default
            self._contour_enabled.setChecked(bool(style.enabled))
            self._levels.setValue(
                style.levels if isinstance(style.levels, int) else 7
            )
            self._contour_width.setValue(float(style.linewidths))
            self._contour_alpha.setValue(float(style.alpha))
            self._labels.setChecked(bool(style.labels))
        except Exception:
            pass
        try:
            from pycsamt.api.mesh import PYCSAMT_MESH as M

            edge = M.review.edge
            self._mesh_edges.setChecked(bool(edge.show))
            self._mesh_width.setValue(float(edge.linewidth))
            self._mesh_alpha.setValue(float(edge.alpha))
            self._mesh_style.setCurrentText(str(edge.linestyle))
        except Exception:
            pass

    def collect(self) -> dict:
        return {
            "contour": {
                "enabled": self._contour_enabled.isChecked(),
                "levels": self._levels.value(),
                "linewidths": self._contour_width.value(),
                "alpha": self._contour_alpha.value(),
                "labels": self._labels.isChecked(),
            },
            "mesh": {
                "review__edge__show": self._mesh_edges.isChecked(),
                "review__edge__linewidth": self._mesh_width.value(),
                "review__edge__alpha": self._mesh_alpha.value(),
                "review__edge__linestyle": self._mesh_style.currentText(),
            },
        }

    def reset(self) -> None:
        try:
            from pycsamt.api.contour import PYCSAMT_CONTOUR
            from pycsamt.api.mesh import PYCSAMT_MESH

            PYCSAMT_CONTOUR.reset()
            PYCSAMT_MESH.reset()
        except Exception:
            pass
        self.populate()
