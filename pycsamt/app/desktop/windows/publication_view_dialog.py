# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PublicationViewDialog — standalone floating window for the MTpy-style
multi-panel publication figure.

Opens on top of the Profile Viewer without overwriting any existing canvas.
The user can export or close independently.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


class PublicationViewDialog(QDialog):
    """Floating publication-quality response figure for a single station.

    Parameters
    ----------
    controller : PlotController
        The active controller that owns the station selection and data.
        ``draw_publication_view`` is called on *controller.fig* so that
        the figure is built lazily only when the dialog opens.
    station_name : str
        Label shown in the window title.
    dark : bool
        Apply dark or light canvas background.
    parent : QWidget, optional
    """

    def __init__(
        self,
        controller,
        station_name: str = "",
        dark: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(
            f"Publication View — {station_name}"
            if station_name
            else "Publication View"
        )
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.resize(1280, 740)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # ── Canvas ────────────────────────────────────────────────────
        self._canvas = MplCanvas(self, toolbar=True)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        layout.addWidget(self._canvas)

        # ── Button row ────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        btn_row.addStretch()

        self._btn_export = QPushButton("⬆  Export…")
        self._btn_export.setToolTip("Save figure to file")
        self._btn_export.clicked.connect(self._on_export)

        self._btn_close = QPushButton("✕  Close")
        self._btn_close.clicked.connect(self.close)

        btn_row.addWidget(self._btn_export)
        btn_row.addWidget(self._btn_close)
        layout.addLayout(btn_row)

        # ── Draw ──────────────────────────────────────────────────────
        controller.draw_publication_view(self._canvas.figure)
        self._canvas.draw()

    # ── Private ───────────────────────────────────────────────────────

    def _on_export(self) -> None:
        try:
            from pycsamt.app.desktop.dialogs.export_dlg import (
                ExportDialog,
            )

            ExportDialog(figure=self._canvas.figure, parent=self).exec()
        except Exception:
            from PySide6.QtWidgets import QFileDialog

            path, _ = QFileDialog.getSaveFileName(
                self,
                "Export figure",
                "publication_view.pdf",
                "PDF (*.pdf);;PNG (*.png);;SVG (*.svg)",
            )
            if path:
                self._canvas.figure.savefig(
                    path, bbox_inches="tight", dpi=200
                )
