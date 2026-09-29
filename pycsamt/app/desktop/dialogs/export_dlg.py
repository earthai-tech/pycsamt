# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ExportDialog — export the current figure to PNG / PDF / SVG / EPS.

Phase 4: QDialog with format selector, DPI spinner, and destination
path browser.  Delegates to ``matplotlib.figure.Figure.savefig``.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

_FORMATS = {
    "PNG  (raster, lossless)": ("png", "*.png"),
    "PDF  (vector)": ("pdf", "*.pdf"),
    "SVG  (vector)": ("svg", "*.svg"),
    "EPS  (vector, legacy)": ("eps", "*.eps"),
    "TIFF (raster, high-res)": ("tiff", "*.tiff"),
}

_DPI_PRESETS = [72, 100, 150, 200, 300, 600]


class ExportDialog(QDialog):
    """
    Modal dialog for exporting a matplotlib Figure to disk.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        The figure to export.
    default_path : str, optional
        Pre-fill the destination path.
    parent : QWidget, optional
    """

    def __init__(
        self,
        figure,
        default_path: str = "",
        figure_factory=None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Export Figure")
        self.setMinimumWidth(480)
        self._figure = figure
        self._figure_factory = figure_factory
        self._build_ui(default_path)

    # ── UI ────────────────────────────────────────────────────────────

    def _build_ui(self, default_path: str) -> None:
        from pycsamt.api.plot import PLOT_CONFIG

        root = QVBoxLayout(self)

        form = QFormLayout()
        form.setSpacing(8)

        # Format
        self._fmt_combo = QComboBox()
        self._fmt_combo.addItems(list(_FORMATS.keys()))
        self._fmt_combo.currentTextChanged.connect(self._update_path_extension)
        form.addRow("Format:", self._fmt_combo)

        # DPI
        self._dpi_spin = QSpinBox()
        self._dpi_spin.setRange(36, 2400)
        self._dpi_spin.setSingleStep(50)
        self._dpi_spin.setValue(max(600, int(PLOT_CONFIG.dpi)))
        self._dpi_spin.setToolTip("Publication preset: at least 600 DPI; adjustable for this export.")
        form.addRow("DPI:", self._dpi_spin)

        # Destination path
        path_row = QHBoxLayout()
        self._path_edit = QLineEdit(
            default_path or str(Path(PLOT_CONFIG.savedir or Path.home()) / "pycsamt_figure.png")
        )
        self._path_edit.setPlaceholderText("Destination file path…")
        btn_browse = QPushButton("Browse…")
        btn_browse.setFixedWidth(70)
        btn_browse.clicked.connect(self._browse)
        path_row.addWidget(self._path_edit)
        path_row.addWidget(btn_browse)
        form.addRow("Save to:", path_row)

        self._transparent = QCheckBox("Transparent background")
        self._transparent.setChecked(bool(PLOT_CONFIG.transparent))
        form.addRow("", self._transparent)
        self._bbox_inches = PLOT_CONFIG.bbox_inches
        preferred = PLOT_CONFIG.resolve_formats()[0]
        for label, (extension, _) in _FORMATS.items():
            if extension == preferred:
                self._fmt_combo.setCurrentText(label)
                break

        root.addLayout(form)

        # OK / Cancel
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_export)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

    # ── Slots ─────────────────────────────────────────────────────────

    def _browse(self) -> None:
        fmt_key = self._fmt_combo.currentText()
        ext_glob = _FORMATS[fmt_key][1]
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Figure",
            self._path_edit.text(),
            f"Figure ({ext_glob});;All files (*)",
        )
        if path:
            self._path_edit.setText(path)

    def _update_path_extension(self, fmt_key: str) -> None:
        """Swap the file extension in the path field when format changes."""
        ext = _FORMATS[fmt_key][0]
        cur = Path(self._path_edit.text())
        self._path_edit.setText(str(cur.with_suffix(f".{ext}")))

    def _on_export(self) -> None:
        fmt_key = self._fmt_combo.currentText()
        ext, _ = _FORMATS[fmt_key]
        dpi = self._dpi_spin.value()
        path = self._path_edit.text().strip()

        if not path:
            QMessageBox.warning(
                self, "Export", "Please specify a destination path."
            )
            return

        try:
            figure = (
                self._figure_factory()
                if self._figure_factory is not None
                else self._figure
            )
            figure.savefig(
                path,
                dpi=dpi,
                format=ext,
                bbox_inches=self._bbox_inches,
                transparent=self._transparent.isChecked(),
                facecolor="none" if self._transparent.isChecked() else "white",
                edgecolor="none",
                pad_inches=0.12,
            )
            if figure is not self._figure:
                from matplotlib.pyplot import close

                close(figure)
            self.accept()
        except Exception as exc:
            QMessageBox.critical(
                self, "Export failed", f"Could not save figure:\n{exc}"
            )
