# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ExportDialog — export a figure to PNG / PDF / SVG / EPS / TIFF.

Given one figure it exports that figure.  Given *sources* (see
:func:`collect_sources`: every figure drawn in the open windows, matplotlib
canvases and PCSF 3-D scenes alike) it adds a **Figure** picker; a 3-D scene
exports as PNG (as shown, camera included) or interactive HTML.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
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

_PLOTLY_FORMATS = {
    "PNG  (as shown, camera included)": ("png", "*.png"),
    "HTML (interactive 3-D)": ("html", "*.html"),
}

_DPI_PRESETS = [72, 100, 150, 200, 300, 600]


# ── figure sources ────────────────────────────────────────────────────────
@dataclass
class FigureSource:
    """One exportable figure found in an open window."""

    label: str
    figure: Any = None  # matplotlib Figure
    plotly: Any = None  # PlotlyView holding a rendered scene
    window: Any = None
    shown: bool = False  # on screen now (not behind another tab)

    @property
    def kind(self) -> str:
        return "plotly" if self.plotly is not None else "mpl"


def _figure_title(fig) -> str:
    sup = getattr(fig, "_suptitle", None)
    if sup is not None and sup.get_text().strip():
        return sup.get_text().strip()
    for ax in fig.axes:
        t = ax.get_title().strip()
        if t:
            return t
    return ""


def _tab_path(widget, stop) -> list[str]:
    """Tab titles between *widget* and its window *stop*."""
    from PySide6.QtWidgets import QStackedWidget, QTabWidget

    out, w = [], widget
    while w is not None and w is not stop:
        parent = w.parentWidget()
        if isinstance(parent, QStackedWidget) and isinstance(
                parent.parentWidget(), QTabWidget):
            tabs = parent.parentWidget()
            i = tabs.indexOf(w)
            if i >= 0 and tabs.tabText(i).strip():
                out.insert(0, tabs.tabText(i).replace("&", "").strip())
        w = parent
    return out


def _label(win, widget, title: str) -> str:
    name = win.windowTitle().strip() or type(win).__name__
    parts = [re.sub(r"^pycsamt\s*[—–-]\s*", "", name, flags=re.I)]
    parts += _tab_path(widget, win)
    if title:
        parts.append(re.sub(r"\$[^$]*\$", "", title).strip() or title)
    text = "  ▸  ".join(dict.fromkeys(p for p in parts if p))
    return text if len(text) <= 90 else text[:87] + "…"


def collect_sources(windows, *, first=None) -> list[FigureSource]:
    """Every non-empty figure drawn in the visible *windows*.

    Canvases are found by walking the widgets (whatever wrapper holds them:
    ``MplCanvas``, a result view, a raw ``FigureCanvasQTAgg``), so no window
    has to register its figures.  Placeholder / message-only figures are
    skipped.  Figures of *first* (the window the user worked in last) come
    first, then the ones on screen.
    """
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

    from pycsamt.app.desktop.controllers.correction_views import (
        figure_blank_reason,
    )
    from pycsamt.app.desktop.widgets.plotly_view import PlotlyView

    out: list[FigureSource] = []
    seen: set[int] = set()
    for win in windows:
        if win is None or not win.isVisible():
            continue
        for canvas in win.findChildren(FigureCanvasQTAgg):
            fig = getattr(canvas, "figure", None)
            if fig is None or id(fig) in seen:
                continue
            try:
                if figure_blank_reason(fig) is not None:
                    continue
            except Exception:
                continue
            seen.add(id(fig))
            out.append(FigureSource(_label(win, canvas, _figure_title(fig)),
                                    figure=fig, window=win,
                                    shown=canvas.isVisible()))
        for view in win.findChildren(PlotlyView):
            if view.plot_ready and getattr(view, "_figure", None) is not None:
                out.append(FigureSource(_label(win, view, "3-D scene"),
                                        plotly=view, window=win,
                                        shown=view.isVisible()))
    out.sort(key=lambda src: (src.window is not first, not src.shown))
    return out


def _slug(text: str) -> str:
    parts = [p.strip() for p in text.split("▸") if p.strip()]
    text = "_".join(dict.fromkeys([parts[0], parts[-1]])) if parts else text
    return (re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_").lower()[:48]
            or "figure")


class ExportDialog(QDialog):
    """
    Modal dialog for exporting a figure to disk.

    Parameters
    ----------
    figure : matplotlib.figure.Figure, optional
        The figure to export (the first of *sources* when omitted).
    default_path : str, optional
        Pre-fill the destination path.
    figure_factory : callable, optional
        Builds the figure to save (e.g. at publication size).
    parent : QWidget, optional
    sources : list of FigureSource, optional
        Figures to choose from (see :func:`collect_sources`).
    """

    def __init__(
        self,
        figure=None,
        default_path: str = "",
        figure_factory=None,
        parent: QWidget | None = None,
        sources: list[FigureSource] | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Export Figure")
        self.setMinimumWidth(520)
        self._sources = list(sources or [])
        if figure is None and self._sources:
            figure = self._sources[0].figure
        self._figure = figure
        self._figure_factory = figure_factory
        self._plotly = None
        self._build_ui(default_path)
        if self._sources:
            self._on_source(0)
        self._update_size()

    # ── UI ────────────────────────────────────────────────────────────

    def _build_ui(self, default_path: str) -> None:
        from pycsamt.api.plot import PLOT_CONFIG

        root = QVBoxLayout(self)

        form = QFormLayout()
        form.setSpacing(8)

        self._source_combo = QComboBox()
        for i, src in enumerate(self._sources):
            self._source_combo.addItem(src.label)
            self._source_combo.setItemData(i, src.label,
                                           Qt.ItemDataRole.ToolTipRole)
        self._source_combo.currentIndexChanged.connect(self._on_source)
        if self._sources:
            form.addRow("Figure:", self._source_combo)

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
        self._dpi_spin.valueChanged.connect(self._update_size)
        form.addRow("DPI:", self._dpi_spin)
        self._size_lbl = QLabel("")
        self._size_lbl.setObjectName("InfoLabel")
        form.addRow("Size:", self._size_lbl)

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
        self._buttons = buttons
        root.addWidget(buttons)

    # ── Slots ─────────────────────────────────────────────────────────

    def _formats(self) -> dict:
        return _PLOTLY_FORMATS if self._plotly is not None else _FORMATS

    def _on_source(self, index: int) -> None:
        if not 0 <= index < len(self._sources):
            return
        src = self._sources[index]
        self._figure, self._plotly = src.figure, src.plotly
        self._figure_factory = None
        previous = self._fmt_combo.currentText()
        self._fmt_combo.blockSignals(True)
        self._fmt_combo.clear()
        self._fmt_combo.addItems(list(self._formats()))
        if previous in self._formats():
            self._fmt_combo.setCurrentText(previous)
        self._fmt_combo.blockSignals(False)
        mpl = self._plotly is None
        self._dpi_spin.setEnabled(mpl)
        self._transparent.setEnabled(mpl)
        cur = Path(self._path_edit.text() or "figure.png")
        ext = self._formats()[self._fmt_combo.currentText()][0]
        self._path_edit.setText(str(cur.with_name(
            f"pycsamt_{_slug(src.label)}.{ext}")))
        self._update_size()

    def _update_size(self, *_a) -> None:
        if self._plotly is not None:
            self._size_lbl.setText("As on screen (PNG at 2× resolution)")
            return
        fig = self._figure
        if fig is None or not hasattr(fig, "get_size_inches"):
            self._size_lbl.setText("")
            return
        w, h = fig.get_size_inches()
        dpi = self._dpi_spin.value()
        self._size_lbl.setText(f"{w:.1f} × {h:.1f} in  →  "
                               f"{int(w * dpi)} × {int(h * dpi)} px")

    def _browse(self) -> None:
        fmt_key = self._fmt_combo.currentText()
        ext_glob = self._formats()[fmt_key][1]
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
        if fmt_key not in self._formats():
            return
        ext = self._formats()[fmt_key][0]
        cur = Path(self._path_edit.text())
        self._path_edit.setText(str(cur.with_suffix(f".{ext}")))

    def _on_export(self) -> None:
        fmt_key = self._fmt_combo.currentText()
        ext, _ = self._formats()[fmt_key]
        dpi = self._dpi_spin.value()
        path = self._path_edit.text().strip()

        if not path:
            QMessageBox.warning(
                self, "Export", "Please specify a destination path."
            )
            return
        if self._plotly is not None:
            self._export_plotly(path, ext)
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

    def _export_plotly(self, path: str, ext: str) -> None:
        if ext == "html":
            try:
                self._plotly.export_html(path)
            except Exception as exc:
                QMessageBox.critical(self, "Export failed",
                                     f"Could not save the scene:\n{exc}")
                return
            self.accept()
            return
        self._buttons.setEnabled(False)
        self._size_lbl.setText("Rendering the scene…")

        def done(saved, error) -> None:
            self._buttons.setEnabled(True)
            if saved:
                self.accept()
            else:
                self._update_size()
                QMessageBox.critical(self, "Export failed",
                                     f"Could not save the scene:\n{error}")

        self._plotly.snapshot_png(path, done)


__all__ = ["ExportDialog", "FigureSource", "collect_sources"]
