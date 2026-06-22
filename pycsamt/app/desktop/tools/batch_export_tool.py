# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
BatchExportDialog — export every open canvas figure to PNG / PDF / SVG.

Strategy
--------
The dialog receives a list of (label, figure) tuples from MainWindow.
It saves each figure to ``<output_dir>/<label>.<ext>`` with the chosen
DPI.  Progress is shown via a QProgressBar and a log widget.

Usage in MainWindow
-------------------
    from pycsamt.app.desktop.tools.batch_export_tool import BatchExportDialog
    figures = self._collect_figures()   # see _collect_figures() in main_window
    dlg = BatchExportDialog(figures, parent=self)
    dlg.exec()
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

_FORMATS = ["PNG", "PDF", "SVG", "EPS", "TIFF"]


# ── Background worker ─────────────────────────────────────────────────────────

class _ExportWorker(QThread):
    progress = Signal(int, int, str)
    done     = Signal(int, str)
    error    = Signal(str)

    def __init__(
        self,
        figures: List[Tuple[str, object]],
        out_dir: Path,
        fmt: str,
        dpi: int,
        transparent: bool,
    ):
        super().__init__()
        self._figures     = figures
        self._out_dir     = out_dir
        self._fmt         = fmt.lower()
        self._dpi         = dpi
        self._transparent = transparent

    def run(self):
        self._out_dir.mkdir(parents=True, exist_ok=True)
        total   = len(self._figures)
        saved   = 0
        for idx, (label, fig) in enumerate(self._figures):
            safe = "".join(c if c.isalnum() or c in "-_ " else "_" for c in label)
            safe = safe.strip().replace(" ", "_")
            out_path = self._out_dir / f"{safe}.{self._fmt}"
            try:
                fig.savefig(
                    str(out_path),
                    dpi=self._dpi,
                    bbox_inches="tight",
                    transparent=self._transparent,
                )
                saved += 1
                self.progress.emit(idx + 1, total, label)
            except Exception as exc:
                self.progress.emit(idx + 1, total, f"{label}  ERROR: {exc}")
        self.done.emit(saved, str(self._out_dir))


# ── Dialog ────────────────────────────────────────────────────────────────────

class BatchExportDialog(QDialog):
    """
    Batch-export all open canvas figures to a folder.

    Parameters
    ----------
    figures : list of (label, Figure)
        Collected from MainWindow via _collect_figures().
    parent : QWidget, optional
    """

    def __init__(
        self,
        figures: List[Tuple[str, object]],
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Batch Plot Export")
        self.setMinimumSize(500, 420)
        self._figures = figures
        self._worker  = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(10)

        # Source summary
        grp_src = QGroupBox("Source")
        form_src = QFormLayout(grp_src)
        n = len(self._figures)
        self._src_lbl = QLabel(f"{n} canvas figure{'s' if n != 1 else ''} found")
        form_src.addRow("Figures:", self._src_lbl)
        if self._figures:
            names = ", ".join(lbl for lbl, _ in self._figures[:5])
            if n > 5:
                names += f", … (+{n - 5} more)"
            form_src.addRow("", QLabel(f"<small>{names}</small>"))
        root.addWidget(grp_src)

        # Output options
        grp_out = QGroupBox("Export Options")
        form_out = QFormLayout(grp_out)
        form_out.setSpacing(8)

        self._fmt_combo = QComboBox()
        self._fmt_combo.addItems(_FORMATS)
        form_out.addRow("Format:", self._fmt_combo)

        self._dpi_spin = QSpinBox()
        self._dpi_spin.setRange(72, 600)
        self._dpi_spin.setValue(150)
        self._dpi_spin.setSuffix(" dpi")
        form_out.addRow("Resolution:", self._dpi_spin)

        self._transp_cb = QCheckBox("Transparent background")
        form_out.addRow("", self._transp_cb)

        dir_row = QHBoxLayout()
        self._dir_edit = QLineEdit(str(Path.home() / "pycsamt_figures"))
        btn_browse = QPushButton("Browse…")
        btn_browse.setFixedWidth(70)
        btn_browse.clicked.connect(self._pick_dir)
        dir_row.addWidget(self._dir_edit)
        dir_row.addWidget(btn_browse)
        form_out.addRow("Output folder:", dir_row)

        root.addWidget(grp_out)

        # Progress
        self._progress = QProgressBar()
        self._progress.setValue(0)
        root.addWidget(self._progress)

        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setFixedHeight(100)
        root.addWidget(self._log)

        # Buttons
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Export All")
        self._run_btn.clicked.connect(self._on_run)
        if not self._figures:
            self._run_btn.setEnabled(False)
        btn_row.addWidget(self._run_btn)
        btn_row.addStretch()
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        btn_row.addWidget(box)
        root.addLayout(btn_row)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _pick_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(
            self, "Select output folder", self._dir_edit.text()
        )
        if d:
            self._dir_edit.setText(d)

    def _on_run(self) -> None:
        out_dir = Path(self._dir_edit.text().strip())
        fmt     = self._fmt_combo.currentText()
        dpi     = self._dpi_spin.value()
        transp  = self._transp_cb.isChecked()

        self._run_btn.setEnabled(False)
        self._progress.setValue(0)
        self._log.clear()
        self._log.append(
            f"Exporting {len(self._figures)} figures → {out_dir}  [{fmt}, {dpi} dpi]"
        )

        self._worker = _ExportWorker(self._figures, out_dir, fmt, dpi, transp)
        self._worker.progress.connect(self._on_progress)
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, cur: int, total: int, label: str) -> None:
        self._progress.setMaximum(total)
        self._progress.setValue(cur)
        self._log.append(f"  [{cur}/{total}]  {label}")

    def _on_done(self, saved: int, folder: str) -> None:
        self._run_btn.setEnabled(True)
        self._log.append(f"\nDone — {saved} figure(s) saved to {folder}")

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._log.append(f"\nError: {msg}")
