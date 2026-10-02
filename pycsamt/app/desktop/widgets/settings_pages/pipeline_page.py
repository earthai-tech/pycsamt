# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Settings page for the processing-pipeline engine (``PYCSAMT_PIPE``).

Covers what matters to a desktop run: where results go, how a failing
step is handled, figure/report output, and the step cache / run history
locations. Console-only knobs (``show_progress``, ``progress_style``,
``repr_width``) are deliberately not exposed: the desktop pipeline shows
its own live progress.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.compact_button import compact_button

from .base_page import SettingsPage

# (value, label shown in the combo)
_ERROR_POLICIES = [
    ("warn", "Warn and continue with the previous data"),
    ("skip", "Skip silently and continue"),
    ("raise", "Stop the run at the first error"),
]
_PLOT_FORMATS = ["png", "pdf", "svg"]
_REPORT_FORMATS = [
    ("html", "HTML report"),
    ("txt", "Plain-text report"),
    ("dashboard", "Dashboard (KPI tiles + charts)"),
]


class _PathEdit(QWidget):
    """Line edit + browse button; an empty value means "use the default"."""

    def __init__(self, placeholder: str, pick_file: bool = False) -> None:
        super().__init__()
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(4)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.setClearButtonEnabled(True)
        browse = QPushButton("…")
        compact_button(browse, 28)
        browse.setToolTip("Browse…")
        browse.clicked.connect(self._browse)
        h.addWidget(self.edit, 1)
        h.addWidget(browse)
        self._pick_file = pick_file

    def _browse(self) -> None:
        if self._pick_file:
            path, _ = QFileDialog.getSaveFileName(
                self, "Run history file", self.edit.text(),
                "JSON Lines (*.jsonl);;All files (*)",
            )
        else:
            path = QFileDialog.getExistingDirectory(
                self, "Select folder", self.edit.text()
            )
        if path:
            self.edit.setText(path)

    def text(self) -> str:
        return self.edit.text().strip()

    def setText(self, value) -> None:  # noqa: N802 (Qt naming)
        self.edit.setText("" if value is None else str(value))


class PipelinePage(SettingsPage):
    """Configure ``PYCSAMT_PIPE`` for desktop and scripted pipeline runs."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        root = QVBoxLayout(self)

        # ── Output ────────────────────────────────────────────────────
        out = QGroupBox("Output  (PYCSAMT_PIPE)")
        f = QFormLayout(out)
        self._output_root = _PathEdit("pipe_results")
        self._output_root.setToolTip(
            "Default results folder when a run does not choose one."
        )
        f.addRow("Results folder:", self._output_root)
        self._processed = QLineEdit()
        self._processed.setToolTip("Sub-folder for the processed EDI files.")
        f.addRow("Processed EDI sub-folder:", self._processed)
        self._plots = QLineEdit()
        self._plots.setToolTip("Sub-folder for QC figures (one per step).")
        f.addRow("Figures sub-folder:", self._plots)
        self._intermediate = QCheckBox("Also save EDI files after every step")
        f.addRow("", self._intermediate)
        root.addWidget(out)

        # ── Errors ────────────────────────────────────────────────────
        err = QGroupBox("When a step fails")
        fe = QFormLayout(err)
        self._on_error = QComboBox()
        for value, label in _ERROR_POLICIES:
            self._on_error.addItem(label, userData=value)
        fe.addRow("Policy:", self._on_error)
        root.addWidget(err)

        # ── Figures & reports ─────────────────────────────────────────
        rep = QGroupBox("Figures and reports")
        fr = QFormLayout(rep)
        self._dpi = QSpinBox()
        self._dpi.setRange(50, 1200)
        self._dpi.setSingleStep(50)
        fr.addRow("Figure resolution (DPI):", self._dpi)
        self._fmt = QComboBox()
        self._fmt.addItems(_PLOT_FORMATS)
        fr.addRow("Figure format:", self._fmt)
        self._reports: dict[str, QCheckBox] = {}
        box = QWidget()
        bv = QVBoxLayout(box)
        bv.setContentsMargins(0, 0, 0, 0)
        bv.setSpacing(2)
        for value, label in _REPORT_FORMATS:
            cb = QCheckBox(label)
            self._reports[value] = cb
            bv.addWidget(cb)
        fr.addRow("Reports:", box)
        root.addWidget(rep)

        # ── Cache & history ───────────────────────────────────────────
        ch = QGroupBox("Step cache and run history")
        fc = QFormLayout(ch)
        self._cache_root = _PathEdit("~/.pycsamt/pipeline_cache  (default)")
        self._cache_root.setToolTip(
            "Where cached step outputs are kept, so an interrupted or "
            "repeated run can resume instead of recomputing."
        )
        fc.addRow("Cache folder:", self._cache_root)
        self._history = _PathEdit(
            "~/.pycsamt/pipeline_history.jsonl  (default)", pick_file=True
        )
        self._history.setToolTip("One line per recorded run, for comparison.")
        fc.addRow("History file:", self._history)
        note = QLabel("Leave a location empty to use the default shown.")
        note.setObjectName("InfoLabel")
        note.setWordWrap(True)
        fc.addRow(note)
        root.addWidget(ch)

        root.addStretch()
        self.populate()

    # ── SettingsPage interface ────────────────────────────────────────

    def populate(self) -> None:
        try:
            from pycsamt.api.pipe import PYCSAMT_PIPE as P
        except Exception:
            return
        self._output_root.setText(P.output_root)
        self._processed.setText(P.processed_subdir)
        self._plots.setText(P.plots_subdir)
        self._intermediate.setChecked(bool(P.save_intermediate))
        idx = self._on_error.findData(P.on_step_error)
        self._on_error.setCurrentIndex(max(idx, 0))
        self._dpi.setValue(int(P.plot_dpi))
        self._fmt.setCurrentText(str(P.plot_fmt))
        formats = set(P.report_formats or ())
        for value, cb in self._reports.items():
            cb.setChecked(value in formats)
        self._cache_root.setText(P.cache_root)
        self._history.setText(P.history_path)

    def collect(self) -> dict:
        return {
            "pipe": {
                "output_root": self._output_root.text() or "pipe_results",
                "processed_subdir": self._processed.text().strip()
                or "processed",
                "plots_subdir": self._plots.text().strip() or "plots",
                "save_intermediate": self._intermediate.isChecked(),
                "on_step_error": self._on_error.currentData(),
                "plot_dpi": self._dpi.value(),
                "plot_fmt": self._fmt.currentText(),
                "report_formats": tuple(
                    v for v, cb in self._reports.items() if cb.isChecked()
                ),
                # empty -> None = the library's lazily resolved default
                "cache_root": self._cache_root.text() or None,
                "history_path": self._history.text() or None,
            }
        }

    def reset(self) -> None:
        try:
            from pycsamt.api.pipe import PYCSAMT_PIPE

            PYCSAMT_PIPE.reset()
        except Exception:
            pass
        self.populate()
