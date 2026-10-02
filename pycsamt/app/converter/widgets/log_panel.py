# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.widgets.log_panel
==========================================

Shared log + progress-bar dock, reused by every page. Mirrors the
``QTextEdit`` + ``QProgressBar`` pairing already used by
``FormatConverterDialog`` in the full desktop app.
"""

from __future__ import annotations

from PySide6.QtWidgets import QProgressBar, QTextEdit, QVBoxLayout, QWidget


class LogPanel(QWidget):
    """A read-only log plus a progress bar."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(4)

        self._progress = QProgressBar()
        self._progress.setValue(0)
        root.addWidget(self._progress)

        self._log = QTextEdit()
        self._log.setObjectName("LogText")
        self._log.setReadOnly(True)
        root.addWidget(self._log)

    def append(self, text: str) -> None:
        self._log.append(text)

    def clear(self) -> None:
        self._log.clear()
        self._progress.setValue(0)

    def set_progress(self, current: int, total: int) -> None:
        self._progress.setMaximum(max(total, 1))
        self._progress.setValue(current)

    def reset_progress(self) -> None:
        self._progress.setMaximum(1)
        self._progress.setValue(0)

    def set_indeterminate(self, active: bool) -> None:
        """Busy/indeterminate mode for jobs with no per-item progress."""
        self._progress.setRange(0, 0 if active else 1)
