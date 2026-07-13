# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LoaderWorker — QThread worker for EDI / AVG / J file loading.

Phase 2: Runs DataController.load() in a background thread so the UI
stays responsive.  Signals:

    progress(int)          0–100 for the progress bar
    finished(object)       Sites collection on success
    error(str)             human-readable error message on failure
"""

from __future__ import annotations

from PySide6.QtCore import QThread, Signal

from pycsamt.app.desktop.controllers.data_controller import (
    DataController,
)


class LoaderWorker(QThread):
    """Background EDI loader thread."""

    progress = Signal(int)
    finished = Signal(object)  # Sites
    error = Signal(str)

    def __init__(self, paths: list[str], parent=None) -> None:
        super().__init__(parent)
        self._paths = list(paths)
        self._ctrl = DataController(progress_callback=self.progress.emit)

    def run(self) -> None:
        try:
            sites = self._ctrl.load(self._paths)
            self.finished.emit(sites)
        except Exception as exc:
            self.error.emit(str(exc))

    @property
    def data_controller(self) -> DataController:
        """Return the DataController after loading completes."""
        return self._ctrl
