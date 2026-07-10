# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PipelineWorker — QThread that executes pipeline steps sequentially,
emitting granular signals so the UI can update in real-time.
"""
from __future__ import annotations

import datetime

from PySide6.QtCore import QThread, Signal


class PipelineWorker(QThread):
    """
    Executes a list of pipeline steps on a background thread.

    Signals
    -------
    step_started(int)           Step index just started.
    step_done(int, object)      Step index, output Sites.
    step_error(int, str)        Step index, error message.
    step_skipped(int)           Step was marked SKIPPED.
    log_line(str)               Timestamped log message.
    progress(int)               Overall 0–100 progress.
    all_done()                  All steps finished (or interrupted).
    """

    step_started  = Signal(int)
    step_done     = Signal(int, object)
    step_error    = Signal(int, str)
    step_skipped  = Signal(int)
    log_line      = Signal(str)
    progress      = Signal(int)
    all_done      = Signal()

    def __init__(
        self,
        controller,
        step_ids: list[int],
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._ctrl     = controller
        self._step_ids = step_ids  # list of step indices to run

    # ── QThread entry-point ───────────────────────────────────────────────────

    def run(self) -> None:
        total = len(self._step_ids)
        for done_count, step_id in enumerate(self._step_ids):
            if self.isInterruptionRequested():
                self._log("Pipeline interrupted by user.")
                break

            step = self._ctrl.steps[step_id]
            from pycsamt.app.desktop.controllers.pipeline_controller import (
                StepStatus,
            )

            # Skip steps explicitly marked as SKIPPED by the user
            if step.status == StepStatus.SKIPPED:
                self._log(f"Step {step_id + 1} ({step.name}) — skipped.")
                self.step_skipped.emit(step_id)
                self.progress.emit(int((done_count + 1) / total * 100))
                continue

            self._log(f"Step {step_id + 1}/{total}: {step.name} …")
            self.step_started.emit(step_id)

            try:
                result = self._ctrl.execute_step(step_id, log_cb=self._log)
                self.step_done.emit(step_id, result)
                self._log(
                    f"  ✓  {step.name} done in {step.elapsed_s:.2f}s "
                    f"— {step.result_info}"
                )
            except Exception as exc:
                self.step_error.emit(step_id, str(exc))
                self._log(f"  ✗  {step.name} FAILED: {exc}")
                if step.abort_on_error:
                    self._log("Aborting pipeline (abort_on_error=True).")
                    break

            self.progress.emit(int((done_count + 1) / total * 100))

        self.all_done.emit()

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _log(self, msg: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_line.emit(f"[{ts}]  {msg}")
