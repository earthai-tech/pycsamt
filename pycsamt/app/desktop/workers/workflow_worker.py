# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
WorkflowWorker — runs planned workflow steps off the GUI thread.

Each step is executed through :meth:`WorkflowController.execute_step` (the
library's own ``Step.transform``). The worker announces every step before
and after it runs — which ``Pipeline.run``'s ``on_step`` hook cannot — and
checks for a stop request between steps.

Signals
-------
step_started(int)   workflow index about to run
step_finished(int)  workflow index finished (DONE or ERROR — see the step)
log_line(str)       human-readable progress line
run_finished(str)   "completed" | "stopped" | "aborted" (error policy) |
                    "blocked" (a step's input was not available)
"""

from __future__ import annotations

import time

from PySide6.QtCore import QThread, Signal


class WorkflowWorker(QThread):
    step_started = Signal(int)
    step_finished = Signal(int)
    log_line = Signal(str)
    run_finished = Signal(str)

    def __init__(self, ctrl, indices: list[int], parent=None) -> None:
        super().__init__(parent)
        self._ctrl = ctrl
        self._indices = list(indices)

    def run(self) -> None:  # noqa: D401 (Qt entry point)
        ctrl = self._ctrl
        t0 = time.perf_counter()
        outcome = "completed"
        for k, idx in enumerate(self._indices, start=1):
            if self.isInterruptionRequested():
                outcome = "stopped"
                break
            sites_in = ctrl.input_for(idx)
            ws = ctrl.steps[idx]
            if sites_in is None:
                self.log_line.emit(
                    f"✗ {ws.label}: no input — run the earlier steps first."
                )
                outcome = "blocked"
                break
            self.step_started.emit(idx)
            self.log_line.emit(f"▶ [{ws.code}] {ws.label} …")
            _, ok, stop = ctrl.execute_step(idx, sites_in, step_idx=k)
            if ok:
                self.log_line.emit(
                    f"✓ [{ws.code}] {ws.label}  {ws.n_in}→{ws.n_out} "
                    f"stations  {ws.elapsed:.2f}s"
                )
            else:
                self.log_line.emit(f"✗ [{ws.code}] {ws.label}  {ws.error}")
            self.step_finished.emit(idx)
            if stop:
                outcome = "aborted"
                break
        ctrl.run_elapsed = time.perf_counter() - t0
        self.run_finished.emit(outcome)


__all__ = ["WorkflowWorker"]
