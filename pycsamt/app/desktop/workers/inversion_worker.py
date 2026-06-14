# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InversionWorker — QThread for running Occam2D inversions in the background.

Phase 5: Wraps ``OccamRunner.run()`` + ``InversionResult`` loading.

Two-stage operation:
  1. Build input files  (``InputBuilder``)
  2. Run Occam2D binary (``OccamRunner.run()``)
  3. Load result        (``InversionResult``)

Signals
-------
stdout_line(str)    Each line printed by the Occam2D subprocess
progress(int)       0–100 (estimated from iteration count)
finished(object)    ``InversionResult`` on completion
error(str)          Human-readable error message on failure
"""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QThread, Signal

logger = logging.getLogger(__name__)


class InversionWorker(QThread):
    """Background thread for an Occam2D inversion run."""

    stdout_line = Signal(str)
    progress    = Signal(int)
    finished    = Signal(object)   # InversionResult
    error       = Signal(str)

    def __init__(
        self,
        workdir: str,
        binary_path: Optional[str] = None,
        max_iter: int = 100,
        target_misfit: float = 1.0,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._workdir       = Path(workdir)
        self._binary_path   = binary_path or None
        self._max_iter      = max_iter
        self._target_misfit = target_misfit
        self._cancelled     = False

    # ── Cancellation ──────────────────────────────────────────────────

    def cancel(self) -> None:
        self._cancelled = True

    # ── Main thread entry ─────────────────────────────────────────────

    def run(self) -> None:
        try:
            self._run_inversion()
        except Exception as exc:
            logger.exception("InversionWorker failed")
            self.error.emit(str(exc))

    def _run_inversion(self) -> None:
        from pycsamt.models.occam2d import OccamRunner, InversionResult

        if not self._workdir.exists():
            raise FileNotFoundError(
                f"Working directory not found: {self._workdir}"
            )

        # Check binary accessibility
        binary = self._binary_path
        if binary and not Path(binary).exists():
            raise FileNotFoundError(
                f"Occam2D binary not found: {binary}\n"
                "Set the correct path in Edit → Preferences → Solvers."
            )

        self.stdout_line.emit(f"Working directory: {self._workdir}")
        self.stdout_line.emit(
            f"Binary: {binary or '(auto-detect / compile)'}"
        )
        self.stdout_line.emit(
            f"Max iterations: {self._max_iter}  "
            f"Target RMS: {self._target_misfit}"
        )
        self.progress.emit(5)

        runner = OccamRunner(
            workdir=self._workdir,
            binary_path=binary,
        )

        # Run with per-iteration stdout capture if supported
        if hasattr(runner, "iter_callback"):
            def _cb(i: int, rms: float) -> None:
                if self._cancelled:
                    raise InterruptedError("Cancelled by user.")
                pct = min(95, int(i / max(self._max_iter, 1) * 90 + 5))
                self.progress.emit(pct)
                self.stdout_line.emit(f"  iter {i:3d}  RMS = {rms:.4f}")
            runner.iter_callback = _cb
        else:
            self.stdout_line.emit("Running Occam2D…  (output after completion)")

        exit_code = runner.run(
            max_iter=self._max_iter,
            target_misfit=self._target_misfit,
        )

        if self._cancelled:
            self.stdout_line.emit("Inversion cancelled.")
            return

        self.progress.emit(95)
        self.stdout_line.emit(f"OccamRunner exit code: {exit_code}")

        result = InversionResult(workdir=str(self._workdir))
        self.stdout_line.emit(result.summary)
        self.progress.emit(100)
        self.finished.emit(result)
