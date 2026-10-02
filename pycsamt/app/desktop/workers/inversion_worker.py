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
from pathlib import Path

from PySide6.QtCore import QThread, Signal

logger = logging.getLogger(__name__)


class InversionWorker(QThread):
    """Background thread for an Occam2D inversion run."""

    stdout_line = Signal(str)
    progress = Signal(int)
    finished = Signal(object)  # InversionResult
    error = Signal(str)

    def __init__(
        self,
        workdir: str,
        binary_path: str | None = None,
        max_iter: int = 100,
        target_misfit: float = 1.0,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._workdir = Path(workdir)
        self._binary_path = binary_path or None
        self._max_iter = max_iter
        self._target_misfit = target_misfit
        self._cancelled = False

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
        from pycsamt.models.occam2d import (
            InversionResult,
            OccamRunner,
        )

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
        self.stdout_line.emit(f"Binary: {binary or '(auto-detect / compile)'}")
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
            self.stdout_line.emit(
                "Running Occam2D…  (output after completion)"
            )

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


class SolverRunWorker(QThread):
    """Run a ModEM or MARE2DEM inversion callable off the GUI thread.

    Same signals as :class:`InversionWorker`, so the Inversion window
    handles every engine identically.  (ModEM and MARE2DEM runs used to go
    through ``InversionWorker``, i.e. *Occam2D's* runner -- they never
    actually ran from the desktop.)  The runners block until the solver
    exits, so progress is reported as start / running / done.
    """

    stdout_line = Signal(str)
    progress = Signal(int)
    finished = Signal(object)
    error = Signal(str)

    def __init__(self, task, *, engine: str, binary: str | None,
                 parent=None) -> None:
        super().__init__(parent)
        self._task = task
        self._engine = engine
        self._binary = binary
        self._cancelled = False

    def cancel(self) -> None:
        # The solver subprocess cannot be interrupted mid-run by the runner
        # API; a cancelled run simply is not reported as finished.
        self._cancelled = True

    def run(self) -> None:
        try:
            self.stdout_line.emit(f"Engine: {self._engine}")
            self.stdout_line.emit(
                f"Binary: {self._binary or '(default name on PATH)'}")
            self.progress.emit(10)
            self.stdout_line.emit(
                f"Running {self._engine} — output is written to the work "
                "directory; this can take a while…")
            result = self._task()
            if self._cancelled:
                self.stdout_line.emit("Run cancelled.")
                return
            self.progress.emit(100)
            self.finished.emit(result)
        except Exception as exc:
            logger.exception("SolverRunWorker failed")
            self.error.emit(f"{type(exc).__name__}: {exc}")


class _SignalReporter:
    """Engine ``RunReporter`` that forwards to the worker's Qt signals."""

    def __init__(self, worker: EngineRunWorker) -> None:
        from pycsamt.app.desktop.controllers.inversion_engines import (
            RunReporter,
        )

        self._base = RunReporter()
        self._w = worker

    # RunReporter interface -------------------------------------------------
    @property
    def lines(self):
        return self._base.lines

    @property
    def history(self):
        return self._base.history

    def log(self, line: str) -> None:
        self._base.log(line)
        self._w.line.emit(line)

    def stage(self, text: str) -> None:
        self._base.stage(text)
        self._w.stage.emit(text)

    def iteration(self, n: int, rms: float, label: str = "") -> None:
        self._base.iteration(n, rms, label)
        self._w.iteration.emit(label, int(n), float(rms))

    def cancel(self) -> None:
        self._base.cancel()

    def cancelled(self) -> bool:
        return self._base.cancelled() or self._w.isInterruptionRequested()


class EngineRunWorker(QThread):
    """Run one Inversion Studio engine task off the GUI thread.

    ``task(reporter)`` comes from
    :meth:`pycsamt.app.desktop.controllers.inversion_engines.Engine.make_task`;
    console lines, accepted iterations and stage changes arrive live, and
    :meth:`cancel` stops the solver process tree (external engines) or the
    Occam1D loop between iterations.

    Signals
    -------
    line(str)                   a console line
    stage(str)                  a new run stage ("Station 3/10 · S07")
    iteration(str, int, float)  label (station for Occam1D), n, RMS
    finished(object)            the task's return value
    error(str)                  failure message
    cancelled()                 the run was stopped by the user
    """

    line = Signal(str)
    stage = Signal(str)
    iteration = Signal(str, int, float)
    finished = Signal(object)
    error = Signal(str)
    cancelled = Signal()

    def __init__(self, task, parent=None) -> None:
        super().__init__(parent)
        self._task = task
        self.reporter = _SignalReporter(self)

    def cancel(self) -> None:
        self.reporter.cancel()
        self.requestInterruption()

    def run(self) -> None:
        try:
            result = self._task(self.reporter)
        except InterruptedError:
            self.cancelled.emit()
            return
        except Exception as exc:
            if self.reporter.cancelled():
                self.cancelled.emit()
                return
            logger.exception("EngineRunWorker failed")
            self.error.emit(f"{type(exc).__name__}: {exc}")
            return
        if self.reporter.cancelled():
            self.cancelled.emit()
        else:
            self.finished.emit(result)
