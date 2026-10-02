# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SolverBuildWorker — dependency checks and install/build runs off the GUI
thread for the Solver Builder window.

It adapts :class:`pycsamt.models.solver_build.Context` callbacks to Qt
signals, so every step (package-manager download, toolchain install,
compilation) streams its output and progress live.

Signals
-------
checks_ready(object)   Readiness of the solver (``mode="check"``)
log_line(str)          one line of tool output
stage(str)             human-readable current stage
progress(int)          0-100, or -1 while the step has no measurable progress
done(bool, str, str)   (ok, message, binary path or "")
"""

from __future__ import annotations

from PySide6.QtCore import QThread, Signal


class SolverBuildWorker(QThread):
    checks_ready = Signal(object)
    log_line = Signal(str)
    stage = Signal(str)
    progress = Signal(int)
    done = Signal(bool, str, str)

    def __init__(self, key: str, *, mode: str, source_dir=None,
                 clean: bool = False, auto_install: bool = True,
                 parent=None) -> None:
        super().__init__(parent)
        self.key = key
        self.mode = mode  # "check" | "install" | "build"
        self.source_dir = source_dir
        self.clean = clean
        self.auto_install = auto_install
        self._stop = False

    def stop(self) -> None:
        self._stop = True

    def run(self) -> None:  # noqa: D401 (Qt entry point)
        from pycsamt.models import solver_build as sb

        if self.mode == "check":
            try:
                self.checks_ready.emit(sb.check(self.key, self.source_dir))
            except Exception as exc:  # report, never crash the thread
                r = sb.Readiness([sb.Check("error", "Dependency check", False,
                                           str(exc), auto=False)])
                self.checks_ready.emit(r)
            return

        ctx = sb.Context(
            log=self.log_line.emit,
            stage=self.stage.emit,
            progress=lambda f: self.progress.emit(
                -1 if f is None else int(round(100 * f))),
            cancelled=lambda: self._stop,
        )
        result: dict = {}
        try:
            if self.mode == "install":
                steps = sb.install_plan(self.key)
            else:
                steps = sb.build_plan(self.key, self.source_dir,
                                      clean=self.clean,
                                      auto_install=self.auto_install,
                                      result=result)
            sb.run_steps(steps, ctx)
        except sb.Cancelled:
            self.done.emit(False, "Stopped by user.", "")
            return
        except Exception as exc:
            self.done.emit(False, str(exc), "")
            return
        if self.mode == "install":
            self.done.emit(True, "Dependencies installed — ready to build.",
                           "")
        else:
            self.done.emit(True, "Build finished.", result.get("binary", ""))


__all__ = ["SolverBuildWorker"]
