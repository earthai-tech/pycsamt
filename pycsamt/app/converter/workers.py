# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.workers
================================

A single, generic background worker shared by every page -- mirrors
``_ConvertWorker`` in :mod:`pycsamt.app.desktop.tools.converter_tool`,
generalized to run any callable from :mod:`pycsamt.app.converter.jobs`
instead of one hard-coded conversion.
"""

from __future__ import annotations

from typing import Any, Callable

from PySide6.QtCore import QThread, Signal


class ConversionWorker(QThread):
    """Run ``fn(*args, **kwargs)`` off the UI thread.

    ``fn`` is one of the plain functions in :mod:`pycsamt.app.converter.jobs`.
    If ``fn`` accepts a ``progress`` keyword, this worker supplies a
    callback that re-emits as the Qt :attr:`progress` signal.

    Signals
    -------
    progress(int, int, str)
        current, total, message -- emitted only if *fn* reports progress.
    done(object)
        the return value of *fn*.
    error(str)
        the exception message, if *fn* raised.
    """

    progress = Signal(int, int, str)
    done = Signal(object)
    error = Signal(str)

    def __init__(
        self,
        fn: Callable[..., Any],
        *args: Any,
        supports_progress: bool = False,
        **kwargs: Any,
    ) -> None:
        super().__init__()
        self._fn = fn
        self._args = args
        self._kwargs = kwargs
        if supports_progress:
            self._kwargs["progress"] = (
                lambda cur, total, msg: self.progress.emit(cur, total, str(msg))
            )

    def run(self) -> None:
        try:
            result = self._fn(*self._args, **self._kwargs)
        except Exception as exc:  # noqa: BLE001 - surface any job failure
            self.error.emit(str(exc))
            return
        self.done.emit(result)
