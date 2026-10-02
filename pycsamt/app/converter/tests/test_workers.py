# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.converter.workers.ConversionWorker.

``ConversionWorker`` is a real ``QThread`` subclass; ``.run()`` is called
directly (synchronously, no real threading spun up) -- the same idiom
already used for ``_ConvertWorker`` in
``pycsamt/app/tests/test_converter_tool.py``. Because the thread object
was never ``.start()``'d, its signal connections resolve to the default
(same-thread) connection type and slots fire synchronously on ``.emit()``.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter.workers import ConversionWorker  # noqa: E402


def test_run_emits_done_with_return_value(qapp):
    worker = ConversionWorker(lambda a, b: a + b, 2, 3)
    done, errors = [], []
    worker.done.connect(done.append)
    worker.error.connect(errors.append)
    worker.run()
    assert done == [5]
    assert errors == []


def test_run_emits_error_on_exception(qapp):
    def _boom():
        raise ValueError("bad input")

    worker = ConversionWorker(_boom)
    done, errors = [], []
    worker.done.connect(done.append)
    worker.error.connect(errors.append)
    worker.run()
    assert done == []
    assert errors == ["bad input"]


def test_kwargs_are_forwarded(qapp):
    worker = ConversionWorker(lambda a, *, mult=1: a * mult, 4, mult=3)
    done = []
    worker.done.connect(done.append)
    worker.run()
    assert done == [12]


def test_supports_progress_injects_progress_callback(qapp):
    calls = []

    def _fn(progress=None):
        progress(1, 3, "a")
        progress(2, 3, "b")
        progress(3, 3, "c")
        return "ok"

    worker = ConversionWorker(_fn, supports_progress=True)
    progress_events = []
    done = []
    worker.progress.connect(lambda *a: progress_events.append(a))
    worker.done.connect(done.append)
    worker.run()
    assert progress_events == [(1, 3, "a"), (2, 3, "b"), (3, 3, "c")]
    assert done == ["ok"]


def test_progress_callback_stringifies_message(qapp):
    def _fn(progress=None):
        progress(1, 1, 42)  # non-str message

    worker = ConversionWorker(_fn, supports_progress=True)
    seen = []
    worker.progress.connect(lambda *a: seen.append(a))
    worker.run()
    assert seen == [(1, 1, "42")]


def test_without_supports_progress_fn_receives_no_progress_kwarg(qapp):
    def _fn(**kwargs):
        assert "progress" not in kwargs
        return "done"

    worker = ConversionWorker(_fn)
    done = []
    worker.done.connect(done.append)
    worker.run()
    assert done == ["done"]


def test_worker_is_a_qthread_subclass(qapp):
    from PySide6.QtCore import QThread

    worker = ConversionWorker(lambda: None)
    assert isinstance(worker, QThread)
