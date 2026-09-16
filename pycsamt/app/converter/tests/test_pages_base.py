# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.converter.pages._base.

``ConverterPage.run_job`` starts a real background ``ConversionWorker``
(QThread) -- it is monkeypatched here for a lightweight fake whose
``.start()`` synchronously replays progress then emits done/error, the
same idiom used for ``_ConvertWorker`` in
``pycsamt/app/tests/test_converter_tool.py``. This exercises
``run_job``'s own orchestration (button enable/disable, log messages,
callbacks) without spinning up a real QThread or an event loop.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QLabel, QPushButton  # noqa: E402

from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages._base import (  # noqa: E402
    ConverterPage,
    mark_primary,
    page_root_layout,
    two_col_form,
)


class _FakeSignal:
    def __init__(self):
        self._fns = []

    def connect(self, fn):
        self._fns.append(fn)

    def emit(self, *a):
        for fn in list(self._fns):
            fn(*a)


def _fake_worker_cls(*, progress_events=(), done_value=None, error_msg=None):
    class _FakeWorker:
        captured = []

        def __init__(self, fn, *args, supports_progress=False, **kwargs):
            self.fn = fn
            self.args = args
            self.kwargs = kwargs
            self.supports_progress = supports_progress
            self.progress = _FakeSignal()
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            _FakeWorker.captured.append(self)

        def start(self):
            for ev in progress_events:
                self.progress.emit(*ev)
            if error_msg is not None:
                self.error.emit(error_msg)
            else:
                self.done.emit(done_value)

    return _FakeWorker


# ── page_root_layout ─────────────────────────────────────────────────────


def test_page_root_layout_wraps_content_in_scroll_area(qapp):
    from PySide6.QtWidgets import QScrollArea, QWidget

    w = QWidget()
    layout = page_root_layout(w)
    scroll = w.findChild(QScrollArea, "PageScrollArea")
    assert scroll is not None
    assert scroll.widgetResizable() is True
    label = QLabel("hi")
    layout.addWidget(label)
    assert label.parentWidget() is not None
    w.close()


# ── two_col_form ─────────────────────────────────────────────────────────


def test_two_col_form_places_fields_in_grid(qapp):
    from PySide6.QtWidgets import QLineEdit

    fields = [QLineEdit() for _ in range(3)]
    grid = two_col_form([(f"Label {i}", f) for i, f in enumerate(fields)])
    # 3 rows of (label, field) -> pairs at (0,0)/(0,1), (0,2)/(0,3), (1,0)/(1,1)
    assert grid.itemAtPosition(0, 1).widget() is fields[0]
    assert grid.itemAtPosition(0, 3).widget() is fields[1]
    assert grid.itemAtPosition(1, 1).widget() is fields[2]
    for f in fields:
        assert f.minimumWidth() >= 90


def test_two_col_form_empty_rows():
    grid = two_col_form([])
    assert grid.count() == 0


# ── mark_primary ─────────────────────────────────────────────────────────


def test_mark_primary_sets_class_property(qapp):
    btn = QPushButton("Go")
    mark_primary(btn)
    assert btn.property("class") == "primary"


# ── ConverterPage.run_job ────────────────────────────────────────────────


def test_run_job_success_path_updates_log_and_reenables_button(qapp, monkeypatch):
    monkeypatch.setattr(
        _base, "ConversionWorker", _fake_worker_cls(done_value={"file": "out.pcsf"})
    )
    page = ConverterPage()
    button = QPushButton("Run")
    button.setEnabled(True)

    page.run_job(
        lambda: None,
        run_button=button,
        success_message=lambda r: f"Wrote {r['file']}",
        busy_message="Working…",
    )

    assert button.isEnabled() is True
    text = page.log._log.toPlainText()
    assert "Working…" in text
    assert "Wrote out.pcsf" in text


def test_run_job_disables_button_while_running(qapp, monkeypatch):
    # A worker whose .start() does NOT immediately emit done/error, so we
    # can observe the button state mid-flight.
    class _PendingWorker:
        def __init__(self, *a, **k):
            self.progress = _FakeSignal()
            self.done = _FakeSignal()
            self.error = _FakeSignal()

        def start(self):
            pass

    monkeypatch.setattr(_base, "ConversionWorker", _PendingWorker)
    page = ConverterPage()
    button = QPushButton("Run")
    page.run_job(lambda: None, run_button=button)
    assert button.isEnabled() is False


def test_run_job_error_path_reports_message_and_reenables_button(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _fake_worker_cls(error_msg="boom"))
    page = ConverterPage()
    button = QPushButton("Run")
    page.run_job(lambda: None, run_button=button)
    assert button.isEnabled() is True
    assert "Error: boom" in page.log._log.toPlainText()


def test_run_job_default_success_message_is_done(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _fake_worker_cls(done_value="anything"))
    page = ConverterPage()
    page.run_job(lambda: None)
    assert "Done." in page.log._log.toPlainText()


def test_run_job_string_success_message(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _fake_worker_cls(done_value="x"))
    page = ConverterPage()
    page.run_job(lambda: None, success_message="Static message")
    assert "Static message" in page.log._log.toPlainText()


def test_run_job_calls_on_success_callback(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _fake_worker_cls(done_value=99))
    page = ConverterPage()
    seen = []
    page.run_job(lambda: None, on_success=seen.append)
    assert seen == [99]


def test_run_job_reports_progress(qapp, monkeypatch):
    monkeypatch.setattr(
        _base,
        "ConversionWorker",
        _fake_worker_cls(progress_events=[(1, 2, "a.edi"), (2, 2, "b.edi")], done_value="ok"),
    )
    page = ConverterPage()
    page.run_job(lambda: None, supports_progress=True)
    text = page.log._log.toPlainText()
    assert "[1/2] a.edi" in text
    assert "[2/2] b.edi" in text


def test_run_job_uses_default_run_button_when_none_passed(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _fake_worker_cls(done_value="ok"))
    page = ConverterPage()
    default_button = QPushButton("Default")
    page._run_button = default_button
    page.run_job(lambda: None)
    assert default_button.isEnabled() is True


def test_run_job_indeterminate_progress_when_not_supported(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _fake_worker_cls(done_value="ok"))
    page = ConverterPage()
    page.run_job(lambda: None, supports_progress=False)
    # After completion, set_indeterminate(False) resets to determinate max=1.
    assert page.log._progress.maximum() == 1
