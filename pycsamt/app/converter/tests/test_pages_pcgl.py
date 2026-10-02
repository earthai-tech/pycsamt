# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.pcgl_page.PcglPage."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.pcgl_page import PcglPage  # noqa: E402


class _FakeSignal:
    def __init__(self):
        self._fns = []

    def connect(self, fn):
        self._fns.append(fn)

    def emit(self, *a):
        for fn in list(self._fns):
            fn(*a)


class _CapturingWorker:
    captured = []

    def __init__(self, fn, *args, supports_progress=False, **kwargs):
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.progress = _FakeSignal()
        self.done = _FakeSignal()
        self.error = _FakeSignal()
        _CapturingWorker.captured.append(self)

    def start(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as exc:  # noqa: BLE001
            self.error.emit(str(exc))
            return
        self.done.emit(result)


@pytest.fixture(autouse=True)
def _clear_captured():
    _CapturingWorker.captured.clear()
    yield
    _CapturingWorker.captured.clear()


@pytest.fixture
def page(qapp, monkeypatch):
    monkeypatch.setattr(_base, "ConversionWorker", _CapturingWorker)
    p = PcglPage()
    yield p
    p.close()


def test_build_requires_source(page, tmp_path):
    page._output.set_path(tmp_path / "out.pcgl.json")
    page._on_build()
    assert "Pick a units CSV" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_requires_output(page, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("name,rho_min,rho_max\nA,1,10\n")
    page._source.set_path(src)
    page._on_build()
    assert "Pick an output" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_runs_job(page, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("name,rho_min,rho_max\nA,1,10\nB,10,100\n")
    dst = tmp_path / "units.pcgl.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._title.setText("Legend Title")
    page._on_build()
    assert dst.exists()
    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.build_pcgl_job
    assert worker.kwargs["title"] == "Legend Title"
    assert "entries: 2" in page.log._log.toPlainText()


def test_build_forwards_document_id_created_by_and_overwrite(page, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("name,rho_min,rho_max\nA,1,10\n")
    dst = tmp_path / "units.pcgl.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._document_id.setText("doc-9")
    page._created_by.setText("Someone")
    page._overwrite.setChecked(True)
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["document_id"] == "doc-9"
    assert worker.kwargs["created_by"] == "Someone"
    assert worker.kwargs["overwrite"] is True


def test_build_error_when_bad_csv(page, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("not,a,valid,csv\n")
    dst = tmp_path / "units.pcgl.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._on_build()
    assert "Error:" in page.log._log.toPlainText()
