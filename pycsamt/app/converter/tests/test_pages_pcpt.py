# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.pcpt_page.PcptPage."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.pcpt_page import PcptPage  # noqa: E402


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
    p = PcptPage()
    yield p
    p.close()


def test_build_requires_source(page, tmp_path):
    page._output.set_path(tmp_path / "out.pcpt.json")
    page._on_build()
    assert "Pick a points CSV/XLSX" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_requires_output(page, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\n")
    page._source.set_path(src)
    page._on_build()
    assert "Pick an output" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_from_csv_runs_job(page, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\nB,4,5,6\n")
    dst = tmp_path / "targets.pcpt.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._on_build()
    assert dst.exists()
    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.build_pcpt_job
    assert "points: 2" in page.log._log.toPlainText()


def test_sheet_blank_becomes_none(page, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\n")
    dst = tmp_path / "targets.pcpt.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["sheet"] is None


def test_sheet_numeric_string_becomes_int(page, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\n")
    dst = tmp_path / "targets.pcpt.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._sheet.setText("2")
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["sheet"] == 2
    assert isinstance(worker.kwargs["sheet"], int)


def test_sheet_name_string_stays_string(page, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\n")
    dst = tmp_path / "targets.pcpt.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._sheet.setText("Sheet1")
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["sheet"] == "Sheet1"


def test_build_forwards_crs_document_id_and_overwrite(page, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\n")
    dst = tmp_path / "targets.pcpt.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._crs.setText("EPSG:4326")
    page._document_id.setText("doc-3")
    page._overwrite.setChecked(True)
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["crs"] == "EPSG:4326"
    assert worker.kwargs["document_id"] == "doc-3"
    assert worker.kwargs["overwrite"] is True
