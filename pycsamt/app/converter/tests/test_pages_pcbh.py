# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.pcbh_page.PcbhPage."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.pcbh_page import PcbhPage  # noqa: E402


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
    p = PcbhPage()
    yield p
    p.close()


_CSV_TEXT = (
    "borehole_id,x,y,z,crs,kind,status,total_depth_md,from_md,to_md,lithology\n"
    "BH-1,100,200,10,EPSG:4326,water,completed,20,0,10,Sand\n"
)


def test_kind_toggle_shows_las_box_only_for_las(page):
    assert not page._las_box.isVisibleTo(page)
    page._kind.setCurrentText("las")
    assert page._las_box.isVisibleTo(page)
    page._kind.setCurrentText("csv")
    assert not page._las_box.isVisibleTo(page)


def test_build_requires_source(page, tmp_path):
    page._output.set_path(tmp_path / "out.pcbh.json")
    page._on_build()
    assert "Pick a source" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_requires_output(page, tmp_path):
    src = tmp_path / "boreholes.csv"
    src.write_text(_CSV_TEXT)
    page._source.set_path(src)
    page._on_build()
    assert "Pick an output" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_from_csv_runs_job(page, tmp_path):
    src = tmp_path / "boreholes.csv"
    src.write_text(_CSV_TEXT)
    dst = tmp_path / "boreholes.pcbh.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._on_build()
    assert dst.exists()
    assert _CapturingWorker.captured[0].fn is jobs.build_pcbh_job
    assert "boreholes: 1" in page.log._log.toPlainText()


def test_las_kind_without_collar_id_blocks_before_job(page, tmp_path):
    fake_las = tmp_path / "hole.las"
    fake_las.write_text("~V\nVERS. 2.0 :\n")
    dst = tmp_path / "hole.pcbh.json"
    page._source.set_path(fake_las)
    page._output.set_path(dst)
    page._kind.setCurrentText("las")
    page._on_build()
    assert "LAS import needs a collar id and CRS" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_las_kind_auto_detected_from_suffix_also_blocks(page, tmp_path):
    fake_las = tmp_path / "hole.las"
    fake_las.write_text("~V\nVERS. 2.0 :\n")
    dst = tmp_path / "hole.pcbh.json"
    page._source.set_path(fake_las)
    page._output.set_path(dst)
    assert page._kind.currentText() == "auto"
    page._on_build()
    assert "LAS import needs a collar id and CRS" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_las_kind_with_collar_and_crs_reaches_job(page, tmp_path):
    fake_las = tmp_path / "hole.las"
    fake_las.write_text("~V\nVERS. 2.0 :\n")
    dst = tmp_path / "hole.pcbh.json"
    page._source.set_path(fake_las)
    page._output.set_path(dst)
    page._kind.setCurrentText("las")
    page._collar_id.setText("BH-1")
    page._crs.setText("EPSG:32650")
    page._on_build()
    # Validation passes and the job is dispatched (it may still fail deep in
    # the LAS parser on this placeholder file -- either way, the guard must
    # not have blocked it).
    assert len(_CapturingWorker.captured) == 1
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["source_kind"] == "las"
    assert worker.kwargs["collar_id"] == "BH-1"
    assert worker.kwargs["crs_horizontal"] == "EPSG:32650"


def test_build_forwards_document_id_and_overwrite(page, tmp_path):
    src = tmp_path / "boreholes.csv"
    src.write_text(_CSV_TEXT)
    dst = tmp_path / "boreholes.pcbh.json"
    page._source.set_path(src)
    page._output.set_path(dst)
    page._document_id.setText("doc-1")
    page._overwrite.setChecked(True)
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["document_id"] == "doc-1"
    assert worker.kwargs["overwrite"] is True
