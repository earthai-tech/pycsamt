# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.inversion_page.InversionPage."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.inversion_page import InversionPage  # noqa: E402


class _FakeSignal:
    def __init__(self):
        self._fns = []

    def connect(self, fn):
        self._fns.append(fn)

    def emit(self, *a):
        for fn in list(self._fns):
            fn(*a)


class _CapturingWorker:
    """Captures the call args run_job would forward to a real worker."""

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
    p = InversionPage()
    yield p
    p.close()


def _npz(tmp_path, name="pred.npz", shape=(4, 5), value=100.0):
    p = tmp_path / name
    np.savez(p, resistivity=np.full(shape, value), x=np.arange(shape[1], dtype=float),
              z=np.arange(shape[0], dtype=float))
    return p


def test_source_detection_updates_hint_label(page, tmp_path):
    p = _npz(tmp_path)
    page._source.set_path(p)
    assert "ai_arrays" in page._detected.text()


def test_source_detection_shows_placeholder_when_no_path(page):
    page._source.set_path("")
    assert "auto-detect" in page._detected.text()


def test_source_detection_handles_bad_path_gracefully(page, tmp_path):
    junk = tmp_path / "notes.txt"
    junk.write_text("not a source")
    page._source.set_path(junk)
    # detect() on an unsupported file should not raise out of the slot;
    # the label reports whatever the classifier returns.
    assert page._detected.text() != ""


def test_convert_requires_source(page, tmp_path):
    page._output_dir.set_path(tmp_path)
    page._on_convert()
    assert "Pick a source" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_convert_requires_output_dir(page, tmp_path):
    page._source.set_path(_npz(tmp_path))
    page._on_convert()
    assert "Pick an output folder" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_convert_runs_job_and_writes_output(page, tmp_path):
    src = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._source.set_path(src)
    page._output_dir.set_path(out_dir)
    page._on_convert()
    assert len(_CapturingWorker.captured) == 1
    assert (out_dir / "pred.pcsf").exists()
    assert "Wrote" in page.log._log.toPlainText()


def test_convert_forwards_advanced_options(page, tmp_path):
    src = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._source.set_path(src)
    page._output_dir.set_path(out_dir)
    page._epsg.setValue(32650)
    page._utm_zone.setText("48N")
    page._created_by.setText("Tester")
    page._description.setText("A description")
    page._overwrite.setChecked(True)
    page._on_convert()

    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.convert_to_pcsf_job
    assert worker.kwargs["epsg"] == 32650
    assert worker.kwargs["utm_zone"] == "48N"
    assert worker.kwargs["created_by"] == "Tester"
    assert worker.kwargs["description"] == "A description"
    assert worker.kwargs["overwrite"] is True


def test_convert_solver_auto_becomes_none(page, tmp_path):
    src = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._source.set_path(src)
    page._output_dir.set_path(out_dir)
    assert page._solver.currentText() == "auto"
    page._on_convert()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["solver"] is None


def test_convert_iteration_zero_becomes_none(page, tmp_path):
    src = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._source.set_path(src)
    page._output_dir.set_path(out_dir)
    page._iteration.setValue(0)
    page._on_convert()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["iteration"] is None


def test_convert_iteration_nonzero_is_forwarded(page, tmp_path):
    src = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._source.set_path(src)
    page._output_dir.set_path(out_dir)
    page._iteration.setValue(7)
    page._on_convert()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["iteration"] == 7


def test_convert_error_reports_message(page, tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    junk = tmp_path / "notes.txt"
    junk.write_text("not a source")
    page._source.set_path(junk)
    page._output_dir.set_path(out_dir)
    page._on_convert()
    assert "Error:" in page.log._log.toPlainText()
