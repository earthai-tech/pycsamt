# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.pcsf_pcsm_page.PcsfPcsmPage."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.pcsf_pcsm_page import PcsfPcsmPage  # noqa: E402


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
    p = PcsfPcsmPage()
    yield p
    p.close()


@pytest.fixture
def pcsf_file(tmp_path):
    npz = tmp_path / "pred.npz"
    np.savez(npz, resistivity=np.full((3, 4), 50.0), x=np.arange(4.0), z=np.arange(3.0))
    dst = tmp_path / "m.pcsf"
    jobs.convert_to_pcsf_job(npz, dst, "pcsf", tmp_path)
    return dst


# ── Transcode ────────────────────────────────────────────────────────────


def test_transcode_requires_source(page, tmp_path):
    page._output.set_path(tmp_path / "out.pcsm")
    page._on_transcode()
    assert "Pick a .pcsf/.pcsm source" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_transcode_requires_output(page, pcsf_file):
    page._source.set_path(pcsf_file)
    page._on_transcode()
    assert "Pick an output file" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_transcode_runs_job(page, pcsf_file, tmp_path):
    dst = tmp_path / "m.pcsm"
    page._source.set_path(pcsf_file)
    page._output.set_path(dst)
    page._on_transcode()
    assert dst.exists()
    assert _CapturingWorker.captured[0].fn is jobs.transcode_pcsf_pcsm_job
    assert f"Wrote {dst}" in page.log._log.toPlainText()


def test_transcode_forwards_log10_view_and_overwrite(page, pcsf_file, tmp_path):
    dst = tmp_path / "m.pcsm"
    page._source.set_path(pcsf_file)
    page._output.set_path(dst)
    page._log10_view.setChecked(True)
    page._overwrite.setChecked(True)
    page._on_transcode()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["log10_view"] is True
    assert worker.kwargs["overwrite"] is True


# ── Validate ─────────────────────────────────────────────────────────────


def test_validate_requires_source(page):
    page._on_validate()
    assert "Pick a .pcsf/.pcsm file" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_validate_runs_job_and_formats_result(page, pcsf_file):
    page._source.set_path(pcsf_file)
    page._on_validate()
    text = page.log._log.toPlainText()
    assert "VALID" in text
    assert "header" in text


def test_format_validation_marks_failures():
    result = {
        "file": "x.pcsf",
        "valid": False,
        "checks": [
            {"step": "header", "ok": True, "detail": "geometry=grid2d"},
            {"step": "load", "ok": False, "detail": "boom"},
        ],
    }
    text = PcsfPcsmPage._format_validation(result)
    assert "INVALID" in text
    assert "✓" in text  # check mark for header
    assert "✗" in text  # cross mark for load
    assert "boom" in text


# ── Info ─────────────────────────────────────────────────────────────────


def test_info_requires_source(page):
    page._on_info()
    assert "Pick a .pcsf/.pcsm file" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_info_runs_job_and_formats_result(page, pcsf_file):
    page._source.set_path(pcsf_file)
    page._on_info()
    text = page.log._log.toPlainText()
    assert str(pcsf_file) in text
    assert "resistivity" in text


def test_format_info_without_resistivity_or_stations():
    result = {
        "file": "x.pcsf",
        "size_bytes": 128,
        "container": "pcsf",
        "geometry_kind": "grid2d",
        "source_backend": "ai_arrays",
        "crs": None,
    }
    text = PcsfPcsmPage._format_info(result)
    assert "x.pcsf" in text
    assert "crs: -" in text
    assert "stations" not in text


def test_format_info_with_all_nan_resistivity():
    result = {
        "file": "x.pcsf",
        "size_bytes": 128,
        "container": "pcsf",
        "geometry_kind": "grid2d",
        "source_backend": "ai_arrays",
        "crs": "EPSG:4326",
        "resistivity": {"shape": [3, 3], "n_nan": 9},
        "stations": {"n": 5, "has_lonlat": True},
    }
    text = PcsfPcsmPage._format_info(result)
    assert "all NaN" in text
    assert "stations: 5" in text
