# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.edi_xml_page.EdiXmlPage."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.edi_xml_page import EdiXmlPage  # noqa: E402

_ROOT = Path(__file__).resolve().parents[4]
_EDI_DIR = _ROOT / "data" / "MT" / "broken-hill" / "edis"


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
    p = EdiXmlPage()
    yield p
    p.close()


def _edi_dir():
    if not _EDI_DIR.exists() or not any(_EDI_DIR.glob("*.edi")):
        pytest.skip(f"No Broken Hill EDI sample data found at {_EDI_DIR}")
    return _EDI_DIR


# ── direction toggle ─────────────────────────────────────────────────────


def test_direction_defaults_to_edi_to_xml(page):
    # Widgets never call .show() in these tests (no display), so isVisible()
    # always reads False regardless of setVisible() -- isVisibleTo(ancestor)
    # reports the "would be visible if shown" flag we actually want here.
    assert page._direction.currentIndex() == 0
    assert page._prefer_spectra.isVisibleTo(page)
    assert not page._on_loss.isVisibleTo(page)
    assert not page._strict.isVisibleTo(page)


def test_direction_toggle_swaps_visible_options(page):
    page._direction.setCurrentIndex(1)
    assert not page._prefer_spectra.isVisibleTo(page)
    assert page._on_loss.isVisibleTo(page)
    assert page._strict.isVisibleTo(page)


# ── validation ───────────────────────────────────────────────────────────


def test_convert_requires_source(page, tmp_path):
    page._output_dir.set_path(tmp_path)
    page._on_convert()
    assert "Pick a source" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_convert_requires_output_dir(page):
    page._source.set_path(_edi_dir())
    page._on_convert()
    assert "Pick an output folder" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


# ── EDI -> XML ───────────────────────────────────────────────────────────


def test_edi_to_xml_runs_job(page, tmp_path):
    edi_dir = _edi_dir()
    out_dir = tmp_path / "xml"
    out_dir.mkdir()
    page._source.set_path(edi_dir)
    page._output_dir.set_path(out_dir)
    page._on_convert()

    assert len(_CapturingWorker.captured) == 1
    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.edi_to_xml_job
    n_edi = len(list(edi_dir.glob("*.edi")))
    assert len(list(out_dir.glob("*.xml"))) == n_edi
    assert f"Wrote {n_edi} EMTF-XML file(s)" in page.log._log.toPlainText()


def test_edi_to_xml_forwards_prefer_spectra_and_overwrite(page, tmp_path):
    edi_dir = _edi_dir()
    out_dir = tmp_path / "xml"
    out_dir.mkdir()
    page._source.set_path(edi_dir)
    page._output_dir.set_path(out_dir)
    page._prefer_spectra.setChecked(False)
    page._overwrite.setChecked(True)
    page._on_convert()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["prefer_spectra"] is False
    assert worker.kwargs["overwrite"] is True


# ── XML -> EDI ───────────────────────────────────────────────────────────


def test_xml_to_edi_runs_job(page, tmp_path):
    edi_dir = _edi_dir()
    xml_dir = tmp_path / "xml"
    jobs.edi_to_xml_job(edi_dir, xml_dir)

    back_dir = tmp_path / "back"
    back_dir.mkdir()
    page._direction.setCurrentIndex(1)
    page._source.set_path(xml_dir)
    page._output_dir.set_path(back_dir)
    page._on_convert()

    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.xml_to_edi_job
    n_xml = len(list(xml_dir.glob("*.xml")))
    assert len(list(back_dir.glob("*.edi"))) == n_xml
    assert f"Wrote {n_xml} EDI file(s)" in page.log._log.toPlainText()


def test_xml_to_edi_forwards_on_loss_and_strict(page, tmp_path):
    edi_dir = _edi_dir()
    xml_dir = tmp_path / "xml"
    jobs.edi_to_xml_job(edi_dir, xml_dir)

    back_dir = tmp_path / "back"
    back_dir.mkdir()
    page._direction.setCurrentIndex(1)
    page._source.set_path(xml_dir)
    page._output_dir.set_path(back_dir)
    page._on_loss.setCurrentText("raise")
    page._strict.setChecked(False)
    page._on_convert()

    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["on_loss"] == "raise"
    assert worker.kwargs["strict"] is False


def test_convert_error_path(page, tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._source.set_path(empty_dir)
    page._output_dir.set_path(out_dir)
    page._on_convert()
    assert "Error:" in page.log._log.toPlainText()
