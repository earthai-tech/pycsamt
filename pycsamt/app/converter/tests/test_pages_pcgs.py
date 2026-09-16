# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.pcgs_page.PcgsPage."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.pcgs_page import PcgsPage  # noqa: E402


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
    p = PcgsPage()
    yield p
    p.close()


def test_build_requires_at_least_one_csv(page, tmp_path):
    page._output.set_path(tmp_path / "s.pcgs.json")
    page._on_build()
    assert "Pick at least one" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_requires_output(page, tmp_path):
    planar = tmp_path / "planar.csv"
    planar.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg\n100,bedding,45,60,135\n"
    )
    page._planar.set_path(planar)
    page._on_build()
    assert "Pick an output" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_build_from_planar_only(page, tmp_path):
    planar = tmp_path / "planar.csv"
    planar.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg\n100,bedding,45,60,135\n"
    )
    dst = tmp_path / "s.pcgs.json"
    page._planar.set_path(planar)
    page._output.set_path(dst)
    page._on_build()
    assert dst.exists()
    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.build_pcgs_job
    assert worker.kwargs["planar_path"] == planar
    assert worker.kwargs["linear_path"] is None
    assert worker.kwargs["faults_path"] is None
    assert "planar: 1" in page.log._log.toPlainText()


def test_build_with_nonexistent_optional_paths_are_dropped(page, tmp_path):
    planar = tmp_path / "planar.csv"
    planar.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg\n100,bedding,45,60,135\n"
    )
    dst = tmp_path / "s.pcgs.json"
    page._planar.set_path(planar)
    # linear/faults fields are typed but point at files that don't exist.
    page._linear.set_path(tmp_path / "does_not_exist.csv")
    page._output.set_path(dst)
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["linear_path"] is None


def test_build_forwards_title_document_id_and_overwrite(page, tmp_path):
    planar = tmp_path / "planar.csv"
    planar.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg\n100,bedding,45,60,135\n"
    )
    dst = tmp_path / "s.pcgs.json"
    page._planar.set_path(planar)
    page._output.set_path(dst)
    page._title.setText("Structure Title")
    page._document_id.setText("doc-7")
    page._overwrite.setChecked(True)
    page._on_build()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["title"] == "Structure Title"
    assert worker.kwargs["document_id"] == "doc-7"
    assert worker.kwargs["overwrite"] is True
