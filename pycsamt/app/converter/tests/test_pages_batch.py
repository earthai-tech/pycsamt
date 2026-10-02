# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.batch_page.BatchPage."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import jobs  # noqa: E402
from pycsamt.app.converter.pages import _base  # noqa: E402
from pycsamt.app.converter.pages.batch_page import BatchPage  # noqa: E402


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
    p = BatchPage()
    yield p
    p.close()


def _npz(tmp_path, name="pred.npz"):
    p = tmp_path / name
    np.savez(p, resistivity=np.full((3, 3), 20.0), x=np.arange(3.0), z=np.arange(3.0))
    return p


# ── queue management ─────────────────────────────────────────────────────


def test_add_paths_populates_table_and_classifies(page, tmp_path):
    npz = _npz(tmp_path)
    page._add_paths([npz])
    assert page._table.rowCount() == 1
    assert page._table.item(0, 0).text() == str(npz)
    assert page._table.item(0, 1).text() == "inversion"
    assert page._table.item(0, 2).text() == "pending"


def test_add_paths_skips_duplicates(page, tmp_path):
    npz = _npz(tmp_path)
    page._add_paths([npz, npz])
    assert page._table.rowCount() == 1
    assert len(page._rows) == 1


def test_add_files_uses_file_dialog(page, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog

    npz = _npz(tmp_path)
    monkeypatch.setattr(
        QFileDialog, "getOpenFileNames", staticmethod(lambda *a, **k: ([str(npz)], ""))
    )
    page._add_files()
    assert page._table.rowCount() == 1


def test_add_folder_finds_edi_and_xml_candidates(page, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog

    folder = tmp_path / "batch"
    folder.mkdir()
    (folder / "a.edi").write_text("data")
    (folder / "b.xml").write_text("<x/>")
    (folder / "c.txt").write_text("ignored")

    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: str(folder))
    )
    page._add_folder()
    assert page._table.rowCount() == 2
    added = {Path(p).name for p in page._rows}
    assert added == {"a.edi", "b.xml"}


def test_add_folder_with_no_candidates_adds_the_folder_itself(page, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog

    folder = tmp_path / "empty_batch"
    folder.mkdir()
    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: str(folder))
    )
    page._add_folder()
    assert page._table.rowCount() == 1
    assert page._rows[0] == folder


def test_add_folder_cancelled_is_a_noop(page, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: "")
    )
    page._add_folder()
    assert page._table.rowCount() == 0


def test_remove_selected(page, tmp_path):
    page._add_paths([_npz(tmp_path, "a.npz"), _npz(tmp_path, "b.npz")])
    page._table.selectRow(0)
    page._remove_selected()
    assert page._table.rowCount() == 1
    assert len(page._rows) == 1


def test_clear_empties_queue(page, tmp_path):
    page._add_paths([_npz(tmp_path, "a.npz"), _npz(tmp_path, "b.npz")])
    page._clear()
    assert page._table.rowCount() == 0
    assert page._rows == []


# ── run all ──────────────────────────────────────────────────────────────


def test_run_all_requires_at_least_one_item(page, tmp_path):
    page._output_dir.set_path(tmp_path)
    page._on_run_all()
    assert "Add at least one file" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_run_all_requires_output_dir(page, tmp_path):
    page._add_paths([_npz(tmp_path)])
    page._on_run_all()
    assert "Pick an output folder" in page.log._log.toPlainText()
    assert _CapturingWorker.captured == []


def test_run_all_runs_batch_job_and_updates_row_status(page, tmp_path):
    npz = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._add_paths([npz])
    page._output_dir.set_path(out_dir)
    page._on_run_all()

    worker = _CapturingWorker.captured[0]
    assert worker.fn is jobs.run_batch_job
    assert page._table.item(0, 2).text() == "done"
    assert (out_dir / "pred.pcsf").exists()
    text = page.log._log.toPlainText()
    assert "1 done, 0 failed, 0 skipped." in text


def test_run_all_records_error_status_without_raising(page, tmp_path):
    npz = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "pred.pcsf").write_bytes(b"placeholder")  # forces overwrite conflict
    page._add_paths([npz])
    page._output_dir.set_path(out_dir)
    page._on_run_all()
    assert page._table.item(0, 2).text() == "error"


def test_run_all_forwards_to_format_and_overwrite(page, tmp_path):
    npz = _npz(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    page._add_paths([npz])
    page._output_dir.set_path(out_dir)
    page._to_format.setCurrentText("pcsm")
    page._overwrite.setChecked(True)
    page._on_run_all()
    worker = _CapturingWorker.captured[0]
    assert worker.kwargs["to_format"] == "pcsm"
    assert worker.kwargs["overwrite"] is True
