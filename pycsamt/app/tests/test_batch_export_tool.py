# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for BatchExportDialog and its background ``_ExportWorker``.

``_ExportWorker`` is a real ``QThread`` whose ``run()`` performs actual
file-saving via ``Figure.savefig`` — it is exercised directly (calling
``.run()`` synchronously, not ``.start()``) so the real save logic is
covered. The dialog's own orchestration is tested separately with a
lightweight fake worker whose ``.start()`` synchronously emits through
plain-callable "signals", mirroring the idiom used in
test_main_window_actions.py / test_recompute_dlg.py / test_inversion_dlg_extra.py.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from PySide6.QtWidgets import QFileDialog

from pycsamt.app.desktop.tools.batch_export_tool import (
    _FORMATS,
    BatchExportDialog,
    _ExportWorker,
)


# ── fixtures ─────────────────────────────────────────────────────────────


@pytest.fixture
def two_figures():
    fig1, ax1 = plt.subplots()
    ax1.plot([1, 2, 3])
    fig2, ax2 = plt.subplots()
    ax2.plot([3, 2, 1])
    figs = [("Station A", fig1), ("Station B", fig2)]
    yield figs
    plt.close(fig1)
    plt.close(fig2)


@pytest.fixture
def dlg(qapp, two_figures):
    d = BatchExportDialog(two_figures, parent=None)
    yield d
    d.close()


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(store):
    """Build a fake _ExportWorker class recording construction args.

    ``store`` is a dict that will be populated with ``kwargs`` and the
    created instance, so the test can inspect what the dialog passed in.
    """

    class _FakeWorker:
        def __init__(self, figures, out_dir, fmt, dpi, transparent):
            self.figures = figures
            self.out_dir = out_dir
            self.fmt = fmt
            self.dpi = dpi
            self.transparent = transparent
            self.progress = _FakeSignal()
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            store["instance"] = self

        def start(self):
            store["started"] = True

    return _FakeWorker


# ── _ExportWorker.run() — real file-saving logic ────────────────────────


@pytest.mark.parametrize("fmt", ["png", "pdf", "svg"])
def test_worker_run_saves_files_for_each_format(two_figures, tmp_path, fmt):
    out_dir = tmp_path / "out"
    worker = _ExportWorker(two_figures, out_dir, fmt, 100, False)

    progress_calls = []
    done_calls = []
    worker.progress.connect(lambda *a: progress_calls.append(a))
    worker.done.connect(lambda *a: done_calls.append(a))

    worker.run()

    assert (out_dir / f"Station_A.{fmt}").exists()
    assert (out_dir / f"Station_B.{fmt}").exists()
    assert len(progress_calls) == 2
    assert progress_calls[0] == (1, 2, "Station A")
    assert progress_calls[1] == (2, 2, "Station B")
    assert done_calls == [(2, str(out_dir))]


def test_worker_run_creates_out_dir_if_missing(two_figures, tmp_path):
    out_dir = tmp_path / "does" / "not" / "exist"
    assert not out_dir.exists()
    worker = _ExportWorker(two_figures, out_dir, "png", 150, False)
    worker.run()
    assert out_dir.exists()
    assert (out_dir / "Station_A.png").exists()


def test_worker_run_sanitizes_label_for_filename(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([1, 2])
    figures = [("Weird/Label:*?<>|Name", fig)]
    worker = _ExportWorker(figures, tmp_path, "png", 100, False)
    worker.run()
    saved = list(tmp_path.glob("*.png"))
    assert len(saved) == 1
    # all disallowed characters replaced with "_", spaces -> "_"
    assert saved[0].name == "Weird_Label______Name.png"
    plt.close(fig)


def test_worker_run_transparent_flag_passed_through(tmp_path, monkeypatch):
    fig, ax = plt.subplots()
    ax.plot([1, 2])
    captured = {}
    orig_savefig = fig.savefig

    def spy_savefig(*a, **k):
        captured.update(k)
        return orig_savefig(*a, **k)

    monkeypatch.setattr(fig, "savefig", spy_savefig)
    worker = _ExportWorker([("Fig1", fig)], tmp_path, "png", 200, True)
    worker.run()
    assert captured["dpi"] == 200
    assert captured["transparent"] is True
    assert captured["bbox_inches"] == "tight"
    plt.close(fig)


def test_worker_run_empty_figures_emits_done_zero(tmp_path):
    worker = _ExportWorker([], tmp_path, "png", 150, False)
    done_calls = []
    worker.progress.connect(lambda *a: done_calls.append(("progress", a)))
    worker.done.connect(lambda *a: done_calls.append(("done", a)))
    worker.run()
    assert done_calls == [("done", (0, str(tmp_path)))]
    assert tmp_path.exists()


def test_worker_run_savefig_error_is_caught_and_reported_via_progress(
    tmp_path, monkeypatch
):
    fig1, ax1 = plt.subplots()
    ax1.plot([1, 2])
    fig2, ax2 = plt.subplots()
    ax2.plot([2, 1])

    def boom(*a, **k):
        raise RuntimeError("disk full")

    monkeypatch.setattr(fig1, "savefig", boom)

    worker = _ExportWorker(
        [("Bad", fig1), ("Good", fig2)], tmp_path, "png", 150, False
    )
    progress_calls = []
    done_calls = []
    worker.progress.connect(lambda *a: progress_calls.append(a))
    worker.done.connect(lambda *a: done_calls.append(a))
    worker.run()

    # worker continues past the failing figure rather than aborting
    assert not (tmp_path / "Bad.png").exists()
    assert (tmp_path / "Good.png").exists()
    assert "ERROR" in progress_calls[0][2]
    assert "disk full" in progress_calls[0][2]
    assert progress_calls[1] == (2, 2, "Good")
    # only the successfully-saved figure is counted
    assert done_calls == [(1, str(tmp_path))]

    plt.close(fig1)
    plt.close(fig2)


# ── BatchExportDialog construction ──────────────────────────────────────


def test_dialog_creates_with_figures(dlg):
    assert dlg.windowTitle() == "Batch Plot Export"
    assert dlg._src_lbl.text() == "2 canvas figures found"
    assert dlg._run_btn.isEnabled()


def test_dialog_format_combo_has_all_formats(dlg):
    items = [dlg._fmt_combo.itemText(i) for i in range(dlg._fmt_combo.count())]
    assert items == _FORMATS


def test_dialog_default_dpi_and_output_folder(dlg):
    assert dlg._dpi_spin.value() == 150
    assert dlg._dir_edit.text() == str(Path.home() / "pycsamt_figures")
    assert not dlg._transp_cb.isChecked()


def test_dialog_empty_figures_disables_run_button(qapp):
    d = BatchExportDialog([], parent=None)
    try:
        assert d._src_lbl.text() == "0 canvas figures found"
        assert not d._run_btn.isEnabled()
    finally:
        d.close()


def test_dialog_single_figure_uses_singular_label(qapp):
    fig, ax = plt.subplots()
    ax.plot([1])
    d = BatchExportDialog([("Only", fig)], parent=None)
    try:
        assert d._src_lbl.text() == "1 canvas figure found"
    finally:
        d.close()
        plt.close(fig)


def test_dialog_more_than_five_figures_shows_overflow_hint(qapp):
    figs = []
    for i in range(7):
        fig, ax = plt.subplots()
        ax.plot([i])
        figs.append((f"Fig{i}", fig))
    d = BatchExportDialog(figs, parent=None)
    try:
        assert d._src_lbl.text() == "7 canvas figures found"
    finally:
        d.close()
        for _, f in figs:
            plt.close(f)


# ── _pick_dir ─────────────────────────────────────────────────────────────


def test_pick_dir_updates_text_when_chosen(dlg, tmp_path, monkeypatch):
    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    dlg._pick_dir()
    assert dlg._dir_edit.text() == str(tmp_path)


def test_pick_dir_cancelled_leaves_text_unchanged(dlg, monkeypatch):
    original = dlg._dir_edit.text()
    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: ""),
    )
    dlg._pick_dir()
    assert dlg._dir_edit.text() == original


# ── _on_run — dialog orchestration with a fake worker ───────────────────


def test_on_run_constructs_worker_with_current_state(dlg, tmp_path, monkeypatch):
    store = {}
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.batch_export_tool._ExportWorker",
        _fake_worker_cls(store),
    )
    dlg._dir_edit.setText(str(tmp_path))
    dlg._fmt_combo.setCurrentText("SVG")
    dlg._dpi_spin.setValue(200)
    dlg._transp_cb.setChecked(True)

    dlg._on_run()

    inst = store["instance"]
    assert inst.figures == dlg._figures
    assert inst.out_dir == Path(tmp_path)
    assert inst.fmt == "SVG"
    assert inst.dpi == 200
    assert inst.transparent is True
    assert store["started"] is True
    assert not dlg._run_btn.isEnabled()
    assert dlg._progress.value() == 0
    assert f"Exporting 2 figures" in dlg._log.toPlainText()


def test_on_run_connects_worker_signals_to_dialog_slots(
    dlg, tmp_path, monkeypatch
):
    store = {}
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.batch_export_tool._ExportWorker",
        _fake_worker_cls(store),
    )
    dlg._dir_edit.setText(str(tmp_path))
    dlg._on_run()

    inst = store["instance"]
    inst.progress.emit(1, 2, "Station A")
    assert "[1/2]" in dlg._log.toPlainText()
    assert dlg._progress.value() == 1

    inst.done.emit(2, str(tmp_path))
    assert "Done" in dlg._log.toPlainText()
    assert dlg._run_btn.isEnabled()


# ── _on_progress ──────────────────────────────────────────────────────────


def test_on_progress_updates_bar_and_log(dlg):
    dlg._on_progress(3, 5, "Station C")
    assert dlg._progress.maximum() == 5
    assert dlg._progress.value() == 3
    assert "[3/5]  Station C" in dlg._log.toPlainText()


# ── _on_done ──────────────────────────────────────────────────────────────


def test_on_done_reenables_button_and_logs_summary(dlg, tmp_path):
    dlg._run_btn.setEnabled(False)
    dlg._on_done(2, str(tmp_path))
    assert dlg._run_btn.isEnabled()
    assert f"Done — 2 figure(s) saved to {tmp_path}" in dlg._log.toPlainText()


# ── _on_error ─────────────────────────────────────────────────────────────


def test_on_error_reenables_button_and_logs_message(dlg):
    dlg._run_btn.setEnabled(False)
    dlg._on_error("boom")
    assert dlg._run_btn.isEnabled()
    assert "Error: boom" in dlg._log.toPlainText()


# ── button wiring ─────────────────────────────────────────────────────────


def test_run_button_triggers_on_run(dlg):
    with mock.patch.object(dlg, "_on_run") as m:
        dlg._run_btn.click()
        m.assert_called_once()


def test_browse_button_triggers_pick_dir(dlg):
    with mock.patch.object(dlg, "_pick_dir") as m:
        btn = dlg.findChild(type(dlg._run_btn), None)
        # locate the "Browse…" button explicitly among children
        from PySide6.QtWidgets import QPushButton

        browse_btns = [
            b
            for b in dlg.findChildren(QPushButton)
            if b.text() == "Browse…"
        ]
        assert len(browse_btns) == 1
        browse_btns[0].click()
        m.assert_called_once()


def test_close_button_rejects_dialog(dlg):
    from PySide6.QtWidgets import QDialogButtonBox, QDialog

    box = dlg.findChild(QDialogButtonBox)
    box.rejected.emit()
    assert dlg.result() == QDialog.DialogCode.Rejected
