# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for RecomputeDialog.

RecomputeWorker is a real QThread that shells out to EDIRecomputer, so it
is replaced with a lightweight stand-in (mirroring the LoaderWorker fake
used in test_main_window_actions.py) whose ``.start()`` synchronously
emits through plain-callable "signals" instead of spinning a real thread.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QDialog

from pycsamt.app.desktop.dialogs.recompute_dlg import (
    RecomputeDialog,
    _icon,
)


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(*, result=None, error=None, is_running=False):
    class _FakeWorker:
        captured_args = []
        captured_kwargs = []

        def __init__(self, source, parent=None, **kw):
            self.source = source
            self.parent = parent
            self.kw = kw
            _FakeWorker.captured_args.append(source)
            _FakeWorker.captured_kwargs.append(kw)
            self.station_done = _FakeSignal()
            self.finished = _FakeSignal()
            self.error = _FakeSignal()
            self._is_running = is_running
            self.interrupted = False
            self.waited = False

        def isRunning(self):
            return self._is_running

        def requestInterruption(self):
            self.interrupted = True

        def wait(self, ms):
            self.waited = True

        def start(self):
            if error is not None:
                self.error.emit(error)
            else:
                self.finished.emit(result)

    return _FakeWorker


@pytest.fixture
def dlg(qapp):
    d = RecomputeDialog(sites=None)
    yield d
    d.close()


@pytest.fixture
def dlg_with_sites(qapp):
    d = RecomputeDialog(sites=["S1", "S2"])
    yield d
    d.close()


def _fake_result(*, records=None, failed=None, output_root="out"):
    return SimpleNamespace(
        records=records or [],
        failed=failed or [],
        output_root=output_root,
    )


# ── construction ──────────────────────────────────────────────────────────


def test_creates(dlg):
    assert dlg.windowTitle() == "Recompute EDI Files"


def test_no_sites_defaults_to_folder_mode(dlg):
    assert dlg._rb_folder.isChecked()
    assert not dlg._rb_loaded.isEnabled()


def test_with_sites_defaults_to_loaded_mode(dlg_with_sites):
    assert dlg_with_sites._rb_loaded.isChecked()
    assert dlg_with_sites._rb_loaded.isEnabled()


def test_folder_edit_enabled_state_follows_radio(dlg_with_sites):
    assert not dlg_with_sites._folder_edit.isEnabled()
    dlg_with_sites._rb_folder.setChecked(True)
    assert dlg_with_sites._folder_edit.isEnabled()


def test_icon_helper_missing_file_returns_empty_icon(qapp):
    icon = _icon("this_icon_definitely_does_not_exist_xyz")
    assert icon.isNull()


# ── browse slots ──────────────────────────────────────────────────────────


def test_browse_source_sets_folder_edit(dlg, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    dlg._browse_source()
    assert dlg._folder_edit.text() == str(tmp_path)


def test_browse_source_cancelled_noop(dlg, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: ""),
    )
    dlg._browse_source()
    assert dlg._folder_edit.text() == ""


def test_browse_output_sets_out_edit(dlg, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    dlg._browse_output()
    assert dlg._out_edit.text() == str(tmp_path)


# ── _resolve_source ───────────────────────────────────────────────────────


def test_resolve_source_uses_loaded_sites(dlg_with_sites):
    assert dlg_with_sites._resolve_source() == ["S1", "S2"]


def test_resolve_source_empty_folder_warns_and_returns_none(dlg):
    dlg._folder_edit.setText("")
    assert dlg._resolve_source() is None


def test_resolve_source_missing_path_warns_and_returns_none(dlg, tmp_path):
    dlg._folder_edit.setText(str(tmp_path / "does_not_exist"))
    assert dlg._resolve_source() is None


def test_resolve_source_valid_folder_returns_path(dlg, tmp_path):
    from pathlib import Path

    dlg._folder_edit.setText(str(tmp_path))
    assert dlg._resolve_source() == Path(tmp_path)


# ── _collect_kwargs ───────────────────────────────────────────────────────


def test_collect_kwargs_defaults(dlg):
    kw = dlg._collect_kwargs()
    assert kw["rotate_angle"] is None  # spin at minimum ("none")
    assert kw["rotate_components"] == ("Z", "Tip")
    assert kw["fmin"] is None
    assert kw["fmax"] is None
    assert kw["fill_missing_values"] is None
    assert kw["output_root"] is None
    assert kw["template"] == "{source_stem}.edi"
    assert kw["preserve_line_dirs"] is True
    assert kw["overwrite"] is False
    assert kw["write"] is True
    assert kw["manifest_csv"] is True
    assert kw["recompute_resphase"] is True


def test_collect_kwargs_rotation_set(dlg):
    dlg._rot_spin.setValue(45.0)
    kw = dlg._collect_kwargs()
    assert kw["rotate_angle"] == pytest.approx(45.0)


@pytest.mark.parametrize(
    "idx,expected",
    [(0, ("Z", "Tip")), (1, ("Z",)), (2, ("Tip",))],
)
def test_collect_kwargs_components(dlg, idx, expected):
    dlg._comp_combo.setCurrentIndex(idx)
    assert dlg._collect_kwargs()["rotate_components"] == expected


def test_collect_kwargs_frequency_range(dlg):
    dlg._fmin_edit.setText("0.1")
    dlg._fmax_edit.setText("100")
    kw = dlg._collect_kwargs()
    assert kw["fmin"] == pytest.approx(0.1)
    assert kw["fmax"] == pytest.approx(100.0)


def test_collect_kwargs_invalid_frequency_is_none(dlg):
    dlg._fmin_edit.setText("not-a-number")
    assert dlg._collect_kwargs()["fmin"] is None


def test_collect_kwargs_fill_missing(dlg):
    dlg._fill_combo.setCurrentText("nan")
    assert dlg._collect_kwargs()["fill_missing_values"] == "nan"


def test_collect_kwargs_output_root_set(dlg):
    dlg._out_edit.setText("my_out")
    assert dlg._collect_kwargs()["output_root"] == "my_out"


def test_collect_kwargs_empty_template_falls_back_to_default(dlg):
    dlg._tpl_edit.setText("   ")
    assert dlg._collect_kwargs()["template"] == "{source_stem}.edi"


# ── _parse_float ──────────────────────────────────────────────────────────


def test_parse_float_empty():
    assert RecomputeDialog._parse_float("") is None
    assert RecomputeDialog._parse_float("   ") is None


def test_parse_float_valid():
    assert RecomputeDialog._parse_float("3.5") == pytest.approx(3.5)


def test_parse_float_invalid():
    assert RecomputeDialog._parse_float("abc") is None


# ── _on_run_clicked ───────────────────────────────────────────────────────


def test_on_run_clicked_no_source_does_not_start_worker(dlg, monkeypatch):
    fake_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.recompute_worker.RecomputeWorker",
        fake_cls,
    )
    dlg._folder_edit.setText("")
    dlg._on_run_clicked()
    assert dlg._worker is None


def test_on_run_clicked_starts_worker_with_loaded_sites(
    dlg_with_sites, monkeypatch
):
    fake_cls = _fake_worker_cls(result=_fake_result())
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.recompute_worker.RecomputeWorker",
        fake_cls,
    )
    dlg_with_sites._on_run_clicked()
    assert dlg_with_sites._worker is not None
    assert dlg_with_sites._worker.source == ["S1", "S2"]


def test_on_run_clicked_cancels_running_worker(dlg):
    fake = _fake_worker_cls(is_running=True)("src")
    dlg._worker = fake
    dlg._on_run_clicked()
    assert fake.interrupted is True
    assert fake.waited is True
    assert dlg._cur_lbl.text() == "Cancelled."
    assert dlg._btn_run.text() == "Recompute"


# ── _on_close_clicked ─────────────────────────────────────────────────────


def test_on_close_clicked_no_worker_rejects(dlg):
    dlg._on_close_clicked()
    assert dlg.result() == QDialog.DialogCode.Rejected


def test_on_close_clicked_running_worker_interrupted(dlg):
    fake = _fake_worker_cls(is_running=True)("src")
    dlg._worker = fake
    dlg._on_close_clicked()
    assert fake.interrupted is True
    assert dlg.result() == QDialog.DialogCode.Rejected


# ── _on_station_done ──────────────────────────────────────────────────────


def test_on_station_done_ok_status(dlg):
    dlg._on_station_done(1, 4, "S001", "ok")
    assert dlg._ok_lbl.text() == "✔  1 ok"
    assert dlg._pbar.value() == 25
    assert dlg._log_list.count() == 1


def test_on_station_done_failed_status_with_message(dlg):
    dlg._on_station_done(2, 4, "S002", "failed", "boom")
    assert dlg._failed_lbl.text() == "✕  1 failed"
    # main row + error-detail row
    assert dlg._log_list.count() == 2


def test_on_station_done_failed_long_message_truncated(dlg):
    long_msg = "x" * 200
    dlg._on_station_done(1, 1, "S003", "failed", long_msg)
    err_text = dlg._log_list.item(1).text()
    assert "…" in err_text
    assert len(err_text) < len(long_msg)


def test_on_station_done_other_status(dlg):
    dlg._on_station_done(1, 2, "S004", "skipped")
    assert dlg._cur_lbl.text() == "Processing:  S004"
    assert dlg._log_list.count() == 1


def test_on_station_done_zero_total_no_division_error(dlg):
    dlg._on_station_done(0, 0, "S005", "ok")
    assert dlg._pbar.value() == 0


# ── _on_finished ──────────────────────────────────────────────────────────


def test_on_finished_all_ok_shows_load_button(dlg):
    result = _fake_result(
        records=[SimpleNamespace(status="ok"), SimpleNamespace(status="ok")]
    )
    dlg._on_finished(result)
    assert not dlg._btn_load.isHidden()
    assert dlg._result is result


def test_on_finished_all_failed_file_exists_hint(dlg):
    result = _fake_result(
        records=[],
        failed=[SimpleNamespace(message="FileExistsError: already there")],
    )
    dlg._on_finished(result)
    items = [dlg._log_list.item(i).text() for i in range(dlg._log_list.count())]
    assert any("already exist" in t for t in items)
    assert dlg._btn_load.isHidden()


def test_on_finished_all_failed_other_reason_hint(dlg):
    result = _fake_result(
        records=[],
        failed=[SimpleNamespace(message="some other parsing error")],
    )
    dlg._on_finished(result)
    items = [dlg._log_list.item(i).text() for i in range(dlg._log_list.count())]
    assert any("First failure" in t for t in items)


def test_on_finished_mixed_ok_and_failed_no_hint(dlg):
    result = _fake_result(
        records=[SimpleNamespace(status="ok")],
        failed=[SimpleNamespace(message="oops")],
    )
    dlg._on_finished(result)
    items = [dlg._log_list.item(i).text() for i in range(dlg._log_list.count())]
    assert not any("hint" in t.lower() for t in items)
    assert not dlg._btn_load.isHidden()


def test_on_finished_updates_current_label(dlg):
    result = _fake_result(records=[SimpleNamespace(status="ok")], failed=[])
    dlg._on_finished(result)
    assert "1 recomputed" in dlg._cur_lbl.text()


# ── _on_error ─────────────────────────────────────────────────────────────


def test_on_error_updates_label_and_log(dlg):
    dlg._on_error("connection lost")
    assert "connection lost" in dlg._cur_lbl.text()
    assert dlg._log_list.count() == 1


# ── _on_load_result ───────────────────────────────────────────────────────


def test_on_load_result_emits_signal_and_accepts(dlg):
    result = _fake_result()
    dlg._result = result
    received = []
    dlg.recompute_committed.connect(received.append)
    dlg._on_load_result()
    assert received == [result]
    assert dlg.result() == QDialog.DialogCode.Accepted


def test_on_load_result_no_result_noop(dlg):
    dlg._result = None
    received = []
    dlg.recompute_committed.connect(received.append)
    dlg._on_load_result()
    assert received == []
    assert dlg.result() != QDialog.DialogCode.Accepted


# ── _reset_progress / _set_running ───────────────────────────────────────


def test_reset_progress_clears_state(dlg):
    dlg._on_station_done(1, 4, "S1", "ok")
    dlg._reset_progress()
    assert dlg._done == 0
    assert dlg._ok == 0
    assert dlg._log_list.count() == 0
    assert dlg._btn_load.isHidden()
    assert dlg._result is None


def test_set_running_true_disables_groups_and_changes_button_text(dlg):
    dlg._set_running(True)
    assert dlg._btn_run.text() == "Stop"
    from PySide6.QtWidgets import QGroupBox

    for gb in dlg.findChildren(QGroupBox):
        if gb.title() in ("Source", "Processing Options", "Output"):
            assert not gb.isEnabled()


def test_set_running_false_reenables_groups(dlg):
    dlg._set_running(True)
    dlg._set_running(False)
    assert dlg._btn_run.text() == "Recompute"
    from PySide6.QtWidgets import QGroupBox

    for gb in dlg.findChildren(QGroupBox):
        if gb.title() in ("Source", "Processing Options", "Output"):
            assert gb.isEnabled()


# ── button wiring ─────────────────────────────────────────────────────────


def test_run_button_triggers_on_run_clicked(dlg):
    with mock.patch.object(dlg, "_on_run_clicked") as m:
        dlg._btn_run.click()
        m.assert_called_once()


def test_close_button_triggers_on_close_clicked(dlg):
    with mock.patch.object(dlg, "_on_close_clicked") as m:
        dlg._btn_close.click()
        m.assert_called_once()
