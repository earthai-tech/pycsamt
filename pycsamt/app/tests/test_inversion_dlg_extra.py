# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Extra tests for InversionWizardDialog, on top of test_inversion_dlg.py
(which covers construction and basic page-0/1 navigation).

This file drives the launch/stop/finish/error flow and the data-page
sites-selection logic. InputBuilder/OccamConfig (pycsamt.models.occam2d)
and InversionWorker are replaced with lightweight stand-ins so no real
Occam2D binary or input-file build is required.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.dialogs.inversion_dlg import InversionWizardDialog
from pycsamt.app.desktop.models.session import SessionState


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls():
    class _FakeWorker:
        instances = []

        def __init__(self, **kw):
            self.kw = kw
            self.stdout_line = _FakeSignal()
            self.progress = _FakeSignal()
            self.finished = _FakeSignal()
            self.error = _FakeSignal()
            self.started = False
            self.cancelled = False
            self.quit_called = False
            _FakeWorker.instances.append(self)

        def start(self):
            self.started = True

        def cancel(self):
            self.cancelled = True

        def quit(self):
            self.quit_called = True

    return _FakeWorker


@pytest.fixture
def fake_sites():
    return SimpleNamespace(
        as_list=lambda: [SimpleNamespace(name="S1"), SimpleNamespace(name="S2")],
        select=mock.Mock(return_value="SELECTED_SUBSET"),
    )


@pytest.fixture
def dlg(qapp):
    d = InversionWizardDialog(sites=None, session=SessionState())
    yield d
    d.close()


@pytest.fixture
def dlg_with_sites(qapp, fake_sites):
    d = InversionWizardDialog(sites=fake_sites, session=SessionState())
    yield d
    d.close()


# ── Page 0: station list population ─────────────────────────────────────


def test_station_list_populated_from_sites(dlg_with_sites):
    items = [
        dlg_with_sites._station_list.item(i).text()
        for i in range(dlg_with_sites._station_list.count())
    ]
    assert items == ["S1", "S2"]


def test_station_list_placeholder_when_no_sites(dlg):
    assert dlg._station_list.count() == 1
    item = dlg._station_list.item(0)
    assert "no stations" in item.text()


# ── _browse_workdir ───────────────────────────────────────────────────────


def test_browse_workdir_sets_text(dlg, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    dlg._browse_workdir()
    assert dlg._workdir_edit.text() == str(tmp_path)


def test_browse_workdir_cancelled_noop(dlg, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: ""),
    )
    dlg._workdir_edit.setText("keepme")
    dlg._browse_workdir()
    assert dlg._workdir_edit.text() == "keepme"


# ── Navigation boundaries ────────────────────────────────────────────────


def test_go_prev_at_first_page_noop(dlg):
    dlg._go_prev()
    assert dlg._stack.currentIndex() == 0


def test_go_next_at_last_page_noop(dlg, tmp_path):
    dlg._workdir_edit.setText(str(tmp_path))
    for _ in range(4):
        dlg._go_next()
    assert dlg._stack.currentIndex() == 3
    dlg._go_next()
    assert dlg._stack.currentIndex() == 3


def test_header_at_last_page_shows_close_hides_next(dlg, tmp_path):
    dlg._workdir_edit.setText(str(tmp_path))
    for _ in range(3):
        dlg._go_next()
    assert dlg._stack.currentIndex() == 3
    assert not dlg._btn_close.isHidden()
    assert dlg._btn_next.isHidden()


def test_validate_page_nonzero_always_true(dlg):
    assert dlg._validate_page(1) is True
    assert dlg._validate_page(2) is True
    assert dlg._validate_page(3) is True


# ── _selected_sites ───────────────────────────────────────────────────────


def test_selected_sites_none_when_no_sites(dlg):
    assert dlg._selected_sites() is None


def test_selected_sites_returns_all_when_no_selection(dlg_with_sites, fake_sites):
    assert dlg_with_sites._selected_sites() is fake_sites


def test_selected_sites_calls_select_on_selection(dlg_with_sites, fake_sites):
    dlg_with_sites._station_list.item(0).setSelected(True)
    result = dlg_with_sites._selected_sites()
    assert result == "SELECTED_SUBSET"
    fake_sites.select.assert_called_once_with(["S1"])


def test_selected_sites_select_raises_falls_back(dlg_with_sites, fake_sites):
    fake_sites.select.side_effect = RuntimeError("boom")
    dlg_with_sites._station_list.item(0).setSelected(True)
    assert dlg_with_sites._selected_sites() is fake_sites


# ── _build_occam_config ───────────────────────────────────────────────────


def test_build_occam_config_reflects_form_values(dlg, monkeypatch):
    fake_cfg_cls = mock.Mock()
    monkeypatch.setattr(
        "pycsamt.models.occam2d.OccamConfig", fake_cfg_cls
    )
    dlg._mode_tm.setChecked(False)
    dlg._n_layers.setValue(40)
    dlg._build_occam_config()
    _, kwargs = fake_cfg_cls.call_args
    assert kwargs["modes"] == ["te"]
    assert kwargs["n_layers"] == 40


# ── _on_launch ────────────────────────────────────────────────────────────


def test_on_launch_no_workdir_logs_error(dlg):
    dlg._workdir_edit.setText("")
    dlg._on_launch()
    assert "no working directory" in dlg._run_log.toPlainText()
    assert dlg._worker is None


def test_on_launch_input_builder_raises_logs_error(dlg, tmp_path, monkeypatch):
    fake_builder_cls = mock.Mock(
        side_effect=RuntimeError("bad input")
    )
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InputBuilder", fake_builder_cls
    )
    monkeypatch.setattr(
        "pycsamt.models.occam2d.OccamConfig", mock.Mock()
    )
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._on_launch()
    assert "InputBuilder error" in dlg._run_log.toPlainText()
    assert dlg._worker is None


def test_on_launch_success_starts_worker(dlg, tmp_path, monkeypatch):
    fake_builder_instance = mock.Mock()
    fake_builder_cls = mock.Mock(return_value=fake_builder_instance)
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InputBuilder", fake_builder_cls
    )
    monkeypatch.setattr(
        "pycsamt.models.occam2d.OccamConfig", mock.Mock()
    )
    fake_worker_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
        fake_worker_cls,
    )
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._on_launch()

    fake_builder_instance.build.assert_called_once()
    assert dlg._worker is not None
    assert dlg._worker.started is True
    assert not dlg._btn_launch.isEnabled()
    assert dlg._btn_stop_inv.isEnabled()
    assert not dlg._run_progress.isHidden()


def test_on_launch_uses_session_binary(dlg, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InputBuilder", mock.Mock(return_value=mock.Mock())
    )
    monkeypatch.setattr("pycsamt.models.occam2d.OccamConfig", mock.Mock())
    fake_worker_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
        fake_worker_cls,
    )
    dlg._session.occam2d_binary = "/opt/occam2d"
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._on_launch()
    assert dlg._worker.kw["binary_path"] == "/opt/occam2d"


def test_on_launch_no_session_binary_is_none(dlg, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InputBuilder", mock.Mock(return_value=mock.Mock())
    )
    monkeypatch.setattr("pycsamt.models.occam2d.OccamConfig", mock.Mock())
    fake_worker_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
        fake_worker_cls,
    )
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._on_launch()
    assert dlg._worker.kw["binary_path"] is None


# ── _on_stop_inversion ────────────────────────────────────────────────────


def test_on_stop_inversion_with_worker(dlg, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InputBuilder", mock.Mock(return_value=mock.Mock())
    )
    monkeypatch.setattr("pycsamt.models.occam2d.OccamConfig", mock.Mock())
    fake_worker_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
        fake_worker_cls,
    )
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._on_launch()
    worker = dlg._worker
    dlg._on_stop_inversion()
    assert worker.cancelled is True
    assert worker.quit_called is True
    assert dlg._btn_launch.isEnabled()
    assert not dlg._btn_stop_inv.isEnabled()
    assert dlg._run_progress.isHidden()
    assert "Stopped." in dlg._run_log.toPlainText()


def test_on_stop_inversion_without_worker_still_toggles_buttons(dlg):
    dlg._on_stop_inversion()  # _worker is None -> `if self._worker:` False branch
    assert dlg._btn_launch.isEnabled()
    assert not dlg._btn_stop_inv.isEnabled()


# ── _on_inversion_finished / _on_inversion_error ─────────────────────────


def test_on_inversion_finished_advances_to_results_page(dlg, tmp_path):
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._go_next()  # -> page 1
    dlg._go_next()  # -> page 2 (Run)
    result = SimpleNamespace(summary="RMS=1.2, converged")
    dlg._on_inversion_finished(result)
    assert dlg._stack.currentIndex() == 3
    assert dlg._result is result
    assert dlg._result_text.toPlainText() == "RMS=1.2, converged"
    assert dlg._btn_open_viewer.isEnabled()
    assert dlg._run_progress.isHidden()
    assert dlg._btn_launch.isEnabled()
    assert not dlg._btn_stop_inv.isEnabled()


def test_on_inversion_finished_result_without_summary_falls_back(dlg, tmp_path):
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._go_next()
    dlg._go_next()
    result = object()  # no .summary attribute
    dlg._on_inversion_finished(result)
    assert dlg._result_text.toPlainText() == "Inversion completed successfully."


def test_on_inversion_error_logs_and_resets_buttons(dlg):
    dlg._btn_launch.setEnabled(False)
    dlg._btn_stop_inv.setEnabled(True)
    dlg._run_progress.setVisible(True)
    dlg._on_inversion_error("solver crashed")
    assert "ERROR: solver crashed" in dlg._run_log.toPlainText()
    assert dlg._btn_launch.isEnabled()
    assert not dlg._btn_stop_inv.isEnabled()
    assert dlg._run_progress.isHidden()
