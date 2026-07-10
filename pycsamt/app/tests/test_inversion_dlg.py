# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for InversionWizardDialog and InversionWorker (Phase 5)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.dialogs.inversion_dlg import (
    InversionWizardDialog,
)
from pycsamt.app.desktop.models.session import SessionState
from pycsamt.app.desktop.workers.inversion_worker import (
    InversionWorker,
)

# ── InversionWizardDialog ─────────────────────────────────────────────────

@pytest.fixture
def dlg(qapp):
    d = InversionWizardDialog(sites=None, session=SessionState())
    yield d
    d.close()


def test_inversion_dlg_creates(qapp):
    d = InversionWizardDialog()
    assert d is not None
    d.close()


def test_initial_page_is_zero(dlg):
    assert dlg._stack.currentIndex() == 0


def test_has_four_pages(dlg):
    assert dlg._stack.count() == 4


def test_header_shows_step_1(dlg):
    assert "1" in dlg._header_lbl.text() or "Data" in dlg._header_lbl.text()


def test_prev_button_disabled_on_first_page(dlg):
    assert not dlg._btn_prev.isEnabled()


def test_next_button_not_hidden_on_first_page(dlg):
    # isVisible() requires parent to be shown; use isHidden() to test own flag
    assert not dlg._btn_next.isHidden()


def test_close_button_hidden_on_first_page(dlg):
    assert not dlg._btn_close.isVisible()


def test_navigation_to_page_1_requires_workdir(dlg):
    # No workdir set → next should stay on page 0
    dlg._workdir_edit.setText("")
    dlg._go_next()
    assert dlg._stack.currentIndex() == 0


def test_navigation_to_page_1_with_workdir(dlg, tmp_path):
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._go_next()
    assert dlg._stack.currentIndex() == 1


def test_go_prev_decrements_page(dlg, tmp_path):
    dlg._workdir_edit.setText(str(tmp_path))
    dlg._go_next()   # page 0 → 1
    dlg._go_prev()   # page 1 → 0
    assert dlg._stack.currentIndex() == 0


def test_config_defaults(dlg):
    assert dlg._max_iter.value() == 100
    assert dlg._target_rms.value() == pytest.approx(1.0)
    assert dlg._initial_rho.value() == pytest.approx(100.0)
    assert dlg._n_layers.value() == 30


def test_mode_checkboxes_default_checked(dlg):
    assert dlg._mode_te.isChecked()
    assert dlg._mode_tm.isChecked()


def test_launch_button_exists(dlg):
    assert dlg._btn_launch is not None


def test_stop_button_disabled_initially(dlg):
    assert not dlg._btn_stop_inv.isEnabled()


def test_get_result_is_none_initially(dlg):
    assert dlg.get_result() is None


# ── InversionWorker ───────────────────────────────────────────────────────

def test_inversion_worker_creates(qapp, tmp_path):
    w = InversionWorker(workdir=str(tmp_path))
    assert w is not None


def test_inversion_worker_cancel(qapp, tmp_path):
    w = InversionWorker(workdir=str(tmp_path))
    w.cancel()
    assert w._cancelled


def test_inversion_worker_missing_workdir_emits_error(qapp, tmp_path):
    errors = []
    missing = str(tmp_path / "no_such_dir")
    w = InversionWorker(workdir=missing)
    w.error.connect(errors.append)
    w.run()
    assert len(errors) == 1
    assert "not found" in errors[0].lower() or "no_such_dir" in errors[0]


def test_inversion_worker_missing_binary_emits_error(qapp, tmp_path):
    errors = []
    w = InversionWorker(
        workdir=str(tmp_path),
        binary_path=str(tmp_path / "fake_binary"),
    )
    w.error.connect(errors.append)
    w.run()
    assert len(errors) == 1


def test_inversion_worker_has_signals(qapp, tmp_path):
    w = InversionWorker(workdir=str(tmp_path))
    assert hasattr(w, "stdout_line")
    assert hasattr(w, "progress")
    assert hasattr(w, "finished")
    assert hasattr(w, "error")


# ── SessionState new fields ───────────────────────────────────────────────

def test_session_new_fields_default():
    s = SessionState()
    assert s.api_key == ""
    assert s.occam2d_binary == ""
    assert s.modem_binary == ""
    assert s.mare2dem_binary == ""
    assert s.inversion_workdir == ""
    assert s.log_level == "WARNING"
    assert s.tile_provider == "OpenStreetMap.Mapnik"
    assert s.max_recent_files == 20


def test_session_new_fields_roundtrip(tmp_path):
    p = tmp_path / "session.json"
    s = SessionState(
        api_key="sk-ant-test",
        occam2d_binary="/opt/Occam2D",
        log_level="DEBUG",
        tile_provider="Esri.WorldTopoMap",
        max_recent_files=15,
    )
    s.save(p)
    s2 = SessionState.load(p)
    assert s2.api_key == "sk-ant-test"
    assert s2.occam2d_binary == "/opt/Occam2D"
    assert s2.log_level == "DEBUG"
    assert s2.tile_provider == "Esri.WorldTopoMap"
    assert s2.max_recent_files == 15
