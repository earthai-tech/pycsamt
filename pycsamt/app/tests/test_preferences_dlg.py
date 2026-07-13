# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for PreferencesDialog (Phase 5)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.dialogs.preferences_dlg import (
    PreferencesDialog,
)
from pycsamt.app.desktop.models.session import SessionState


@pytest.fixture
def session():
    return SessionState()


@pytest.fixture
def dlg(qapp, session):
    d = PreferencesDialog(session=session)
    yield d
    d.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_prefs_dlg_creates(qapp, session):
    d = PreferencesDialog(session=session)
    assert d is not None
    d.close()


def test_has_four_tabs(dlg):
    assert dlg._tabs.count() == 4


def test_tab_names(dlg):
    labels = [dlg._tabs.tabText(i) for i in range(4)]
    assert "General" in labels
    assert "Solvers" in labels
    assert "AI / LLM" in labels
    assert "Advanced" in labels


# ── General tab ───────────────────────────────────────────────────────────


def test_theme_combo_preselects_session_theme(qapp):
    sess = SessionState(theme="light")
    d = PreferencesDialog(session=sess)
    assert d._theme_combo.currentText() == "light"
    d.close()


def test_max_recent_files_default(dlg):
    assert dlg._max_recent_spin.value() == 20


# ── Solvers tab ───────────────────────────────────────────────────────────


def test_occam2d_edit_prepopulated(qapp):
    sess = SessionState(occam2d_binary="/usr/bin/Occam2D")
    d = PreferencesDialog(session=sess)
    assert d._occam2d_edit.text() == "/usr/bin/Occam2D"
    d.close()


def test_workdir_edit_exists(dlg):
    assert dlg._workdir_edit is not None


# ── AI / LLM tab ─────────────────────────────────────────────────────────


def test_api_key_echo_is_password(dlg):
    from PySide6.QtWidgets import QLineEdit

    assert dlg._api_key_edit.echoMode() == QLineEdit.EchoMode.Password


def test_api_key_prepopulated(qapp):
    sess = SessionState(api_key="sk-ant-test123")
    d = PreferencesDialog(session=sess)
    assert d._api_key_edit.text() == "sk-ant-test123"
    d.close()


def test_provider_change_updates_model_combo(dlg):
    dlg._provider_combo.setCurrentText("openai")
    model_items = [
        dlg._model_combo.itemText(i) for i in range(dlg._model_combo.count())
    ]
    assert any("gpt" in m for m in model_items)

    dlg._provider_combo.setCurrentText("claude")
    model_items = [
        dlg._model_combo.itemText(i) for i in range(dlg._model_combo.count())
    ]
    assert any("claude" in m for m in model_items)


# ── Advanced tab ──────────────────────────────────────────────────────────


def test_log_level_default(qapp):
    sess = SessionState(log_level="WARNING")
    d = PreferencesDialog(session=sess)
    assert d._log_level_combo.currentText() == "WARNING"
    d.close()


def test_tile_provider_preselected(qapp):
    sess = SessionState(tile_provider="Esri.WorldTopoMap")
    d = PreferencesDialog(session=sess)
    assert d._tile_combo.currentText() == "Esri.WorldTopoMap"
    d.close()


# ── Accepted → session updated ────────────────────────────────────────────


def test_accept_updates_session_theme(qapp):
    sess = SessionState(theme="dark")
    d = PreferencesDialog(session=sess)
    d._theme_combo.setCurrentText("light")
    d._on_accepted()
    assert sess.theme == "light"
    d.close()


def test_accept_updates_session_api_key(qapp):
    sess = SessionState()
    d = PreferencesDialog(session=sess)
    d._api_key_edit.setText("new-key-xyz")
    d._on_accepted()
    assert sess.api_key == "new-key-xyz"
    d.close()


def test_accept_updates_solver_paths(qapp):
    sess = SessionState()
    d = PreferencesDialog(session=sess)
    d._occam2d_edit.setText("/opt/occam/Occam2D")
    d._on_accepted()
    assert sess.occam2d_binary == "/opt/occam/Occam2D"
    d.close()


def test_accept_updates_log_level(qapp):
    sess = SessionState()
    d = PreferencesDialog(session=sess)
    d._log_level_combo.setCurrentText("DEBUG")
    d._on_accepted()
    assert sess.log_level == "DEBUG"
    d.close()
