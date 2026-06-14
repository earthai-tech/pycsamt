# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for RunAgentDialog (Phase 4)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.dialogs.run_agent_dlg import RunAgentDialog, _ParamPage
from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY, agent_names


@pytest.fixture
def dlg(qapp):
    d = RunAgentDialog()
    yield d
    d.close()


# ── Construction ──────────────────────────────────────────────────────────

def test_run_agent_dlg_creates(qapp):
    d = RunAgentDialog()
    assert d is not None
    d.close()


def test_list_populated(dlg):
    # Should have at least as many items as agents (plus group headers)
    assert dlg._list.count() >= len(AGENT_REGISTRY)


def test_stack_has_one_page_per_agent(dlg):
    assert dlg._stack.count() == len(AGENT_REGISTRY)


def test_api_key_field_exists(dlg):
    assert dlg._api_key_edit is not None


def test_api_key_prepopulated(qapp):
    d = RunAgentDialog(api_key="test-key-123")
    assert d._api_key_edit.text() == "test-key-123"
    d.close()


def test_api_key_echo_mode_password(dlg):
    from PySide6.QtWidgets import QLineEdit
    assert dlg._api_key_edit.echoMode() == QLineEdit.EchoMode.Password


# ── _ParamPage ────────────────────────────────────────────────────────────

def test_param_page_creates_for_all_agents(qapp):
    for name in agent_names():
        entry = AGENT_REGISTRY[name]
        page = _ParamPage(entry)
        assert page is not None
        page.close()


def test_param_page_current_params_returns_dict(qapp):
    entry = AGENT_REGISTRY["QC Analysis"]
    page = _ParamPage(entry)
    params = page.current_params()
    assert isinstance(params, dict)
    page.close()


def test_param_page_defaults_match_registry(qapp):
    entry = AGENT_REGISTRY["QC Analysis"]
    page = _ParamPage(entry)
    params = page.current_params()
    for p_name, spec in entry["params"].items():
        assert p_name in params
        expected = spec["default"]
        # Floats may differ by spinbox precision — allow close comparison
        if isinstance(expected, float):
            assert abs(params[p_name] - expected) < 0.001
        else:
            assert params[p_name] == expected
    page.close()


def test_param_page_no_params(qapp):
    entry = AGENT_REGISTRY["QC Quicklook"]
    page = _ParamPage(entry)
    assert page.current_params() == {}
    page.close()


# ── Agent selection ───────────────────────────────────────────────────────

def test_first_item_selected_on_open(dlg):
    # First selectable item should have the first agent selected
    assert dlg._stack.currentIndex() >= 0


def test_agent_selection_switches_page(dlg):
    # Find the list item for "QC Quicklook"
    for i in range(dlg._list.count()):
        item = dlg._list.item(i)
        from PySide6.QtCore import Qt
        if item.data(Qt.ItemDataRole.UserRole) == "QC Quicklook":
            dlg._list.setCurrentRow(i)
            current_page = dlg._stack.currentWidget()
            assert current_page is dlg._pages.get("QC Quicklook")
            break
