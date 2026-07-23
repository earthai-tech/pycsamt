# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for NoDataDialog."""

from __future__ import annotations

from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QDialog

from pycsamt.app.desktop.dialogs.no_data_dialog import NoDataDialog


def test_creates_with_tool_name(qapp):
    dlg = NoDataDialog("Strike Analyzer")
    assert dlg.windowTitle() == "No Survey Data Loaded"
    assert dlg._tool_name == "Strike Analyzer"
    dlg.close()


def test_creates_without_tool_name(qapp):
    dlg = NoDataDialog()
    assert dlg._tool_name == ""
    dlg.close()


def test_fixed_width(qapp):
    dlg = NoDataDialog()
    assert dlg.maximumWidth() == 420
    dlg.close()


def test_load_button_accepts(qapp):
    dlg = NoDataDialog("X")
    dlg._load_btn.click()
    assert dlg.result() == QDialog.DialogCode.Accepted
    dlg.close()


def test_cancel_button_rejects(qapp):
    from PySide6.QtWidgets import QPushButton

    dlg = NoDataDialog("X")
    buttons = dlg.findChildren(QPushButton)
    cancel_btn = next(b for b in buttons if b.text() == "Cancel")
    cancel_btn.click()
    assert dlg.result() == QDialog.DialogCode.Rejected
    dlg.close()


def test_require_returns_true_when_accepted(qapp):
    with mock.patch.object(
        NoDataDialog, "exec", return_value=QDialog.DialogCode.Accepted
    ):
        assert NoDataDialog.require(None, "Y") is True


def test_require_returns_false_when_rejected(qapp):
    with mock.patch.object(
        NoDataDialog, "exec", return_value=QDialog.DialogCode.Rejected
    ):
        assert NoDataDialog.require(None, "Y") is False


def test_require_default_tool_name(qapp):
    with mock.patch.object(
        NoDataDialog, "exec", return_value=QDialog.DialogCode.Rejected
    ):
        assert NoDataDialog.require(None) is False
