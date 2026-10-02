# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.dialogs.license_expired_dialog."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.dialogs.license_expired_dialog import (
    LicenseExpiredDialog,
)
from pycsamt.app.desktop.licensing.interfaces import LicenseStatus


class _FakeManager:
    def __init__(self, status, activate_result=LicenseStatus.LICENSED):
        self._status = status
        self._activate_result = activate_result
        self.activated_with: list[str] = []

    def status(self):
        return self._status

    def trial_state(self):
        return None

    def activate(self, key):
        self.activated_with.append(key)
        return self._activate_result


def test_trial_expired_headline(qapp):
    dlg = LicenseExpiredDialog(_FakeManager(LicenseStatus.TRIAL_EXPIRED))
    assert "trial has ended" in dlg._title_lbl.text()


def test_invalid_status_headline(qapp):
    dlg = LicenseExpiredDialog(_FakeManager(LicenseStatus.INVALID))
    assert "could not be verified" in dlg._title_lbl.text()
    assert dlg.windowTitle().endswith("License required")


def test_unknown_status_gets_generic_headline(qapp):
    dlg = LicenseExpiredDialog(_FakeManager(LicenseStatus.UNKNOWN))
    assert "required to continue" in dlg._title_lbl.text()


def test_exit_button_rejects_dialog(qapp):
    dlg = LicenseExpiredDialog(_FakeManager(LicenseStatus.TRIAL_EXPIRED))
    dlg._exit_btn.click()
    assert dlg.result() == dlg.DialogCode.Rejected


def test_successful_activation_accepts_dialog(qapp):
    mgr = _FakeManager(
        LicenseStatus.TRIAL_EXPIRED, activate_result=LicenseStatus.LICENSED
    )
    dlg = LicenseExpiredDialog(mgr)
    dlg._license_page._key_edit.setText("REAL-KEY")
    dlg._license_page._on_activate()
    assert mgr.activated_with == ["REAL-KEY"]
    assert dlg.result() == dlg.DialogCode.Accepted


def test_failed_activation_does_not_close_dialog(qapp):
    mgr = _FakeManager(
        LicenseStatus.TRIAL_EXPIRED, activate_result=LicenseStatus.INVALID
    )
    dlg = LicenseExpiredDialog(mgr)
    dlg._license_page._key_edit.setText("BAD-KEY")
    dlg._license_page._on_activate()
    assert dlg.result() != dlg.DialogCode.Accepted


def test_buy_button_opens_docs_url(qapp, monkeypatch):
    opened = []
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.license_expired_dialog.QDesktopServices.openUrl",
        lambda url: opened.append(url.toString()),
    )
    dlg = LicenseExpiredDialog(_FakeManager(LicenseStatus.TRIAL_EXPIRED))
    dlg._buy_btn.click()
    assert opened
