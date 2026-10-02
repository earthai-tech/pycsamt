# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for LicensePage (pycsamt.app.desktop.widgets.license_page).

Strategy
--------
* Default (no ``license_manager`` given) exercises the real
  ``NullLicenseManager`` — the current, honest "free trial for now"
  product state (unlimited trial, every key rejected).
* A fake ``LicenseManager``-shaped object drives the other status
  branches (``TRIAL_EXPIRED``/``LICENSED``/``INVALID``/``UNKNOWN``) that
  ``NullLicenseManager`` can never itself produce.
"""

from __future__ import annotations

import datetime as dt

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.licensing.interfaces import LicenseStatus, TrialState
from pycsamt.app.desktop.licensing.null_manager import NullLicenseManager
from pycsamt.app.desktop.widgets.license_page import (
    LicensePage,
    status_badge_html,
)


class _FakeManager:
    def __init__(self, status, trial=None, activate_result=None):
        self._status = status
        self._trial = trial
        self._activate_result = (
            activate_result if activate_result is not None else status
        )
        self.activated_with: list[str] = []

    def status(self):
        return self._status

    def trial_state(self):
        return self._trial

    def activate(self, license_key):
        self.activated_with.append(license_key)
        self._status = self._activate_result
        return self._activate_result


# ── status_badge_html ────────────────────────────────────────────────────────


def test_status_badge_html_covers_every_status():
    for status in LicenseStatus:
        html = status_badge_html(status)
        assert "<span" in html
        assert "#" in html  # a colour is present


def test_status_badge_html_unknown_status_falls_back_gracefully():
    html = status_badge_html(LicenseStatus.UNKNOWN)
    assert "Unknown" in html


# ── Default (NullLicenseManager) ─────────────────────────────────────────────


def test_defaults_to_null_license_manager(qapp):
    page = LicensePage()
    assert isinstance(page._manager, NullLicenseManager)


def test_null_manager_shows_trial_active(qapp):
    page = LicensePage()
    assert "Trial active" in page._status_lbl.text()
    assert page._trial_lbl.text() == ""  # NullLicenseManager.trial_state() is None
    assert page._btn_activate.isEnabled()
    assert page._key_edit.isEnabled()


def test_null_manager_activate_always_reports_invalid(qapp):
    page = LicensePage()
    page._key_edit.setText("SOME-KEY")
    page._on_activate()
    assert "could not be activated" in page._activate_msg.text()


def test_activate_with_empty_key_shows_hint_without_calling_manager(qapp):
    calls = []
    page = LicensePage(license_manager=_FakeManager(LicenseStatus.TRIAL_ACTIVE))
    page._manager.activate = lambda k: calls.append(k)
    page._key_edit.setText("")
    page._on_activate()
    assert "Enter a license key" in page._activate_msg.text()
    assert calls == []


# ── Fake manager: every status branch ────────────────────────────────────────


def test_trial_active_with_days_remaining(qapp):
    trial = TrialState(
        started_at=dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1),
        expires_at=dt.datetime.now(dt.timezone.utc) + dt.timedelta(days=6),
    )
    page = LicensePage(
        license_manager=_FakeManager(LicenseStatus.TRIAL_ACTIVE, trial=trial)
    )
    assert "remaining in your trial" in page._trial_lbl.text()


def test_trial_expired_shows_ended_message(qapp):
    trial = TrialState(
        started_at=dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=40),
        expires_at=dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=10),
    )
    page = LicensePage(
        license_manager=_FakeManager(LicenseStatus.TRIAL_EXPIRED, trial=trial)
    )
    assert page._trial_lbl.text() == "Your trial has ended."


def test_licensed_status_disables_key_entry_and_activate(qapp):
    page = LicensePage(license_manager=_FakeManager(LicenseStatus.LICENSED))
    assert "Thank you for licensing" in page._trial_lbl.text()
    assert not page._btn_activate.isEnabled()
    assert not page._key_edit.isEnabled()


def test_successful_activation_clears_key_and_shows_thanks(qapp):
    mgr = _FakeManager(
        LicenseStatus.TRIAL_ACTIVE, activate_result=LicenseStatus.LICENSED
    )
    page = LicensePage(license_manager=mgr)
    page._key_edit.setText("REAL-KEY")
    page._on_activate()
    assert mgr.activated_with == ["REAL-KEY"]
    assert "activated" in page._activate_msg.text().lower()
    assert page._key_edit.text() == ""
    # refresh() ran again after activation -- status now reflects LICENSED
    assert not page._btn_activate.isEnabled()


def test_refresh_can_be_called_repeatedly(qapp):
    page = LicensePage()
    page.refresh()
    page.refresh()
    assert "Trial active" in page._status_lbl.text()
