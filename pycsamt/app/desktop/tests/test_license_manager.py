# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.licensing.manager.OfflineLicenseManager.

Uses a throwaway Ed25519 keypair (monkeypatched over ``signing.
PUBLIC_KEY_HEX``, same technique as test_signing.py) plus a dict-backed
fake settings object -- no QSettings/Qt, no real filesystem/registry
writes.
"""

from __future__ import annotations

import datetime as dt

import pytest

pytest.importorskip("cryptography", reason="cryptography required")

from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import (
    Ed25519PrivateKey,
)

from pycsamt.app.desktop.licensing import manager as manager_mod
from pycsamt.app.desktop.licensing import signing
from pycsamt.app.desktop.licensing.interfaces import LicenseStatus
from pycsamt.app.desktop.licensing.manager import (
    OfflineLicenseManager,
    get_default_manager,
)


class _FakeSettings:
    def __init__(self):
        self._data: dict[str, str] = {}

    def value(self, key, default="", type=None):  # noqa: A002
        return self._data.get(key, default)

    def setValue(self, key, value):
        self._data[key] = value

    def sync(self):
        pass


@pytest.fixture()
def signed_key(monkeypatch):
    priv = Ed25519PrivateKey.generate()
    pub = priv.public_key()
    raw = pub.public_bytes(
        encoding=serialization.Encoding.Raw,
        format=serialization.PublicFormat.Raw,
    )
    monkeypatch.setattr(signing, "PUBLIC_KEY_HEX", raw.hex())

    def _make(**kwargs):
        payload_bytes = signing.build_signing_payload(**kwargs)
        sig = priv.sign(payload_bytes)
        return signing.encode_license_key(payload_bytes, sig)

    return _make


def test_no_stored_key_reports_trial_active_on_first_use():
    mgr = OfflineLicenseManager(_FakeSettings())
    assert mgr.status() == LicenseStatus.TRIAL_ACTIVE
    assert mgr.trial_state() is not None


def test_activate_with_valid_key_persists_and_reports_licensed(signed_key):
    settings = _FakeSettings()
    mgr = OfflineLicenseManager(settings)
    key = signed_key(customer="Acme")

    result = mgr.activate(key)
    assert result == LicenseStatus.LICENSED
    assert mgr.status() == LicenseStatus.LICENSED
    assert mgr.trial_state() is None
    assert settings.value("license/key") == key


def test_activate_with_garbage_key_reports_invalid_and_does_not_persist():
    settings = _FakeSettings()
    mgr = OfflineLicenseManager(settings)

    result = mgr.activate("not-a-real-key")
    assert result == LicenseStatus.INVALID
    assert settings.value("license/key") == ""
    assert mgr.status() == LicenseStatus.TRIAL_ACTIVE  # falls back to trial


def test_activate_with_already_expired_key_reports_invalid(signed_key):
    settings = _FakeSettings()
    mgr = OfflineLicenseManager(settings)
    past = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)
    key = signed_key(customer="Acme", expires_at=past)

    result = mgr.activate(key)
    assert result == LicenseStatus.INVALID
    assert settings.value("license/key") == ""


def test_stored_expired_key_reports_invalid_on_status_check(signed_key):
    settings = _FakeSettings()
    # Sign a key that is valid *now* so activate() accepts it, then move
    # time forward by monkeypatching a pre-expired one directly into
    # storage to exercise the "previously valid, now expired" read path.
    past = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)
    key = signed_key(customer="Acme", expires_at=past)
    settings.setValue("license/key", key)

    mgr = OfflineLicenseManager(settings)
    assert mgr.status() == LicenseStatus.INVALID


def test_stored_corrupted_key_reports_invalid_not_a_crash():
    settings = _FakeSettings()
    settings.setValue("license/key", "garbage.garbage")
    mgr = OfflineLicenseManager(settings)
    assert mgr.status() == LicenseStatus.INVALID


def test_get_default_manager_returns_the_same_instance(monkeypatch):
    # Fake out OfflineLicenseManager itself so this never touches a real
    # QSettings (the real default writes into the actual OS user-config
    # store under branding.ORG_NAME/APP_NAME -- not something a test run
    # should do as a side effect).
    built = []

    class _FakeOfflineLicenseManager:
        def __init__(self):
            built.append(self)

    monkeypatch.setattr(manager_mod, "OfflineLicenseManager", _FakeOfflineLicenseManager)
    manager_mod._default_manager = None
    try:
        first = get_default_manager()
        second = get_default_manager()
        assert first is second
        assert len(built) == 1
    finally:
        manager_mod._default_manager = None
