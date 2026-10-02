# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.licensing.fingerprint."""

from __future__ import annotations

from pycsamt.app.desktop.licensing import fingerprint as fp_mod
from pycsamt.app.desktop.licensing.fingerprint import machine_fingerprint


def test_machine_fingerprint_is_stable_across_calls():
    assert machine_fingerprint() == machine_fingerprint()


def test_machine_fingerprint_is_a_32_char_hex_string():
    fp = machine_fingerprint()
    assert len(fp) == 32
    int(fp, 16)  # raises ValueError if not valid hex


def test_fingerprint_unaffected_by_mac_address_changing(monkeypatch):
    """Regression test for a real bug: the original implementation hashed
    in ``uuid.getnode()`` (the MAC address) directly, which silently
    changed on the same physical machine after nothing more than using
    WSL2 in between two launches (it spins up a virtual network adapter,
    which changed which MAC Windows' enumeration returned first). Since
    the trial tracker is fail-closed, that looked exactly like tampering
    and killed a legitimate trial. When a stable OS-level id is available,
    the fingerprint must not move just because the MAC address does.
    """
    monkeypatch.setattr(fp_mod, "_stable_machine_id", lambda: "fixed-machine-id")
    monkeypatch.setattr(fp_mod.uuid, "getnode", lambda: 0x001122334455)
    first = machine_fingerprint()
    monkeypatch.setattr(fp_mod.uuid, "getnode", lambda: 0xAABBCCDDEEFF)
    second = machine_fingerprint()
    assert first == second


def test_fingerprint_falls_back_to_mac_when_no_stable_id(monkeypatch):
    """When no OS-level id can be read at all, the old MAC-based scheme is
    still the fallback -- better than crashing or a constant fingerprint
    that can't distinguish machines."""
    monkeypatch.setattr(fp_mod, "_stable_machine_id", lambda: None)
    monkeypatch.setattr(fp_mod.uuid, "getnode", lambda: 0x001122334455)
    first = machine_fingerprint()
    monkeypatch.setattr(fp_mod.uuid, "getnode", lambda: 0xAABBCCDDEEFF)
    second = machine_fingerprint()
    assert first != second


def test_stable_machine_id_never_raises(monkeypatch):
    """Whatever the platform/permissions, callers rely on this returning
    ``None`` rather than propagating an exception."""
    assert fp_mod._stable_machine_id() is None or isinstance(
        fp_mod._stable_machine_id(), str
    )
