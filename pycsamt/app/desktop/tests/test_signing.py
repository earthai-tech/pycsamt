# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.licensing.signing.

Uses a throwaway Ed25519 keypair generated in-process (NOT the real
embedded ``PUBLIC_KEY_HEX``) to exercise the happy path deterministically,
plus real-key smoke tests to confirm the actually-embedded public key is
internally consistent with itself.
"""

from __future__ import annotations

import datetime as dt

import pytest

pytest.importorskip("cryptography", reason="cryptography required")

from cryptography.hazmat.primitives.asymmetric.ed25519 import (
    Ed25519PrivateKey,
)

from pycsamt.app.desktop.licensing import signing


def _sign_with(private_key, **payload_kwargs):
    payload_bytes = signing.build_signing_payload(**payload_kwargs)
    signature = private_key.sign(payload_bytes)
    return signing.encode_license_key(payload_bytes, signature)


@pytest.fixture()
def keypair():
    priv = Ed25519PrivateKey.generate()
    return priv, priv.public_key()


def _patched_public_key(monkeypatch, public_key):
    from cryptography.hazmat.primitives import serialization

    raw = public_key.public_bytes(
        encoding=serialization.Encoding.Raw,
        format=serialization.PublicFormat.Raw,
    )
    monkeypatch.setattr(signing, "PUBLIC_KEY_HEX", raw.hex())


def test_verify_and_parse_accepts_a_validly_signed_key(monkeypatch, keypair):
    priv, pub = keypair
    _patched_public_key(monkeypatch, pub)
    key = _sign_with(priv, customer="Acme", seats=3)

    lic = signing.verify_and_parse(key)
    assert lic.customer == "Acme"
    assert lic.seats == 3
    assert lic.expires_at is None
    assert not lic.is_expired


def test_verify_and_parse_reports_expiry(monkeypatch, keypair):
    priv, pub = keypair
    _patched_public_key(monkeypatch, pub)
    past = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)
    key = _sign_with(priv, customer="Acme", expires_at=past)

    lic = signing.verify_and_parse(key)
    assert lic.is_expired


def test_verify_and_parse_rejects_tampered_payload(monkeypatch, keypair):
    priv, pub = keypair
    _patched_public_key(monkeypatch, pub)
    key = _sign_with(priv, customer="Acme")
    payload_b64, sig_b64 = key.split(".")
    tampered = signing._b64url_encode(
        signing._b64url_decode(payload_b64).replace(b"Acme", b"Evil")
    )
    with pytest.raises(signing.LicenseError):
        signing.verify_and_parse(f"{tampered}.{sig_b64}")


def test_verify_and_parse_rejects_signature_from_a_different_key(monkeypatch, keypair):
    _priv, pub = keypair
    other_priv = Ed25519PrivateKey.generate()
    _patched_public_key(monkeypatch, pub)
    key = _sign_with(other_priv, customer="Acme")
    with pytest.raises(signing.LicenseError):
        signing.verify_and_parse(key)


@pytest.mark.parametrize(
    "bad_key",
    ["", "not-a-license-key", "only.one.part.too.many", "abc"],
)
def test_verify_and_parse_rejects_malformed_keys(bad_key):
    with pytest.raises(signing.LicenseError):
        signing.verify_and_parse(bad_key)


def test_embedded_public_key_hex_is_a_valid_32_byte_ed25519_key():
    # The matching private half is intentionally not in this repository
    # (see tools/issue_license.py) -- this only confirms the embedded
    # constant round-trips through the cryptography API without raising.
    from cryptography.hazmat.primitives.asymmetric.ed25519 import (
        Ed25519PublicKey,
    )

    raw = bytes.fromhex(signing.PUBLIC_KEY_HEX)
    assert len(raw) == 32
    Ed25519PublicKey.from_public_bytes(raw)  # raises if malformed
