# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Offline Ed25519 license-key verification.

Product decision (PYCSAMT-DESKTOP-V2.6-MODERNIZATION-PLAN.md, Phase 11 /
Section 4): a license "key" is a compact, JWT-like string --
``<base64url(json payload)>.<base64url(Ed25519 signature)>`` -- signed
offline by earthai-tech with a private key that never ships with the app,
and verified here against the matching public key, which *is* embedded
below. No network call is ever made to verify a key.

Issuing keys is a separate, deliberately-not-shipped concern -- see
``tools/issue_license.py`` in this package, which needs the private key
(kept outside the repository) to run at all.
"""

from __future__ import annotations

import base64
import datetime as _dt
import json
from dataclasses import dataclass
from typing import Optional, Tuple

# ``cryptography`` is imported lazily, inside verify_and_parse(), rather
# than at module scope: this module is reached by importing the
# ``licensing`` package itself (via manager.py), and that package must
# stay importable even in an environment missing the ``desktop`` extra's
# dependencies -- see the package docstring's "no broken license check can
# crash the import path" guarantee.

# Public half of the pycsamt-desktop license-signing keypair (generated
# 2026-09-23 for Phase 11). The matching private key is held by
# earthai-tech outside this repository -- see the ``tools/issue_license.py``
# docstring for how it is used to issue real keys.
PUBLIC_KEY_HEX = (
    "d16b9956ca773741aea0929f1451f8ff5eb46dca82bc7d7cc19d4a8f9ee0cae1"
)


class LicenseError(ValueError):
    """A license key is malformed, its signature is invalid, or it is
    otherwise unusable."""


@dataclass(frozen=True)
class LicenseData:
    """A verified license's payload."""

    customer: str
    seats: int
    issued_at: _dt.datetime
    expires_at: Optional[_dt.datetime]  # None == perpetual license
    features: Tuple[str, ...] = ()

    @property
    def is_expired(self) -> bool:
        if self.expires_at is None:
            return False
        return _dt.datetime.now(_dt.timezone.utc) >= self.expires_at


def _b64url_encode(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).rstrip(b"=").decode("ascii")


def _b64url_decode(text: str) -> bytes:
    padding = "=" * (-len(text) % 4)
    return base64.urlsafe_b64decode(text + padding)


def build_signing_payload(
    customer: str,
    seats: int = 1,
    issued_at: Optional[_dt.datetime] = None,
    expires_at: Optional[_dt.datetime] = None,
    features: Tuple[str, ...] = (),
) -> bytes:
    """The exact bytes that get signed / verified for a license.

    Shared between :func:`verify_and_parse` and the issuing tool so the two
    can never silently drift into incompatible payload shapes.
    """
    issued_at = issued_at or _dt.datetime.now(_dt.timezone.utc)
    payload = {
        "customer": customer,
        "seats": int(seats),
        "issued": issued_at.isoformat(),
        "expires": expires_at.isoformat() if expires_at else None,
        "features": list(features),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )


def encode_license_key(payload_bytes: bytes, signature: bytes) -> str:
    """Assemble the public compact license-key string from its parts."""
    return f"{_b64url_encode(payload_bytes)}.{_b64url_encode(signature)}"


def verify_and_parse(license_key: str) -> LicenseData:
    """Verify *license_key*'s signature and return its parsed payload.

    Raises :class:`LicenseError` for anything not usable: wrong shape, bad
    base64, bad signature, or a malformed payload. Never raises for an
    *expired* license -- callers check :attr:`LicenseData.is_expired`
    themselves, since "expired but genuine" and "forged" are different
    situations worth distinguishing.
    """
    key = (license_key or "").strip()
    parts = key.split(".")
    if len(parts) != 2:
        raise LicenseError("Malformed license key.")
    payload_b64, sig_b64 = parts

    try:
        payload_bytes = _b64url_decode(payload_b64)
        signature = _b64url_decode(sig_b64)
    except Exception as exc:  # noqa: BLE001 - any decode failure is fatal
        raise LicenseError("Malformed license key.") from exc

    try:
        from cryptography.exceptions import InvalidSignature
        from cryptography.hazmat.primitives.asymmetric.ed25519 import (
            Ed25519PublicKey,
        )
    except ImportError as exc:
        raise LicenseError(
            "The 'cryptography' package is required to verify a license "
            "key (part of the 'desktop' extra)."
        ) from exc

    public_key = Ed25519PublicKey.from_public_bytes(
        bytes.fromhex(PUBLIC_KEY_HEX)
    )
    try:
        public_key.verify(signature, payload_bytes)
    except InvalidSignature as exc:
        raise LicenseError("Invalid license signature.") from exc

    try:
        payload = json.loads(payload_bytes.decode("utf-8"))
        customer = str(payload["customer"])
        seats = int(payload.get("seats", 1))
        issued_at = _dt.datetime.fromisoformat(payload["issued"])
        expires_raw = payload.get("expires")
        expires_at = (
            _dt.datetime.fromisoformat(expires_raw) if expires_raw else None
        )
        features = tuple(payload.get("features") or ())
    except (KeyError, ValueError, TypeError) as exc:
        raise LicenseError("Malformed license payload.") from exc

    return LicenseData(
        customer=customer,
        seats=seats,
        issued_at=issued_at,
        expires_at=expires_at,
        features=features,
    )
