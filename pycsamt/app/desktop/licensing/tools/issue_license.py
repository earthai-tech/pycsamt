# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
issue_license.py -- earthai-tech's own license-signing CLI.

Signs a license payload with the *private* half of the Ed25519 keypair
whose public half is embedded in ``pycsamt.app.desktop.licensing.signing``.
That private key is **not** in this repository. Supply it via the
``PYCSAMT_LICENSE_PRIVATE_KEY`` environment variable (64 hex chars, 32
raw bytes) -- never pass it as a bare CLI argument (shell history / process
list exposure) and never commit it anywhere.

Usage
-----
    PYCSAMT_LICENSE_PRIVATE_KEY=<hex> python -m \\
        pycsamt.app.desktop.licensing.tools.issue_license \\
        --customer "Acme Geophysics" --seats 5 --days 365

    # Perpetual license (no --days / --expires):
    PYCSAMT_LICENSE_PRIVATE_KEY=<hex> python -m \\
        pycsamt.app.desktop.licensing.tools.issue_license \\
        --customer "Acme Geophysics" --seats 1

Prints the compact license-key string to stdout -- paste it into
Settings ▸ License in the app, or the trial-expiry dialog's key field.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys

from cryptography.hazmat.primitives.asymmetric.ed25519 import (
    Ed25519PrivateKey,
)

from ..signing import build_signing_payload, encode_license_key

_ENV_VAR = "PYCSAMT_LICENSE_PRIVATE_KEY"


def _load_private_key() -> Ed25519PrivateKey:
    hex_key = os.environ.get(_ENV_VAR, "")
    if not hex_key:
        raise SystemExit(
            f"Set {_ENV_VAR} to the signing private key (hex) before "
            "running this tool. It is intentionally not stored in this "
            "repository."
        )
    try:
        raw = bytes.fromhex(hex_key)
    except ValueError as exc:
        raise SystemExit(f"{_ENV_VAR} is not valid hex.") from exc
    return Ed25519PrivateKey.from_private_bytes(raw)


def issue(
    customer: str,
    seats: int = 1,
    days: int | None = None,
    features: tuple[str, ...] = (),
    private_key: Ed25519PrivateKey | None = None,
) -> str:
    """Build and sign one license key. Returns the compact key string."""
    private_key = private_key or _load_private_key()
    expires_at = None
    if days is not None:
        expires_at = _dt.datetime.now(_dt.timezone.utc) + _dt.timedelta(
            days=days
        )
    payload_bytes = build_signing_payload(
        customer=customer, seats=seats, expires_at=expires_at, features=features
    )
    signature = private_key.sign(payload_bytes)
    return encode_license_key(payload_bytes, signature)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--customer", required=True, help="Customer name/id.")
    parser.add_argument("--seats", type=int, default=1)
    parser.add_argument(
        "--days",
        type=int,
        default=None,
        help="Validity in days from now. Omit for a perpetual license.",
    )
    parser.add_argument(
        "--feature",
        action="append",
        dest="features",
        default=[],
        help="Repeatable feature flag to embed in the license.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    key = issue(
        customer=args.customer,
        seats=args.seats,
        days=args.days,
        features=tuple(args.features),
    )
    print(key)
    return 0


if __name__ == "__main__":
    sys.exit(main())
