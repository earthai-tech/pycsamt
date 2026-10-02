# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Lightweight local machine fingerprint.

Used only to bind a locally-recorded trial start date to "this machine" so
that copying a settings file between machines doesn't silently transplant
an already-running (or already-expired) trial -- see ``trial.py``. This is
not a hardware-attestation or anti-fraud mechanism; it is a deliberately
simple hash of a stable-per-machine identifier, acknowledged as crackable
by a determined user (see PYCSAMT-DESKTOP-V2.6-MODERNIZATION-PLAN.md,
Phase 11 / Section 4).

**Why not just ``uuid.getnode()`` (MAC address)?** The original
implementation did exactly that, and it broke a real, legitimate trial in
testing: the fingerprint changed on the same physical machine between two
launches with nothing "wrong" happening in between -- just using WSL2,
which creates a virtual network adapter, which changed which MAC address
Windows' adapter enumeration handed back first. The same can happen from
connecting/disconnecting a VPN, Docker Desktop starting, USB tethering, or
even a laptop switching between Wi-Fi and Ethernet. Since ``trial.py`` is
fail-closed by design, a fingerprint mismatch is indistinguishable from
tampering and permanently kills the trial -- a real usability bug, not the
"crackable by a determined user" tradeoff Phase 11 knowingly accepted.
This module now prefers a genuine OS-level machine identifier instead,
which doesn't move when network hardware comes and goes.
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
import uuid


def _stable_machine_id() -> str | None:
    """A real OS-level machine identifier, when one can be read.

    Unlike a MAC address, none of these change when a network adapter is
    added or removed:

    - Windows: the registry's ``MachineGuid`` (generated once at OS
      install, readable by any authenticated user without elevation).
    - Linux: ``/etc/machine-id`` (systemd standard; falls back to the
      older ``/var/lib/dbus/machine-id`` location), world-readable.
    - macOS: ``IOPlatformUUID``, read via ``ioreg`` (no elevation needed).

    Returns ``None`` (never raises) if the platform is unrecognized or the
    identifier can't be read -- callers fall back to the old MAC-based
    scheme in that case, which is strictly better than crashing.
    """
    system = platform.system()
    try:
        if system == "Windows":
            import winreg

            with winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE,
                r"SOFTWARE\Microsoft\Cryptography",
                0,
                winreg.KEY_READ | winreg.KEY_WOW64_64KEY,
            ) as key:
                value, _ = winreg.QueryValueEx(key, "MachineGuid")
                return str(value).strip() or None

        if system == "Linux":
            for path in ("/etc/machine-id", "/var/lib/dbus/machine-id"):
                try:
                    with open(path, encoding="utf-8") as fh:
                        value = fh.read().strip()
                except OSError:
                    continue
                if value:
                    return value
            return None

        if system == "Darwin":
            out = subprocess.run(
                ["ioreg", "-rd1", "-c", "IOPlatformExpertDevice"],
                capture_output=True,
                text=True,
                timeout=5,
                check=False,
            ).stdout
            for line in out.splitlines():
                if "IOPlatformUUID" in line:
                    parts = line.split('"')
                    if len(parts) >= 4:
                        return parts[-2].strip() or None
            return None
    except Exception:
        return None
    return None


def machine_fingerprint() -> str:
    """Return a short, stable-per-machine hash.

    Prefers a real OS-level machine identifier (see
    :func:`_stable_machine_id`) over a MAC address -- that identifier
    doesn't move when network hardware comes and goes, unlike
    ``uuid.getnode()`` (see the module docstring for why that matters in
    practice, not just in theory). Falls back to hostname + OS + arch +
    ``uuid.getnode()`` only when no OS-level id could be read at all
    (locked-down permissions, an unrecognized platform, ...). Truncated to
    32 hex chars; this only needs to distinguish machines from each other,
    not serve as a cryptographic key.
    """
    stable_id = _stable_machine_id()
    if stable_id:
        raw = "|".join([platform.system(), platform.machine(), stable_id])
    else:
        raw = "|".join(
            [
                platform.node(),
                platform.system(),
                platform.machine(),
                str(uuid.getnode()),
            ]
        )
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:32]
