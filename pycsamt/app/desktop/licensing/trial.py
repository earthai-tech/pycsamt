# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TrialTracker -- local, QSettings-backed 30-day trial bookkeeping.

Product decision (PYCSAMT-DESKTOP-V2.6-MODERNIZATION-PLAN.md, Phase 11 /
Section 4): 30-day trial, timed from first launch, bound to a lightweight
machine-fingerprint hash so the trial survives an app reinstall (QSettings
persists in the OS user-config location, not the install directory) and a
casual attempt to hand-edit the stored start date. This is explicitly *not*
hardened DRM -- a determined user who reads this open-source file can patch
the check itself. That is an accepted tradeoff for a first paid release, not
an oversight.

Tamper handling is fail-closed by design: if the stored record doesn't pass
its integrity check (edited timestamp, mismatched fingerprint, corrupted
value), the trial is treated as already expired rather than silently
granting a fresh one. A genuinely first-ever launch (no stored record at
all) is the only path that starts a new trial.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import hmac
from typing import Optional

from .fingerprint import machine_fingerprint
from .interfaces import TrialState

TRIAL_DAYS = 30

# Casual-tamper resistance only -- visible to anyone reading this
# open-source file, not a real secret. See module docstring.
_GUARD_SECRET = b"pycsamt-desktop-trial-guard-v1"

_KEY_STARTED = "trial/started_at"
_KEY_FINGERPRINT = "trial/fingerprint"
_KEY_GUARD = "trial/guard"


def _guard_digest(started_iso: str, fingerprint: str) -> str:
    msg = f"{started_iso}|{fingerprint}".encode("utf-8")
    return hmac.new(_GUARD_SECRET, msg, hashlib.sha256).hexdigest()


class TrialTracker:
    """Reads/writes trial bookkeeping through a ``QSettings``-like object.

    Parameters
    ----------
    settings : object, optional
        Anything exposing Qt's ``QSettings`` ``value``/``setValue``/``sync``
        methods. Defaults to a real ``QSettings(branding.ORG_NAME,
        branding.APP_NAME)`` -- imported lazily so this module stays
        importable without Qt (tests pass a plain fake).
    """

    def __init__(self, settings=None) -> None:
        if settings is None:
            from PySide6.QtCore import QSettings

            from .. import branding

            settings = QSettings(branding.ORG_NAME, branding.APP_NAME)
        self._settings = settings

    def state(self) -> TrialState:
        started = self._load_or_start()
        expires = started + _dt.timedelta(days=TRIAL_DAYS)
        return TrialState(started_at=started, expires_at=expires)

    # ── internals ────────────────────────────────────────────────────

    def _load_or_start(self) -> _dt.datetime:
        started_iso = self._settings.value(_KEY_STARTED, "", type=str) or ""
        if started_iso:
            return self._validate_or_expire(started_iso)
        return self._start_new_trial()

    def _validate_or_expire(self, started_iso: str) -> _dt.datetime:
        stored_fp = self._settings.value(_KEY_FINGERPRINT, "", type=str) or ""
        stored_guard = self._settings.value(_KEY_GUARD, "", type=str) or ""
        expected_guard = _guard_digest(started_iso, stored_fp)
        current_fp = machine_fingerprint()

        intact = hmac.compare_digest(stored_guard, expected_guard)
        same_machine = stored_fp == current_fp
        if intact and same_machine:
            try:
                started = _dt.datetime.fromisoformat(started_iso)
            except ValueError:
                return self._already_expired()
            if started.tzinfo is None:
                started = started.replace(tzinfo=_dt.timezone.utc)
            return started

        # Tampered timestamp, corrupted record, or a settings file copied
        # in from a different machine -- fail closed.
        return self._already_expired()

    def _already_expired(self) -> _dt.datetime:
        return _dt.datetime.now(_dt.timezone.utc) - _dt.timedelta(
            days=TRIAL_DAYS + 1
        )

    def _start_new_trial(self) -> _dt.datetime:
        now = _dt.datetime.now(_dt.timezone.utc)
        fp = machine_fingerprint()
        now_iso = now.isoformat()
        self._settings.setValue(_KEY_STARTED, now_iso)
        self._settings.setValue(_KEY_FINGERPRINT, fp)
        self._settings.setValue(_KEY_GUARD, _guard_digest(now_iso, fp))
        sync = getattr(self._settings, "sync", None)
        if callable(sync):
            sync()
        return now
