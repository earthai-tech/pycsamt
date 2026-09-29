# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
OfflineLicenseManager -- the real ``LicenseManager`` implementation.

Combines a locally persisted, Ed25519-signed license key (``signing.py``)
with ``TrialTracker``'s local 30-day trial bookkeeping. Every check is
offline; nothing here ever makes a network call. This is the concrete
class :class:`~pycsamt.app.desktop.licensing.interfaces.LicenseManager`
callers (splash screen, About dialog, Settings ▸ License tab) should be
handed once Phase 11 lands -- see :func:`get_default_manager`.
"""

from __future__ import annotations

from typing import Optional

from .interfaces import LicenseStatus, TrialState
from .signing import LicenseError, verify_and_parse
from .trial import TrialTracker

_KEY_LICENSE = "license/key"


class OfflineLicenseManager:
    """Offline, signed-license-file ``LicenseManager`` implementation."""

    def __init__(self, settings=None) -> None:
        if settings is None:
            from PySide6.QtCore import QSettings

            from .. import branding

            settings = QSettings(branding.ORG_NAME, branding.APP_NAME)
        self._settings = settings
        self._trial = TrialTracker(settings)

    # ── LicenseManager protocol ─────────────────────────────────────

    def status(self) -> LicenseStatus:
        raw_key = self._settings.value(_KEY_LICENSE, "", type=str) or ""
        if raw_key:
            lic = self._verify(raw_key)
            if lic is None:
                return LicenseStatus.INVALID
            return (
                LicenseStatus.INVALID if lic.is_expired else LicenseStatus.LICENSED
            )
        trial = self._trial.state()
        return (
            LicenseStatus.TRIAL_EXPIRED
            if trial.is_expired
            else LicenseStatus.TRIAL_ACTIVE
        )

    def trial_state(self) -> Optional[TrialState]:
        raw_key = self._settings.value(_KEY_LICENSE, "", type=str) or ""
        if raw_key and self._verify(raw_key) is not None:
            return None  # a real license is active; no trial to report
        return self._trial.state()

    def activate(self, license_key: str) -> LicenseStatus:
        lic = self._verify(license_key)
        if lic is None or lic.is_expired:
            return LicenseStatus.INVALID
        self._settings.setValue(_KEY_LICENSE, (license_key or "").strip())
        sync = getattr(self._settings, "sync", None)
        if callable(sync):
            sync()
        return LicenseStatus.LICENSED

    # ── internals ────────────────────────────────────────────────────

    @staticmethod
    def _verify(license_key: str):
        try:
            return verify_and_parse(license_key)
        except LicenseError:
            return None


_default_manager: Optional[OfflineLicenseManager] = None


def get_default_manager() -> OfflineLicenseManager:
    """The one ``OfflineLicenseManager`` the running app should share.

    All callers (About dialog, Preferences ▸ License tab, the future
    startup gate) should go through this rather than constructing their
    own, so activating a key in one place is immediately reflected
    everywhere else without re-reading ``QSettings`` from scratch.
    """
    global _default_manager
    if _default_manager is None:
        _default_manager = OfflineLicenseManager()
    return _default_manager
