# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
NullLicenseManager — always-trial stand-in.

Grants an unlimited trial and rejects every key -- no persistence, no
verification, no expiry. Phase 11 landed the real implementation,
:class:`~pycsamt.app.desktop.licensing.manager.OfflineLicenseManager`
(offline Ed25519-signed license + local 30-day trial tracker); the running
app should go through
:func:`~pycsamt.app.desktop.licensing.manager.get_default_manager` rather
than this class. ``NullLicenseManager`` remains useful for tests and any
call site that deliberately wants an unconditional trial with zero
persistence side effects.
"""

from __future__ import annotations

from typing import Optional

from .interfaces import LicenseStatus, TrialState


class NullLicenseManager:
    """No-op ``LicenseManager``: always reports an active, non-expiring trial."""

    def status(self) -> LicenseStatus:
        return LicenseStatus.TRIAL_ACTIVE

    def trial_state(self) -> Optional[TrialState]:
        return None

    def activate(self, license_key: str) -> LicenseStatus:
        return LicenseStatus.INVALID
