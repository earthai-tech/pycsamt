# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Interfaces for the licensing/trial subsystem.

Nothing here performs real verification. ``LicenseManager`` is a
``Protocol`` so Phase 11 of ``PYCSAMT-DESKTOP-V2.6-MODERNIZATION-PLAN.md``
can supply a concrete implementation (offline signed-license file is the
currently recommended approach -- see that plan's Section 4) without every
caller needing to import a concrete class today.
"""

from __future__ import annotations

import datetime as _dt
import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Protocol, runtime_checkable


class LicenseStatus(Enum):
    """Coarse state the UI branches on (splash gate, About panel, Settings)."""

    TRIAL_ACTIVE = "trial_active"
    TRIAL_EXPIRED = "trial_expired"
    LICENSED = "licensed"
    INVALID = "invalid"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class TrialState:
    """Local trial bookkeeping.

    ``started_at`` is set once, on first launch, and never rewritten.
    ``expires_at`` is derived at creation time (``started_at`` + trial
    length) rather than recomputed from a stored length, so a later change
    to the default trial length cannot silently extend an already-running
    trial.
    """

    started_at: _dt.datetime
    expires_at: _dt.datetime

    @property
    def is_expired(self) -> bool:
        return _dt.datetime.now(_dt.timezone.utc) >= self.expires_at

    @property
    def days_remaining(self) -> int:
        # Count part-days as a day: ``timedelta.days`` truncates, so a trial
        # with 2 d 23 h left read "2 days" and a fresh 30-day trial "29".
        delta = self.expires_at - _dt.datetime.now(_dt.timezone.utc)
        return max(0, math.ceil(delta.total_seconds() / 86400.0))


@runtime_checkable
class LicenseManager(Protocol):
    """What the splash screen / About dialog / Settings license page need.

    Phase 11 supplies a concrete implementation. Code written before Phase
    11 lands (splash screen, About dialog, Settings license page) should
    depend on this Protocol and default to ``NullLicenseManager`` rather
    than skip the license check entirely, so the real implementation drops
    in later without call-site changes.
    """

    def status(self) -> LicenseStatus:
        """Current license/trial status."""
        ...

    def trial_state(self) -> Optional[TrialState]:
        """Local trial bookkeeping, or ``None`` once a paid license is active."""
        ...

    def activate(self, license_key: str) -> LicenseStatus:
        """Validate and persist a license key; return the resulting status."""
        ...
