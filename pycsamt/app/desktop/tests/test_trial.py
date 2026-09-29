# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.licensing.trial.TrialTracker.

Uses a plain dict-backed fake in place of QSettings (this module only
calls ``value``/``setValue``/``sync``), so these tests need no Qt/qapp
fixture at all and exercise the exact persistence contract TrialTracker
depends on.
"""

from __future__ import annotations

import datetime as dt

from pycsamt.app.desktop.licensing.fingerprint import machine_fingerprint
from pycsamt.app.desktop.licensing.trial import (
    TRIAL_DAYS,
    TrialTracker,
    _guard_digest,
)


class _FakeSettings:
    def __init__(self):
        self._data: dict[str, str] = {}
        self.synced = False

    def value(self, key, default="", type=None):  # noqa: A002 - match QSettings
        return self._data.get(key, default)

    def setValue(self, key, value):
        self._data[key] = value

    def sync(self):
        self.synced = True


def test_first_launch_starts_a_new_trial_and_persists_it():
    settings = _FakeSettings()
    tracker = TrialTracker(settings)
    state = tracker.state()

    assert not state.is_expired
    assert state.days_remaining == TRIAL_DAYS - 1 or state.days_remaining == TRIAL_DAYS
    assert settings.synced
    assert settings.value("trial/started_at")
    assert settings.value("trial/fingerprint") == machine_fingerprint()


def test_second_call_reuses_the_same_start_date():
    settings = _FakeSettings()
    tracker = TrialTracker(settings)
    first = tracker.state()
    second = TrialTracker(settings).state()
    assert first.started_at == second.started_at


def test_trial_expires_after_trial_days():
    settings = _FakeSettings()
    started = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=TRIAL_DAYS + 1)
    fp = machine_fingerprint()
    started_iso = started.isoformat()
    settings.setValue("trial/started_at", started_iso)
    settings.setValue("trial/fingerprint", fp)
    settings.setValue("trial/guard", _guard_digest(started_iso, fp))

    state = TrialTracker(settings).state()
    assert state.is_expired
    assert state.days_remaining == 0


def test_hand_edited_timestamp_fails_closed_to_expired():
    """Editing only the timestamp (not the guard) must not extend a trial."""
    settings = _FakeSettings()
    # A legitimate, freshly-started trial...
    TrialTracker(settings).state()
    # ...tampered with: push the start date back without recomputing the
    # guard digest (exactly what a casual ini/registry edit would do).
    settings.setValue(
        "trial/started_at",
        (dt.datetime.now(dt.timezone.utc) + dt.timedelta(days=100)).isoformat(),
    )

    state = TrialTracker(settings).state()
    assert state.is_expired


def test_settings_copied_from_a_different_machine_fails_closed():
    settings = _FakeSettings()
    started_iso = dt.datetime.now(dt.timezone.utc).isoformat()
    other_fp = "0" * 32
    settings.setValue("trial/started_at", started_iso)
    settings.setValue("trial/fingerprint", other_fp)
    settings.setValue("trial/guard", _guard_digest(started_iso, other_fp))

    state = TrialTracker(settings).state()
    assert state.is_expired


def test_corrupted_timestamp_fails_closed():
    settings = _FakeSettings()
    fp = machine_fingerprint()
    settings.setValue("trial/started_at", "not-a-real-timestamp")
    settings.setValue("trial/fingerprint", fp)
    settings.setValue("trial/guard", _guard_digest("not-a-real-timestamp", fp))

    state = TrialTracker(settings).state()
    assert state.is_expired
