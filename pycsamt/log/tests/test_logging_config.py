# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""The bundled YAML logging config must load, and a failed log
rollover must never crash logging."""

from __future__ import annotations

import logging
import os

import pytest

from pycsamt.log.logger import (
    SafeRotatingFileHandler,
    configure_logging,
)


@pytest.fixture(autouse=True)
def _restore_logging():
    """Re-point logging back at the real data home afterwards, so a
    tmp-path reconfigure here doesn't leave the rest of the suite
    logging into a deleted directory."""
    yield
    for h in logging.getLogger().handlers[:]:
        try:
            h.close()
        except Exception:  # noqa: BLE001
            pass
    configure_logging(force=True)


def test_yaml_config_applies_with_custom_handler(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_DATA", str(tmp_path))
    configure_logging(force=True)
    root = logging.getLogger()
    # dictConfig succeeded (not the basicConfig single-StreamHandler fallback)
    assert [type(h).__name__ for h in root.handlers] == [
        "SafeRotatingFileHandler"
    ] * 3


def test_werkzeug_request_spam_is_kept_off_the_info_log(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_DATA", str(tmp_path))
    configure_logging(force=True)
    wz = logging.getLogger("werkzeug")
    assert wz.level == logging.WARNING
    assert wz.propagate is False
    # its INFO per-request lines never reach a file handler
    assert not wz.isEnabledFor(logging.INFO)


def test_failed_rollover_is_swallowed(tmp_path):
    path = tmp_path / "roll.log"
    h = SafeRotatingFileHandler(path, maxBytes=120, backupCount=1)

    def _boom():
        raise PermissionError(32, "in use by another process")

    h.rotate = lambda *a, **k: _boom()  # force every rollover to fail
    lg = logging.getLogger("pycsamt._rolltest")
    lg.handlers = [h]
    lg.propagate = False
    lg.setLevel(logging.INFO)
    for i in range(40):  # would raise from doRollover without the guard
        lg.info("padding-line %03d", i)
    h.close()
    assert path.exists()  # still writing to a usable stream
    assert os.path.getsize(path) > 120  # grew past maxBytes, rollover skipped
