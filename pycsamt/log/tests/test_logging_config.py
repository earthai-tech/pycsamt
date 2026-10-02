# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""The bundled YAML logging config must load, and a failed log
rollover must never crash logging."""

from __future__ import annotations

import logging
import os
import sys

import pytest
import yaml

import pycsamt.log._config as log_config
import pycsamt.log.logger as logger_module
from pycsamt.log.logger import (
    SafeRotatingFileHandler,
    configure_logging,
    enable_console_logging,
    get_data_home,
    get_logger,
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
        logging.getLogger().removeHandler(h)
    logger_module._CONFIGURED = False


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


def test_init_logging_path_selection(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(
        log_config,
        "configure_logging",
        lambda **kwargs: calls.append(kwargs),
    )

    explicit = tmp_path / "explicit.yml"
    log_config.init_logging(str(explicit), use_default=True)
    assert calls[-1] == {
        "config_path": str(explicit),
        "use_default": True,
    }

    env_path = tmp_path / "environment.yml"
    monkeypatch.setenv("PYCSAMT_LOG_CONFIG", str(env_path))
    log_config.init_logging("  ")
    assert calls[-1]["config_path"] == str(env_path)

    monkeypatch.setenv("PYCSAMT_LOG_CONFIG", "  ")
    log_config.init_logging()
    assert calls[-1]["config_path"].endswith("p.configlog.yml")


def test_get_data_home_precedence_and_creation(tmp_path, monkeypatch):
    explicit = tmp_path / "explicit"
    assert get_data_home(str(explicit)) == str(explicit)
    assert explicit.is_dir()

    environment = tmp_path / "environment"
    monkeypatch.setenv("PYCSAMT_DATA", str(environment))
    assert get_data_home() == str(environment)
    assert environment.is_dir()


def test_get_data_home_warns_when_creation_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(
        logger_module.os,
        "makedirs",
        lambda *args, **kwargs: (_ for _ in ()).throw(OSError("denied")),
    )
    with pytest.warns(UserWarning, match="Could not create"):
        assert get_data_home(str(tmp_path / "blocked"))


def test_configure_logging_is_idempotent(monkeypatch):
    logger_module._CONFIGURED = True
    called = False

    def unexpected(*args, **kwargs):
        nonlocal called
        called = True

    monkeypatch.setattr(logger_module.logging, "basicConfig", unexpected)
    configure_logging(use_default=True)
    assert called is False


def test_configure_logging_docs_and_default_modes(monkeypatch):
    calls = []
    monkeypatch.setattr(
        logger_module.logging,
        "basicConfig",
        lambda **kwargs: calls.append(kwargs),
    )
    monkeypatch.setenv("PYCSAMT_DOCS_BUILD", "1")
    configure_logging(force=True)
    assert calls[-1]["level"] == logging.WARNING

    monkeypatch.delenv("PYCSAMT_DOCS_BUILD")
    configure_logging(use_default=True, force=True)
    assert calls[-1]["level"] == logging.INFO
    assert "datefmt" in calls[-1]


def test_missing_and_invalid_yaml_fall_back(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        logger_module.logging,
        "basicConfig",
        lambda **kwargs: calls.append(kwargs),
    )
    configure_logging(config_path=str(tmp_path / "missing.yml"), force=True)
    assert calls[-1]["level"] == logging.INFO

    invalid = tmp_path / "invalid.yml"
    invalid.write_text("handlers: [")
    configure_logging(config_path=str(invalid), force=True)
    assert calls[-1]["level"] == logging.INFO


def test_yaml_rewrites_relative_but_not_absolute_log_paths(
    tmp_path, monkeypatch
):
    absolute = tmp_path / "absolute.log"
    config = {
        "version": 1,
        "handlers": {
            "relative": {
                "class": "logging.FileHandler",
                "filename": "relative.log",
            },
            "absolute": {
                "class": "logging.FileHandler",
                "filename": str(absolute),
            },
            "null": {"class": "logging.NullHandler"},
        },
        "root": {
            "level": "INFO",
            "handlers": ["relative", "absolute", "null"],
        },
    }
    path = tmp_path / "logging.yml"
    path.write_text(yaml.safe_dump(config))
    data_home = tmp_path / "data"
    monkeypatch.setenv("PYCSAMT_DATA", str(data_home))

    configure_logging(config_path=str(path), force=True)

    filenames = {
        h.baseFilename
        for h in logging.getLogger().handlers
        if isinstance(h, logging.FileHandler)
    }
    assert str(data_home / "logs" / "relative.log") in filenames
    assert str(absolute) in filenames


def test_enable_console_logging_adds_then_updates_handler(monkeypatch):
    logger = logging.getLogger("pycsamt")
    previous = logger.handlers[:]
    non_console = logging.NullHandler()
    logger.handlers = [non_console]
    try:
        enable_console_logging(logging.DEBUG)
        assert len(logger.handlers) == 2
        handler = logger.handlers[1]
        assert handler.stream is sys.stdout
        assert handler.level == logging.DEBUG

        enable_console_logging(logging.ERROR)
        assert logger.handlers == [non_console, handler]
        assert handler.level == logging.ERROR
    finally:
        logger.handlers = previous


def test_get_logger_named_and_default():
    assert get_logger("custom.name").name == "custom.name"
    assert get_logger().name == logger_module.__name__


def test_rollover_recovers_closed_stream(tmp_path, monkeypatch):
    handler = SafeRotatingFileHandler(tmp_path / "recover.log")
    handler.close()
    monkeypatch.setattr(
        logging.handlers.RotatingFileHandler,
        "doRollover",
        lambda self: (_ for _ in ()).throw(OSError("locked")),
    )
    handler.doRollover()
    assert handler.stream is not None
    handler.close()


def test_rollover_failure_keeps_existing_stream(tmp_path, monkeypatch):
    handler = SafeRotatingFileHandler(tmp_path / "existing.log")
    stream = handler.stream
    monkeypatch.setattr(
        logging.handlers.RotatingFileHandler,
        "doRollover",
        lambda self: (_ for _ in ()).throw(OSError("locked")),
    )
    handler.doRollover()
    assert handler.stream is stream
    handler.close()


def test_rollover_swallows_reopen_failure(tmp_path, monkeypatch):
    handler = SafeRotatingFileHandler(tmp_path / "recover.log")
    handler.close()
    monkeypatch.setattr(
        logging.handlers.RotatingFileHandler,
        "doRollover",
        lambda self: (_ for _ in ()).throw(OSError("locked")),
    )
    monkeypatch.setattr(
        handler,
        "_open",
        lambda: (_ for _ in ()).throw(OSError("still locked")),
    )
    handler.doRollover()
    assert handler.stream is None


def test_module_import_falls_back_to_basic_config_on_unexpected_error(
    monkeypatch,
):
    """The auto-configure-on-import guard at the bottom of logger.py
    must never let a broken configure_logging() crash the import.

    dictConfig is patched (not configure_logging itself) because
    importlib.reload() re-executes the module top-to-bottom: patching
    configure_logging as a module attribute would just be overwritten
    by the module's own `def configure_logging` before the bottom
    auto-run block ever calls it.
    """
    import importlib

    logger_module._CONFIGURED = False
    monkeypatch.setattr(
        logging.config,
        "dictConfig",
        lambda cfg: (_ for _ in ()).throw(RuntimeError("boom")),
    )
    try:
        importlib.reload(logger_module)
    finally:
        monkeypatch.undo()
        logger_module._CONFIGURED = False
        importlib.reload(logger_module)
