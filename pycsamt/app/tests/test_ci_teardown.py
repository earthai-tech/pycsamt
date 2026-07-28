from __future__ import annotations

import sys
from types import SimpleNamespace

import conftest as root_conftest


def test_root_sessionfinish_terminates_interface_run(monkeypatch):
    calls = []
    session = SimpleNamespace(
        config=SimpleNamespace(args=["pycsamt/app/tests", "pycsamt/map/tests"])
    )
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_sessionfinish(session, 7)

    assert calls == [7]


def test_root_sessionfinish_ignores_non_interface_run(monkeypatch):
    calls = []
    session = SimpleNamespace(config=SimpleNamespace(args=["pycsamt/core/tests"]))
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_sessionfinish(session, 0)

    assert calls == []


def test_root_terminal_summary_terminates_after_reports(monkeypatch):
    calls = []
    config = SimpleNamespace(args=["pycsamt/app/tests"])
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_terminal_summary(None, 3, config)

    assert calls == [3]
