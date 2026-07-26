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
    session = SimpleNamespace(
        config=SimpleNamespace(args=["pycsamt/core/tests"])
    )
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_sessionfinish(session, 0)

    assert calls == []
