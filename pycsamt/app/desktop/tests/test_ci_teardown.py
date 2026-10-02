from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

import conftest as root_conftest


@pytest.fixture(autouse=True)
def _not_an_xdist_worker(monkeypatch):
    # These tests call the hooks directly but may themselves run inside an
    # xdist worker, where the root hooks deliberately defer the exit.
    monkeypatch.setattr(root_conftest, "_IS_XDIST_WORKER", False)


def test_root_sessionfinish_defers_exit_in_xdist_worker(monkeypatch):
    # A worker must not exit before xdist sends ``workerfinished``; it
    # leaves from pytest_unconfigure instead.
    calls = []
    session = SimpleNamespace(
        config=SimpleNamespace(args=["pycsamt/app/desktop/tests"])
    )
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(root_conftest, "_IS_XDIST_WORKER", True)
    monkeypatch.setattr(root_conftest, "_worker_exit_status", None)
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_sessionfinish(session, 5)
    assert calls == []
    root_conftest.pytest_unconfigure(session.config)
    assert calls == [5]


def test_root_sessionfinish_terminates_interface_run(monkeypatch):
    calls = []
    session = SimpleNamespace(
        config=SimpleNamespace(
            args=["pycsamt/app/desktop/tests", "pycsamt/map/tests"]
        )
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
    config = SimpleNamespace(args=["pycsamt/app/desktop/tests"])
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_terminal_summary(None, 3, config)

    assert calls == [3]


def test_root_sessionfinish_ignores_mapview_and_web_only_run(monkeypatch):
    """mapview/web are pure Dash suites; neither imports PySide6 anywhere,
    so a run touching only those two directories must never trigger the
    Shiboken-teardown workaround (unlike desktop/agent_master/converter)."""
    calls = []
    session = SimpleNamespace(
        config=SimpleNamespace(
            args=["pycsamt/app/mapview/tests", "pycsamt/app/web/tests"]
        )
    )
    monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
    monkeypatch.setattr(
        root_conftest, "_terminate_process", lambda code: calls.append(code)
    )

    root_conftest.pytest_sessionfinish(session, 0)

    assert calls == []


def test_root_sessionfinish_terminates_agent_master_or_converter_alone(
    monkeypatch,
):
    """Each Qt test directory triggers the workaround on its own, not only
    when bundled alongside pycsamt/app/desktop/tests in one invocation."""
    for qt_dir in (
        "pycsamt/app/agent_master/tests",
        "pycsamt/app/converter/tests",
    ):
        calls = []
        session = SimpleNamespace(config=SimpleNamespace(args=[qt_dir]))
        monkeypatch.setitem(sys.modules, "PySide6", SimpleNamespace())
        monkeypatch.setattr(
            root_conftest,
            "_terminate_process",
            lambda code: calls.append(code),
        )

        root_conftest.pytest_sessionfinish(session, 0)

        assert calls == [0], f"expected termination for {qt_dir!r}"
