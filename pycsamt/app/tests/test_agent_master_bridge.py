# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.agent_master_bridge.

No Qt dependency — this module is pure subprocess/threading/urllib
plumbing, so these tests run in any environment. The module-level
``_PROCESS`` singleton is saved/restored around every test so state does
not leak between tests or into a real desktop-app session.
"""

from __future__ import annotations

import os
import subprocess
from unittest import mock

import pytest

from pycsamt.app.desktop import agent_master_bridge as bridge


@pytest.fixture(autouse=True)
def reset_process_singleton():
    saved = bridge._PROCESS
    bridge._PROCESS = None
    yield
    bridge._PROCESS = saved


# ── agent_master_url ─────────────────────────────────────────────────────────


def test_agent_master_url_defaults():
    assert bridge.agent_master_url() == "http://127.0.0.1:8765"


def test_agent_master_url_custom_host_port():
    assert bridge.agent_master_url("0.0.0.0", 9000) == "http://0.0.0.0:9000"


def test_agent_master_url_coerces_port_to_int():
    assert bridge.agent_master_url("localhost", "1234") == "http://localhost:1234"


# ── is_agent_master_running ──────────────────────────────────────────────────


def test_is_agent_master_running_true_on_2xx():
    resp = mock.MagicMock()
    resp.status = 200
    resp.__enter__.return_value = resp
    with mock.patch.object(bridge, "urlopen", return_value=resp):
        assert bridge.is_agent_master_running() is True


def test_is_agent_master_running_true_on_4xx():
    resp = mock.MagicMock()
    resp.status = 404
    resp.__enter__.return_value = resp
    with mock.patch.object(bridge, "urlopen", return_value=resp):
        assert bridge.is_agent_master_running() is True


def test_is_agent_master_running_false_on_5xx():
    resp = mock.MagicMock()
    resp.status = 500
    resp.__enter__.return_value = resp
    with mock.patch.object(bridge, "urlopen", return_value=resp):
        assert bridge.is_agent_master_running() is False


@pytest.mark.parametrize(
    "exc",
    [OSError("refused"), TimeoutError("timeout"), ValueError("bad url")],
)
def test_is_agent_master_running_false_on_exceptions(exc):
    with mock.patch.object(bridge, "urlopen", side_effect=exc):
        assert bridge.is_agent_master_running() is False


def test_is_agent_master_running_passes_host_port_timeout():
    resp = mock.MagicMock()
    resp.status = 200
    resp.__enter__.return_value = resp
    with mock.patch.object(bridge, "urlopen", return_value=resp) as m:
        bridge.is_agent_master_running("1.2.3.4", 9999, timeout=1.5)
        args, kwargs = m.call_args
        assert args[0] == "http://1.2.3.4:9999"
        assert kwargs["timeout"] == 1.5


# ── launch_agent_master ──────────────────────────────────────────────────────


def test_launch_already_running_opens_browser(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: True)
    opened = []
    monkeypatch.setattr(bridge.webbrowser, "open", lambda url: opened.append(url))
    result = bridge.launch_agent_master(open_browser=True)
    assert result.started is False
    assert result.url == bridge.agent_master_url()
    assert opened == [bridge.agent_master_url()]


def test_launch_already_running_no_browser(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: True)
    opened = []
    monkeypatch.setattr(bridge.webbrowser, "open", lambda url: opened.append(url))
    result = bridge.launch_agent_master(open_browser=False)
    assert result.started is False
    assert opened == []


def test_launch_existing_process_still_alive_reuses_it(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    fake_proc = mock.MagicMock()
    fake_proc.poll.return_value = None  # still running
    bridge._PROCESS = fake_proc

    fake_thread_cls = mock.MagicMock()
    monkeypatch.setattr(bridge.threading, "Thread", fake_thread_cls)

    result = bridge.launch_agent_master(open_browser=True)
    assert result.started is False
    assert result.process is fake_proc
    fake_thread_cls.assert_called_once()
    _, kwargs = fake_thread_cls.call_args
    assert kwargs["target"] is bridge._open_when_ready
    assert kwargs["daemon"] is True
    fake_thread_cls.return_value.start.assert_called_once()


def test_launch_existing_process_alive_no_browser_skips_thread(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    fake_proc = mock.MagicMock()
    fake_proc.poll.return_value = None
    bridge._PROCESS = fake_proc

    fake_thread_cls = mock.MagicMock()
    monkeypatch.setattr(bridge.threading, "Thread", fake_thread_cls)

    result = bridge.launch_agent_master(open_browser=False)
    assert result.started is False
    fake_thread_cls.assert_not_called()


def test_launch_spawns_new_process(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    bridge._PROCESS = None

    fake_proc = mock.MagicMock()
    fake_popen = mock.MagicMock(return_value=fake_proc)
    monkeypatch.setattr(bridge.subprocess, "Popen", fake_popen)

    fake_thread_cls = mock.MagicMock()
    monkeypatch.setattr(bridge.threading, "Thread", fake_thread_cls)

    result = bridge.launch_agent_master(host="127.0.0.1", port=8765, open_browser=True)

    assert result.started is True
    assert result.process is fake_proc
    assert bridge._PROCESS is fake_proc

    cmd, kwargs = fake_popen.call_args[0][0], fake_popen.call_args[1]
    assert cmd[0] == __import__("sys").executable
    assert "pycsamt.app.agent_master" in cmd
    assert "--host" in cmd and "127.0.0.1" in cmd
    assert "--port" in cmd and "8765" in cmd
    assert "--no-browser" in cmd
    assert kwargs["stdout"] == subprocess.DEVNULL
    assert kwargs["env"]["PYTHONDONTWRITEBYTECODE"] == "1"
    fake_thread_cls.assert_called_once()


def test_launch_spawns_new_process_sets_windows_creationflags(monkeypatch):
    process_os_name = os.name
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    bridge._PROCESS = None
    fake_popen = mock.MagicMock(return_value=mock.MagicMock())
    monkeypatch.setattr(bridge.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(bridge.threading, "Thread", mock.MagicMock())
    # ``bridge.os`` normally references the process-wide ``os`` module.
    # Patching ``bridge.os.name`` would therefore also change ``os.name``
    # inside pathlib and pytest, making Linux try to instantiate WindowsPath.
    fake_os = mock.Mock(wraps=bridge.os)
    fake_os.name = "nt"
    fake_os.environ = bridge.os.environ
    monkeypatch.setattr(bridge, "os", fake_os)
    monkeypatch.setattr(
        bridge.subprocess, "CREATE_NO_WINDOW", 0x08000000, raising=False
    )
    monkeypatch.setattr(
        bridge.subprocess,
        "CREATE_NEW_PROCESS_GROUP",
        0x00000200,
        raising=False,
    )

    bridge.launch_agent_master(open_browser=False)

    assert os.name == process_os_name
    kwargs = fake_popen.call_args[1]
    assert "creationflags" in kwargs
    assert "start_new_session" not in kwargs


def test_launch_spawns_new_process_sets_posix_start_new_session(monkeypatch):
    process_os_name = os.name
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    bridge._PROCESS = None
    fake_popen = mock.MagicMock(return_value=mock.MagicMock())
    monkeypatch.setattr(bridge.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(bridge.threading, "Thread", mock.MagicMock())
    fake_os = mock.Mock(wraps=bridge.os)
    fake_os.name = "posix"
    fake_os.environ = bridge.os.environ
    monkeypatch.setattr(bridge, "os", fake_os)

    bridge.launch_agent_master(open_browser=False)

    assert os.name == process_os_name
    kwargs = fake_popen.call_args[1]
    assert kwargs.get("start_new_session") is True
    assert "creationflags" not in kwargs


def test_launch_no_browser_skips_open_when_ready_thread(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    bridge._PROCESS = None
    monkeypatch.setattr(
        bridge.subprocess,
        "Popen",
        mock.MagicMock(return_value=mock.MagicMock()),
    )
    fake_thread_cls = mock.MagicMock()
    monkeypatch.setattr(bridge.threading, "Thread", fake_thread_cls)

    bridge.launch_agent_master(open_browser=False)
    fake_thread_cls.assert_not_called()


# ── _open_when_ready ──────────────────────────────────────────────────────────


def test_open_when_ready_opens_once_server_responds(monkeypatch):
    calls = {"n": 0}

    def _fake_running(h, p):
        calls["n"] += 1
        return calls["n"] >= 3

    monkeypatch.setattr(bridge, "is_agent_master_running", _fake_running)
    monkeypatch.setattr(bridge.time, "sleep", lambda s: None)
    opened = []
    monkeypatch.setattr(bridge.webbrowser, "open", lambda url: opened.append(url))

    bridge._open_when_ready("h", 1, "http://h:1", attempts=10, interval=0)

    assert opened == ["http://h:1"]
    assert calls["n"] == 3


def test_open_when_ready_opens_anyway_after_exhausting_attempts(monkeypatch):
    monkeypatch.setattr(bridge, "is_agent_master_running", lambda h, p: False)
    monkeypatch.setattr(bridge.time, "sleep", lambda s: None)
    opened = []
    monkeypatch.setattr(bridge.webbrowser, "open", lambda url: opened.append(url))

    bridge._open_when_ready("h", 1, "http://h:1", attempts=3, interval=0)

    assert opened == ["http://h:1"]
