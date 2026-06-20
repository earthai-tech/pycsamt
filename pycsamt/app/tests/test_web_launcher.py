# -*- coding: utf-8 -*-
"""Tests for the ``pycsamt-web`` launcher."""

from __future__ import annotations

def test_browser_url_uses_loopback_for_public_bind():
    from pycsamt.app.web.app import _browser_url

    assert _browser_url("0.0.0.0", 8050) == "http://127.0.0.1:8050"
    assert _browser_url("127.0.0.1", 8051) == "http://127.0.0.1:8051"


def test_find_free_port_returns_requested_port_when_available(monkeypatch):
    from pycsamt.app.web.app import _find_free_port

    class FakeSocket:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def setsockopt(self, *args):
            return None

        def bind(self, addr):
            self.addr = addr

        def getsockname(self):
            return ("127.0.0.1", 49000)

    monkeypatch.setattr("socket.socket", lambda *args, **kwargs: FakeSocket())

    port = _find_free_port("127.0.0.1", 0)
    assert isinstance(port, int)
    assert port > 0

    port = _find_free_port("127.0.0.1", 8123)
    assert port == 8123


def test_find_free_port_moves_when_requested_port_is_busy(monkeypatch):
    from pycsamt.app.web import app as web_app

    class BusySocket:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def setsockopt(self, *args):
            return None

        def bind(self, addr):
            raise OSError("busy")

    monkeypatch.setattr("socket.socket", lambda *args, **kwargs: BusySocket())
    monkeypatch.setattr(web_app, "_allocate_free_port", lambda host: 49001)

    fallback = web_app._find_free_port("127.0.0.1", 8050)

    assert fallback == 49001


def test_main_parses_advanced_launcher_options(monkeypatch):
    from pycsamt.app.web import app as web_app

    calls = {}

    class FakeDash:
        def run(self, **kwargs):
            calls["run"] = kwargs

    def fake_create_app(debug=False):
        calls["debug"] = debug
        return FakeDash()

    monkeypatch.setattr(web_app, "create_app", fake_create_app)
    monkeypatch.setattr(web_app, "_find_free_port", lambda host, port: port)
    monkeypatch.setattr(
        web_app,
        "_open_browser_later",
        lambda url, delay: calls.update(browser=(url, delay)),
    )

    web_app.main(
        argv=[
            "--host",
            "127.0.0.1",
            "--port",
            "8123",
            "--debug",
            "--browser-delay",
            "0.25",
        ]
    )

    assert calls["debug"] is True
    assert calls["browser"] == ("http://127.0.0.1:8123", 0.25)
    assert calls["run"] == {
        "host": "127.0.0.1",
        "port": 8123,
        "debug": True,
        "use_reloader": False,
    }


def test_main_no_browser_does_not_schedule_browser(monkeypatch):
    from pycsamt.app.web import app as web_app

    calls = {}

    class FakeDash:
        def run(self, **kwargs):
            calls["run"] = kwargs

    monkeypatch.setattr(web_app, "create_app", lambda debug=False: FakeDash())
    monkeypatch.setattr(web_app, "_find_free_port", lambda host, port: port)
    monkeypatch.setattr(
        web_app,
        "_open_browser_later",
        lambda url, delay: calls.update(browser=(url, delay)),
    )

    web_app.main(argv=["--port", "8124", "--no-browser"])

    assert "browser" not in calls
    assert calls["run"]["port"] == 8124


def test_module_main_forwards_command_line_arguments(monkeypatch):
    from pycsamt.app.web import __main__ as web_main

    calls = {}

    def fake_main(argv=None):
        calls["argv"] = argv

    monkeypatch.setattr("pycsamt.app.web.app.main", fake_main)
    web_main.main(["--port", "8125", "--no-browser"])

    assert calls["argv"] == ["--port", "8125", "--no-browser"]
