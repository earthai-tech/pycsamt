# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for python -m pycsamt.app.agent_master (CLI entry point)."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestParse:
    def test_defaults(self, monkeypatch):
        import sys

        from pycsamt.app.agent_master.__main__ import _parse

        monkeypatch.setattr(sys, "argv", ["pycsamt-agent-master"])
        args = _parse()
        assert args.host == "127.0.0.1"
        assert args.port == 8765
        assert args.debug is False
        assert args.no_browser is False

    def test_custom_flags(self, monkeypatch):
        import sys

        from pycsamt.app.agent_master.__main__ import _parse

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "pycsamt-agent-master",
                "--host",
                "0.0.0.0",
                "--port",
                "9001",
                "--debug",
                "--no-browser",
            ],
        )
        args = _parse()
        assert args.host == "0.0.0.0"
        assert args.port == 9001
        assert args.debug is True
        assert args.no_browser is True


class TestMain:
    def test_main_calls_launch_with_parsed_args(self, monkeypatch):
        import sys

        import pycsamt.app.agent_master.__main__ as main_mod

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "pycsamt-agent-master",
                "--port",
                "9002",
                "--no-browser",
            ],
        )
        recorded = {}
        monkeypatch.setattr(
            "pycsamt.app.agent_master.app.launch",
            lambda **kwargs: recorded.update(kwargs),
        )
        assert main_mod.main() == 0
        assert recorded["host"] == "127.0.0.1"
        assert recorded["port"] == 9002
        assert recorded["open_browser"] is False
        assert recorded["debug"] is False
