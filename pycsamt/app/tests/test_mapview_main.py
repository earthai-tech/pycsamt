# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for python -m pycsamt.app.mapview (CLI entry point)."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestParse:
    def test_defaults(self, monkeypatch):
        import sys

        from pycsamt.app.mapview.__main__ import _parse

        monkeypatch.setattr(sys, "argv", ["pycsamt-mapview"])
        args = _parse()
        assert args.host == "127.0.0.1"
        assert args.port == 8770
        assert args.data is None
        assert args.debug is False
        assert args.no_browser is False

    def test_custom_flags(self, monkeypatch):
        import sys

        from pycsamt.app.mapview.__main__ import _parse

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "pycsamt-mapview",
                "--host",
                "0.0.0.0",
                "--port",
                "9000",
                "--data",
                "/some/folder",
                "--debug",
                "--no-browser",
            ],
        )
        args = _parse()
        assert args.host == "0.0.0.0"
        assert args.port == 9000
        assert args.data == "/some/folder"
        assert args.debug is True
        assert args.no_browser is True


class TestMain:
    def test_main_launches_without_data(self, monkeypatch):
        import sys

        import pycsamt.app.mapview.__main__ as main_mod

        monkeypatch.setattr(
            sys, "argv", ["pycsamt-mapview", "--no-browser"]
        )
        recorded = {}

        def fake_launch(**kwargs):
            recorded.update(kwargs)

        monkeypatch.setattr(
            "pycsamt.app.mapview.app.launch", fake_launch
        )
        assert main_mod.main() == 0
        assert recorded["host"] == "127.0.0.1"
        assert recorded["port"] == 8770
        assert recorded["open_browser"] is False
        assert recorded["view"] is None

    def test_main_preloads_data_folder(self, monkeypatch, tmp_path):
        import sys

        import pycsamt.app.mapview.__main__ as main_mod

        monkeypatch.setattr(
            sys,
            "argv",
            ["pycsamt-mapview", "--data", str(tmp_path), "--no-browser"],
        )
        recorded = {}
        monkeypatch.setattr(
            "pycsamt.app.mapview.app.launch",
            lambda **kwargs: recorded.update(kwargs),
        )
        monkeypatch.setattr(
            "pycsamt.map.MapView.from_folder",
            classmethod(lambda cls, path, **kw: "fake-view"),
        )
        assert main_mod.main() == 0
        assert recorded["view"] == "fake-view"
