# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the pyCSAMT Map View Dash app factory (no real server)."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestAppFactory:
    def test_create_app_returns_dash_instance(self):
        import dash

        from pycsamt.app.mapview.app import create_app

        app = create_app()
        assert isinstance(app, dash.Dash)

    def test_app_has_layout(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        assert app.layout is not None

    def test_app_title(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        assert app.title == "pyCSAMT — Map View"

    def test_app_callbacks_registered(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        assert len(app.callback_map) > 0

    def test_debug_flag_accepted(self):
        from pycsamt.app.mapview.app import create_app

        # `debug` only takes effect via Dash.run(); create_app just
        # needs to accept the flag without raising.
        app = create_app(debug=True)
        assert app is not None

    def test_icon_route_registered(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        rules = [str(r) for r in app.server.url_map.iter_rules()]
        assert any("mv-icons" in r for r in rules)


class TestLaunch:
    def test_launch_seeds_view_and_runs_app(self, monkeypatch):
        from pycsamt.app.mapview import app as app_mod

        recorded = {}

        class _FakeApp:
            def run(self, **kwargs):
                recorded["run_kwargs"] = kwargs

        monkeypatch.setattr(
            app_mod, "create_app", lambda debug=False: _FakeApp()
        )

        seeded = {}

        def fake_set_seed(view):
            seeded["view"] = view

        monkeypatch.setattr(
            "pycsamt.app.mapview.cache.set_seed", fake_set_seed
        )

        sentinel = object()
        app_mod.launch(open_browser=False, view=sentinel)

        assert seeded["view"] is sentinel
        assert recorded["run_kwargs"]["host"] == "127.0.0.1"
        assert recorded["run_kwargs"]["port"] == 8770

    def test_launch_without_view_does_not_seed(self, monkeypatch):
        from pycsamt.app.mapview import app as app_mod

        class _FakeApp:
            def run(self, **kwargs):
                pass

        monkeypatch.setattr(
            app_mod, "create_app", lambda debug=False: _FakeApp()
        )

        called = {"n": 0}
        monkeypatch.setattr(
            "pycsamt.app.mapview.cache.set_seed",
            lambda v: called.__setitem__("n", called["n"] + 1),
        )

        app_mod.launch(open_browser=False, view=None)
        assert called["n"] == 0

    def test_launch_opens_browser_thread(self, monkeypatch):
        from pycsamt.app.mapview import app as app_mod

        class _FakeApp:
            def run(self, **kwargs):
                pass

        monkeypatch.setattr(
            app_mod, "create_app", lambda debug=False: _FakeApp()
        )

        started = {"n": 0}

        class _FakeThread:
            def __init__(self, target=None, daemon=None):
                self._target = target

            def start(self):
                started["n"] += 1

        monkeypatch.setattr("threading.Thread", _FakeThread)

        app_mod.launch(open_browser=True)
        assert started["n"] == 1
