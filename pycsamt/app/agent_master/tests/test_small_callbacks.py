# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the small agent_master callback modules: help, plus, session,
splash — mostly registration/wiring smoke checks plus the few pure helpers
each module exposes."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestHelp:
    def test_register_help_is_callable(self):
        from pycsamt.app.agent_master.callbacks.help import register_help

        assert callable(register_help)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.MODAL_HELP in cb_outputs
        assert IDs.INPUT in cb_outputs


class TestPlus:
    def test_register_plus_is_callable(self):
        from pycsamt.app.agent_master.callbacks.plus import register_plus

        assert callable(register_plus)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.CANVAS_EDI in cb_outputs
        assert IDs.CANVAS_SETTINGS in cb_outputs


class TestSession:
    def test_session_path_default_is_timestamped(self, monkeypatch, tmp_path):
        import pycsamt.app.agent_master.callbacks.session as session_mod

        monkeypatch.setattr(session_mod, "_SESSION_DIR", tmp_path)
        path = session_mod._session_path()
        assert path.parent == tmp_path
        assert path.name.startswith("session_")
        assert path.name.endswith(".json")
        assert tmp_path.is_dir()

    def test_session_path_honors_explicit_name(self, monkeypatch, tmp_path):
        import pycsamt.app.agent_master.callbacks.session as session_mod

        monkeypatch.setattr(session_mod, "_SESSION_DIR", tmp_path)
        path = session_mod._session_path("custom.json")
        assert path == tmp_path / "custom.json"

    def test_register_session_is_callable(self):
        from pycsamt.app.agent_master.callbacks.session import (
            register_session,
        )

        assert callable(register_session)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.SAVE_STATUS in cb_outputs


class TestSplash:
    def test_splash_cards_reference_valid_ids(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.callbacks.splash import _SPLASH_CARDS

        assert IDs.SPLASH_CARD_LOAD in _SPLASH_CARDS
        assert IDs.SPLASH_CARD_CHAT in _SPLASH_CARDS
        assert IDs.SPLASH_CARD_AI in _SPLASH_CARDS
        assert IDs.SPLASH_CARD_REPORT in _SPLASH_CARDS
        assert len(_SPLASH_CARDS) == 4

    def test_register_splash_is_callable(self):
        from pycsamt.app.agent_master.callbacks.splash import (
            register_splash,
        )

        assert callable(register_splash)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.CANVAS_EDI in cb_outputs
        assert IDs.SPLASH_OVERLAY in cb_outputs
