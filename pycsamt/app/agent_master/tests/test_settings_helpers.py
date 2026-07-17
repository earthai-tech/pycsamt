# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pure helpers in pycsamt.app.agent_master.callbacks.settings."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")


class TestEnvKey:
    def test_reads_env_var_for_provider(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _env_key

        monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-abc")
        assert _env_key("claude") == "sk-abc"

    def test_missing_env_var_returns_empty(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _env_key

        monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
        assert _env_key("claude") == ""

    def test_offline_provider_has_no_env(self):
        from pycsamt.app.agent_master.callbacks.settings import _env_key

        assert _env_key("offline") == ""


class TestStoredKey:
    def test_prefers_saved_key_over_env(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _stored_key

        monkeypatch.setenv("ANTHROPIC_API_KEY", "from-env")
        cfg = {"key_claude": "from-disk"}
        assert _stored_key("claude", cfg) == "from-disk"

    def test_falls_back_to_env_when_unsaved(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _stored_key

        monkeypatch.setenv("ANTHROPIC_API_KEY", "from-env")
        assert _stored_key("claude", {}) == "from-env"


class TestSourceFor:
    def test_empty_typed_with_env_present(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _source_for

        monkeypatch.setenv("ANTHROPIC_API_KEY", "e")
        assert _source_for("claude", "", {}) == "env"

    def test_empty_typed_without_env(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _source_for

        monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
        assert _source_for("claude", "", {}) == "none"

    def test_typed_matches_saved(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _source_for

        monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
        cfg = {"key_claude": "saved-key"}
        assert _source_for("claude", "saved-key", cfg) == "saved"

    def test_typed_matches_env_and_not_saved(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _source_for

        monkeypatch.setenv("ANTHROPIC_API_KEY", "env-key")
        assert _source_for("claude", "env-key", {}) == "env"

    def test_typed_diverges_is_unsaved(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks.settings import _source_for

        monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
        cfg = {"key_claude": "saved-key"}
        assert _source_for("claude", "new-typed-key", cfg) == "unsaved"


class TestStatusBadge:
    def test_offline_provider_is_zero_cost(self):
        from pycsamt.app.agent_master.callbacks.settings import (
            _status_badge,
        )

        badge = _status_badge("offline", "none")
        assert "Zero cost" in str(badge)

    def test_saved_source_badge(self):
        from pycsamt.app.agent_master.callbacks.settings import (
            _status_badge,
        )

        assert "Key saved" in str(_status_badge("claude", "saved"))

    def test_env_source_badge(self):
        from pycsamt.app.agent_master.callbacks.settings import (
            _status_badge,
        )

        assert "From environment" in str(_status_badge("claude", "env"))

    def test_unsaved_source_badge(self):
        from pycsamt.app.agent_master.callbacks.settings import (
            _status_badge,
        )

        assert "Unsaved" in str(_status_badge("claude", "unsaved"))

    def test_no_key_badge(self):
        from pycsamt.app.agent_master.callbacks.settings import (
            _status_badge,
        )

        assert "No key" in str(_status_badge("claude", "none"))


class TestHint:
    def test_contains_env_var_fallback(self):
        from pycsamt.app.agent_master.callbacks.settings import _hint

        text = str(_hint("claude"))
        assert "ANTHROPIC_API_KEY" in text
        assert "agent_master.json" in text

    def test_unknown_provider_shows_dash(self):
        from pycsamt.app.agent_master.callbacks.settings import _hint

        text = str(_hint("mystery"))
        assert "—" in text


class TestBadge:
    def test_badge_contains_text_and_icon(self):
        from pycsamt.app.agent_master.callbacks.settings import _badge

        b = _badge("Hello", "--tag-ok", "bi-check")
        text = str(b)
        assert "Hello" in text
        assert "bi-check" in text


class TestRegisterSettings:
    def test_register_settings_is_callable(self):
        from pycsamt.app.agent_master.callbacks.settings import (
            register_settings,
        )

        assert callable(register_settings)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_SETTINGS in cb_outputs
        assert IDs.CANVAS_SETTINGS in cb_outputs
        assert IDs.PROVIDER_PANEL in cb_outputs
