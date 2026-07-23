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


class TestLoadSaveCfg:
    def test_load_cfg_missing_file_returns_empty(self, tmp_path, monkeypatch):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        monkeypatch.setattr(st_mod, "_CFG_FILE", tmp_path / "nope.json")
        assert st_mod._load_cfg() == {}

    def test_load_cfg_corrupt_file_returns_empty(self, tmp_path, monkeypatch):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        f = tmp_path / "bad.json"
        f.write_text("{not json")
        monkeypatch.setattr(st_mod, "_CFG_FILE", f)
        assert st_mod._load_cfg() == {}

    def test_save_then_load_roundtrip(self, tmp_path, monkeypatch):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_dir = tmp_path / "cfgdir"
        monkeypatch.setattr(st_mod, "_CFG_DIR", cfg_dir)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_dir / "am.json")
        st_mod._save_cfg({"provider": "claude"})
        assert st_mod._load_cfg() == {"provider": "claude"}


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == input_id
        and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


class TestInitSettings:
    def _fn(self, agent_app):
        return _find(agent_app, "am-root", "am-store-settings.data")

    def test_returns_no_update_when_no_cfg(
        self, agent_app, monkeypatch, tmp_path
    ):
        from dash import no_update

        from pycsamt.app.agent_master.callbacks import settings as st_mod

        monkeypatch.setattr(st_mod, "_CFG_FILE", tmp_path / "nope.json")
        fn = self._fn(agent_app)
        assert fn("am-root") is no_update

    def test_returns_cfg_when_present(self, agent_app, monkeypatch, tmp_path):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_file = tmp_path / "am.json"
        monkeypatch.setattr(st_mod, "_CFG_DIR", tmp_path)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_file)
        st_mod._save_cfg({"provider": "claude"})
        fn = self._fn(agent_app)
        assert fn("am-root") == {"provider": "claude"}


class TestToggleSettings:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_SETTINGS, IDs.CANVAS_SETTINGS)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, False)

    def test_flips_open_state(self, agent_app):
        fn = self._fn(agent_app)
        assert fn(1, False) is True
        assert fn(1, True) is False


class TestSyncProviderPanel:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.ACTIVE_PROVIDER, IDs.PROVIDER_PANEL)

    def test_offline_provider_shows_note(self, agent_app):
        fn = self._fn(agent_app)
        result = fn("offline", {})
        panel_style, note_style = result[0], result[1]
        assert panel_style == {"display": "none"}
        assert note_style == {"display": "block"}
        assert "Zero cost" in str(result[6])

    def test_llm_provider_uses_draft_key_over_disk(
        self, agent_app, monkeypatch, tmp_path
    ):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_file = tmp_path / "am.json"
        monkeypatch.setattr(st_mod, "_CFG_DIR", tmp_path)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_file)
        st_mod._save_cfg({"key_claude": "disk-key"})

        fn = self._fn(agent_app)
        result = fn("claude", {"key_claude": "draft-key"})
        panel_style, note_style, key, options, model = result[:5]
        assert panel_style == {"display": "block"}
        assert key == "draft-key"
        assert isinstance(options, list)

    def test_llm_provider_falls_back_to_stored_key(
        self, agent_app, monkeypatch, tmp_path
    ):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_file = tmp_path / "am.json"
        monkeypatch.setattr(st_mod, "_CFG_DIR", tmp_path)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_file)
        st_mod._save_cfg({"key_claude": "disk-key"})

        fn = self._fn(agent_app)
        result = fn("claude", {})
        assert result[2] == "disk-key"


class TestStashDraft:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.KEY_INPUT, IDs.STORE_KEY_DRAFTS)

    def test_prevent_update_for_offline_provider(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn("some-key", "some-model", "offline", {})

    def test_stashes_key_and_model_draft(
        self, agent_app, monkeypatch, tmp_path
    ):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        monkeypatch.setattr(st_mod, "_CFG_FILE", tmp_path / "nope.json")
        fn = self._fn(agent_app)
        drafts, badge = fn("typed-key", "opus", "claude", {})
        assert drafts == {
            "key_claude": "typed-key",
            "model_claude": "opus",
        }
        assert "Unsaved" in str(badge)


class TestRevealKey:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_KEY_REVEAL, IDs.KEY_INPUT)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_odd_clicks_reveal_text(self, agent_app):
        fn = self._fn(agent_app)
        input_type, _icon = fn(1)
        assert input_type == "text"

    def test_even_clicks_hide_password(self, agent_app):
        fn = self._fn(agent_app)
        input_type, _icon = fn(2)
        assert input_type == "password"


class TestLoadPrefs:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.CANVAS_SETTINGS, IDs.ACTIVE_PROVIDER)

    def test_prevent_update_when_closed(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(False)

    def test_returns_defaults_without_cfg(
        self, agent_app, monkeypatch, tmp_path
    ):
        from pycsamt.app.agent_master.callbacks import settings as st_mod
        from pycsamt.app.agent_master._providers import OFFLINE

        monkeypatch.setattr(st_mod, "_CFG_FILE", tmp_path / "nope.json")
        fn = self._fn(agent_app)
        provider, export_fmt, out_dir, line_reg = fn(True)
        assert provider == OFFLINE
        assert export_fmt == "png"
        assert out_dir == ""
        assert line_reg == ""

    def test_returns_saved_prefs(self, agent_app, monkeypatch, tmp_path):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_file = tmp_path / "am.json"
        monkeypatch.setattr(st_mod, "_CFG_DIR", tmp_path)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_file)
        st_mod._save_cfg(
            {
                "provider": "claude",
                "export_fmt": "svg",
                "output_dir": "/out",
                "line_registry": "reg",
            }
        )
        fn = self._fn(agent_app)
        assert fn(True) == ("claude", "svg", "/out", "reg")


class TestSaveSettings:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_SAVE_KEYS, IDs.STORE_SETTINGS)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, "offline", "", None, "png", "", "", {})

    def test_saves_llm_provider_key_and_model(
        self, agent_app, monkeypatch, tmp_path
    ):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_file = tmp_path / "am.json"
        monkeypatch.setattr(st_mod, "_CFG_DIR", tmp_path)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_file)

        fn = self._fn(agent_app)
        cfg, status, badge = fn(
            1, "claude", "sk-typed", "opus", "pdf", "/out", "reg", {}
        )
        assert cfg["key_claude"] == "sk-typed"
        assert cfg["model_claude"] == "opus"
        assert cfg["provider"] == "claude"
        assert cfg["export_fmt"] == "pdf"
        assert "Saved" in str(status)
        # persisted to disk too
        assert st_mod._load_cfg()["key_claude"] == "sk-typed"

    def test_drafts_from_other_providers_are_merged(
        self, agent_app, monkeypatch, tmp_path
    ):
        from pycsamt.app.agent_master.callbacks import settings as st_mod

        cfg_file = tmp_path / "am.json"
        monkeypatch.setattr(st_mod, "_CFG_DIR", tmp_path)
        monkeypatch.setattr(st_mod, "_CFG_FILE", cfg_file)

        fn = self._fn(agent_app)
        drafts = {"key_openai": "other-provider-key"}
        cfg, _status, _badge = fn(
            1, "offline", "", None, "png", "", "", drafts
        )
        assert cfg["key_openai"] == "other-provider-key"
        assert cfg["provider"] == "offline"


class TestToggleThemeSettings:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_THEME, IDs.STORE_THEME)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, "dark")

    def test_light_to_dark(self, agent_app):
        fn = self._fn(agent_app)
        theme, icon_cls = fn(1, "light")
        assert theme == "dark"
        assert "sun" in icon_cls

    def test_dark_to_light(self, agent_app):
        fn = self._fn(agent_app)
        theme, icon_cls = fn(1, "dark")
        assert theme == "light"
        assert "moon" in icon_cls
