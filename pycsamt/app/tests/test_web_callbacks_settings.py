# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.settings (Settings off-canvas panel).

Pure dict/JSON manipulation over the SETTINGS_STORE -- no matplotlib, no
Sites data required. Exercised directly via the callback functions found
in ``web_app.callback_map`` (no Dash server needed), following the shared
pattern used across the other ``test_web_callbacks_*`` files.
"""

from __future__ import annotations

import base64
import json

import pytest
from dash import no_update

from pycsamt.app.web.callbacks.settings import (
    _DEFAULTS,
    _KEY_PREFIX,
    _MODEL_DEFAULT,
    _MODEL_OPTIONS,
    _status_icon,
)
from pycsamt.app.web.layout import IDs

# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(
        k for k in web_app.callback_map if all(s in k for s in substrings)
    )
    return _unwrap(web_app.callback_map[key])


def _cb_where(web_app, include, exclude=()):
    """Find a callback by required/forbidden output-key substrings.

    Needed because settings.py registers several callbacks sharing
    overlapping output ids (e.g. multiple ones touch the provider-btn
    pattern-matching Output or SETTINGS_FEEDBACK), so a plain
    ``_cb_multi`` (first match wins) can pick the wrong one.
    """
    for k, v in web_app.callback_map.items():
        if all(s in k for s in include) and not any(s in k for s in exclude):
            return _unwrap(v)
    raise AssertionError(f"no callback found: include={include} exclude={exclude}")


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


# ── 1. _toggle (offcanvas open/close) ────────────────────────────────────────


class TestToggleOffcanvas:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.SETTINGS_CANVAS}.is_open")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, False) is no_update

    def test_toggles_open(self, web_app):
        assert self._fn(web_app)(1, False) is True

    def test_toggles_closed(self, web_app):
        assert self._fn(web_app)(1, True) is False


# ── 2. _populate (fill form fields when the offcanvas opens) ────────────────


class TestPopulate:
    def _fn(self, web_app):
        return _cb_multi(web_app, "settings-interp.value", "settings-ai-key.value")

    def test_closed_no_update(self, web_app):
        out = self._fn(web_app)(False, {"interp": "nearest"})
        assert out == (no_update,) * 7

    def test_open_no_stored_uses_defaults(self, web_app):
        out = self._fn(web_app)(True, None)
        assert out == (
            _DEFAULTS["interp"],
            _DEFAULTS["skindepth"],
            _DEFAULTS["dem"],
            _DEFAULTS["cmap"],
            _DEFAULTS["ai_provider"],
            _DEFAULTS["ai_key"],
            _DEFAULTS["ai_model"],
        )

    def test_open_with_partial_stored_merges_defaults(self, web_app):
        stored = {"interp": "nearest", "ai_key": "sk-ant-abc"}
        out = self._fn(web_app)(True, stored)
        assert out[0] == "nearest"
        assert out[5] == "sk-ant-abc"
        # Untouched fields fall back to defaults.
        assert out[1] == _DEFAULTS["skindepth"]
        assert out[6] == _DEFAULTS["ai_model"]


# ── 3. _persist (form fields -> SETTINGS_STORE) ──────────────────────────────


class TestPersist:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.SETTINGS_STORE}.data")

    def test_all_values_provided(self, web_app):
        out = self._fn(web_app)(
            "nearest", "niblett", "srtm", "plasma", "openai", "sk-abc", "gpt-4o"
        )
        assert out == {
            "interp": "nearest",
            "skindepth": "niblett",
            "dem": "srtm",
            "cmap": "plasma",
            "ai_provider": "openai",
            "ai_key": "sk-abc",
            "ai_model": "gpt-4o",
        }

    def test_falsy_values_fall_back_to_defaults(self, web_app):
        out = self._fn(web_app)(None, "", None, "", None, None, "")
        assert out == {
            "interp": _DEFAULTS["interp"],
            "skindepth": _DEFAULTS["skindepth"],
            "dem": _DEFAULTS["dem"],
            "cmap": _DEFAULTS["cmap"],
            "ai_provider": _DEFAULTS["ai_provider"],
            "ai_key": "",
            "ai_model": _DEFAULTS["ai_model"],
        }


# ── 4. _pick_provider (provider pill click) ──────────────────────────────────


class TestPickProvider:
    def _fn(self, web_app):
        return _cb_where(
            web_app,
            include=["settings-ai-provider-store.data", "settings-provider-btn"],
        )

    def test_no_trigger_all_no_update(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.settings as settings_mod

        monkeypatch.setattr(
            settings_mod, "ctx", type("C", (), {"triggered_id": None})()
        )
        out = self._fn(web_app)([None, None, None, None], "claude")
        assert out == (no_update, no_update, no_update, no_update)

    def test_non_dict_trigger_all_no_update(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.settings as settings_mod

        monkeypatch.setattr(
            settings_mod, "ctx", type("C", (), {"triggered_id": "settings-interp"})()
        )
        out = self._fn(web_app)([1], "claude")
        assert out == (no_update, no_update, no_update, no_update)

    def test_picking_openai_updates_store_classes_options_model(
        self, monkeypatch, web_app
    ):
        import pycsamt.app.web.callbacks.settings as settings_mod

        monkeypatch.setattr(
            settings_mod,
            "ctx",
            type(
                "C",
                (),
                {"triggered_id": {"type": "settings-provider-btn", "index": "openai"}},
            )(),
        )
        provider, classes, opts, model = self._fn(web_app)(
            [None, 1, None, None], "claude"
        )
        assert provider == "openai"
        assert classes == [
            "settings-provider-btn",
            "settings-provider-btn active",
            "settings-provider-btn",
            "settings-provider-btn",
        ]
        assert opts == _MODEL_OPTIONS["openai"]
        assert model == _MODEL_DEFAULT["openai"]

    def test_unknown_provider_falls_back_to_claude_options(
        self, monkeypatch, web_app
    ):
        import pycsamt.app.web.callbacks.settings as settings_mod

        monkeypatch.setattr(
            settings_mod,
            "ctx",
            type(
                "C",
                (),
                {
                    "triggered_id": {
                        "type": "settings-provider-btn",
                        "index": "mystery-llm",
                    }
                },
            )(),
        )
        provider, classes, opts, model = self._fn(web_app)([1], "claude")
        assert provider == "mystery-llm"
        # No entry in _MODEL_OPTIONS/_MODEL_DEFAULT for "mystery-llm" ->
        # falls back to claude's option list, and the model default falls
        # back to that list's first option value.
        assert opts == _MODEL_OPTIONS["claude"]
        assert model == _MODEL_OPTIONS["claude"][0]["value"]
        # None of the four fixed pill buttons matches "mystery-llm", so
        # every className stays inactive.
        assert classes == ["settings-provider-btn"] * 4


# ── 5. _sync_provider_ui (restore provider pills on offcanvas open) ─────────


class TestSyncProviderUi:
    def _fn(self, web_app):
        return _cb_where(
            web_app,
            include=["settings-provider-btn", "settings-ai-model.options"],
            exclude=["settings-ai-provider-store.data"],
        )

    def test_none_provider_defaults_to_claude(self, web_app):
        classes, opts = self._fn(web_app)(None)
        assert classes == [
            "settings-provider-btn active",
            "settings-provider-btn",
            "settings-provider-btn",
            "settings-provider-btn",
        ]
        assert opts == _MODEL_OPTIONS["claude"]

    def test_gemini_provider(self, web_app):
        classes, opts = self._fn(web_app)("gemini")
        assert classes == [
            "settings-provider-btn",
            "settings-provider-btn",
            "settings-provider-btn active",
            "settings-provider-btn",
        ]
        assert opts == _MODEL_OPTIONS["gemini"]

    def test_unknown_provider_falls_back_to_claude_options_all_inactive(
        self, web_app
    ):
        classes, opts = self._fn(web_app)("mystery-llm")
        assert classes == ["settings-provider-btn"] * 4
        assert opts == _MODEL_OPTIONS["claude"]


# ── 6. _test_connection (API-key format validation) ──────────────────────────


class TestTestConnection:
    def _fn(self, web_app):
        return _cb_multi(
            web_app, "settings-ai-status.children", "settings-ai-status.className"
        )

    def test_no_clicks_no_update(self, web_app):
        out = self._fn(web_app)(None, "sk-ant-abc", "claude")
        assert out == (no_update, no_update)

    def test_empty_key_shows_error(self, web_app):
        children, cls = self._fn(web_app)(1, "   ", "claude")
        assert "No API key entered" in children[1]
        assert cls == "settings-ai-status settings-ai-err"

    def test_bad_prefix_shows_warning(self, web_app):
        children, cls = self._fn(web_app)(1, "totally-wrong-format", "claude")
        assert "Key format invalid for Claude" in children[1]
        assert cls == "settings-ai-status settings-ai-warn"

    def test_valid_claude_prefix_accepted(self, web_app):
        children, cls = self._fn(web_app)(1, "sk-ant-abcdefghijklmno", "claude")
        assert "Key accepted" in children[1]
        assert cls == "settings-ai-status settings-ai-ok"

    def test_gemini_has_no_prefix_requirement_any_key_accepted(self, web_app):
        # _KEY_PREFIX["gemini"] == () -> the "prefixes and not any(...)"
        # guard short-circuits on the empty tuple, so *any* non-empty key
        # is accepted for gemini regardless of shape.
        assert _KEY_PREFIX["gemini"] == ()
        children, cls = self._fn(web_app)(1, "anything-goes-here", "gemini")
        assert "Key accepted" in children[1]
        assert cls == "settings-ai-status settings-ai-ok"

    def test_default_provider_is_claude_when_none(self, web_app):
        children, cls = self._fn(web_app)(1, "bad-format", None)
        assert "Key format invalid for Claude" in children[1]

    def test_masked_key_shows_first_8_chars_then_dots(self, web_app):
        key = "sk-ant-" + "x" * 30
        children, _cls = self._fn(web_app)(1, key, "claude")
        text = children[1]
        assert key[:8] in text
        assert "•" in text

    def test_short_key_under_8_chars_masks_to_empty_suffix(self, web_app):
        # len(key) - 8 is negative -> "•" * negative == "" (no dots).
        key = "sk-"  # 3 chars, no prefix issue since it doesn't match claude
        # Use deepseek (also requires "sk-" prefix) so this short key passes
        # the format check and reaches the masking branch.
        children, cls = self._fn(web_app)(1, key, "deepseek")
        assert cls == "settings-ai-status settings-ai-ok"
        assert key in children[1]


# ── 7. _reset (reset all settings to defaults) ────────────────────────────────


class TestReset:
    def _fn(self, web_app):
        return _cb_where(web_app, include=["settings-ai-key.value@"])

    def test_no_clicks_all_no_update(self, web_app):
        out = self._fn(web_app)(None)
        assert out == (no_update,) * 9

    def test_reset_returns_defaults_and_feedback(self, web_app):
        out = self._fn(web_app)(1)
        assert out[0] == _DEFAULTS["interp"]
        assert out[1] == _DEFAULTS["skindepth"]
        assert out[2] == _DEFAULTS["dem"]
        assert out[3] == _DEFAULTS["cmap"]
        assert out[4] == _DEFAULTS["ai_provider"]
        assert out[5] == _DEFAULTS["ai_key"]
        assert out[6] == _DEFAULTS["ai_model"]
        assert out[7] == _DEFAULTS
        assert "All settings reset to defaults." in out[8][1]


# ── 8. _save (export settings profile as downloadable JSON) ──────────────────


class TestSave:
    def _fn(self, web_app):
        return _cb_where(web_app, include=["settings-download.data"])

    def test_no_clicks_no_update(self, web_app):
        out = self._fn(web_app)(None, {"interp": "nearest"})
        assert out == (no_update, no_update)

    def test_save_excludes_ai_key_and_uses_defaults_for_missing(self, web_app):
        stored = {"interp": "nearest", "ai_key": "sk-ant-should-not-export"}
        download, feedback = self._fn(web_app)(1, stored)
        assert download["filename"] == "pycsamt_settings.json"
        cfg = json.loads(download["content"])
        assert "ai_key" not in cfg
        assert cfg["interp"] == "nearest"
        assert cfg["skindepth"] == _DEFAULTS["skindepth"]
        assert "Profile saved (API key excluded)." in feedback[1]

    def test_save_with_no_stored_data_uses_all_defaults(self, web_app):
        download, _feedback = self._fn(web_app)(1, None)
        cfg = json.loads(download["content"])
        assert cfg["interp"] == _DEFAULTS["interp"]


# ── 9. _load (import settings profile from an uploaded JSON file) ───────────


class TestLoad:
    def _fn(self, web_app):
        return _cb_where(
            web_app,
            include=["settings-feedback.children@", "settings-store.data@"],
        )

    def _b64_contents(self, cfg: dict) -> str:
        encoded = base64.b64encode(json.dumps(cfg).encode()).decode()
        return f"data:application/json;base64,{encoded}"

    def test_no_contents_all_no_update(self, web_app):
        out = self._fn(web_app)(None, "profile.json")
        assert out == (no_update,) * 8

    def test_valid_profile_loaded(self, web_app):
        cfg = {"interp": "nearest", "cmap": "plasma", "ai_provider": "openai"}
        contents = self._b64_contents(cfg)
        out = self._fn(web_app)(contents, "myprofile.json")
        assert out[0] == "nearest"
        assert out[3] == "plasma"
        assert out[4] == "openai"
        assert out[6]["interp"] == "nearest"
        assert "Profile loaded from myprofile.json." in out[7][1]

    def test_malformed_base64_shows_parse_error(self, web_app):
        out = self._fn(web_app)("data:application/json;base64,not-valid-b64!!", "bad.json")
        assert out[:7] == (no_update,) * 7
        assert "Could not parse bad.json" in out[7][1]

    def test_valid_base64_invalid_json_shows_parse_error(self, web_app):
        encoded = base64.b64encode(b"not json at all").decode()
        contents = f"data:application/json;base64,{encoded}"
        out = self._fn(web_app)(contents, "notjson.json")
        assert out[:7] == (no_update,) * 7
        assert "Could not parse notjson.json" in out[7][1]


# ── _status_icon ──────────────────────────────────────────────────────────────


class TestStatusIcon:
    @pytest.mark.parametrize(
        "variant,color",
        [("ok", "var(--green)"), ("warn", "var(--yellow)"), ("danger", "var(--red)")],
    )
    def test_known_variants_map_to_color(self, variant, color):
        icon = _status_icon("bi-check-circle-fill", variant)
        assert icon.className == "bi bi-check-circle-fill me-1"
        assert icon.style == {"color": color}

    def test_unknown_variant_falls_back_to_inherit(self):
        icon = _status_icon("bi-x", "some-other-variant")
        assert icon.style == {"color": "inherit"}
