# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/settings.py — Settings off-canvas panel.

Responsibilities
----------------
- Open/close the Settings offcanvas.
- Persist all settings (computation + AI) in SETTINGS_STORE (localStorage).
- Provider pill button active state + model options update.
- API key show/hide toggle (clientside).
- Test connection (format validation).
- Save / Load / Reset profile JSON.
"""

from __future__ import annotations

import base64
import json

from dash import (
    ALL,
    Input,
    Output,
    State,
    no_update,
)
from dash import (
    callback_context as ctx,
)

from pycsamt.app.web.layout import IDs

# ── Defaults ─────────────────────────────────────────────────────────────────
_DEFAULTS = {
    "interp": "bilinear",
    "skindepth": "cagniard",
    "dem": "open-elevation",
    "cmap": "viridis",
    "ai_provider": "claude",
    "ai_key": "",
    "ai_model": "claude-sonnet-4-6",
}

_MODEL_OPTIONS = {
    "claude": [
        {
            "label": "claude-sonnet-4-6  (default)",
            "value": "claude-sonnet-4-6",
        },
        {
            "label": "claude-opus-4-8",
            "value": "claude-opus-4-8",
        },
        {
            "label": "claude-haiku-4-5",
            "value": "claude-haiku-4-5-20251001",
        },
    ],
    "openai": [
        {
            "label": "gpt-4o  (default)",
            "value": "gpt-4o",
        },
        {"label": "gpt-4o-mini", "value": "gpt-4o-mini"},
        {"label": "gpt-4-turbo", "value": "gpt-4-turbo"},
        {"label": "gpt-3.5-turbo", "value": "gpt-3.5-turbo"},
    ],
    "gemini": [
        {
            "label": "gemini-2.0-flash  (default)",
            "value": "gemini-2.0-flash",
        },
        {"label": "gemini-1.5-pro", "value": "gemini-1.5-pro"},
        {"label": "gemini-1.5-flash", "value": "gemini-1.5-flash"},
    ],
    "deepseek": [
        {
            "label": "deepseek-chat  (default)",
            "value": "deepseek-chat",
        },
        {
            "label": "deepseek-reasoner",
            "value": "deepseek-reasoner",
        },
    ],
}

_MODEL_DEFAULT = {
    "claude": "claude-sonnet-4-6",
    "openai": "gpt-4o",
    "gemini": "gemini-2.0-flash",
    "deepseek": "deepseek-chat",
}

_KEY_PREFIX = {
    "claude": ("sk-ant-",),
    "openai": ("sk-",),
    "gemini": (),
    "deepseek": ("sk-",),
}


def register_settings(app) -> None:

    # ── 1. Toggle offcanvas ───────────────────────────────────────────────────
    @app.callback(
        Output(IDs.SETTINGS_CANVAS, "is_open"),
        Input(IDs.BTN_SETTINGS, "n_clicks"),
        State(IDs.SETTINGS_CANVAS, "is_open"),
        prevent_initial_call=True,
    )
    def _toggle(n, is_open):
        return not is_open if n else no_update

    # ── 2. Populate form on open ──────────────────────────────────────────────
    @app.callback(
        Output("settings-interp", "value"),
        Output("settings-skindepth", "value"),
        Output("settings-dem", "value"),
        Output("settings-cmap", "value"),
        Output(IDs.SETTINGS_AI_PROVIDER, "data"),
        Output(IDs.SETTINGS_AI_KEY, "value"),
        Output(IDs.SETTINGS_AI_MODEL, "value"),
        Input(IDs.SETTINGS_CANVAS, "is_open"),
        State(IDs.SETTINGS_STORE, "data"),
        prevent_initial_call=True,
    )
    def _populate(is_open, stored):
        if not is_open:
            return (no_update,) * 7
        cfg = {**_DEFAULTS, **(stored or {})}
        return (
            cfg["interp"],
            cfg["skindepth"],
            cfg["dem"],
            cfg["cmap"],
            cfg["ai_provider"],
            cfg["ai_key"],
            cfg["ai_model"],
        )

    # ── 3. Persist all settings to store ─────────────────────────────────────
    @app.callback(
        Output(IDs.SETTINGS_STORE, "data"),
        Input("settings-interp", "value"),
        Input("settings-skindepth", "value"),
        Input("settings-dem", "value"),
        Input("settings-cmap", "value"),
        Input(IDs.SETTINGS_AI_PROVIDER, "data"),
        Input(IDs.SETTINGS_AI_KEY, "value"),
        Input(IDs.SETTINGS_AI_MODEL, "value"),
        prevent_initial_call=True,
    )
    def _persist(interp, skindepth, dem, cmap, provider, key, model):
        return {
            "interp": interp or _DEFAULTS["interp"],
            "skindepth": skindepth or _DEFAULTS["skindepth"],
            "dem": dem or _DEFAULTS["dem"],
            "cmap": cmap or _DEFAULTS["cmap"],
            "ai_provider": provider or _DEFAULTS["ai_provider"],
            "ai_key": key or "",
            "ai_model": model or _DEFAULTS["ai_model"],
        }

    # ── 4. Provider pill → store + button classNames + model options ──────────
    @app.callback(
        Output(IDs.SETTINGS_AI_PROVIDER, "data", allow_duplicate=True),
        Output({"type": "settings-provider-btn", "index": ALL}, "className"),
        Output(IDs.SETTINGS_AI_MODEL, "options"),
        Output(IDs.SETTINGS_AI_MODEL, "value", allow_duplicate=True),
        Input({"type": "settings-provider-btn", "index": ALL}, "n_clicks"),
        State(IDs.SETTINGS_AI_PROVIDER, "data"),
        prevent_initial_call=True,
    )
    def _pick_provider(n_clicks, current):
        triggered = ctx.triggered_id
        if not triggered or not isinstance(triggered, dict):
            return no_update, no_update, no_update, no_update
        provider = triggered["index"]
        providers = ["claude", "openai", "gemini", "deepseek"]
        classes = [
            "settings-provider-btn active" if p == provider else "settings-provider-btn"
            for p in providers
        ]
        opts = _MODEL_OPTIONS.get(provider, _MODEL_OPTIONS["claude"])
        model = _MODEL_DEFAULT.get(provider, opts[0]["value"])
        return provider, classes, opts, model

    # ── 5. Restore provider button state on open ──────────────────────────────
    @app.callback(
        Output(
            {"type": "settings-provider-btn", "index": ALL},
            "className",
            allow_duplicate=True,
        ),
        Output(
            IDs.SETTINGS_AI_MODEL,
            "options",
            allow_duplicate=True,
        ),
        Input(IDs.SETTINGS_AI_PROVIDER, "data"),
        prevent_initial_call=True,
    )
    def _sync_provider_ui(provider):
        provider = provider or "claude"
        providers = ["claude", "openai", "gemini", "deepseek"]
        classes = [
            "settings-provider-btn active" if p == provider else "settings-provider-btn"
            for p in providers
        ]
        return (
            classes,
            _MODEL_OPTIONS.get(provider, _MODEL_OPTIONS["claude"]),
        )

    # ── 6. Show / hide API key (clientside) ───────────────────────────────────
    app.clientside_callback(
        """
        function(n, cur_type) {
            const new_type = (cur_type === 'password') ? 'text' : 'password';
            const icon = (new_type === 'text') ? 'bi bi-eye-slash' : 'bi bi-eye';
            return [new_type, icon];
        }
        """,
        Output(IDs.SETTINGS_AI_KEY, "type"),
        Output("settings-ai-eye-icon", "className"),
        Input(IDs.SETTINGS_AI_KEY_SHOW, "n_clicks"),
        State(IDs.SETTINGS_AI_KEY, "type"),
        prevent_initial_call=True,
    )

    # ── 7. Test connection (key format validation) ────────────────────────────
    @app.callback(
        Output(IDs.SETTINGS_AI_STATUS, "children"),
        Output(IDs.SETTINGS_AI_STATUS, "className"),
        Input(IDs.SETTINGS_AI_TEST, "n_clicks"),
        State(IDs.SETTINGS_AI_KEY, "value"),
        State(IDs.SETTINGS_AI_PROVIDER, "data"),
        prevent_initial_call=True,
    )
    def _test_connection(n, key, provider):
        if not n:
            return no_update, no_update
        key = (key or "").strip()
        if not key:
            return (
                [
                    _status_icon("bi-x-circle-fill", "danger"),
                    " No API key entered",
                ],
                "settings-ai-status settings-ai-err",
            )
        prefixes = _KEY_PREFIX.get(provider or "claude", ())
        if prefixes and not any(key.startswith(p) for p in prefixes):
            prov_name = (provider or "claude").capitalize()
            return (
                [
                    _status_icon("bi-exclamation-triangle-fill", "warn"),
                    f" Key format invalid for {prov_name}",
                ],
                "settings-ai-status settings-ai-warn",
            )
        masked = key[:8] + "•" * min(len(key) - 8, 20)
        return (
            [
                _status_icon("bi-check-circle-fill", "ok"),
                f" Key accepted · {masked}",
            ],
            "settings-ai-status settings-ai-ok",
        )

    # ── 8. Reset to defaults ──────────────────────────────────────────────────
    @app.callback(
        Output("settings-interp", "value", allow_duplicate=True),
        Output("settings-skindepth", "value", allow_duplicate=True),
        Output("settings-dem", "value", allow_duplicate=True),
        Output("settings-cmap", "value", allow_duplicate=True),
        Output(IDs.SETTINGS_AI_PROVIDER, "data", allow_duplicate=True),
        Output(IDs.SETTINGS_AI_KEY, "value", allow_duplicate=True),
        Output(IDs.SETTINGS_AI_MODEL, "value", allow_duplicate=True),
        Output(IDs.SETTINGS_STORE, "data", allow_duplicate=True),
        Output(IDs.SETTINGS_FEEDBACK, "children"),
        Input(IDs.BTN_SETTINGS_RESET, "n_clicks"),
        prevent_initial_call=True,
    )
    def _reset(n):
        if not n:
            return (no_update,) * 9
        d = _DEFAULTS
        return (
            d["interp"],
            d["skindepth"],
            d["dem"],
            d["cmap"],
            d["ai_provider"],
            d["ai_key"],
            d["ai_model"],
            d,
            [
                _status_icon("bi-check-circle-fill", "ok"),
                " All settings reset to defaults.",
            ],
        )

    # ── 9. Save profile JSON ──────────────────────────────────────────────────
    app.layout.children  # safe-guard: ensure layout is initialised

    @app.callback(
        Output("settings-download", "data"),
        Output(IDs.SETTINGS_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SETTINGS_SAVE, "n_clicks"),
        State(IDs.SETTINGS_STORE, "data"),
        prevent_initial_call=True,
    )
    def _save(n, stored):
        if not n:
            return no_update, no_update
        cfg = {**_DEFAULTS, **(stored or {})}
        cfg.pop("ai_key", None)  # never export key
        return (
            dict(
                content=json.dumps(cfg, indent=2),
                filename="pycsamt_settings.json",
            ),
            [
                _status_icon("bi-check-circle-fill", "ok"),
                " Profile saved (API key excluded).",
            ],
        )

    # ── 10. Load profile from uploaded JSON ───────────────────────────────────
    @app.callback(
        Output("settings-interp", "value", allow_duplicate=True),
        Output("settings-skindepth", "value", allow_duplicate=True),
        Output("settings-dem", "value", allow_duplicate=True),
        Output("settings-cmap", "value", allow_duplicate=True),
        Output(IDs.SETTINGS_AI_PROVIDER, "data", allow_duplicate=True),
        Output(IDs.SETTINGS_AI_MODEL, "value", allow_duplicate=True),
        Output(IDs.SETTINGS_STORE, "data", allow_duplicate=True),
        Output(IDs.SETTINGS_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.SETTINGS_UPLOAD, "contents"),
        State(IDs.SETTINGS_UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def _load(contents, filename):
        if not contents:
            return (no_update,) * 8
        try:
            _, encoded = contents.split(",", 1)
            cfg = {
                **_DEFAULTS,
                **json.loads(base64.b64decode(encoded).decode()),
            }
            return (
                cfg["interp"],
                cfg["skindepth"],
                cfg["dem"],
                cfg["cmap"],
                cfg["ai_provider"],
                cfg["ai_model"],
                cfg,
                [
                    _status_icon("bi-check-circle-fill", "ok"),
                    f" Profile loaded from {filename}.",
                ],
            )
        except Exception as exc:
            return (no_update,) * 7 + (
                [
                    _status_icon("bi-x-circle-fill", "danger"),
                    f" Could not parse {filename}: {exc}",
                ],
            )


def _status_icon(icon_cls: str, variant: str) -> object:
    from dash import html

    color = {
        "ok": "var(--green)",
        "warn": "var(--yellow)",
        "danger": "var(--red)",
    }.get(variant, "inherit")
    return html.I(className=f"bi {icon_cls} me-1", style={"color": color})
