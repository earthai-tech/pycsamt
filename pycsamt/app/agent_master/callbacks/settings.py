# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Settings callbacks — API keys, theme, export."""

from __future__ import annotations

import json
import os
from pathlib import Path

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs
from .._providers import (
    OFFLINE,
    default_model,
    env_var,
    is_llm,
    label_for,
    models_for,
)

_CFG_DIR = Path.home() / ".config" / "pycsamt"
_CFG_FILE = _CFG_DIR / "agent_master.json"


def _load_cfg() -> dict:
    if _CFG_FILE.exists():
        try:
            return json.loads(_CFG_FILE.read_text())
        except Exception:
            pass
    return {}


def _save_cfg(cfg: dict) -> None:
    _CFG_DIR.mkdir(parents=True, exist_ok=True)
    _CFG_FILE.write_text(json.dumps(cfg, indent=2))


# ── credential resolution ──────────────────────────────────────────────────────


def _env_key(provider: str) -> str:
    var = env_var(provider)
    return (os.environ.get(var, "") if var else "").strip()


def _stored_key(provider: str, cfg: dict) -> str:
    """The key to prefill: the saved one, else the environment fallback."""
    return (cfg.get(f"key_{provider}") or "").strip() or _env_key(provider)


def _source_for(provider: str, typed: str | None, cfg: dict) -> str:
    """Where the key in the box came from: saved / env / unsaved / none.

    One helper for every badge, so a key typed but not yet saved reads
    "Unsaved" the moment it diverges from what is on disk — rather than
    still claiming "Key saved".
    """
    typed = (typed or "").strip()
    saved = (cfg.get(f"key_{provider}") or "").strip()
    env = _env_key(provider)
    if not typed:
        return "env" if env else "none"
    if typed == saved:
        return "saved"
    if typed == env and not saved:
        return "env"
    return "unsaved"


def _badge(text: str, color: str, icon: str) -> html.Span:
    return html.Span(
        [html.I(className=f"bi {icon} me-1"), text],
        className="am-prov-badge",
        style={"color": f"var({color})"},
    )


def _status_badge(provider: str, source: str) -> html.Span:
    """Where the active key came from — so the user never guesses."""
    if not is_llm(provider):
        return _badge("Zero cost", "--fg-muted", "bi-shield-check")
    if source == "saved":
        return _badge("Key saved", "--tag-ok", "bi-check-circle-fill")
    if source == "env":
        return _badge("From environment", "--blue", "bi-terminal")
    if source == "unsaved":
        return _badge("Unsaved", "--peach", "bi-pencil-fill")
    return _badge("No key", "--fg-muted", "bi-exclamation-circle")


def _hint(provider: str) -> html.Span:
    var = env_var(provider)
    return html.Span(
        [
            "Stored locally in ",
            html.Code("~/.config/pycsamt/agent_master.json"),
            ". Leave empty to fall back to ",
            html.Code(var or "—"),
            ".",
        ]
    )


def register_settings(app) -> None:

    # On page load: hydrate STORE_SETTINGS from disk
    # so _run_agent always has provider / key even
    # before the user opens the settings panel.
    @app.callback(
        Output(IDs.STORE_SETTINGS, "data"),
        Input("am-root", "id"),
        prevent_initial_call=False,
    )
    def init_settings(_):
        cfg = _load_cfg()
        return cfg if cfg else no_update

    # Open settings offcanvas
    @app.callback(
        Output(IDs.CANVAS_SETTINGS, "is_open"),
        Input(IDs.BTN_SETTINGS, "n_clicks"),
        State(IDs.CANVAS_SETTINGS, "is_open"),
        prevent_initial_call=True,
    )
    def toggle_settings(n, is_open):
        if not n:
            raise PreventUpdate
        return not is_open

    # Selecting a provider swaps in *its* credentials: the key box + model
    # list appear, prefilled from unsaved drafts → saved config → env var.
    # Offline shows the zero-cost note instead. Runs on load too, so the
    # default "offline" selection starts correct.
    @app.callback(
        Output(IDs.PROVIDER_PANEL, "style"),
        Output(IDs.OFFLINE_NOTE, "style"),
        Output(IDs.KEY_INPUT, "value"),
        Output(IDs.MODEL_SELECT, "options"),
        Output(IDs.MODEL_SELECT, "value"),
        Output(IDs.PROVIDER_HINT, "children"),
        Output(IDs.PROVIDER_BADGE, "children"),
        Input(IDs.ACTIVE_PROVIDER, "value"),
        State(IDs.STORE_KEY_DRAFTS, "data"),
        prevent_initial_call=False,
    )
    def sync_provider_panel(provider, drafts):
        hide = {"display": "none"}
        show = {"display": "block"}
        if not is_llm(provider):
            return (
                hide,
                show,
                "",
                [],
                None,
                "",
                _status_badge(OFFLINE, "none"),
            )

        cfg = _load_cfg()
        drafts = drafts or {}
        # A key typed this session (not yet saved) wins over disk/env.
        draft_key = drafts.get(f"key_{provider}")
        key = (
            draft_key if draft_key is not None else _stored_key(provider, cfg)
        )

        models = models_for(provider)
        model = (
            drafts.get(f"model_{provider}")
            or cfg.get(f"model_{provider}")
            or default_model(provider)
        )
        if model not in models and models:
            model = models[0]
        return (
            show,
            hide,
            key,
            [{"label": m, "value": m} for m in models],
            model,
            _hint(provider),
            _status_badge(provider, _source_for(provider, key, cfg)),
        )

    # Keep unsaved edits per provider so flipping the dropdown (or closing
    # the panel) never silently discards a key you just typed — and keep the
    # badge honest while typing ("Unsaved" the moment it diverges from disk).
    @app.callback(
        Output(IDs.STORE_KEY_DRAFTS, "data"),
        Output(
            IDs.PROVIDER_BADGE,
            "children",
            allow_duplicate=True,
        ),
        Input(IDs.KEY_INPUT, "value"),
        Input(IDs.MODEL_SELECT, "value"),
        State(IDs.ACTIVE_PROVIDER, "value"),
        State(IDs.STORE_KEY_DRAFTS, "data"),
        prevent_initial_call=True,
    )
    def stash_draft(key, model, provider, drafts):
        if not is_llm(provider):
            raise PreventUpdate
        drafts = dict(drafts or {})
        drafts[f"key_{provider}"] = key or ""
        if model:
            drafts[f"model_{provider}"] = model
        badge = _status_badge(
            provider, _source_for(provider, key, _load_cfg())
        )
        return drafts, badge

    # Show / hide the key (a password field you can never read back is a
    # good way to save the wrong key twice).
    @app.callback(
        Output(IDs.KEY_INPUT, "type"),
        Output(IDs.BTN_KEY_REVEAL, "children"),
        Input(IDs.BTN_KEY_REVEAL, "n_clicks"),
        prevent_initial_call=True,
    )
    def reveal_key(n):
        if not n:
            raise PreventUpdate
        shown = bool(n % 2)
        icon = "bi bi-eye-slash" if shown else "bi bi-eye"
        return ("text" if shown else "password"), html.I(className=icon)

    # Opening the panel restores the saved provider + preferences; the key
    # and model then follow from sync_provider_panel.
    @app.callback(
        Output(IDs.ACTIVE_PROVIDER, "value"),
        Output(IDs.EXPORT_FORMAT, "value"),
        Output(IDs.OUTPUT_DIR, "value"),
        Output(IDs.LINE_REGISTRY, "value"),
        Input(IDs.CANVAS_SETTINGS, "is_open"),
        prevent_initial_call=True,
    )
    def load_prefs(is_open):
        if not is_open:
            raise PreventUpdate
        cfg = _load_cfg()
        return (
            cfg.get("provider", OFFLINE),
            cfg.get("export_fmt", "png"),
            cfg.get("output_dir", ""),
            cfg.get("line_registry", ""),
        )

    # Save: merge this session's drafts over the stored config so switching
    # providers never wipes another provider's key, then persist.
    @app.callback(
        Output(
            IDs.STORE_SETTINGS,
            "data",
            allow_duplicate=True,
        ),
        Output(IDs.KEYS_STATUS, "children"),
        Output(
            IDs.PROVIDER_BADGE,
            "children",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_SAVE_KEYS, "n_clicks"),
        State(IDs.ACTIVE_PROVIDER, "value"),
        State(IDs.KEY_INPUT, "value"),
        State(IDs.MODEL_SELECT, "value"),
        State(IDs.EXPORT_FORMAT, "value"),
        State(IDs.OUTPUT_DIR, "value"),
        State(IDs.LINE_REGISTRY, "value"),
        State(IDs.STORE_KEY_DRAFTS, "data"),
        prevent_initial_call=True,
    )
    def save_settings(
        n,
        provider,
        key,
        model,
        export_fmt,
        output_dir,
        line_registry,
        drafts,
    ):
        if not n:
            raise PreventUpdate
        cfg = _load_cfg()
        # Drafts first (other providers edited this session), then the live
        # inputs for the provider currently on screen.
        for field, value in (drafts or {}).items():
            if value:
                cfg[field] = value
        if is_llm(provider):
            cfg[f"key_{provider}"] = (key or "").strip()
            cfg[f"model_{provider}"] = model or default_model(provider) or ""
        cfg["provider"] = provider or OFFLINE
        cfg["export_fmt"] = export_fmt or "png"
        cfg["output_dir"] = output_dir or ""
        cfg["line_registry"] = line_registry or ""
        _save_cfg(cfg)

        source = (
            _source_for(provider, key, cfg) if is_llm(provider) else "none"
        )
        summary = f"Saved — {label_for(provider)}"
        if source == "saved":
            summary += f" · {cfg[f'model_{provider}']}"
        status = html.Span(
            [
                html.I(
                    className="bi bi-check-circle-fill me-1",
                    style={"color": "var(--tag-ok)"},
                ),
                summary,
            ],
            style={"color": "var(--tag-ok)"},
        )
        return dict(cfg), status, _status_badge(provider, source)

    # Theme toggle
    @app.callback(
        Output(IDs.STORE_THEME, "data"),
        Output("am-theme-icon", "className"),
        Input(IDs.BTN_THEME, "n_clicks"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def toggle_theme(n, current):
        if not n:
            raise PreventUpdate
        new = "dark" if current == "light" else "light"
        icon = "bi bi-sun-fill" if new == "dark" else "bi bi-moon-stars"
        return new, icon
