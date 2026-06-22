# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Settings callbacks — API keys, theme, export."""

from __future__ import annotations

import json
import os
from pathlib import Path

from dash import Input, Output, State
from dash import html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs

_CFG_DIR = Path.home() / ".config" / "pycsamt"
_CFG_FILE = _CFG_DIR / "agent_master.json"

_DEF_MODELS = {
    "claude":   "claude-sonnet-4-6",
    "openai":   "gpt-4o",
    "gemini":   "gemini-2.0-flash",
    "deepseek": "deepseek-chat",
}


def _load_cfg() -> dict:
    if _CFG_FILE.exists():
        try:
            return json.loads(
                _CFG_FILE.read_text()
            )
        except Exception:
            pass
    return {}


def _save_cfg(cfg: dict) -> None:
    _CFG_DIR.mkdir(parents=True, exist_ok=True)
    _CFG_FILE.write_text(
        json.dumps(cfg, indent=2)
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

    # Populate keys + models from disk on open
    @app.callback(
        Output(IDs.KEY_CLAUDE,   "value"),
        Output(IDs.KEY_OPENAI,   "value"),
        Output(IDs.KEY_GEMINI,   "value"),
        Output(IDs.KEY_DEEPSEEK, "value"),
        Output(IDs.MODEL_CLAUDE,   "value"),
        Output(IDs.MODEL_OPENAI,   "value"),
        Output(IDs.MODEL_GEMINI,   "value"),
        Output(IDs.MODEL_DEEPSEEK, "value"),
        Output(IDs.ACTIVE_PROVIDER, "value"),
        Output(IDs.EXPORT_FORMAT, "value"),
        Output(IDs.OUTPUT_DIR, "value"),
        Output(IDs.LINE_REGISTRY, "value"),
        Input(IDs.CANVAS_SETTINGS, "is_open"),
        prevent_initial_call=True,
    )
    def load_keys(is_open):
        if not is_open:
            raise PreventUpdate
        cfg = _load_cfg()
        return (
            cfg.get("key_claude")
            or os.environ.get(
                "ANTHROPIC_API_KEY", ""
            ),
            cfg.get("key_openai")
            or os.environ.get(
                "OPENAI_API_KEY", ""
            ),
            cfg.get("key_gemini")
            or os.environ.get(
                "GOOGLE_API_KEY", ""
            ),
            cfg.get("key_deepseek")
            or os.environ.get(
                "DEEPSEEK_API_KEY", ""
            ),
            cfg.get(
                "model_claude",
                _DEF_MODELS["claude"],
            ),
            cfg.get(
                "model_openai",
                _DEF_MODELS["openai"],
            ),
            cfg.get(
                "model_gemini",
                _DEF_MODELS["gemini"],
            ),
            cfg.get(
                "model_deepseek",
                _DEF_MODELS["deepseek"],
            ),
            cfg.get("provider", "offline"),
            cfg.get("export_fmt", "png"),
            cfg.get("output_dir", ""),
            cfg.get("line_registry", ""),
        )

    # Save keys + models button
    @app.callback(
        Output(
            IDs.STORE_SETTINGS, "data",
            allow_duplicate=True,
        ),
        Output(IDs.KEYS_STATUS, "children"),
        Input(IDs.BTN_SAVE_KEYS, "n_clicks"),
        State(IDs.KEY_CLAUDE,   "value"),
        State(IDs.KEY_OPENAI,   "value"),
        State(IDs.KEY_GEMINI,   "value"),
        State(IDs.KEY_DEEPSEEK, "value"),
        State(IDs.MODEL_CLAUDE,   "value"),
        State(IDs.MODEL_OPENAI,   "value"),
        State(IDs.MODEL_GEMINI,   "value"),
        State(IDs.MODEL_DEEPSEEK, "value"),
        State(IDs.ACTIVE_PROVIDER, "value"),
        State(IDs.EXPORT_FORMAT, "value"),
        State(IDs.OUTPUT_DIR, "value"),
        State(IDs.LINE_REGISTRY, "value"),
        prevent_initial_call=True,
    )
    def save_keys(
        n,
        k_claude, k_openai,
        k_gemini, k_deepseek,
        m_claude, m_openai,
        m_gemini, m_deepseek,
        provider, export_fmt,
        output_dir,
        line_registry,
    ):
        if not n:
            raise PreventUpdate
        cfg = {
            "key_claude":     k_claude or "",
            "key_openai":     k_openai or "",
            "key_gemini":     k_gemini or "",
            "key_deepseek":   k_deepseek or "",
            "model_claude":   (
                m_claude
                or _DEF_MODELS["claude"]
            ),
            "model_openai":   (
                m_openai
                or _DEF_MODELS["openai"]
            ),
            "model_gemini":   (
                m_gemini
                or _DEF_MODELS["gemini"]
            ),
            "model_deepseek": (
                m_deepseek
                or _DEF_MODELS["deepseek"]
            ),
            "provider":       provider or "offline",
            "export_fmt":     export_fmt or "png",
            "output_dir":     output_dir or "",
            "line_registry":  line_registry or "",
        }
        _save_cfg(cfg)
        store = {k: v for k, v in cfg.items()}
        status = html.Span(
            [
                html.I(
                    className=(
                        "bi bi-check-circle-fill"
                        " me-1"
                    ),
                    style={
                        "color": "var(--tag-ok)"
                    },
                ),
                "Keys saved.",
            ],
            style={"color": "var(--tag-ok)"},
        )
        return store, status

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
        new = (
            "dark" if current == "light"
            else "light"
        )
        icon = (
            "bi bi-sun-fill"
            if new == "dark"
            else "bi bi-moon-stars"
        )
        return new, icon
