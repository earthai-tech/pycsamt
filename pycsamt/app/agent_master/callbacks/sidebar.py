# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Sidebar toggle, history, and figure callbacks."""

from __future__ import annotations

from datetime import datetime

from dash import ALL, Input, Output, State
from dash import ctx, dcc, html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs
from ..layout import _chat_welcome


# ── helpers ────────────────────────────────────────

def _restore_bubble(msg: dict) -> html.Div:
    role = msg.get("role", "user")
    text = msg.get("content", "")
    ts   = msg.get("ts", "")
    if role == "user":
        return html.Div(
            [
                html.Div(
                    html.Div(
                        [
                            html.Div(
                                text,
                                className=(
                                    "am-bubble user"
                                ),
                            ),
                            html.Div(
                                ts,
                                className="am-ts",
                            ),
                        ]
                    ),
                    style={"maxWidth": "100%"},
                ),
                html.Div(
                    html.I(
                        className=(
                            "bi bi-person-fill"
                        )
                    ),
                    className="am-avatar user",
                ),
            ],
            className="am-msg-row user",
        )
    body = (
        dcc.Markdown(text, className="am-md", link_target="_blank")
        if (text or "").strip()
        else html.P("(no content)")
    )
    return html.Div(
        [
            html.Div(
                html.I(
                    className="bi bi-robot"
                ),
                className="am-avatar agent",
            ),
            html.Div(
                [
                    html.Div(body),
                    html.Div(
                        ts,
                        className="am-ts",
                    ),
                ],
                className="am-bubble agent am-bubble-flat",
            ),
        ],
        className="am-msg-row",
    )


def register_sidebar(app) -> None:

    # ── toggle open / collapsed ───────────────
    app.clientside_callback(
        """
        function(n1, n2, cls) {
            if (!n1 && !n2) {
                return window.dash_clientside.no_update;
            }
            var c = cls
                || "am-sidebar am-sidebar-open";
            if (c.includes(
                    "am-sidebar-collapsed"
            )) {
                return c.replace(
                    "am-sidebar-collapsed",
                    "am-sidebar-open"
                );
            }
            return c.replace(
                "am-sidebar-open",
                "am-sidebar-collapsed"
            );
        }
        """,
        Output(IDs.SIDEBAR, "className"),
        Input(IDs.BTN_SIDEBAR, "n_clicks"),
        Input(
            "am-sidebar-close", "n_clicks"
        ),
        State(IDs.SIDEBAR, "className"),
        prevent_initial_call=True,
    )

    # ── segmented switcher (Chat / Session) ──
    app.clientside_callback(
        """
        function(nC, nS) {
            var ctx = window.dash_clientside.callback_context;
            var trig = (ctx && ctx.triggered && ctx.triggered.length)
                ? ctx.triggered[0].prop_id : "";
            var tab = (trig.indexOf("am-seg-session") === 0)
                ? "session" : "chat";
            var off = "am-tab-panel";
            var on = "am-tab-panel am-tab-active";
            var seg = "am-seg-btn";
            var segA = "am-seg-btn am-seg-active";
            return [
                tab === "chat" ? on : off,
                tab === "session" ? on : off,
                tab === "chat" ? segA : seg,
                tab === "session" ? segA : seg,
                tab
            ];
        }
        """,
        Output(IDs.PANEL_CHAT, "className"),
        Output(IDs.PANEL_SESSION, "className"),
        Output(IDs.SEG_CHAT, "className"),
        Output(IDs.SEG_SESSION, "className"),
        Output(IDs.STORE_SB_TAB, "data"),
        Input(IDs.SEG_CHAT, "n_clicks"),
        Input(IDs.SEG_SESSION, "n_clicks"),
        prevent_initial_call=True,
    )

    # ── "open trace" (↗ on a message) → open the sidebar on the Session
    # tab, where Recent runs lives. Reuses the segmented-switch outputs
    # (allow_duplicate) so it stays a pure clientside toggle. The guard
    # ignores the all-zero firing that happens when a new bubble adds a
    # button to the pattern-matching group.
    app.clientside_callback(
        """
        function(clicks){
            var nu = window.dash_clientside.no_update;
            var any = (clicks || []).some(function(c){ return c; });
            if(!any){ return [nu, nu, nu, nu, nu, nu]; }
            return [
                "am-sidebar am-sidebar-open",
                "am-tab-panel",
                "am-tab-panel am-tab-active",
                "am-seg-btn",
                "am-seg-btn am-seg-active",
                "session"
            ];
        }
        """,
        Output(IDs.SIDEBAR, "className", allow_duplicate=True),
        Output(IDs.PANEL_CHAT, "className", allow_duplicate=True),
        Output(IDs.PANEL_SESSION, "className", allow_duplicate=True),
        Output(IDs.SEG_CHAT, "className", allow_duplicate=True),
        Output(IDs.SEG_SESSION, "className", allow_duplicate=True),
        Output(IDs.STORE_SB_TAB, "data", allow_duplicate=True),
        Input({"type": "am-trace-open", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )

    # ── new chat ──────────────────────────────
    @app.callback(
        Output(
            IDs.STORE_MESSAGES,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_EDI,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_FIGS,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_JOB,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.CHAT_WINDOW,
            "children",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_NEW_CHAT, "n_clicks"),
        Input(IDs.BTN_NEW_SESSION, "n_clicks"),
        prevent_initial_call=True,
    )
    def _new_chat(n_chat, n_session):
        if not n_chat and not n_session:
            raise PreventUpdate
        # clear the assistant session so a new chat starts context-free
        try:
            from .chat import _reset_session
            _reset_session()
        except Exception:  # noqa: BLE001
            pass
        return [], {}, {}, {}, [_chat_welcome()]

    # ── auto-save to history ──────────────────
    # Fires whenever STORE_MESSAGES changes.
    # Only saves when last message is assistant.
    # Skips if messages are an exact restore of
    # an existing entry (no new content).
    @app.callback(
        Output(
            IDs.STORE_HISTORY,
            "data",
        ),
        Input(IDs.STORE_MESSAGES, "data"),
        State(IDs.STORE_HISTORY, "data"),
        prevent_initial_call=True,
    )
    def _auto_save(messages, history):
        if not messages or len(messages) < 2:
            raise PreventUpdate
        last = messages[-1]
        if last.get("role") != "assistant":
            raise PreventUpdate
        history = list(history or [])
        # Skip if messages are identical to any
        # stored entry — this is a pure restore.
        for h in history:
            if h.get("messages") == messages:
                raise PreventUpdate
        preview = ""
        for m in messages:
            if m.get("role") == "user":
                preview = (
                    (m.get("content") or "")[:60]
                )
                break
        entry = {
            "ts": datetime.now().strftime(
                "%Y-%m-%d %H:%M"
            ),
            "preview": preview or "(empty)",
            "messages": messages,
        }
        # Replace any existing entry with the
        # same preview (across whole history list)
        # and move it to the top.
        for j, h in enumerate(history):
            if h.get("preview") == preview:
                history.pop(j)
                break
        history.insert(0, entry)
        return history[:20]

    # ── explicit save (Save button) ───────────
    @app.callback(
        Output(
            IDs.STORE_HISTORY,
            "data",
            allow_duplicate=True,
        ),
        Input(
            IDs.BTN_SAVE_SESSION, "n_clicks"
        ),
        State(IDs.STORE_MESSAGES, "data"),
        State(IDs.STORE_HISTORY, "data"),
        prevent_initial_call=True,
    )
    def _save_history(n, messages, history):
        if not n or not messages:
            raise PreventUpdate
        history = list(history or [])
        preview = ""
        for m in messages:
            if m.get("role") == "user":
                preview = (
                    (m.get("content") or "")[:60]
                )
                break
        entry = {
            "ts": datetime.now().strftime(
                "%Y-%m-%d %H:%M"
            ),
            "preview": preview or "(empty)",
            "messages": messages,
        }
        # Replace any existing same-preview entry
        for j, h in enumerate(history):
            if h.get("preview") == preview:
                history.pop(j)
                break
        history.insert(0, entry)
        return history[:20]

    # ── restore session ───────────────────────
    @app.callback(
        Output(
            IDs.STORE_MESSAGES,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.CHAT_WINDOW,
            "children",
            allow_duplicate=True,
        ),
        Input(
            {"type": "am-hist-item",
             "index": ALL},
            "n_clicks",
        ),
        State(IDs.STORE_HISTORY, "data"),
        prevent_initial_call=True,
    )
    def _restore_session(clicks, history):
        if not any(clicks) or not history:
            raise PreventUpdate
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict):
            raise PreventUpdate
        idx = triggered.get("index", 0)
        if idx >= len(history):
            raise PreventUpdate
        entry = history[idx]
        messages = entry.get("messages", [])
        if not messages:
            raise PreventUpdate
        banner = html.Div(
            [
                html.I(
                    className=(
                        "bi bi-clock-history me-2"
                    )
                ),
                f"Session from "
                f"{entry.get('ts', '')}",
                html.Span(
                    " — figures not restored.",
                    style={"opacity": ".65"},
                ),
            ],
            className="am-restore-banner",
        )
        children = [banner] + [
            _restore_bubble(m)
            for m in messages
        ]
        return messages, children

    # ── delete history entry ──────────────────
    @app.callback(
        Output(
            IDs.STORE_HISTORY,
            "data",
            allow_duplicate=True,
        ),
        Input(
            {"type": "am-hist-del",
             "index": ALL},
            "n_clicks",
        ),
        State(IDs.STORE_HISTORY, "data"),
        prevent_initial_call=True,
    )
    def _delete_entry(clicks, history):
        if not any(clicks) or not history:
            raise PreventUpdate
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict):
            raise PreventUpdate
        idx = triggered.get("index", 0)
        history = list(history or [])
        if idx < len(history):
            history.pop(idx)
        return history

    # ── render history list ───────────────────
    @app.callback(
        Output(
            IDs.SIDEBAR_HISTORY, "children"
        ),
        Input(IDs.STORE_HISTORY, "data"),
    )
    def _render_history(history):
        if not history:
            return html.Div(
                "No saved sessions yet.",
                className="am-sidebar-empty",
            )
        items = []
        for i, entry in enumerate(history):
            items.append(
                html.Div(
                    [
                        html.Button(
                            [
                                html.Div(
                                    entry.get(
                                        "ts", ""
                                    ),
                                    className=(
                                        "am-hitem-ts"
                                    ),
                                ),
                                html.Div(
                                    entry.get(
                                        "preview",
                                        "...",
                                    ),
                                    className=(
                                        "am-hitem-"
                                        "preview"
                                    ),
                                ),
                            ],
                            className=(
                                "am-history-item"
                            ),
                            id={
                                "type": (
                                    "am-hist-item"
                                ),
                                "index": i,
                            },
                            n_clicks=0,
                        ),
                        html.Button(
                            html.I(
                                className="bi bi-x"
                            ),
                            className=(
                                "am-hist-del-btn"
                            ),
                            id={
                                "type": (
                                    "am-hist-del"
                                ),
                                "index": i,
                            },
                            title="Remove",
                            n_clicks=0,
                        ),
                    ],
                    className="am-hist-row",
                )
            )
        return items

    # ── render figure thumbnails ──────────────
    @app.callback(
        Output(
            IDs.SIDEBAR_FIGS, "children"
        ),
        Input(IDs.STORE_FIGS, "data"),
    )
    def _render_figs(figs):
        if not figs:
            return html.Div(
                "No figures yet.",
                className="am-sidebar-empty",
            )
        return [
            html.Div(
                html.Button(
                    html.Img(
                        src=(
                            "data:image/png;"
                            f"base64,{info['b64']}"
                        ),
                        className=(
                            "am-sb-fig-thumb"
                        ),
                    ),
                    id={
                        "type": "am-fig-open",
                        "key": key,
                    },
                    className="am-sb-fig-btn",
                    title=info.get(
                        "title", key
                    ),
                    n_clicks=0,
                ),
                className="am-sb-fig-item",
            )
            for key, info in figs.items()
        ]
