# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Start on the survey handed over by the desktop (``?handoff=<token>``).

The desktop writes its in-memory survey (edits and processing included)
with :func:`pycsamt.app.agent_master._handoff.write_handoff` and opens
Agent Master on ``/?handoff=<token>``.  This callback loads that session
into ``STORE_EDI`` exactly as the Load dialog would, skips the welcome
splash and says where the data came from.  Without a token Agent Master
starts from scratch as before.
"""

from __future__ import annotations

from dash import Input, Output, html
from dash.exceptions import PreventUpdate

from .._handoff import read_handoff, token_from_search
from .._ids import IDs


def handoff_note(session: dict) -> html.Div:
    n = int(session.get("n_edi", 0))
    lines = len(session.get("groups") or {})
    state = ("edited in the desktop" if session.get("edited")
             else "as loaded in the desktop")
    names = ", ".join(list(session.get("groups") or {})[:6])
    return html.Div(
        [
            html.I(className="bi bi-pc-display-horizontal me-2"),
            html.Span([
                html.B(f"{n} station{'s' if n != 1 else ''} · "
                       f"{lines} line{'s' if lines != 1 else ''}"),
                f" received from the pyCSAMT desktop ({state})",
                html.Br(),
                html.Small(names + (" …" if lines > 6 else "")),
            ]),
        ],
        className="am-handoff-note",
        style={"display": "flex", "alignItems": "flex-start",
               "gap": "4px", "margin": "10px auto", "maxWidth": "560px",
               "padding": "10px 14px", "borderRadius": "10px",
               "border": "1px solid rgba(28,126,214,.35)",
               "background": "rgba(28,126,214,.08)", "textAlign": "left"},
    )


def register_handoff(app) -> None:
    @app.callback(
        Output(IDs.STORE_EDI, "data", allow_duplicate=True),
        Output(IDs.EDI_BADGE, "className", allow_duplicate=True),
        Output(IDs.EDI_BADGE_TEXT, "children", allow_duplicate=True),
        Output(IDs.SPLASH_OVERLAY, "className", allow_duplicate=True),
        Output(IDs.WELCOME_NOTE, "children"),
        Input(IDs.URL, "search"),
        prevent_initial_call="initial_duplicate",
    )
    def load_handoff(search):
        session = read_handoff(token_from_search(search))
        if session is None:
            raise PreventUpdate  # no hand-off: the usual welcome
        n = int(session["n_edi"])
        lines = len(session["groups"])
        badge = (f"{n} station{'s' if n != 1 else ''} · {lines} line(s) · "
                 "desktop")
        return (session, "am-edi-badge visible", badge,
                "wlc-overlay wlc-gone", handoff_note(session))
