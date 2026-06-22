# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Line-selection modal callbacks.

When the user references an ambiguous line
(e.g. 'line 1') and multiple survey lines
are loaded, this modal presents a checklist
of available lines so the user can pick
which one(s) to process.
"""

from __future__ import annotations

import threading

from dash import ALL, Input, Output, State
from dash import ctx, html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs
from .chat import (
    _new_job,
    _run_agent,
    _thinking_bubble,
    _ts,
    _waiting_bubble,
    _quick_workflow,
    _WF_LABELS,
)


def register_line_sel(app) -> None:

    # 1. Open modal when STORE_PENDING has
    #    "disambiguation": "lines"
    @app.callback(
        Output(
            IDs.MODAL_LINE_SEL, "is_open"
        ),
        Output(
            IDs.LINE_SEL_LIST, "options"
        ),
        Output(
            IDs.LINE_SEL_LIST, "value"
        ),
        Input(
            IDs.STORE_PENDING, "data"
        ),
        prevent_initial_call=True,
    )
    def open_line_modal(pending):
        if (
            not pending
            or pending.get("disambiguation")
            != "lines"
        ):
            raise PreventUpdate
        groups = pending.get("groups", {})
        names = list(groups.keys())
        opts = [
            {"label": nm, "value": nm}
            for nm in names
        ]
        return True, opts, []

    # 2. Run selected / Run all — start job
    @app.callback(
        Output(
            IDs.CHAT_WINDOW,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_JOB,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.INTERVAL_POLL,
            "disabled",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_MESSAGES,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_PENDING,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.MODAL_LINE_SEL,
            "is_open",
            allow_duplicate=True,
        ),
        Output(
            IDs.LINE_SEL_STATUS, "children"
        ),
        Input(
            IDs.BTN_LINE_RUN_SEL, "n_clicks"
        ),
        Input(
            IDs.BTN_LINE_RUN_ALL, "n_clicks"
        ),
        State(IDs.LINE_SEL_LIST, "value"),
        State(IDs.STORE_PENDING, "data"),
        State(IDs.STORE_EDI, "data"),
        State(IDs.STORE_SETTINGS, "data"),
        State(
            IDs.STORE_INV_CONFIG, "data"
        ),
        State(
            IDs.STORE_MESSAGES, "data"
        ),
        State(IDs.CHAT_WINDOW, "children"),
        prevent_initial_call=True,
    )
    def run_with_lines(
        n_sel, n_all,
        sel_lines,
        pending,
        edi_store,
        settings,
        inv_config,
        stored_msgs,
        current_msgs,
    ):
        triggered = ctx.triggered_id
        if not triggered:
            raise PreventUpdate

        text = (pending or {}).get("text", "")
        groups = (pending or {}).get(
            "groups", {}
        )

        if triggered == IDs.BTN_LINE_RUN_ALL:
            sel = list(groups.keys())
        else:
            sel = sel_lines or []
            if not sel:
                return (
                    no_update, no_update,
                    no_update, no_update,
                    no_update, True,
                    html.Span(
                        "Select at least one"
                        " line first.",
                        style={
                            "color": "orange"
                        },
                    ),
                )

        msgs = list(current_msgs or [])
        new_stored = list(stored_msgs or [])

        # Replace line-waiting bubble if present
        _replaced = False
        for i, c in enumerate(msgs):
            cid = (
                c.get("props", {}).get("id", "")
                if isinstance(c, dict)
                else ""
            )
            if cid == "am-line-waiting-bubble":
                _replaced = True
                # Check if param modal needed
                wf = _quick_workflow(text)
                if wf:
                    msgs[i] = _waiting_bubble(wf)
                    new_pend = {
                        "workflow": wf,
                        "text": text,
                        "selected_lines": sel,
                    }
                    return (
                        msgs, no_update,
                        no_update, no_update,
                        new_pend, False,
                        no_update,
                    )
                else:
                    msgs[i] = _thinking_bubble(
                        [{
                            "label": (
                                f"Loading"
                                f" {len(sel)}"
                                " line(s)..."
                            ),
                            "status": "running",
                        }]
                    )
                break

        if not _replaced:
            wf = _quick_workflow(text)
            if wf:
                msgs.append(
                    _waiting_bubble(wf)
                )
                new_pend = {
                    "workflow": wf,
                    "text": text,
                    "selected_lines": sel,
                }
                return (
                    msgs, no_update,
                    no_update, no_update,
                    new_pend, False,
                    no_update,
                )
            msgs.append(
                _thinking_bubble(
                    [{
                        "label": (
                            f"Loading {len(sel)}"
                            " line(s)..."
                        ),
                        "status": "running",
                    }]
                )
            )

        _edi_use = dict(edi_store or {})
        _edi_use["selected_lines"] = sel
        jid = _new_job()
        threading.Thread(
            target=_run_agent,
            args=(
                jid, text, _edi_use,
                settings or {},
                inv_config or {},
            ),
            daemon=True,
        ).start()

        return (
            msgs,
            {"jid": jid},
            False,
            new_stored,
            {},
            False,
            no_update,
        )
