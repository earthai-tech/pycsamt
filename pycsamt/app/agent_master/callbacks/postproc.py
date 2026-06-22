# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Post-correction action modal callbacks.

After static-shift, denoising, or QC
workflows finish, the post-proc modal lets
the user either apply corrected data to the
active session or export to a custom folder.
"""

from __future__ import annotations

import time
from pathlib import Path

from dash import Input, Output, State
from dash import html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs


def _export(sites, dest: Path) -> list:
    """Write corrected sites as EDI files."""
    from pycsamt.site.export import write_sites
    return write_sites(
        sites, dest, exist_ok=True
    )


def register_postproc(app) -> None:

    # 1. Open modal when STORE_POSTPROC is set
    @app.callback(
        Output(
            IDs.MODAL_POSTPROC, "is_open"
        ),
        Output(
            IDs.POSTPROC_MSG, "children"
        ),
        Input(
            IDs.STORE_POSTPROC, "data"
        ),
        prevent_initial_call=True,
    )
    def open_modal(data):
        if not (data or {}).get("jid"):
            raise PreventUpdate
        wtype = data.get(
            "workflow", "processing"
        )
        return True, (
            f"The {wtype} correction finished."
            " What would you like to do"
            " with the corrected data?"
        )

    # 2. "Apply to session" — export corrected
    #    EDIs to output_dir, update STORE_EDI.
    @app.callback(
        Output(
            IDs.STORE_EDI,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.POSTPROC_STATUS, "children"
        ),
        Output(
            IDs.MODAL_POSTPROC,
            "is_open",
            allow_duplicate=True,
        ),
        Input(
            IDs.BTN_POSTPROC_APPLY, "n_clicks"
        ),
        State(
            IDs.STORE_POSTPROC, "data"
        ),
        State(
            IDs.STORE_SETTINGS, "data"
        ),
        State(IDs.STORE_EDI, "data"),
        prevent_initial_call=True,
    )
    def apply_to_session(
        n, postproc, settings, edi_store
    ):
        if not n:
            raise PreventUpdate
        from .chat import _CORR_CACHE
        jid = (postproc or {}).get("jid")
        wtype = (postproc or {}).get(
            "workflow", "corrected"
        )
        base = (
            (postproc or {}).get(
                "output_dir", ""
            )
            or (settings or {}).get(
                "output_dir", ""
            )
            or "pycsamt_workflow_output"
        ).strip()
        sites = _CORR_CACHE.pop(jid, None)
        if sites is None:
            return (
                no_update,
                html.Span(
                    "Corrected data not in "
                    "memory — re-run the "
                    "workflow.",
                    style={"color": "red"},
                ),
                True,
            )
        dest = (
            Path(base)
            / f"{wtype}_{int(time.time())}"
        )
        try:
            paths = _export(sites, dest)
            new_edi = dict(edi_store or {})
            new_edi["path"] = str(dest)
            msg = html.Span(
                f"{len(paths)} EDI(s) saved"
                f" to {dest}",
                style={
                    "color": "var(--tag-ok)"
                },
            )
            return new_edi, msg, False
        except Exception as exc:
            return (
                no_update,
                html.Span(
                    str(exc),
                    style={"color": "red"},
                ),
                True,
            )

    # 3. "Export to folder" — toggle collapse
    @app.callback(
        Output(
            IDs.POSTPROC_COLLAPSE, "is_open"
        ),
        Input(
            IDs.BTN_POSTPROC_EXPORT, "n_clicks"
        ),
        prevent_initial_call=True,
    )
    def show_export_path(n):
        return bool(n and n % 2 == 1)

    # 4. Confirm export to custom folder path
    @app.callback(
        Output(
            IDs.POSTPROC_STATUS,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.MODAL_POSTPROC,
            "is_open",
            allow_duplicate=True,
        ),
        Output(
            IDs.POSTPROC_COLLAPSE,
            "is_open",
            allow_duplicate=True,
        ),
        Input(
            IDs.BTN_POSTPROC_OK, "n_clicks"
        ),
        State(IDs.POSTPROC_PATH, "value"),
        State(
            IDs.STORE_POSTPROC, "data"
        ),
        prevent_initial_call=True,
    )
    def confirm_export(n, path_val, postproc):
        if not n:
            raise PreventUpdate
        if not (path_val or "").strip():
            return (
                html.Span(
                    "Enter a folder path "
                    "first.",
                    style={"color": "orange"},
                ),
                True,
                True,
            )
        from .chat import _CORR_CACHE
        jid = (postproc or {}).get("jid")
        sites = _CORR_CACHE.pop(jid, None)
        if sites is None:
            return (
                html.Span(
                    "Corrected data not in "
                    "memory — re-run workflow.",
                    style={"color": "red"},
                ),
                True,
                True,
            )
        dest = Path(path_val.strip())
        try:
            paths = _export(sites, dest)
            return (
                html.Span(
                    f"Exported {len(paths)}"
                    f" EDI(s) to {dest}",
                    style={
                        "color": "var(--tag-ok)"
                    },
                ),
                False,
                False,
            )
        except Exception as exc:
            return (
                html.Span(
                    str(exc),
                    style={"color": "red"},
                ),
                True,
                True,
            )
