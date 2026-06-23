# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""EDI load + line detection callbacks."""

from __future__ import annotations

import base64
import os
import re
import tempfile
from pathlib import Path

from dash import ALL, Input, Output, State
from dash import ctx, html, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .._ids import IDs

_LINE_COLORS = [
    "#1c7ed6", "#2f9e44", "#e67700",
    "#7048e8", "#c92a2a", "#0c8599",
    "#a61e4d", "#5c7cfa", "#f76707",
    "#37b24d",
]
_SHOW = {"display": "block"}
_HIDE = {"display": "none"}


def _color(i: int) -> str:
    return _LINE_COLORS[i % len(_LINE_COLORS)]


# ── Line detection helpers ─────────────────────────

def _detect_from_ids(
    station_ids: list[str],
) -> dict[str, list[str]]:
    r"""
    Auto-detect line groups from station IDs.

    Prefix before first '-' or '_';
    numeric prefix gets an L-prefix
    (e.g. '22-001' -> L22).
    """
    groups: dict[str, list[str]] = {}
    for sid in station_ids:
        m = re.match(
            r"^([A-Za-z]*\d+)", str(sid)
        )
        prefix = (
            m.group(1) if m else str(sid)[:4]
        )
        if prefix.isdigit():
            prefix = f"L{prefix}"
        groups.setdefault(
            prefix, []
        ).append(sid)
    return groups


def _folder_to_lines(
    edi_dir: str,
) -> dict[str, list[str]]:
    """Return {subfolder_name: [edi_files]}."""
    p = Path(edi_dir)
    groups: dict[str, list[str]] = {}
    edis = sorted(p.rglob("*.edi"))
    if not edis:
        edis = sorted(p.rglob("*.EDI"))
    for edi in edis:
        line = edi.parent.name
        groups.setdefault(line, []).append(
            str(edi)
        )
    # single flat folder: use parent name
    if len(groups) == 1:
        only_key = list(groups.keys())[0]
        if only_key == p.name:
            pass  # already correct
    if not groups and edis:
        groups["Default"] = [
            str(e) for e in edis
        ]
    return groups


def _detect_lines_to_files(
    edi_dir: str,
) -> dict[str, list[str]]:
    """Return ``{detected_line: [edi_file_paths]}`` grouped by the station-ID
    prefix (same rule as :func:`_detect_from_ids`, e.g. ``22-001`` -> ``L22``),
    but mapping to real file paths so the groups can be loaded/filtered."""
    p = Path(edi_dir)
    edis = sorted(p.rglob("*.edi")) or sorted(p.rglob("*.EDI"))
    groups: dict[str, list[str]] = {}
    for edi in edis:
        sid = edi.stem
        m = re.match(r"^([A-Za-z]*\d+)", sid)
        prefix = m.group(1) if m else sid[:4]
        if prefix.isdigit():
            prefix = f"L{prefix}"
        groups.setdefault(prefix, []).append(str(edi))
    return groups


def _build_lines_panel(
    groups: dict[str, list[str]],
    mode: str = "folder",
    editable: bool = False,
) -> list:
    if not groups:
        return [
            html.Div(
                "No EDI files found.",
                style={
                    "fontSize": "12px",
                    "color": "var(--fg-muted)",
                },
            )
        ]
    rows = []
    for i, (name, files) in enumerate(
        groups.items()
    ):
        count = len(files)
        name_el = (
            dbc.Input(
                value=name,
                id={
                    "type": "am-line-rename",
                    "index": i,
                },
                size="sm",
                style={"flex": "1"},
                debounce=True,
            )
            if editable
            else html.Span(
                name, className="am-line-name"
            )
        )
        rows.append(
            html.Div(
                [
                    html.Div(
                        className="am-line-dot",
                        style={
                            "background": _color(i)
                        },
                    ),
                    name_el,
                    html.Span(
                        f"{count} EDI",
                        className="am-line-count",
                    ),
                ],
                className="am-line-row",
                id={
                    "type": "am-line-row",
                    "index": i,
                },
            )
        )
    return [
        html.Div(
            f"{len(groups)} line(s) detected",
            style={
                "fontSize": "11px",
                "color": "var(--fg-muted)",
                "marginBottom": "8px",
            },
        ),
        *rows,
    ]


# ── Error helper ───────────────────────────────────

def _err_div(msg: str) -> html.Div:
    return html.Div(
        [
            html.I(
                className=(
                    "bi bi-exclamation-triangle"
                    " me-1"
                ),
                style={"color": "var(--tag-err)"},
            ),
            html.Span(
                msg,
                style={
                    "fontSize": "12px",
                    "color": "var(--tag-err)",
                },
            ),
        ],
        className="d-flex align-items-center",
    )


# ── Callbacks ──────────────────────────────────────

def register_edi(app) -> None:

    # Open / close EDI canvas
    @app.callback(
        Output(IDs.CANVAS_EDI, "is_open"),
        [
            Input(IDs.BTN_LOAD_EDI, "n_clicks"),
            Input(IDs.BTN_ATTACH, "n_clicks"),
        ],
        State(IDs.CANVAS_EDI, "is_open"),
        prevent_initial_call=True,
    )
    def toggle_edi_canvas(n1, n2, is_open):
        return not is_open

    # JS folder store → decode + write temp → set path
    @app.callback(
        Output(
            IDs.EDI_PATH_INPUT, "value",
            allow_duplicate=True,
        ),
        Input(IDs.FOLDER_STORE, "data"),
        prevent_initial_call=True,
    )
    def handle_folder_store(store_data):
        if not store_data:
            raise PreventUpdate
        filenames = store_data.get(
            "filenames", []
        )
        contents = store_data.get("contents", [])
        if not filenames:
            raise PreventUpdate

        tmp = tempfile.mkdtemp(
            prefix="am_edi_folder_"
        )
        for fname, content in zip(
            filenames, contents
        ):
            # fname is like "L18PLT/18-001A.edi"
            parts = fname.replace("\\", "/").split(
                "/"
            )
            dest = Path(tmp)
            for part in parts[:-1]:
                dest = dest / part
                dest.mkdir(exist_ok=True)
            out = dest / parts[-1]
            _, b64 = content.split(",", 1)
            out.write_bytes(
                base64.b64decode(b64)
            )
        return tmp

    # Mode → show/hide detect / rename buttons
    @app.callback(
        Output(IDs.BTN_DETECT_LINES, "style"),
        Output(IDs.BTN_APPLY_RENAME, "style"),
        Input(IDs.LINES_MODE, "value"),
        prevent_initial_call=True,
    )
    def toggle_mode_btns(mode):
        det = _SHOW if mode == "auto" else _HIDE
        ren = _SHOW if mode == "edit" else _HIDE
        return det, ren

    # Path input → scan folder + show lines
    @app.callback(
        Output(IDs.LINES_PANEL, "children"),
        Output(IDs.EDI_LOAD_STATUS, "children"),
        Input(IDs.EDI_PATH_INPUT, "value"),
        Input(IDs.BTN_DETECT_LINES, "n_clicks"),
        State(IDs.LINES_MODE, "value"),
        prevent_initial_call=True,
    )
    def update_lines_panel(path, _n, mode):
        if not path:
            raise PreventUpdate
        p = Path(path.strip())
        if not p.exists():
            return no_update, _err_div(
                f"Path not found: {p}"
            )
        if p.is_file():
            groups = {"Default": [str(p)]}
        else:
            groups = _folder_to_lines(str(p))
        if (
            mode == "auto"
            and ctx.triggered_id
            == IDs.BTN_DETECT_LINES
        ):
            all_ids = [
                Path(f).stem
                for files in groups.values()
                for f in files
            ]
            groups = _detect_from_ids(all_ids)
        panel = _build_lines_panel(
            groups,
            mode=mode,
            editable=(mode == "edit"),
        )
        n_edi = sum(
            len(v) for v in groups.values()
        )
        status = html.Div(
            [
                html.I(
                    className=(
                        "bi bi-check-circle me-1"
                    ),
                    style={"color": "var(--tag-ok)"},
                ),
                html.Span(
                    f"{n_edi} EDI file(s) in "
                    f"{len(groups)} line(s)",
                    style={"fontSize": "12px"},
                ),
            ],
            className="d-flex align-items-center",
        )
        return panel, status

    # Upload → decode to temp folder → set path
    @app.callback(
        Output(IDs.EDI_PATH_INPUT, "value"),
        Input(IDs.UPLOAD_EDI, "contents"),
        State(IDs.UPLOAD_EDI, "filename"),
        prevent_initial_call=True,
    )
    def handle_upload(contents, filenames):
        if not contents:
            raise PreventUpdate
        if isinstance(contents, str):
            contents = [contents]
        if isinstance(filenames, str):
            filenames = [filenames]
        filenames = filenames or []
        tmp = tempfile.mkdtemp(
            prefix="am_edi_upload_"
        )
        for content, fname in zip(
            contents, filenames
        ):
            _, b64 = content.split(",", 1)
            data = base64.b64decode(b64)
            out = Path(tmp) / fname
            out.write_bytes(data)
        return tmp

    # Confirm load → write to STORE_EDI
    @app.callback(
        Output(IDs.STORE_EDI, "data"),
        Output(IDs.EDI_BADGE, "className"),
        Output(IDs.EDI_BADGE_TEXT, "children"),
        Output(
            IDs.CANVAS_EDI, "is_open",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_LOAD_CONFIRM, "n_clicks"),
        State(IDs.EDI_PATH_INPUT, "value"),
        State(IDs.LINES_MODE, "value"),
        State({"type": "am-line-rename", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def confirm_load(n, path, mode, rename_vals):
        if not n or not path:
            raise PreventUpdate
        p = Path(path.strip())
        if not p.exists():
            raise PreventUpdate
        # Persist groups for the *mode the user picked* — not always the
        # folder layout. "auto" stores the station-ID line names (L18…),
        # "edit" applies the renamed labels, "folder" uses subfolder names.
        if p.is_file():
            groups = {"Default": [str(p)]}
        elif mode == "auto":
            groups = _detect_lines_to_files(str(p)) or _folder_to_lines(str(p))
        elif mode == "edit":
            base = _folder_to_lines(str(p))
            names = rename_vals or []
            groups = {}
            for i, (orig, files) in enumerate(base.items()):
                nm = names[i] if (i < len(names) and names[i]) else orig
                groups[str(nm)] = list(files)
        else:  # folder (default)
            groups = _folder_to_lines(str(p))
        n_edi = sum(
            len(v) for v in groups.values()
        )
        store = {
            "path": str(p),
            "groups": {
                k: list(v)
                for k, v in groups.items()
            },
            "n_edi": n_edi,
            "mode": mode,
        }
        badge_text = (
            f"{n_edi} EDI · "
            f"{len(groups)} line(s)"
        )
        return (
            store,
            "am-edi-badge visible",
            badge_text,
            False,
        )
