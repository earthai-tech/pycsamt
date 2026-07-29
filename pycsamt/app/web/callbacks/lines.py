# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/lines.py — survey-line management panel.

Three detection modes:
  folder     — folder/subfolder name is the line name  (default)
  auto       — prefix before first '-'/'_' in station ID; numeric → L-prefix
  edit       — user types new names; Apply commits them to station_records
"""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import (
    ALL,
    Input,
    Output,
    State,
    ctx,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.web.layout import IDs

# ── Mode → hint text ─────────────────────────────────────────────────────────

_MODE_HINTS = {
    "folder": (
        "Lines come from subfolder names. Works automatically with "
        "multi-folder uploads — each subfolder = one line."
    ),
    "auto": (
        "Parses the prefix before the first '-' or '_' in each station ID. "
        "Purely numeric prefixes get an 'L' prepended "
        "(e.g. '18-001A' → L18, '22-345A' → L22). "
        "Click Detect Lines to apply."
    ),
    "edit": (
        "Type new names in the inputs that appear next to each line, "
        "then click Apply Renames to commit."
    ),
}


def register_lines(app) -> None:
    _register_init(app)
    _register_open(app)
    _register_mode_ui(app)
    _register_auto_detect(app)
    _register_apply_renames(app)
    _register_panel(app)
    _register_toggle(app)
    _register_all_none(app)


# ── 1. Initialise / update STORE_ACTIVE_LINES when data loads ────────────────


def _register_init(app) -> None:
    @app.callback(
        Output(IDs.STORE_ACTIVE_LINES, "data"),
        Input(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State("load-mode-store", "data"),
        prevent_initial_call=True,
    )
    def init_active_lines(store_data, prev_active_store, load_mode):
        if not store_data:
            return {"active": [], "all": []}

        records = store_data.get("station_records", [])
        lines: list[str] = []
        seen: set = set()
        for r in records:
            ln = r.get("Line", "")
            if ln and ln != "—" and ln not in seen:
                lines.append(str(ln))
                seen.add(ln)

        # Append mode: preserve existing active/inactive choices.
        # Only newly added lines are defaulted to active.
        if load_mode == "append" and prev_active_store and prev_active_store.get("all"):
            prev_active = set(prev_active_store.get("active", []))
            prev_all = set(prev_active_store.get("all", []))
            new_lines = [ln for ln in lines if ln not in prev_all]
            active = [ln for ln in lines if ln in prev_active or ln in new_lines]
            return {"active": active, "all": lines}

        # Replace mode (or first load): all lines default to active
        return {"active": lines, "all": lines}


# ── 2. Toggle offcanvas ───────────────────────────────────────────────────────


def _register_open(app) -> None:
    @app.callback(
        Output(IDs.LINES_OFFCANVAS, "is_open"),
        Input(IDs.BTN_LINES_OPEN, "n_clicks"),
        State(IDs.LINES_OFFCANVAS, "is_open"),
        prevent_initial_call=True,
    )
    def toggle_lines_offcanvas(n_clicks, is_open):
        if n_clicks:
            return not is_open
        return no_update


# ── 3. Mode → button visibility + hint text ───────────────────────────────────


def _register_mode_ui(app) -> None:
    @app.callback(
        Output(IDs.BTN_LINES_DETECT, "style"),
        Output(IDs.BTN_LINES_RENAME_APPLY, "style"),
        Output("lines-mode-hint", "children"),
        Input(IDs.LINES_DETECT_MODE, "value"),
    )
    def sync_mode_ui(mode):
        mode = mode or "folder"
        show = {"display": "block"}
        hide = {"display": "none"}
        return (
            show if mode == "auto" else hide,
            show if mode == "edit" else hide,
            _MODE_HINTS.get(mode, ""),
        )


# ── 4. Auto-detect: parse station IDs → update STORE_DATA ─────────────────────


def _register_auto_detect(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Input(IDs.BTN_LINES_DETECT, "n_clicks"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def detect_lines_auto(n_clicks, store_data):
        if not n_clicks or not store_data:
            raise PreventUpdate

        from pycsamt.site.lines import (
            detect_lines_from_station_ids,
        )

        records = store_data.get("station_records", [])
        if not records:
            raise PreventUpdate

        station_ids = [r.get("ID", "") for r in records if r.get("ID")]
        groups = detect_lines_from_station_ids(station_ids)

        # Build station_id → new_line map
        id_to_line: dict[str, str] = {}
        for line_name, ids in groups.items():
            for sid in ids:
                id_to_line[sid] = line_name

        # Update records in-place
        updated = []
        for rec in records:
            new_rec = dict(rec)
            sid = new_rec.get("ID", "")
            if sid in id_to_line:
                new_rec["Line"] = id_to_line[sid]
            updated.append(new_rec)

        # Rebuild line_counts summary
        line_counts: dict[str, int] = {}
        for rec in updated:
            ln = rec.get("Line", "")
            if ln and ln != "—":
                line_counts[ln] = line_counts.get(ln, 0) + 1

        new_store = dict(store_data)
        new_store["station_records"] = updated
        new_store["n_lines"] = len(line_counts)
        new_store["line_counts"] = line_counts
        return new_store


# ── 5. Apply renames: text inputs → update STORE_DATA ─────────────────────────


def _register_apply_renames(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Input(IDs.BTN_LINES_RENAME_APPLY, "n_clicks"),
        State({"type": "line-rename-input", "index": ALL}, "value"),
        State({"type": "line-rename-input", "index": ALL}, "id"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def apply_renames(n_clicks, new_names, input_ids, store_data):
        if not n_clicks or not store_data or not input_ids:
            raise PreventUpdate

        from pycsamt.site.lines import apply_line_renames

        # Build rename map: original_key → new_name
        renames: dict[str, str] = {}
        for item_id, new_name in zip(input_ids, new_names):
            original = item_id.get("index", "")
            name = (new_name or "").strip()
            if name and name != original:
                renames[original] = name

        if not renames:
            raise PreventUpdate

        records = store_data.get("station_records", [])
        updated = apply_line_renames(records, renames)

        # Rebuild line_counts
        line_counts: dict[str, int] = {}
        for rec in updated:
            ln = rec.get("Line", "")
            if ln and ln != "—":
                line_counts[ln] = line_counts.get(ln, 0) + 1

        new_store = dict(store_data)
        new_store["station_records"] = updated
        new_store["n_lines"] = len(line_counts)
        new_store["line_counts"] = line_counts
        return new_store


# ── 6. Rebuild panel body ─────────────────────────────────────────────────────


def _register_panel(app) -> None:
    @app.callback(
        Output(IDs.LINES_PANEL_BODY, "children"),
        Input(IDs.LINES_OFFCANVAS, "is_open"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.LINES_DETECT_MODE, "value"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def build_lines_panel(is_open, store_data, mode, active_store, theme):
        if not is_open:
            return no_update

        if not store_data:
            return html.P(
                "Load survey data first.",
                className="text-muted small px-2 pt-2",
            )

        records = store_data.get("station_records", [])
        line_counts: dict[str, int] = {}
        for r in records:
            ln = r.get("Line", "")
            if ln and ln != "—":
                line_counts[str(ln)] = line_counts.get(str(ln), 0) + 1

        if not line_counts:
            return _empty_hint(mode)

        active_set = set((active_store or {}).get("active", list(line_counts.keys())))

        try:
            from pycsamt.app.web.utils import _LINE_COLORS

            palette = _LINE_COLORS
        except ImportError:
            palette = [
                "#89b4fa",
                "#a6e3a1",
                "#fab387",
                "#f38ba8",
                "#cba6f7",
                "#89dceb",
                "#f9e2af",
                "#b4befe",
            ]

        edit_mode = (mode or "folder") == "edit"
        rows = []
        for idx, (line_name, count) in enumerate(line_counts.items()):
            color = palette[idx % len(palette)]
            is_active = line_name in active_set
            rows.append(_build_row(line_name, count, color, is_active, edit_mode))

        children = [html.Div(rows, className="lines-list")]

        if edit_mode:
            children.append(
                html.Div(
                    f"{len(line_counts)} line{'s' if len(line_counts) != 1 else ''}  "
                    "— edit names above, then click Apply Renames.",
                    className="text-muted",
                    style={
                        "fontSize": "10px",
                        "marginTop": "6px",
                        "padding": "0 4px",
                    },
                )
            )
        elif (mode or "folder") == "auto":
            n = len(line_counts)
            children.append(
                html.Div(
                    f"{n} line{'s' if n != 1 else ''} detected from station ID prefixes.",
                    className="text-muted",
                    style={
                        "fontSize": "10px",
                        "marginTop": "6px",
                        "padding": "0 4px",
                    },
                )
            )

        return children


def _empty_hint(mode: str) -> html.P:
    msgs = {
        "folder": "No named lines found. Try 'Auto-detect from ID' if you uploaded "
        "all EDIs into a single folder.",
        "auto": "No station records loaded yet. Load data, then click Detect Lines.",
        "edit": "No lines to rename. Load data first.",
    }
    return html.P(
        msgs.get(mode or "folder", "No lines found."),
        className="text-muted small px-2 pt-2",
    )


def _build_row(
    line_name: str,
    count: int,
    color: str,
    is_active: bool,
    edit_mode: bool,
) -> html.Div:
    dot = html.Div(
        style={
            "width": "10px",
            "height": "10px",
            "borderRadius": "50%",
            "backgroundColor": color,
            "flexShrink": "0",
        },
        className="lines-color-dot",
    )

    toggle = dbc.Switch(
        id={"type": "line-switch", "index": line_name},
        value=is_active,
        className="lines-toggle ms-auto",
        label="",
        style={"marginBottom": "0"},
    )
    count_badge = html.Span(f"{count} sta.", className="lines-count")

    if edit_mode:
        name_widget = dbc.Input(
            id={"type": "line-rename-input", "index": line_name},
            value=line_name,
            size="sm",
            debounce=False,
            className="lines-rename-input",
            style={
                "maxWidth": "130px",
                "fontSize": "12px",
                "height": "26px",
                "padding": "2px 6px",
            },
        )
    else:
        name_widget = html.Span(line_name, className="lines-name")

    return html.Div(
        [dot, name_widget, count_badge, toggle],
        className="lines-row",
    )


# ── 7. Sync STORE_ACTIVE_LINES when any switch is toggled ─────────────────────


def _register_toggle(app) -> None:
    @app.callback(
        Output(IDs.STORE_ACTIVE_LINES, "data", allow_duplicate=True),
        Input({"type": "line-switch", "index": ALL}, "value"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def sync_active_lines(switch_vals, active_store):
        triggered = ctx.inputs_list[0]
        if not triggered:
            return no_update

        all_lines = (active_store or {}).get("all", [])
        active: list[str] = []
        for item in triggered:
            if item.get("value"):
                active.append(item["id"]["index"])

        return {"active": active, "all": all_lines}


# ── 8. All / None quick-select buttons ────────────────────────────────────────


def _register_all_none(app) -> None:
    @app.callback(
        Output(
            {"type": "line-switch", "index": ALL},
            "value",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_LINES_ALL, "n_clicks"),
        Input(IDs.BTN_LINES_NONE, "n_clicks"),
        State({"type": "line-switch", "index": ALL}, "id"),
        prevent_initial_call=True,
    )
    def all_or_none(n_all, n_none, switch_ids):
        trigger = ctx.triggered_id
        if not switch_ids:
            return no_update
        n = len(switch_ids)
        if trigger == IDs.BTN_LINES_ALL:
            return [True] * n
        if trigger == IDs.BTN_LINES_NONE:
            return [False] * n
        return no_update
