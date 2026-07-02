# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/data.py — modal toggle, EDI loading, station table update.

Two load triggers are handled in a single callback:
  • dcc.Upload (drag-and-drop / file picker) → decodes base64, writes tmpdir
  • BTN_LOAD_CONFIRM (server-path text input) → scans filesystem

Note: load_data is intentionally synchronous so it works even when diskcache
is not installed.  Install diskcache for background agent execution:
    pip install diskcache
"""

from __future__ import annotations

from dash import ALL, clientside_callback, ctx, dcc, html, Input, Output, State, no_update

from pycsamt.app.desktop.controllers.data_controller import DataController
from pycsamt.app.web.cache import cache_set, cache_merge_sites
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    decode_upload_to_tempdir,
    decode_folder_upload_to_tempdir,
    find_edi_files_by_line,
)

_MODE_REPLACE = "replace"
_MODE_APPEND  = "append"


def register_data(app) -> None:
    _register_modal_toggle(app)
    _register_load_mode(app)
    _register_upload_file_manager(app)
    _register_folder_filter(app)
    _register_load_btn_state(app)
    _register_drop_zone_swap(app)
    _register_preflight_preview(app)
    _register_load_progress(app)
    _register_load_data(app)
    _register_table_update(app)
    _register_station_list(app)
    _register_sta_load_bar(app)


# ──────────────────────────────────────────────────────────────────────────────

def _register_modal_toggle(app) -> None:
    @app.callback(
        Output(IDs.MODAL_LOAD,          "is_open"),
        Output("load-mode-store",       "data"),
        Output("load-mode-btn-replace", "className"),
        Output("load-mode-btn-append",  "className"),
        Output("load-modal-title",      "children"),
        Output("load-mode-hint-text",   "children"),
        Input(IDs.BTN_LOAD,             "n_clicks"),
        Input(IDs.BTN_ADD_LINE,         "n_clicks"),
        Input("modal-close-btn",        "n_clicks"),
        Input(IDs.WELCOME_CTA,          "n_clicks"),
        Input(IDs.WELCOME_CARD_LOAD,    "n_clicks"),
        Input(IDs.WELCOME_CARD_PIPE,    "n_clicks"),
        Input(IDs.WELCOME_CARD_AGENTS,  "n_clicks"),
        Input(IDs.WELCOME_CARD_VIZ,     "n_clicks"),
        State(IDs.MODAL_LOAD,           "is_open"),
        prevent_initial_call=True,
    )
    def toggle_modal(_nb, _nadd, _close, _cta, _cl, _cp, _ca, _cv, is_open):
        from dash import html as _html
        tid = ctx.triggered_id
        if tid == "modal-close-btn":
            return False, no_update, no_update, no_update, no_update, no_update
        mode = _MODE_APPEND if tid == IDs.BTN_ADD_LINE else _MODE_REPLACE
        return (
            True,
            mode,
            "load-mode-btn active" if mode == _MODE_REPLACE else "load-mode-btn",
            "load-mode-btn active" if mode == _MODE_APPEND  else "load-mode-btn",
            _modal_title(mode),
            _mode_hint(mode),
        )


def _register_load_mode(app) -> None:
    """In-modal toggle buttons update the mode store independently."""
    from dash import clientside_callback as _csc
    _csc(
        """
        function(nr, na, cur) {
            const tid = window.dash_clientside.callback_context.triggered.map(
                t => t.prop_id.split('.')[0]
            )[0];
            if (!tid) return window.dash_clientside.no_update;
            const mode = tid === 'load-mode-btn-append' ? 'append' : 'replace';
            return [
                mode,
                mode === 'replace' ? 'load-mode-btn active' : 'load-mode-btn',
                mode === 'append'  ? 'load-mode-btn active' : 'load-mode-btn',
            ];
        }
        """,
        Output("load-mode-store",       "data",      allow_duplicate=True),
        Output("load-mode-btn-replace", "className", allow_duplicate=True),
        Output("load-mode-btn-append",  "className", allow_duplicate=True),
        Input("load-mode-btn-replace",  "n_clicks"),
        Input("load-mode-btn-append",   "n_clicks"),
        State("load-mode-store",        "data"),
        prevent_initial_call=True,
    )


def _modal_title(mode: str):
    from dash import html as _html
    if mode == _MODE_APPEND:
        return [_html.I(className="bi bi-layer-forward me-2"), "Add Lines to Survey"]
    return [_html.I(className="bi bi-cloud-upload-fill me-2"), "Load Survey Data"]


def _mode_hint(mode: str) -> str:
    if mode == _MODE_APPEND:
        return ("Files will be merged into the current survey. "
                "Existing stations are kept; new ones are appended. "
                "New lines default to Active in the line manager.")
    return "Existing survey data will be replaced."


def _normalise_names(names) -> list[str]:
    if not names:
        return []
    if isinstance(names, str):
        return [names]
    return [str(n) for n in names if n]


def _sanitize_upload_name(name: str, fallback: str) -> str:
    """Return a safe display/load filename while preserving relative folders."""
    raw = str(name or "").strip().replace("\\", "/")
    if not raw:
        raw = fallback
    parts = []
    for part in raw.split("/"):
        clean = "".join(
            ch if ch.isalnum() or ch in "._- " else "_"
            for ch in part.strip()
        ).strip(" .")
        if clean:
            parts.append(clean)
    return "/".join(parts) or fallback


def _upload_entries(contents, filenames, *, source: str = "upload") -> list[dict]:
    names = _normalise_names(filenames)
    if not contents:
        return []
    if isinstance(contents, str):
        contents = [contents]
    entries = []
    for i, content in enumerate(contents):
        original = names[i] if i < len(names) else f"file_{i + 1}.edi"
        entries.append({
            "id": f"upload-{i}",
            "source": source,
            "original": original,
            "filename": _sanitize_upload_name(original, f"file_{i + 1}.edi"),
            "content": content,
        })
    return entries


def _entry_names(entries: list[dict]) -> list[str]:
    return [str(e.get("filename") or e.get("original") or "") for e in entries]


def _ext_counts(names: list[str]) -> dict[str, int]:
    counts = {"edi": 0, "avg": 0, "j": 0, "other": 0}
    for name in names:
        ext = name.rsplit(".", 1)[-1].lower() if "." in name else ""
        if ext in counts:
            counts[ext] += 1
        else:
            counts["other"] += 1
    return counts


def _infer_line_counts(names: list[str]) -> dict[str, int]:
    """Infer line names from relative paths produced by folder upload."""
    line_counts: dict[str, int] = {}
    for name in names:
        line = _line_from_path(name)
        line_counts[line] = line_counts.get(line, 0) + 1
    return line_counts


def _line_from_path(name: str) -> str:
    parts = [p for p in str(name).replace("\\", "/").split("/") if p]
    if len(parts) >= 3:
        return parts[1]
    if len(parts) == 2:
        return parts[0]
    return "flat files"


def _entry_line(entry: dict) -> str:
    name = entry.get("filename") or entry.get("original") or ""
    return _line_from_path(name)


def _filtered_upload_entries(entries, selected_lines):
    entries = list(entries or [])
    selected = {str(v) for v in (selected_lines or []) if str(v).strip()}
    if not selected:
        return entries
    return [entry for entry in entries if _entry_line(entry) in selected]


def _preflight_children(names: list[str], *, mode: str, source: str):
    if not names:
        return html.Div([
            html.I(className="bi bi-info-circle me-1"),
            "Choose a folder or drop files to preview stations and line grouping.",
        ], className="load-preflight-empty")

    counts = _ext_counts(names)
    line_counts = _infer_line_counts(names)
    shown_lines = list(line_counts.items())[:6]
    hidden = max(0, len(line_counts) - len(shown_lines))
    fmt_badges = [
        html.Span(f"{label.upper()} {count}", className="load-preflight-badge")
        for label, count in counts.items()
        if count
    ]
    mode_label = "Add lines" if mode == _MODE_APPEND else "Replace survey"

    return html.Div([
        html.Div([
            html.Div([
                html.Div(str(len(names)), className="load-preflight-stat-value"),
                html.Div("recognized files", className="load-preflight-stat-label"),
            ], className="load-preflight-stat"),
            html.Div([
                html.Div(str(len(line_counts)), className="load-preflight-stat-value"),
                html.Div("detected line groups", className="load-preflight-stat-label"),
            ], className="load-preflight-stat"),
            html.Div([
                html.Div(mode_label, className="load-preflight-stat-value small"),
                html.Div(f"{source} source", className="load-preflight-stat-label"),
            ], className="load-preflight-stat"),
        ], className="load-preflight-stats"),
        html.Div(fmt_badges, className="load-preflight-badges"),
        html.Div([
            html.Div([
                html.Span(str(line), className="load-line-name"),
                html.Span(f"{count} file{'s' if count != 1 else ''}",
                          className="load-line-count"),
            ], className="load-line-row")
            for line, count in shown_lines
        ], className="load-line-preview"),
        html.Div(
            f"+ {hidden} more line group{'s' if hidden != 1 else ''}",
            className="load-preflight-more",
            style={"display": "block" if hidden else "none"},
        ),
    ], className="load-preflight-ready")


def _detected_summary(names: list[str], *, mode: str, source: str):
    if not names:
        return [
            html.I(className="bi bi-info-circle me-1"),
            "No data selected yet.",
        ]
    counts = _ext_counts(names)
    line_counts = _infer_line_counts(names)
    fmt = " · ".join(
        f"{label.upper()} {count}"
        for label, count in counts.items()
        if count and label != "other"
    ) or f"{len(names)} files"
    mode_label = "add lines" if mode == _MODE_APPEND else "replace survey"
    source_label = "folder" if source == "folder" else "files"
    return [
        html.I(className="bi bi-check-circle-fill me-1"),
        f"{len(names)} file{'s' if len(names) != 1 else ''} · "
        f"{len(line_counts)} line group{'s' if len(line_counts) != 1 else ''} · "
        f"{fmt} · {source_label} · {mode_label}",
    ]


def _register_upload_file_manager(app) -> None:
    @app.callback(
        Output("load-upload-selection", "data"),
        Input(IDs.UPLOAD, "contents"),
        Input(IDs.FOLDER_UPLOAD_STORE, "data"),
        State(IDs.UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def capture_upload_selection(contents, folder_store, filenames):
        trigger = ctx.triggered_id
        if trigger == IDs.FOLDER_UPLOAD_STORE and folder_store:
            return _upload_entries(
                folder_store.get("contents", []),
                folder_store.get("filenames", []),
                source="folder",
            )
        if trigger == IDs.UPLOAD:
            return _upload_entries(contents, filenames, source="upload")
        return []

    @app.callback(
        Output("load-upload-selection", "data", allow_duplicate=True),
        Input({"type": "load-file-remove", "index": ALL}, "n_clicks"),
        Input({"type": "load-file-clear", "index": ALL}, "n_clicks"),
        Input({"type": "load-file-name", "index": ALL}, "value"),
        State("load-upload-selection", "data"),
        State({"type": "load-file-name", "index": ALL}, "id"),
        prevent_initial_call=True,
    )
    def edit_upload_selection(_remove_clicks, _clear_clicks,
                              name_values, entries, name_ids):
        # Guard: Dash pattern-matching callbacks fire when new matching
        # components are dynamically added to the DOM (even with
        # prevent_initial_call=True).  When the file-manager rows are first
        # rendered, the remove/clear buttons arrive with n_clicks=0 and the
        # name inputs arrive with their initial filename values — none of these
        # represent real user actions.  We skip them by checking the triggered
        # value: n_clicks=0 means "just rendered", not "clicked".
        triggered_value = ctx.triggered[0].get("value") if ctx.triggered else None
        trigger = ctx.triggered_id
        entries = list(entries or [])

        if isinstance(trigger, dict) and trigger.get("type") == "load-file-clear":
            if not triggered_value:   # n_clicks=0 → initial render, not a click
                return no_update
            return []

        if isinstance(trigger, dict) and trigger.get("type") == "load-file-remove":
            if not triggered_value:   # n_clicks=0 → initial render, not a click
                return no_update
            remove_id = trigger.get("index")
            return [e for e in entries if e.get("id") != remove_id]

        if isinstance(trigger, dict) and trigger.get("type") == "load-file-name":
            if not entries:
                return no_update
            value_by_id = {
                ident.get("index"): value
                for ident, value in zip(name_ids or [], name_values or [])
            }
            updated = []
            any_changed = False
            for i, entry in enumerate(entries):
                eid = entry.get("id") or f"upload-{i}"
                new_name = value_by_id.get(eid, entry.get("filename"))
                sanitized = _sanitize_upload_name(
                    new_name,
                    entry.get("original") or f"file_{i + 1}.edi",
                )
                updated_entry = dict(entry)
                if sanitized != entry.get("filename"):
                    updated_entry["filename"] = sanitized
                    any_changed = True
                updated.append(updated_entry)
            # Return no_update when filenames are unchanged to prevent the
            # cascading re-render loop (render → new inputs → edit → render …).
            return updated if any_changed else no_update

        return entries

    @app.callback(
        Output("load-upload-file-manager", "children"),
        Output(IDs.UPLOAD_FILELIST, "children", allow_duplicate=True),
        Input("load-upload-selection", "data"),
        prevent_initial_call=True,
    )
    def render_upload_manager(entries):
        entries = list(entries or [])
        if not entries:
            return "", ""

        rows = []
        for entry in entries:
            eid = entry.get("id")
            filename = entry.get("filename") or entry.get("original") or "file.edi"
            original = entry.get("original") or filename
            ext = filename.rsplit(".", 1)[-1].upper() if "." in filename else "FILE"
            rows.append(html.Div([
                html.Div(ext, className="load-file-ext"),
                html.Div([
                    dcc.Input(
                        id={"type": "load-file-name", "index": eid},
                        value=filename,
                        debounce=True,
                        className="load-file-name-input",
                    ),
                    html.Div(f"original: {original}", className="load-file-original"),
                ], className="load-file-main"),
                html.Button(
                    html.I(className="bi bi-trash"),
                    id={"type": "load-file-remove", "index": eid},
                    className="load-file-remove",
                    title="Remove from this import",
                    n_clicks=0,
                ),
            ], className="load-file-row"))

        manager = html.Div([
            html.Div([
                html.Span(f"{len(entries)} selected file"
                          f"{'s' if len(entries) != 1 else ''}",
                          className="load-file-manager-title"),
                html.Button(
                    [html.I(className="bi bi-x-circle me-1"), "Clear"],
                    id={"type": "load-file-clear", "index": "all"},
                    className="load-file-clear",
                    n_clicks=0,
                ),
            ], className="load-file-manager-head"),
            html.Div(rows, className="load-file-list"),
            html.Div(
                "Rename affects the import name only; source files are not modified.",
                className="load-file-note",
            ),
        ])
        return manager, f"✓ {len(entries)} file(s) ready. Rename or remove before loading."


def _register_folder_filter(app) -> None:
    @app.callback(
        Output(IDs.LOAD_FOLDER_FILTER, "options"),
        Output(IDs.LOAD_FOLDER_FILTER, "value"),
        Output(IDs.LOAD_FOLDER_FILTER_WRAP, "style"),
        Input("load-upload-selection", "data"),
        prevent_initial_call=True,
    )
    def populate(entries):
        entries = list(entries or [])
        if not entries:
            return [], [], {"display": "none"}
        names = _entry_names(entries)
        counts = _infer_line_counts(names)
        if len(counts) <= 1:
            return [], [], {"display": "none"}
        options = [
            {
                "label": f"{line} ({count} file{'s' if count != 1 else ''})",
                "value": line,
            }
            for line, count in sorted(counts.items())
        ]
        return options, [opt["value"] for opt in options], {"display": "block"}


def _register_load_btn_state(app) -> None:
    """Disable/enable the Load Survey button based on whether usable files are selected."""
    clientside_callback(
        """
        function(entries, selected) {
            var has = entries && entries.some(function(e){ return e && e.content; });
            if (has && selected && selected.length) {
                function lineOf(name) {
                    var parts = String(name || '').replace(/\\\\/g, '/')
                        .split('/').filter(Boolean);
                    if (parts.length >= 3) return parts[1];
                    if (parts.length === 2) return parts[0];
                    return 'flat files';
                }
                has = entries.some(function(e) {
                    var name = (e && (e.filename || e.original)) || '';
                    return e && e.content && selected.indexOf(lineOf(name)) >= 0;
                });
            }
            return !has;  // disabled=true when no files
        }
        """,
        Output(IDs.BTN_LOAD_CONFIRM, "disabled"),
        Input("load-upload-selection", "data"),
        Input(IDs.LOAD_FOLDER_FILTER, "value"),
    )


def _register_drop_zone_swap(app) -> None:
    """
    Swap the drag zone out for the file list once files are selected.

    When load-upload-file-manager has content → hide the drop zone and update
    the section title.  When it's cleared → restore the drop zone.
    This keeps the modal compact: the file list occupies exactly the space
    the drag box used to occupy, with no double-stacking.
    """
    clientside_callback(
        """
        function(manager_children) {
            var has_files = (
                manager_children !== null &&
                manager_children !== undefined &&
                manager_children !== '' &&
                !(Array.isArray(manager_children) && manager_children.length === 0)
            );
            var drop_style  = has_files ? {display: 'none'} : {display: 'block'};
            var title_text  = has_files ? 'Selected files' : 'Or drag & drop files / folders';
            return [drop_style, title_text];
        }
        """,
        Output("upload-drop-wrap",     "style"),
        Output("upload-section-title", "children"),
        Input("load-upload-file-manager", "children"),
    )


def _register_preflight_preview(app) -> None:
    @app.callback(
        Output("load-preflight-preview", "children"),
        Output("load-detected-summary", "children"),
        Output("load-source-selection", "data"),
        Input("load-upload-selection", "data"),
        Input(IDs.LOAD_FOLDER_FILTER, "value"),
        Input("load-mode-store", "data"),
        State("load-source-selection", "data"),
    )
    def update_preflight(upload_entries, selected_lines, mode, current_source):
        source = current_source or "none"
        if upload_entries:
            filtered = _filtered_upload_entries(upload_entries, selected_lines)
            source = str(upload_entries[0].get("source") or "upload")
            names = _entry_names(filtered)
            source_label = "folder" if source == "folder" else "file"
            return (
                _preflight_children(
                    names,
                    mode=mode or _MODE_REPLACE,
                    source=source_label,
                ),
                _detected_summary(
                    names,
                    mode=mode or _MODE_REPLACE,
                    source=source,
                ),
                source,
            )
        return (
            _preflight_children([], mode=mode or _MODE_REPLACE, source="none"),
            _detected_summary([], mode=mode or _MODE_REPLACE, source="none"),
            "none",
        )


def _register_load_progress(app) -> None:
    """Clientside show/hide of the animated progress bar during EDI loading."""

    # Show the progress bar as soon as Load Survey is clicked.
    clientside_callback(
        """
        function(n_clicks, entries, selected) {
            if (!n_clicks) return [window.dash_clientside.no_update,
                                   window.dash_clientside.no_update,
                                   window.dash_clientside.no_update];

            function lineOf(name) {
                var parts = String(name || '').replace(/\\\\/g, '/')
                    .split('/').filter(Boolean);
                if (parts.length >= 3) return parts[1];
                if (parts.length === 2) return parts[0];
                return 'flat files';
            }
            var filtered = entries || [];
            if (selected && selected.length) {
                filtered = filtered.filter(function(e) {
                    var name = (e && (e.filename || e.original)) || '';
                    return selected.indexOf(lineOf(name)) >= 0;
                });
            }
            var n = filtered.filter(function(e){ return e && e.content; }).length;
            var label = n > 0 ? ('Parsing ' + n + ' file' + (n !== 1 ? 's' : '') + '…') : 'Processing stations…';

            // Restart the CSS fill animation by cloning the fill node
            requestAnimationFrame(function() {
                var fill = document.getElementById('load-prog-fill');
                if (fill) {
                    var clone = fill.cloneNode(true);
                    clone.classList.remove('load-prog-done');
                    fill.parentNode.replaceChild(clone, fill);
                }
            });

            return [{'display': 'block'}, label, ''];
        }
        """,
        Output("load-progress-wrap",   "style"),
        Output("load-prog-label",      "children"),
        Output("load-prog-sublabel",   "children"),
        Input(IDs.BTN_LOAD_CONFIRM,    "n_clicks"),
        State("load-upload-selection", "data"),
        State(IDs.LOAD_FOLDER_FILTER,  "value"),
        prevent_initial_call=True,
    )

    # When LOAD_FEEDBACK receives a result: complete the progress bar animation
    # and fade it out after 900 ms via direct DOM manipulation.
    clientside_callback(
        """
        function(feedback) {
            if (!feedback) return window.dash_clientside.no_update;
            var fill = document.getElementById('load-prog-fill');
            if (fill) fill.classList.add('load-prog-done');
            setTimeout(function() {
                var wrap = document.getElementById('load-progress-wrap');
                if (wrap) wrap.style.display = 'none';
            }, 900);
            return window.dash_clientside.no_update;
        }
        """,
        Output("load-progress-wrap",   "style",    allow_duplicate=True),
        Input(IDs.LOAD_FEEDBACK,       "children"),
        prevent_initial_call=True,
    )


def _register_load_data(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA,           "data"),
        Output(IDs.LOAD_FEEDBACK,        "children"),
        Output(IDs.MODAL_LOAD,           "is_open",  allow_duplicate=True),
        Output(IDs.UPLOAD_FILELIST,      "children"),
        Output("load-upload-selection",  "data",     allow_duplicate=True),
        Input(IDs.BTN_LOAD_CONFIRM,      "n_clicks"),
        State("load-upload-selection",   "data"),
        State(IDs.LOAD_FOLDER_FILTER,    "value"),
        State(IDs.SESSION_ID,            "data"),
        State("load-mode-store",         "data"),
        State(IDs.STORE_DATA,            "data"),
        prevent_initial_call=True,
    )
    def load_data(_n_clicks, upload_entries, selected_lines,
                  session_id, load_mode,
                  existing_store):
        _SKIP = (no_update, no_update, no_update, no_update, no_update)
        mode  = load_mode or _MODE_REPLACE

        if not session_id:
            return no_update, "⚠ Session not initialised — refresh the page.", no_update, no_update, no_update

        # ── Derive source directly from entries — avoids race with the
        #    load-source-selection store (whose value may lag behind the
        #    upload-selection store when the user clicks quickly).
        upload_entries = _filtered_upload_entries(upload_entries, selected_lines)
        usable = [e for e in (upload_entries or []) if e.get("content")]
        if usable:
            selected_source = str(usable[0].get("source") or "upload")
        else:
            selected_source = "none"

        if selected_source == "none":
            return (no_update,
                    "⚠ Choose a survey folder or drop files first.",
                    no_update, no_update, no_update)

        import shutil

        if selected_source == "folder":
            filenames = [e.get("filename") for e in usable]
            contents  = [e.get("content")  for e in usable]

            if not filenames:
                return no_update, "⚠ No EDI / AVG / J files found in selected folder.", no_update, no_update, no_update

            paths, tmpdir, line_map = decode_folder_upload_to_tempdir(contents, filenames)
            if not paths:
                return no_update, "⚠ No recognised files could be decoded.", no_update, no_update, no_update

            path_to_line = {p: line for line, plist in line_map.items() for p in plist}
            n_lines      = len(line_map)

            try:
                ctrl        = DataController()
                sites       = ctrl.load(paths, path_to_line=path_to_line)
                df          = ctrl.dataframe
                new_records = df.to_dict("records")
                line_counts = {line: len(plist) for line, plist in line_map.items()}

                if mode == _MODE_APPEND and existing_store:
                    store, sites = _merge_store(
                        existing_store, new_records, sites, session_id,
                        new_dir="[browsed]", n_new_lines=n_lines,
                        new_line_counts=line_counts,
                    )
                    n_added  = len(new_records)
                    n_total  = store["n_stations"]
                    feedback = (
                        f"✓ Added {n_added} station{'s' if n_added != 1 else ''} "
                        f"from {n_lines} line{'s' if n_lines != 1 else ''} "
                        f"({n_total} stations total)."
                    )
                else:
                    cache_set(session_id, sites)
                    shown        = list(line_counts.items())[:5]
                    line_summary = ", ".join(f"{ln}: {cnt}" for ln, cnt in shown)
                    if n_lines > 5:
                        line_summary += f" … +{n_lines - 5} more"
                    store = {
                        "station_records": new_records,
                        "data_dir":        "[browsed]",
                        "n_stations":      len(new_records),
                        "n_lines":         n_lines,
                        "line_counts":     line_counts,
                    }
                    n_s = len(new_records)
                    feedback = (
                        f"✓ Loaded {n_s} station{'s' if n_s != 1 else ''} "
                        f"from {n_lines} line{'s' if n_lines != 1 else ''}"
                        + (f" ({line_summary})" if n_lines > 1 else "")
                        + "."
                    )

                # Close modal and clear selection so next open starts fresh
                return store, feedback, False, no_update, []
            except Exception as exc:
                return no_update, f"✗ Error: {exc}", no_update, no_update, no_update
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        # Individual file upload branch
        upload_contents  = [e.get("content")  for e in usable]
        upload_filenames = [e.get("filename") for e in usable]
        paths, tmpdir = decode_upload_to_tempdir(upload_contents, upload_filenames)
        if not paths:
            return (no_update,
                    "⚠ No recognised files in upload (.edi / .avg / .j).",
                    no_update, no_update, no_update)

        names   = upload_filenames if isinstance(upload_filenames, list) else [upload_filenames]
        preview = ", ".join(names[:5]) + (f" … and {len(names) - 5} more" if len(names) > 5 else "")

        try:
            ctrl        = DataController()
            sites       = ctrl.load(paths)
            df          = ctrl.dataframe
            new_records = df.to_dict("records")

            if mode == _MODE_APPEND and existing_store:
                store, sites = _merge_store(
                    existing_store, new_records, sites, session_id,
                    new_dir="[uploaded]", n_new_lines=1,
                    new_line_counts={"uploaded": len(new_records)},
                )
                verb = f"✓ Added {len(new_records)} station(s) to survey ({store['n_stations']} total)."
            else:
                cache_set(session_id, sites)
                store = {
                    "station_records": new_records,
                    "data_dir":        "[uploaded]",
                    "n_stations":      len(new_records),
                    "n_lines":         1,
                    "line_counts":     {"uploaded": len(new_records)},
                }
                verb = f"✓ Loaded {len(new_records)} station(s) from {len(paths)} file(s)."

            # Close modal and clear selection
            return store, verb, False, f"✓ {len(paths)} file(s): {preview}", []
        except Exception as exc:
            return no_update, f"✗ Error loading files: {exc}", no_update, "", no_update
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


def _merge_store(existing_store, new_records, new_sites, session_id,
                 new_dir, n_new_lines, new_line_counts):
    """
    Merge *new_records* and *new_sites* into the existing survey store.

    Returns (merged_store_dict, merged_sites_object).
    """
    from collections import Counter

    existing_records = existing_store.get("station_records", [])
    existing_dir     = existing_store.get("data_dir", "")

    # Deduplicate by station ID — new data wins on collision
    by_id = {r["ID"]: r for r in existing_records}
    for r in new_records:
        by_id[r["ID"]] = r
    merged_records = list(by_id.values())

    # Recompute line counts from merged records
    merged_line_counts: dict = dict(Counter(
        r.get("Line", "") for r in merged_records if r.get("Line")
    ))
    merged_n_lines = len(merged_line_counts)

    merged_dir = existing_dir if existing_dir == new_dir else f"{existing_dir} + {new_dir}"

    merged_store = {
        "station_records": merged_records,
        "data_dir":        merged_dir,
        "n_stations":      len(merged_records),
        "n_lines":         merged_n_lines,
        "line_counts":     merged_line_counts,
    }

    merged_sites = cache_merge_sites(session_id, new_sites)
    return merged_store, merged_sites


def _register_table_update(app) -> None:
    @app.callback(
        Output(IDs.STATION_TABLE,   "data"),
        Output(IDs.STATION_SUMMARY, "children"),
        Output(IDs.STATUS_LABEL,    "children"),
        Output(IDs.KPI_STATIONS,    "children"),
        Output(IDs.KPI_FREQ,        "children"),
        Output(IDs.KPI_PROFILES,    "children"),
        Output(IDs.KPI_SURVEY,      "children"),
        Input(IDs.STORE_DATA,       "data"),
    )
    def update_table(store_data):
        if not store_data:
            return [], "No stations loaded", "No data loaded", "—", "—", "—", "—"

        import os
        import statistics

        records = store_data.get("station_records", [])
        n       = store_data.get("n_stations", 0)
        label   = f"{n} station{'s' if n != 1 else ''} loaded"

        kpi_stations = str(n) if n else "—"

        freq_vals = [r.get("N_freq", 0) for r in records if r.get("N_freq")]
        kpi_freq  = str(int(statistics.median(freq_vals))) if freq_vals else "—"

        # Use stored n_lines when available; fall back to distinct lat clusters
        n_lines = store_data.get("n_lines", 0)
        if n_lines:
            kpi_profiles = str(n_lines)
        else:
            lats = sorted({round(r.get("Latitude", 0), 3) for r in records if r.get("Latitude")})
            kpi_profiles = str(len(lats)) if lats else "—"

        data_dir = store_data.get("data_dir", "")
        if data_dir == "[uploaded]":
            kpi_survey = "upload"
        elif data_dir:
            kpi_survey = os.path.basename(data_dir.rstrip("/\\")) or data_dir[:10]
        else:
            kpi_survey = "—"

        return records, label, f"✓ {label}", kpi_stations, kpi_freq, kpi_profiles, kpi_survey


# ── Station list (visual two-line rows + line filter pills) ───────────────────

_PALETTE_DARK = [
    "#89b4fa", "#a6e3a1", "#fab387", "#cba6f7",
    "#94e2d5", "#f5c2e7", "#f9e2af", "#f38ba8", "#89dceb", "#b4befe",
]
_PALETTE_LIGHT = [
    "#1e66f5", "#40a02b", "#fe640b", "#8839ef",
    "#179299", "#ea76cb", "#df8e1d", "#d20f39", "#04a5e5", "#7287fd",
]


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _register_station_list(app) -> None:

    # ── 1. render pills + rows ────────────────────────────────────────────────
    @app.callback(
        Output(IDs.STATION_LINE_PILLS, "children"),
        Output(IDs.STATION_LIST_BODY,  "children"),
        Input(IDs.STORE_DATA,          "data"),
        Input(IDs.STORE_ACTIVE_LINES,  "data"),
        Input(IDs.STORE_SELECTION,     "data"),
        Input(IDs.STATION_LINE_FILTER, "data"),
        Input(IDs.STORE_THEME,         "data"),
    )
    def render_station_list(store_data, active_lines, selection, line_filter, theme):
        dark    = (theme or "dark") == "dark"
        palette = _PALETTE_DARK if dark else _PALETTE_LIGHT

        _empty = html.Div(
            "No stations loaded",
            style={"color": "var(--subtext0)", "fontSize": "11px",
                   "padding": "14px 12px"},
        )
        if not store_data:
            return [], [_empty]

        records = store_data.get("station_records", [])
        if not records:
            return [], [_empty]

        active_set = set((active_lines or {}).get("active", []))
        all_lines  = (active_lines or {}).get("all", [])

        # line → colour index (stable order from all_lines)
        line_col: dict[str, str] = {
            ln: palette[i % len(palette)]
            for i, ln in enumerate(all_lines)
        }

        # ── pills ────────────────────────────────────────────────────────────
        all_cls = "sta-line-pill sta-pill-all" + (
            " pill-filtered" if line_filter is None else "")
        pills: list = [
            html.Span(
                "All",
                id={"type": "sta-line-pill", "index": "__all__"},
                className=all_cls,
                n_clicks=0,
            )
        ]
        for ln in all_lines:
            col        = line_col.get(ln, palette[0])
            r, g, b    = _hex_to_rgb(col)
            is_active  = ln in active_set
            is_focused = line_filter == ln
            cls = (
                "sta-line-pill"
                + ("" if is_active else " pill-inactive")
                + (" pill-filtered" if is_focused else "")
            )
            pills.append(html.Span(
                ["● ", ln],
                id={"type": "sta-line-pill", "index": ln},
                className=cls,
                n_clicks=0,
                style={
                    "color":       col,
                    "background":  f"rgba({r},{g},{b},0.13)",
                    "borderColor": f"rgba({r},{g},{b},0.40)",
                },
            ))

        # ── station rows ─────────────────────────────────────────────────────
        selected_id = (selection or {}).get("station_id")
        rows: list  = []

        for rec in records:
            sid  = rec.get("ID",   "")
            line = rec.get("Line", "")

            if line_filter and line_filter != "__all__" and line != line_filter:
                continue

            col     = line_col.get(line, "#6c7086")
            r, g, b = _hex_to_rgb(col) if col.startswith("#") else (108, 112, 134)
            is_sel  = sid == selected_id

            # metadata: use explicit ± format for coords, comma-sep elev
            parts: list[str] = []
            lat  = rec.get("Latitude")
            lon  = rec.get("Longitude")
            elev = rec.get("Elevation")
            nf   = rec.get("N_freq")
            try:
                if lat  is not None: parts.append(f"{float(lat):+.4f}°")
                if lon  is not None: parts.append(f"{float(lon):+.4f}°")
                if elev is not None: parts.append(f"{int(elev):,} m")
                if nf   is not None: parts.append(f"{int(nf)} f")
            except (TypeError, ValueError):
                pass
            meta = " · ".join(parts) or "—"

            has_tipper = bool(rec.get("Tipper", False))

            rows.append(html.Div(
                [
                    html.Div(
                        [
                            html.Span(
                                line or "—",
                                className="sta-line-badge",
                                style={
                                    "color":      col,
                                    "background": f"rgba({r},{g},{b},0.15)",
                                },
                            ),
                            html.Span(sid, className="sta-row-id"),
                            html.Span(
                                className="sta-tip-dot",
                                title="Has tipper",
                                style={"display": "inline-block"
                                       if has_tipper else "none"},
                            ),
                        ],
                        className="sta-row-top",
                    ),
                    html.Div(meta, className="sta-row-meta"),
                ],
                id={"type": "sta-row", "index": sid},
                className="sta-row row-selected" if is_sel else "sta-row",
                n_clicks=0,
                style={"borderLeftColor": col},
            ))

        if not rows:
            rows = [html.Div(
                "No stations match the current filter",
                style={"color": "var(--subtext0)", "fontSize": "11px",
                       "padding": "12px"},
            )]

        return pills, rows

    # ── 2. pill click → toggle local line filter ──────────────────────────────
    @app.callback(
        Output(IDs.STATION_LINE_FILTER, "data"),
        Input({"type": "sta-line-pill", "index": ALL}, "n_clicks"),
        State(IDs.STATION_LINE_FILTER, "data"),
        prevent_initial_call=True,
    )
    def toggle_line_filter(_, current):
        tid = ctx.triggered_id
        if tid is None:
            return no_update
        line = tid.get("index") if isinstance(tid, dict) else None
        if not line or line == "__all__":
            return None
        return None if current == line else line


def _register_sta_load_bar(app) -> None:
    """
    Clientside-only loading strip inside the station panel.

    Shows an infinite animated gradient bar below the panel header as soon as
    STORE_DATA changes (the server hasn't started building the list yet), and
    hides it the moment STATION_LIST_BODY receives its updated children.
    """

    # Activate bar + lemniscate spinner when new data lands in STORE_DATA
    clientside_callback(
        """
        function(store_data) {
            if (!store_data || !store_data.n_stations) {
                return ['sta-load-bar sta-load-bar--hidden',
                        'sta-loading-spinner sta-loading-spinner--hidden'];
            }
            var n = store_data.n_stations || 0;
            var msg = document.getElementById('sta-loader-msg');
            if (msg) msg.textContent = 'Building ' + n + ' station' + (n !== 1 ? 's' : '') + '…';
            return ['sta-load-bar',
                    'sta-loading-spinner'];
        }
        """,
        Output("sta-load-bar",          "className"),
        Output("sta-loading-spinner",   "className"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )

    # Deactivate both once the list body has been populated
    clientside_callback(
        """
        function(children) {
            return ['sta-load-bar sta-load-bar--hidden',
                    'sta-loading-spinner sta-loading-spinner--hidden'];
        }
        """,
        Output("sta-load-bar",          "className", allow_duplicate=True),
        Output("sta-loading-spinner",   "className", allow_duplicate=True),
        Input(IDs.STATION_LIST_BODY, "children"),
        prevent_initial_call=True,
    )
