# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Multi-line EDI loading — modeled on ``app/web`` with replace/append,
spinner, file manager, preflight, and an animated progress bar.

The decoded ``{line: [paths]}`` grouping feeds
:meth:`pycsamt.map.MapView.from_lines`; the heavy view is cached
server-side while only station records reach the dcc store.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from dash import (
    ALL, ctx, dcc, html, no_update,
    clientside_callback, Input, Output, State,
)

from .._ids import IDs
from ..cache import get_view, set_view
from .._render import merge_views, store_from_view

_MODE_REPLACE = "replace"
_MODE_APPEND = "append"
_EXTS = {".edi", ".avg", ".j"}


def register_load(app) -> None:
    _register_modal_toggle(app)
    _register_mode_buttons(app)
    _register_capture_selection(app)
    _register_edit_selection(app)
    _register_file_manager(app)
    _register_drop_swap(app)
    _register_preflight(app)
    _register_btn_state(app)
    _register_progress(app)
    _register_confirm_load(app)


# ── modal open / mode ──────────────────────────────────


def _register_modal_toggle(app) -> None:
    @app.callback(
        Output(IDs.MODAL_LOAD, "is_open"),
        Output(IDs.LOAD_MODE_STORE, "data"),
        Output(IDs.MODE_BTN_REPLACE, "className"),
        Output(IDs.MODE_BTN_APPEND, "className"),
        Output(IDs.MODAL_TITLE, "children"),
        Output(IDs.MODE_HINT, "children"),
        Input(IDs.BTN_LOAD, "n_clicks"),
        Input(IDs.BTN_ADD_LINE, "n_clicks"),
        prevent_initial_call=True,
    )
    def open_modal(_n_load, _n_add):
        append = ctx.triggered_id == IDs.BTN_ADD_LINE
        mode = _MODE_APPEND if append else _MODE_REPLACE
        return (
            True,
            mode,
            "mv-mode-btn" if append else "mv-mode-btn active",
            "mv-mode-btn active" if append else "mv-mode-btn",
            _modal_title(mode),
            _mode_hint(mode),
        )


def _register_mode_buttons(app) -> None:
    clientside_callback(
        """
        function(nr, na, cur) {
            const trig = window.dash_clientside.callback_context.triggered;
            if (!trig.length) return window.dash_clientside.no_update;
            const tid = trig[0].prop_id.split('.')[0];
            const mode = tid.indexOf('append') >= 0 ? 'append' : 'replace';
            return [
                mode,
                mode === 'replace' ? 'mv-mode-btn active' : 'mv-mode-btn',
                mode === 'append'  ? 'mv-mode-btn active' : 'mv-mode-btn',
            ];
        }
        """,
        Output(IDs.LOAD_MODE_STORE, "data", allow_duplicate=True),
        Output(IDs.MODE_BTN_REPLACE, "className", allow_duplicate=True),
        Output(IDs.MODE_BTN_APPEND, "className", allow_duplicate=True),
        Input(IDs.MODE_BTN_REPLACE, "n_clicks"),
        Input(IDs.MODE_BTN_APPEND, "n_clicks"),
        State(IDs.LOAD_MODE_STORE, "data"),
        prevent_initial_call=True,
    )


def _modal_title(mode):
    if mode == _MODE_APPEND:
        return [html.I(className="bi bi-layer-forward me-2"),
                "Add Lines to Survey"]
    return [html.I(className="bi bi-cloud-upload-fill me-2"),
            "Load Survey Lines"]


def _mode_hint(mode):
    if mode == _MODE_APPEND:
        return ("New lines are merged into the current survey; existing "
                "stations are kept and new ones appended.")
    return "Existing survey data will be replaced."


# ── selection capture ──────────────────────────────────


def _sanitize(name, fallback):
    raw = str(name or "").strip().replace("\\", "/")
    if not raw:
        return fallback
    parts = []
    for part in raw.split("/"):
        clean = "".join(
            ch if ch.isalnum() or ch in "._- " else "_"
            for ch in part.strip()
        ).strip(" .")
        if clean:
            parts.append(clean)
    return "/".join(parts) or fallback


def _entries(contents, filenames, *, source):
    if not contents:
        return []
    if isinstance(contents, str):
        contents = [contents]
    if isinstance(filenames, str):
        filenames = [filenames]
    filenames = filenames or []
    out = []
    for i, content in enumerate(contents):
        original = filenames[i] if i < len(filenames) else f"file_{i+1}.edi"
        if Path(original).suffix.lower() not in _EXTS:
            continue
        out.append({
            "id": f"f-{i}",
            "source": source,
            "original": original,
            "filename": _sanitize(original, f"file_{i+1}.edi"),
            "content": content,
        })
    return out


def _register_capture_selection(app) -> None:
    @app.callback(
        Output(IDs.UPLOAD_SELECTION, "data"),
        Input(IDs.UPLOAD, "contents"),
        Input(IDs.FOLDER_STORE, "data"),
        State(IDs.UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def capture(contents, folder_store, filenames):
        trig = ctx.triggered_id
        if trig == IDs.FOLDER_STORE and folder_store:
            return _entries(
                folder_store.get("contents", []),
                folder_store.get("filenames", []),
                source="folder",
            )
        if trig == IDs.UPLOAD:
            return _entries(contents, filenames, source="upload")
        return []


def _register_edit_selection(app) -> None:
    @app.callback(
        Output(IDs.UPLOAD_SELECTION, "data", allow_duplicate=True),
        Input({"type": "mv-file-remove", "index": ALL}, "n_clicks"),
        Input({"type": "mv-file-clear", "index": ALL}, "n_clicks"),
        State(IDs.UPLOAD_SELECTION, "data"),
        prevent_initial_call=True,
    )
    def edit(_rm, _clr, entries):
        trig = ctx.triggered_id
        value = ctx.triggered[0].get("value") if ctx.triggered else None
        entries = list(entries or [])
        if isinstance(trig, dict) and trig.get("type") == "mv-file-clear":
            return [] if value else no_update
        if isinstance(trig, dict) and trig.get("type") == "mv-file-remove":
            if not value:
                return no_update
            return [e for e in entries if e.get("id") != trig.get("index")]
        return no_update


# ── file manager + preflight + summary ─────────────────


def _register_file_manager(app) -> None:
    @app.callback(
        Output(IDs.FILE_MANAGER, "children"),
        Input(IDs.UPLOAD_SELECTION, "data"),
        prevent_initial_call=True,
    )
    def render_manager(entries):
        entries = list(entries or [])
        if not entries:
            return ""
        rows = []
        for e in entries:
            fname = e.get("filename") or e.get("original") or "file.edi"
            ext = fname.rsplit(".", 1)[-1].upper() if "." in fname else "FILE"
            rows.append(html.Div([
                html.Span(ext, className="mv-file-ext"),
                html.Span(fname, className="mv-file-name"),
                html.Button(
                    html.I(className="bi bi-trash"),
                    id={"type": "mv-file-remove", "index": e.get("id")},
                    className="mv-file-remove",
                    n_clicks=0,
                    title="Remove from import",
                ),
            ], className="mv-file-row"))
        return html.Div([
            html.Div([
                html.Span(f"{len(entries)} file(s) selected",
                          className="mv-file-head-title"),
                html.Button(
                    [html.I(className="bi bi-x-circle me-1"), "Clear"],
                    id={"type": "mv-file-clear", "index": "all"},
                    className="mv-file-clear",
                    n_clicks=0,
                ),
            ], className="mv-file-head"),
            html.Div(rows, className="mv-file-list"),
        ])


def _infer_lines(names):
    counts = {}
    for name in names:
        parts = [p for p in str(name).replace("\\", "/").split("/") if p]
        if len(parts) >= 3:
            line = parts[1]
        elif len(parts) == 2:
            line = parts[0]
        else:
            line = "flat files"
        counts[line] = counts.get(line, 0) + 1
    return counts


def _register_preflight(app) -> None:
    @app.callback(
        Output(IDs.PREFLIGHT, "children"),
        Output(IDs.DETECTED_SUMMARY, "children"),
        Input(IDs.UPLOAD_SELECTION, "data"),
        Input(IDs.LOAD_MODE_STORE, "data"),
        prevent_initial_call=True,
    )
    def preflight(entries, mode):
        entries = list(entries or [])
        if not entries:
            return "", [html.I(className="bi bi-info-circle me-1"),
                        "No data selected yet."]
        names = [e.get("original") or e.get("filename") for e in entries]
        line_counts = _infer_lines(names)
        source = entries[0].get("source", "upload")
        mode_lbl = "Add lines" if mode == _MODE_APPEND else "Replace survey"
        stats = html.Div([
            _stat(str(len(names)), "files"),
            _stat(str(len(line_counts)), "line groups"),
            _stat(mode_lbl, f"{source} source"),
        ], className="mv-preflight-stats")
        rows = html.Div([
            html.Div([
                html.Span(str(line), className="mv-pf-line"),
                html.Span(f"{cnt} file{'s' if cnt != 1 else ''}",
                          className="mv-pf-count"),
            ], className="mv-pf-row")
            for line, cnt in list(line_counts.items())[:8]
        ], className="mv-pf-rows")
        summary = [
            html.I(className="bi bi-check-circle-fill me-1"),
            f"{len(names)} file(s) · {len(line_counts)} line group(s)",
        ]
        return html.Div([stats, rows], className="mv-preflight-ready"), summary


def _stat(value, label):
    return html.Div([
        html.Div(value, className="mv-stat-value"),
        html.Div(label, className="mv-stat-label"),
    ], className="mv-stat")


def _register_drop_swap(app) -> None:
    clientside_callback(
        """
        function(manager) {
            var has = manager && !(Array.isArray(manager) && !manager.length)
                      && manager !== '';
            return [has ? {display:'none'} : {display:'block'},
                    has ? 'Selected files' : 'Or drag & drop files / folders'];
        }
        """,
        Output(IDs.DROP_WRAP, "style"),
        Output(IDs.DROP_TITLE, "children"),
        Input(IDs.FILE_MANAGER, "children"),
        prevent_initial_call=True,
    )


def _register_btn_state(app) -> None:
    clientside_callback(
        """
        function(entries) {
            var ok = entries && entries.some(function(e){ return e && e.content; });
            return !ok;
        }
        """,
        Output(IDs.BTN_LOAD_CONFIRM, "disabled"),
        Input(IDs.UPLOAD_SELECTION, "data"),
    )


# ── progress bar ───────────────────────────────────────


def _register_progress(app) -> None:
    clientside_callback(
        """
        function(n, entries) {
            if (!n) return [window.dash_clientside.no_update,
                            window.dash_clientside.no_update];
            var c = entries ? entries.filter(function(e){return e && e.content;}).length : 0;
            requestAnimationFrame(function() {
                var fill = document.getElementById('mv-prog-fill');
                if (fill) {
                    var clone = fill.cloneNode(true);
                    clone.classList.remove('mv-prog-done');
                    fill.parentNode.replaceChild(clone, fill);
                }
            });
            return [{'display':'block'},
                    'Parsing ' + c + ' file' + (c !== 1 ? 's' : '') + '…'];
        }
        """,
        Output(IDs.PROGRESS_WRAP, "style"),
        Output(IDs.PROG_LABEL, "children"),
        Input(IDs.BTN_LOAD_CONFIRM, "n_clicks"),
        State(IDs.UPLOAD_SELECTION, "data"),
        prevent_initial_call=True,
    )

    clientside_callback(
        """
        function(feedback) {
            if (!feedback) return window.dash_clientside.no_update;
            var fill = document.getElementById('mv-prog-fill');
            if (fill) fill.classList.add('mv-prog-done');
            setTimeout(function() {
                var w = document.getElementById('mv-progress-wrap');
                if (w) w.style.display = 'none';
            }, 800);
            return window.dash_clientside.no_update;
        }
        """,
        Output(IDs.PROGRESS_WRAP, "style", allow_duplicate=True),
        Input(IDs.LOAD_FEEDBACK, "children"),
        prevent_initial_call=True,
    )


# ── confirm load → MapView ─────────────────────────────


def _register_confirm_load(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Output(IDs.LOAD_FEEDBACK, "children"),
        Output(IDs.MODAL_LOAD, "is_open", allow_duplicate=True),
        Output(IDs.UPLOAD_SELECTION, "data", allow_duplicate=True),
        Output(IDs.DATA_BADGE_TEXT, "children", allow_duplicate=True),
        Output(IDs.DATA_BADGE, "className", allow_duplicate=True),
        Input(IDs.BTN_LOAD_CONFIRM, "n_clicks"),
        State(IDs.UPLOAD_SELECTION, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.LOAD_MODE_STORE, "data"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def confirm(_n, entries, session_id, mode, existing, theme):
        from pycsamt.map import MapView

        usable = [e for e in (entries or []) if e.get("content")]
        if not usable:
            return (no_update, "⚠ Choose a folder or drop files first.",
                    no_update, no_update, no_update, no_update)
        if not session_id:
            return (no_update, "⚠ Session not initialised — refresh.",
                    no_update, no_update, no_update, no_update)

        source = usable[0].get("source", "upload")
        theme = theme or "light"
        tmpdir = None
        try:
            line_map, tmpdir = _decode(usable, source)
            if not line_map:
                return (no_update, "⚠ No recognised EDI/AVG/J files.",
                        no_update, no_update, no_update, no_update)
            view = MapView.from_lines(line_map, theme=theme)
            if view.n_stations == 0:
                return (no_update, "⚠ No stations could be parsed.",
                        no_update, no_update, no_update, no_update)

            if mode == _MODE_APPEND:
                old = get_view(session_id)
                if old is not None:
                    view = merge_views(old, view)
            set_view(session_id, view)
            store = store_from_view(view)
            n_s, n_l = store["n_stations"], store["n_lines"]
            feedback = (f"✓ Loaded {n_s} station(s) from {n_l} line(s).")
            badge = f"{n_s} stations · {n_l} line(s)"
            return store, feedback, False, [], badge, "mv-data-badge visible"
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            return (no_update, f"✗ Error: {exc}",
                    no_update, no_update, no_update, no_update)
        finally:
            if tmpdir:
                shutil.rmtree(tmpdir, ignore_errors=True)


def _decode(usable, source):
    """Return ``({line: [paths]}, tmpdir)`` for the staged files."""
    from pycsamt.app.web.utils import (
        decode_folder_upload_to_tempdir,
        decode_upload_to_tempdir,
    )

    contents = [e.get("content") for e in usable]
    filenames = [e.get("filename") for e in usable]
    if source == "folder":
        _paths, tmpdir, line_map = decode_folder_upload_to_tempdir(
            contents, filenames
        )
        return line_map, tmpdir
    paths, tmpdir = decode_upload_to_tempdir(contents, filenames)
    return ({"uploaded": list(paths)} if paths else {}), tmpdir
