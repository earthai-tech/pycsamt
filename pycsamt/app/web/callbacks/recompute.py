# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/recompute.py — EDI recompute workflow (offcanvas panel).

Flow
----
1. Navbar "⟳ Recompute" button → opens offcanvas.
2. Source radio: "current" uses cached Sites; "path" loads EDIs from disk first.
3. User clicks "Recompute & Load":
   - If STORE_DATA.recomputed is True  → shows inline confirm card (no job yet).
   - Otherwise → starts background job immediately.
4. "Yes, recompute again" button → always starts job (bypasses confirm).
5. "Cancel" → hides confirm card.
6. dcc.Interval polls the job:
   - Phase-aware ticker: "Scanning folder…" → "⟳ HBH07 · 7/128" → "✓ Done".
   - Spin icon: rotating while running, check-circle when done, x-circle on error.
7. On done: STORE_DATA updated with recomputed:True, nav badge appears.
"""

from __future__ import annotations

import os
import threading
import time
from typing import Any

from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    html,
    no_update,
)

from pycsamt.app.web.cache import cache_get, cache_set
from pycsamt.app.web.layout import IDs

# ── Per-session job registry ──────────────────────────────────────────────────
_JOBS: dict[str, dict[str, Any]] = {}
_LOCK = threading.Lock()

_MAX_LOG  = 50
_LOG_SHOWN = 14

# Icon classNames for the spin indicator
_ICON_SPIN  = "bi bi-arrow-repeat recompute-spin me-1"
_ICON_OK    = "bi bi-check-circle-fill me-1 rcp-icon-ok"
_ICON_ERR   = "bi bi-x-circle-fill me-1 rcp-icon-err"
_ICON_WAIT  = "bi bi-arrow-repeat me-1"


# ── Background worker ─────────────────────────────────────────────────────────

def _worker(session_id: str, sites_or_path, opts: dict) -> None:
    """
    Two-phase worker.  ``sites_or_path`` is either a loaded Sites object
    (source=="current") or a folder path string (source=="path").
    """
    from pycsamt.site.base import Sites
    from pycsamt.site.recompute import EDIRecomputer

    try:
        # ── Phase 1: load from folder if needed ──────────────────────────
        if isinstance(sites_or_path, str):
            with _LOCK:
                _JOBS[session_id]["phase"] = "loading"

            edi_paths = _collect_edis(sites_or_path)

            with _LOCK:
                _JOBS[session_id]["found"] = len(edi_paths)

            if not edi_paths:
                with _LOCK:
                    _JOBS[session_id]["status"]  = "error"
                    _JOBS[session_id]["err_msg"] = (
                        f"No EDI files found in: {sites_or_path}"
                    )
                return

            sites = Sites(edi_paths)
        else:
            sites = sites_or_path

        # ── Phase 2: recompute ────────────────────────────────────────────
        edis  = sites.as_list()
        total = len(edis)

        with _LOCK:
            _JOBS[session_id].update({
                "phase":   "recomputing",
                "total":   total,
                "status":  "running",
                "t_start": time.time(),
            })

        rotate_angle = opts.get("rotate_angle") or None
        comps    = opts.get("comps") or ["Z", "Tip"]
        fmin     = opts.get("fmin") or None
        fmax     = opts.get("fmax") or None
        fill     = opts.get("fill") or None
        resphase = "on" in (opts.get("resphase") or [])

        recomputer = EDIRecomputer(
            rotate_angle        = rotate_angle,
            rotate_components   = tuple(comps),
            fmin                = fmin,
            fmax                = fmax,
            fill_missing_values = fill,
            recompute_resphase  = resphase,
            write               = False,
            copy                = True,
            strict              = False,
            progress            = False,
        )

        result_edis = []
        failed = 0

        for i, edi in enumerate(edis):
            sname = _safe_name(edi, i)
            try:
                out = recomputer.recompute_one(edi)
                result_edis.append(out)
                entry = f"✓ {sname}"
            except Exception as exc:
                result_edis.append(edi)
                failed += 1
                entry = f"✗ {sname}: {str(exc)[:70]}"

            with _LOCK:
                job = _JOBS[session_id]
                job["done"]    = i + 1
                job["current"] = sname
                job["failed"]  = failed
                job["log"].append(entry)
                if len(job["log"]) > _MAX_LOG:
                    job["log"] = job["log"][-_MAX_LOG:]

        result_sites = Sites(result_edis)
        with _LOCK:
            _JOBS[session_id]["sites"]  = result_sites
            _JOBS[session_id]["status"] = "done"

    except Exception as exc:
        with _LOCK:
            _JOBS[session_id]["status"]  = "error"
            _JOBS[session_id]["err_msg"] = str(exc)


def _collect_edis(folder: str) -> list[str]:
    paths = []
    for root, _dirs, files in os.walk(folder):
        for f in files:
            if f.lower().endswith(".edi"):
                paths.append(os.path.join(root, f))
    return sorted(paths)


def _safe_name(edi, idx: int) -> str:
    try:
        from pycsamt.site.utils import station_name as _sn
        return _sn(edi) or f"site_{idx+1:04d}"
    except Exception:
        return f"site_{idx+1:04d}"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _sites_to_store(sites, old_store: dict | None = None) -> dict:
    from pycsamt.app.desktop.controllers.data_controller import (
        DataController,
    )

    ctrl = DataController()
    ctrl._sites = sites
    if old_store and old_store.get("station_records"):
        records = old_store["station_records"]
        ctrl._station_to_line = {r["ID"]: r.get("Line", "—") for r in records}
    ctrl._df = ctrl._build_dataframe()
    df = ctrl.dataframe

    n_lines     = df["Line"].nunique() if "Line" in df.columns else 1
    line_counts = df.groupby("Line").size().to_dict() if "Line" in df.columns else {}

    return {
        "station_records": df.to_dict("records"),
        "data_dir":        "[recomputed]",
        "n_stations":      len(df),
        "n_lines":         n_lines,
        "line_counts":     line_counts,
        "recomputed":      True,   # flag consumed by nav badge + confirm gate
    }


def _log_chips(entries: list[str]) -> list:
    shown = entries[-_LOG_SHOWN:]
    chips = []
    for e in reversed(shown):
        ok = e.startswith("✓")
        chips.append(html.Span(
            e,
            className="recompute-log-chip "
                      + ("recompute-log-ok" if ok else "recompute-log-err"),
        ))
    return chips


# ── Callback registration ─────────────────────────────────────────────────────

def register_recompute(app) -> None:
    _register_open_canvas(app)
    _register_path_section(app)
    _register_nav_badge(app)
    _register_run(app)
    _register_cancel_confirm(app)
    _register_poll(app)


def _register_open_canvas(app) -> None:
    @app.callback(
        Output(IDs.RECOMPUTE_CANVAS, "is_open"),
        Input(IDs.BTN_RECOMPUTE_OPEN, "n_clicks"),
        State(IDs.RECOMPUTE_CANVAS,   "is_open"),
        prevent_initial_call=True,
    )
    def toggle_canvas(n_clicks, is_open):
        if n_clicks:
            return not is_open
        return no_update


def _register_path_section(app) -> None:
    clientside_callback(
        """
        function(source) {
            return source === 'path'
                ? {display: 'block'}
                : {display: 'none'};
        }
        """,
        Output(IDs.RECOMPUTE_PATH_SEC, "style"),
        Input(IDs.RECOMPUTE_SOURCE,    "value"),
    )


def _register_nav_badge(app) -> None:
    """Show/hide the navbar '✓ recomputed' pill from STORE_DATA flag."""
    clientside_callback(
        """
        function(data) {
            if (data && data.recomputed) {
                var n = data.n_stations || '';
                return {display: 'inline-flex'};
            }
            return {display: 'none'};
        }
        """,
        Output(IDs.RECOMPUTE_NAV_BADGE, "style"),
        Input(IDs.STORE_DATA, "data"),
    )


def _register_cancel_confirm(app) -> None:
    @app.callback(
        Output(IDs.RECOMPUTE_CONFIRM_SEC, "style"),
        Input(IDs.BTN_RECOMPUTE_CANCEL,   "n_clicks"),
        prevent_initial_call=True,
    )
    def cancel_confirm(n_clicks):
        if n_clicks:
            return {"display": "none"}
        return no_update


def _register_run(app) -> None:
    # 12 outputs ──────────────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.RECOMPUTE_PROG_SEC,    "style"),
        Output(IDs.RECOMPUTE_BADGE,       "children"),
        Output(IDs.RECOMPUTE_BADGE,       "color"),
        Output(IDs.RECOMPUTE_TICKER,      "children"),
        Output(IDs.RECOMPUTE_PROGRESS,    "value"),
        Output(IDs.RECOMPUTE_PROGRESS,    "animated"),
        Output(IDs.RECOMPUTE_PROGRESS,    "striped"),
        Output(IDs.RECOMPUTE_LOG,         "children"),
        Output(IDs.RECOMPUTE_INTERVAL,    "disabled"),
        Output(IDs.RECOMPUTE_FEEDBACK,    "children"),
        Output(IDs.RECOMPUTE_CONFIRM_SEC, "style",     allow_duplicate=True),
        Output(IDs.RECOMPUTE_SPIN_ICON,   "className"),
        Input(IDs.BTN_RECOMPUTE,          "n_clicks"),
        Input(IDs.BTN_RECOMPUTE_YES,      "n_clicks"),
        State(IDs.RECOMPUTE_SOURCE,       "value"),
        State(IDs.RECOMPUTE_PATH,         "value"),
        State(IDs.RECOMPUTE_RESPHASE,     "value"),
        State(IDs.RECOMPUTE_ROTATE_ANG,   "value"),
        State(IDs.RECOMPUTE_COMPS,        "value"),
        State(IDs.RECOMPUTE_FMIN,         "value"),
        State(IDs.RECOMPUTE_FMAX,         "value"),
        State(IDs.RECOMPUTE_FILL,         "value"),
        State(IDs.SESSION_ID,             "data"),
        State(IDs.STORE_DATA,             "data"),
        prevent_initial_call=True,
    )
    def start_recompute(n_main, n_yes, source, path, resphase, rotate_ang,
                        comps, fmin, fmax, fill, session_id, store_data):
        _SHOW     = {"display": "block"}
        _HIDE     = {"display": "none"}
        _SKIP     = (no_update,) * 12

        triggered = ctx.triggered_id

        if not triggered:
            return _SKIP

        if not session_id:
            return (_HIDE, "No session", "danger", "", 0, False, False, [],
                    True, "⚠ Session not initialised — refresh the page.",
                    _HIDE, _ICON_WAIT)

        source = source or "current"

        # ── Confirm gate — only on the main button, not on "Yes" ─────────
        already_done = bool(store_data and store_data.get("recomputed"))
        if triggered == IDs.BTN_RECOMPUTE and already_done:
            return (_HIDE, no_update, no_update, no_update, no_update,
                    no_update, no_update, no_update, True, "",
                    _SHOW, _ICON_WAIT)

        # ── Validate data source ──────────────────────────────────────────
        if source == "current":
            sites = cache_get(session_id)
            if sites is None:
                if not store_data:
                    msg = ("⚠ No data loaded — use Load Data first, "
                           "or switch source to 'Load from folder path'.")
                else:
                    msg = ("⚠ Site cache expired — reload your EDI files "
                           "or switch to 'Load from folder path'.")
                return (_HIDE, "No data", "warning", "", 0, False, False,
                        [], True, msg, _HIDE, _ICON_WAIT)
            sites_or_path = sites

        else:  # "path"
            folder = (path or "").strip()
            if not folder:
                return (_HIDE, "No path", "warning", "", 0, False, False, [],
                        True, "⚠ Enter a folder path.", _HIDE, _ICON_WAIT)
            if not os.path.isdir(folder):
                return (_HIDE, "Not found", "danger", "", 0, False, False, [],
                        True, f"⚠ Folder not found: {folder}", _HIDE, _ICON_WAIT)
            sites_or_path = folder

        # ── Start job ─────────────────────────────────────────────────────
        with _LOCK:
            _JOBS[session_id] = {
                "status":  "starting",
                "phase":   "loading" if source == "path" else "recomputing",
                "total":   0,
                "done":    0,
                "failed":  0,
                "found":   0,
                "current": "",
                "log":     [],
                "sites":   None,
                "err_msg": "",
                "t_start": time.time(),
            }

        opts = {
            "resphase":     resphase,
            "rotate_angle": float(rotate_ang) if rotate_ang is not None else None,
            "comps":        comps or ["Z", "Tip"],
            "fmin":         float(fmin) if fmin is not None else None,
            "fmax":         float(fmax) if fmax is not None else None,
            "fill":         fill or None,
        }

        t = threading.Thread(
            target=_worker, args=(session_id, sites_or_path, opts), daemon=True
        )
        t.start()

        init_ticker = (
            "Scanning folder for EDI files…"
            if source == "path" else "Initialising recompute…"
        )

        return (
            _SHOW,            # show progress section
            "Starting…",      # badge text
            "warning",        # badge colour
            init_ticker,      # ticker
            0,                # progress value
            True,             # animated
            True,             # striped
            [],               # log (empty)
            False,            # enable interval
            "",               # clear feedback
            _HIDE,            # hide confirm card
            _ICON_SPIN,       # start spinning icon
        )


def _register_poll(app) -> None:
    @app.callback(
        Output(IDs.RECOMPUTE_TICKER,    "children",   allow_duplicate=True),
        Output(IDs.RECOMPUTE_PROGRESS,  "value",      allow_duplicate=True),
        Output(IDs.RECOMPUTE_PROGRESS,  "animated",   allow_duplicate=True),
        Output(IDs.RECOMPUTE_PROGRESS,  "striped",    allow_duplicate=True),
        Output(IDs.RECOMPUTE_LOG,       "children",   allow_duplicate=True),
        Output(IDs.RECOMPUTE_BADGE,     "children",   allow_duplicate=True),
        Output(IDs.RECOMPUTE_BADGE,     "color",      allow_duplicate=True),
        Output(IDs.RECOMPUTE_INTERVAL,  "disabled",   allow_duplicate=True),
        Output(IDs.STORE_DATA,          "data",       allow_duplicate=True),
        Output(IDs.RECOMPUTE_FEEDBACK,  "children",   allow_duplicate=True),
        Output(IDs.RECOMPUTE_SPIN_ICON, "className",  allow_duplicate=True),
        Input(IDs.RECOMPUTE_INTERVAL,   "n_intervals"),
        State(IDs.SESSION_ID,           "data"),
        State(IDs.STORE_DATA,           "data"),
        prevent_initial_call=True,
    )
    def poll_recompute(n_intervals, session_id, store_data):
        _skip = (no_update,) * 11

        if not session_id:
            return _skip

        with _LOCK:
            job = _JOBS.get(session_id)
            if job is None:
                return _skip
            status  = job["status"]
            phase   = job.get("phase", "recomputing")
            total   = job["total"]
            done    = job["done"]
            failed  = job["failed"]
            found   = job.get("found", 0)
            current = job["current"]
            log     = list(job["log"])
            err_msg = job["err_msg"]
            t_start = job["t_start"]
            sites   = job.get("sites")

        pct     = int(done / total * 100) if total > 0 else 0
        elapsed = time.time() - t_start

        # ── Running ───────────────────────────────────────────────────────
        if status in ("starting", "running"):
            if phase == "loading":
                ticker    = (f"Scanning folder…  {found} EDI files found"
                             if found else "Scanning folder for EDI files…")
                badge_txt = "Loading"
                pct       = 0
            else:
                if total > 0:
                    eta_s  = (elapsed / done * (total - done)) if done > 0 else 0
                    eta    = f"  ·  ETA {eta_s:.0f} s" if eta_s > 0 else ""
                    ticker = (
                        f"⟳  {current}  ·  {done} / {total}"
                        + (f"  ({failed} ✗)" if failed else "")
                        + eta
                    )
                else:
                    ticker = "Initialising recompute…"
                badge_txt = "Running" if status == "running" else "Starting"

            return (
                ticker, pct, True, True,
                _log_chips(log),
                badge_txt, "warning",
                False,        # keep interval alive
                no_update,    # don't touch store yet
                no_update,
                _ICON_SPIN,   # keep spinning
            )

        # ── Done ──────────────────────────────────────────────────────────
        if status == "done" and sites is not None:
            cache_set(session_id, sites)

            try:
                new_store = _sites_to_store(sites, old_store=store_data)
            except Exception:
                new_store = store_data

            badge_txt = (
                f"✓ {total} stations"
                + (f"  ·  {failed} failed" if failed else "")
            )
            ticker = (
                f"✓ {total} station{'s' if total != 1 else ''} recomputed"
                + (f"  ·  {failed} failed" if failed else "")
                + f"  ·  {elapsed:.1f} s"
            )
            feedback = (
                f"✓ {total} station(s) recomputed in {elapsed:.1f} s"
                + (f"  ·  {failed} failed (originals kept)" if failed else "")
                + "  —  app now uses recomputed data."
            )

            with _LOCK:
                _JOBS[session_id]["sites"] = None

            return (
                ticker, 100, False, False,
                _log_chips(log),
                badge_txt, "success",
                True,         # disable interval
                new_store,
                feedback,
                _ICON_OK,     # swap to check-circle, no spin
            )

        # ── Error ─────────────────────────────────────────────────────────
        if status == "error":
            ticker = f"✗ {err_msg[:80]}"
            return (
                ticker, pct, False, False,
                _log_chips(log),
                "✗ Error", "danger",
                True,         # disable interval
                no_update,
                f"✗ Error: {err_msg}",
                _ICON_ERR,    # swap to x-circle, no spin
            )

        return _skip
