# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/inversion.py — MT inversion launch + coloured solver log."""

from __future__ import annotations

import io
import re
import sys
import traceback

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import empty_src, fig_to_src


# ── log-line colour rules (first match wins) ──────────────────────────────────
_LOG_RULES: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^\s*(={3,}|-{3,}|\*{3,}|#{2,})"),            "log-head"),
    (re.compile(r"iter|iteration|step\s+\d", re.I),             "log-iter"),
    (re.compile(r"converge|misfit\s*<|rms.*<|done|finish", re.I), "log-conv"),
    (re.compile(r"warn|caution", re.I),                          "log-warn"),
    (re.compile(r"error|fail|traceback|exception", re.I),        "log-err"),
    (re.compile(r"load|read|write|sav|export", re.I),            "log-info"),
]


def _classify(line: str) -> str:
    for pattern, cls in _LOG_RULES:
        if pattern.search(line):
            return cls
    return ""


def _coloured_log(lines: list[str]) -> list[html.Span]:
    spans = []
    for line in lines:
        cls = _classify(line)
        if cls:
            spans.append(html.Span(line + "\n", className=cls))
        else:
            spans.append(html.Span(line + "\n"))
    return spans


def register_inversion(app) -> None:

    @app.callback(
        Output(IDs.INV_LOG,      "children"),
        Output(IDs.IMG_INV,      "src"),
        Output(IDs.INV_SPINNER,  "children"),
        Output(IDs.TOAST_ERROR,  "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,   "children", allow_duplicate=True),
        Input(IDs.BTN_INV_RUN,  "n_clicks"),
        State(IDs.INV_SOLVER,   "value"),
        State(IDs.INV_DIM,      "value"),
        State(IDs.INV_DATA_DIR, "value"),
        State(IDs.SESSION_ID,   "data"),
        prevent_initial_call=True,
    )
    def launch_inversion(n_clicks, solver, dim, data_dir, session_id):
        if not n_clicks:
            raise PreventUpdate

        log_lines: list[str] = []
        src = empty_src(dark=True)

        def _log(msg: str) -> None:
            log_lines.append(msg)

        solver  = (solver  or "occam2d").lower()
        dim     = (dim     or "2d").lower()

        _log("=" * 60)
        _log(f"  pycsamt Inversion  —  solver={solver}  dim={dim}")
        _log("=" * 60)

        try:
            sites = cache_get(session_id)

            # ── AI / Neural solvers ───────────────────────────────────────
            if solver in ("em1d", "em2d", "gcn3d"):
                _log(f"Loading AI solver: {solver.upper()} …")
                import pycsamt.ai as ai_mod

                cls_map = {
                    "em1d":  ("nets", "EMInverter1D"),
                    "em2d":  ("nets", "EMInverter2D"),
                    "gcn3d": ("nets", "GCNInverter3D"),
                }
                pkg, cls_name = cls_map[solver]
                pkg_mod = getattr(ai_mod, pkg, None)
                if pkg_mod is None:
                    raise ImportError(f"pycsamt.ai.{pkg} not available.")
                inverter_cls = getattr(pkg_mod, cls_name)
                inverter = inverter_cls()
                _log(f"Instantiated {cls_name}.")

                if sites is None:
                    _log("INFO  No survey data in session — running in demo mode.")
                    _log("      Load EDI files to invert real data.")
                    _log("-" * 40)
                    _log("Iteration   1 / 50   misfit=12.45")
                    _log("Iteration  10 / 50   misfit= 4.12")
                    _log("Iteration  25 / 50   misfit= 1.83")
                    _log("Converged:  RMS misfit = 0.97  (threshold=1.05)")
                else:
                    buf = io.StringIO()
                    old_stdout = sys.stdout
                    sys.stdout = buf
                    try:
                        result = inverter.fit(sites)
                    finally:
                        sys.stdout = old_stdout
                    for ln in buf.getvalue().splitlines():
                        _log(ln)
                    _log("Done.")

                    # Try to produce a result plot
                    try:
                        import matplotlib
                        matplotlib.use("Agg")
                        import matplotlib.pyplot as plt
                        from pycsamt.app.web.utils import apply_web_dark_theme
                        apply_web_dark_theme()
                        fig = plt.figure(figsize=(10, 4))
                        ax  = fig.add_subplot(111)
                        if hasattr(result, "plot"):
                            result.plot(ax=ax)
                        src = fig_to_src(fig)
                    except Exception:
                        pass

            # ── Traditional solvers (external binaries) ───────────────────
            else:
                _log(f"Launching external solver: {solver} …")
                if sites is None:
                    _log("WARNING  No survey data loaded — cannot write input files.")
                    _log("         Load EDI data first, then re-run.")
                else:
                    from pycsamt.inversion.backends.occam2d import Occam2DBackend
                    backend_map = {"occam2d": Occam2DBackend}
                    if solver not in backend_map:
                        _log(f"INFO  Backend for '{solver}' not yet wired.")
                        _log("      Supported: occam2d")
                    else:
                        backend = backend_map[solver](sites=sites)
                        buf = io.StringIO()
                        old_stdout = sys.stdout
                        sys.stdout = buf
                        try:
                            result = backend.run()
                        finally:
                            sys.stdout = old_stdout
                        for ln in buf.getvalue().splitlines():
                            _log(ln)
                        _log(f"Solver finished.  Status: {getattr(result, 'status', 'done')}")

        except Exception as exc:
            _log(f"ERROR  {exc}")
            _log(traceback.format_exc()[-600:])
            spans = _coloured_log(log_lines)
            return spans, src, "", True, str(exc)

        _log("=" * 60)
        spans = _coloured_log(log_lines)
        return spans, src, "", False, ""
