# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Post-correction action modal callbacks.

After static-shift, denoising, or QC
workflows finish, the post-proc modal lets
the user either apply corrected data to the
active session or export to a custom folder.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs

logger = logging.getLogger(__name__)


def _export(sites, dest: Path) -> list:
    """Write corrected sites as EDI files."""
    from pycsamt.site.export import write_sites

    return write_sites(sites, dest, exist_ok=True)


def _count_valid_impedance(sites) -> int:
    """Count stations in a Sites object that carry parseable Z."""
    from pycsamt.agents.loader import _quality_scan

    rows, _ = _quality_scan(sites)
    return sum(1 for r in rows if r.get("has_z"))


def _validate_export(dest: Path) -> tuple[int, int, str]:
    """Re-load an exported folder and count stations with valid Z.

    Returns ``(n_valid_impedance, n_total, detail)``. ``detail`` carries a
    short diagnostic when something went wrong (an exception, or files that
    loaded but lacked Z) so the failure isn't silently reported as "0 valid"
    — which previously made a validation crash look identical to genuinely
    corrupt data. Used to refuse repointing the session at a folder that
    would later fail every workflow with "No stations with valid impedance".
    """
    try:
        from pycsamt.emtools._core import ensure_sites

        n_files = len(list(Path(dest).rglob("*.edi")))
        s = ensure_sites(str(dest), recursive=True, verbose=0)
        from pycsamt.agents.loader import _quality_scan

        rows, _ = _quality_scan(s)
        n_valid = sum(1 for r in rows if r.get("has_z"))
        detail = (
            ""
            if n_valid
            else f"{n_files} .edi file(s) on disk, {len(rows)} loaded, "
            f"0 with a parseable Z block"
        )
        return n_valid, len(rows), detail
    except Exception as exc:  # noqa: BLE001
        logger.exception("Export validation failed for %s", dest)
        return -1, 0, f"validation error: {type(exc).__name__}: {exc}"


def register_postproc(app) -> None:

    # 1. Open modal when STORE_POSTPROC is set
    @app.callback(
        Output(IDs.MODAL_POSTPROC, "is_open"),
        Output(IDs.POSTPROC_MSG, "children"),
        Input(IDs.STORE_POSTPROC, "data"),
        prevent_initial_call=True,
    )
    def open_modal(data):
        if not (data or {}).get("jid"):
            raise PreventUpdate
        wtype = data.get("workflow", "processing")
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
        Output(IDs.POSTPROC_STATUS, "children"),
        Output(
            IDs.MODAL_POSTPROC,
            "is_open",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_POSTPROC_APPLY, "n_clicks"),
        State(IDs.STORE_POSTPROC, "data"),
        State(IDs.STORE_SETTINGS, "data"),
        State(IDs.STORE_EDI, "data"),
        prevent_initial_call=True,
    )
    def apply_to_session(n, postproc, settings, edi_store):
        if not n:
            raise PreventUpdate
        from .chat import _CORR_CACHE

        jid = (postproc or {}).get("jid")
        wtype = (postproc or {}).get("workflow", "corrected")
        base = (
            (postproc or {}).get("output_dir", "")
            or (settings or {}).get("output_dir", "")
            or "pycsamt_workflow_output"
        ).strip()
        # Peek (don't pop): if the export turns out invalid we keep the
        # corrected data in memory so the user can retry / export instead.
        sites = _CORR_CACHE.get(jid)
        if sites is None:
            return (
                no_update,
                html.Span(
                    "Corrected data not in memory — re-run the workflow.",
                    style={"color": "red"},
                ),
                True,
            )
        # Absolute path so the re-load validation and the later workflow
        # loader never depend on the process CWD.
        dest = Path(base).expanduser().resolve() / f"{wtype}_{int(time.time())}"
        # Sanity-check the in-memory corrected data first (source of truth).
        try:
            n_mem = _count_valid_impedance(sites)
        except Exception:  # noqa: BLE001
            n_mem = -1
        try:
            paths = _export(sites, dest)
        except Exception as exc:
            return (
                no_update,
                html.Span(
                    str(exc),
                    style={"color": "red"},
                ),
                True,
            )

        # Validate before repointing the session. Repointing at a folder
        # whose EDIs carry no impedance silently corrupts every later
        # workflow ("No stations with valid impedance data found"). If the
        # export is bad, leave STORE_EDI untouched and tell the user *why*.
        n_valid, n_total, detail = _validate_export(dest)
        if n_valid <= 0:
            return (
                no_update,
                html.Span(
                    f"Could not apply corrected data — session left "
                    f"unchanged (kept in memory). {len(paths)} file(s) "
                    f"written to {dest}; in-memory stations with valid "
                    f"Z: {n_mem}. Reason: {detail or 'unknown'}.",
                    style={"color": "red"},
                ),
                True,
            )

        _CORR_CACHE.pop(jid, None)
        new_edi = dict(edi_store or {})
        new_edi["path"] = str(dest)
        # Drop any stale in-memory sites / groups so the loader re-reads
        # the corrected directory rather than the original upload.
        new_edi.pop("sites", None)
        new_edi.pop("groups", None)
        new_edi.pop("selected_lines", None)
        msg = html.Span(
            f"{n_valid}/{n_total} corrected EDI(s) applied " f"to session → {dest}",
            style={"color": "var(--tag-ok)"},
        )
        return new_edi, msg, False

    # 3. "Export to folder" — toggle collapse
    @app.callback(
        Output(IDs.POSTPROC_COLLAPSE, "is_open"),
        Input(IDs.BTN_POSTPROC_EXPORT, "n_clicks"),
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
        Input(IDs.BTN_POSTPROC_OK, "n_clicks"),
        State(IDs.POSTPROC_PATH, "value"),
        State(IDs.STORE_POSTPROC, "data"),
        prevent_initial_call=True,
    )
    def confirm_export(n, path_val, postproc):
        if not n:
            raise PreventUpdate
        if not (path_val or "").strip():
            return (
                html.Span(
                    "Enter a folder path first.",
                    style={"color": "orange"},
                ),
                True,
                True,
            )
        from .chat import _CORR_CACHE

        jid = (postproc or {}).get("jid")
        sites = _CORR_CACHE.get(jid)
        if sites is None:
            return (
                html.Span(
                    "Corrected data not in memory — re-run workflow.",
                    style={"color": "red"},
                ),
                True,
                True,
            )
        dest = Path(path_val.strip()).expanduser().resolve()
        try:
            n_mem = _count_valid_impedance(sites)
        except Exception:  # noqa: BLE001
            n_mem = -1
        try:
            paths = _export(sites, dest)
        except Exception as exc:
            return (
                html.Span(
                    str(exc),
                    style={"color": "red"},
                ),
                True,
                True,
            )
        n_valid, n_total, detail = _validate_export(dest)
        if n_valid <= 0:
            return (
                html.Span(
                    f"Wrote {len(paths)} file(s) to {dest}, but the "
                    f"re-load check did not confirm valid impedance "
                    f"(in-memory valid stations: {n_mem}). "
                    f"Reason: {detail or 'unknown'}. Corrected data "
                    f"kept in memory.",
                    style={"color": "red"},
                ),
                True,
                True,
            )
        _CORR_CACHE.pop(jid, None)
        return (
            html.Span(
                f"Exported {n_valid}/{n_total} corrected EDI(s) to {dest}",
                style={"color": "var(--tag-ok)"},
            ),
            False,
            False,
        )
