# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Import a backend-neutral PCSF/PCSM inversion result into the map-view
session.

Two independent sources can turn up ``.pcsf``/``.pcsm``/``.pcsm.gz``
candidate files, both landing in :data:`IDs.INV_CANDIDATES_STORE`:

- **Browse folder** (``assets/pcsf_loader.js``, modeled on
  ``assets/edi_loader.js``): a ``webkitdirectory`` input is injected
  over the button, the folder is scanned for matching extensions
  client-side, and every match found is staged into
  :data:`IDs.INV_FOLDER_STORE`.
- **Drop / browse files** (:data:`IDs.INV_PCSF_UPLOAD`, a plain
  :class:`dash.dcc.Upload`, ``multiple=True``) — needs no custom JS.

``_register_capture_candidates`` merges whichever source last fired
into one candidate list. ``_register_classify_candidates`` then peeks
each candidate's PCSF ``geometry.kind`` (via
:func:`pycsamt.format.peek_kind` — no full array load) and renders a
picker: every one of PCSF's four geometry kinds is selectable —
:func:`pycsamt.map.inversion.load_pcsf_lines` builds a real
per-station curtain for all of them, including a native ModEM
``grid3d`` volume (nearest-cell sampling, same registration approach
:func:`pycsamt.models.modem.section.station_curtain` uses for a live
ModEM folder) and a MARE2DEM ``mesh_unstructured`` mesh
(point-location on its real triangulation, since an unstructured mesh
has no rectilinear index to look a column up by). Only a candidate
this peek could not even read (a corrupt/truncated upload) is
disabled. Auto-selects when exactly one candidate is importable, and
labels the rest instead of letting the user pick one and hit a raw
import error.
``_register_resolve_pick`` turns the chosen candidate into
:data:`IDs.INV_RESOLVED_STORE`, and ``_register_confirm`` decodes that
one file and loads it through :meth:`pycsamt.map.MapView.from_pcsf` —
into the same session store as EDI loading, so the 3-D fence/depth
views work identically regardless of data source.

A raw ModEM/Occam2D/MARE2DEM result folder is not accepted directly
here any more — convert it to PCSF/PCSM first (see
``pycsamt.format.adapters``).
:meth:`pycsamt.map.MapView.from_inversion_results` still exists for
script/notebook use; it is just no longer wired to this UI.
"""

from __future__ import annotations

import base64
import os
import shutil
import tempfile
from pathlib import Path

from dash import Input, Output, State, ctx, no_update

from .._ids import IDs
from .._render import merge_views, store_from_view
from ..cache import get_view, set_view

_MODE_APPEND = "append"
_PCSF_EXTS = (".pcsf", ".pcsm", ".pcsm.gz")
_IMPORTABLE_KINDS = {"grid2d", "multiline", "grid3d", "mesh_unstructured"}


def register_inversion_import(app) -> None:
    _register_capture_candidates(app)
    _register_classify_candidates(app)
    _register_resolve_pick(app)
    _register_confirm(app)


def _is_pcsf_name(name: str) -> bool:
    return str(name).lower().endswith(_PCSF_EXTS)


def _decode_to_tempfile(filename: str, content: str) -> str:
    """Write one staged base64 file to a real temp file, preserving its
    (possibly double, e.g. ``.pcsm.gz``) extension so format-dispatch
    by suffix works the same as it does for a file already on disk."""
    tmpdir = tempfile.mkdtemp(prefix="pycsamt_pcsf_upload_")
    _, b64 = content.split(",", 1)
    raw = base64.b64decode(b64)
    path = os.path.join(tmpdir, Path(filename).name)
    with open(path, "wb") as fh:
        fh.write(raw)
    return path


def _peek_kind_safe(filename: str, content: str) -> str | None:
    """Best-effort :func:`pycsamt.format.peek_kind` for one staged
    candidate; ``None`` on any failure (corrupt/unreadable upload) --
    the picker shows it as unreadable rather than crashing the tab."""
    from pycsamt.format import peek_kind

    tmpdir = None
    try:
        path = _decode_to_tempfile(filename, content)
        tmpdir = os.path.dirname(path)
        return peek_kind(path)
    except Exception:  # noqa: BLE001 - a bad candidate must not break the list
        return None
    finally:
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ── capture: merge folder-scan + drop/browse into one candidate list ──


def _register_capture_candidates(app) -> None:
    @app.callback(
        Output(IDs.INV_CANDIDATES_STORE, "data"),
        Input(IDs.INV_PCSF_UPLOAD, "contents"),
        Input(IDs.INV_FOLDER_STORE, "data"),
        State(IDs.INV_PCSF_UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def capture(upload_contents, folder_data, upload_filenames):
        trig = ctx.triggered_id
        if trig == IDs.INV_FOLDER_STORE:
            data = folder_data or {}
            filenames = list(data.get("filenames") or [])
            contents = list(data.get("contents") or [])
        else:
            contents = upload_contents
            filenames = upload_filenames
            if contents is None:
                filenames, contents = [], []
            elif isinstance(contents, str):
                filenames, contents = [filenames], [contents]
            else:
                filenames = list(filenames or [])
                contents = list(contents or [])
        pairs = [
            (n, c) for n, c in zip(filenames, contents) if _is_pcsf_name(n)
        ]
        return {
            "filenames": [n for n, _ in pairs],
            "contents": [c for _, c in pairs],
        }


# ── classify: peek kind, render the picker ─────────────────────────


def _register_classify_candidates(app) -> None:
    @app.callback(
        Output(IDs.INV_CANDIDATE_PICKER, "options"),
        Output(IDs.INV_CANDIDATE_PICKER, "value"),
        Output(IDs.INV_CANDIDATE_WRAP, "style"),
        Output(IDs.INV_FILE_COUNT, "children"),
        Output(IDs.INV_BROWSE_STATUS, "children"),
        Input(IDs.INV_CANDIDATES_STORE, "data"),
        prevent_initial_call=True,
    )
    def classify(candidates):
        filenames = (candidates or {}).get("filenames") or []
        contents = (candidates or {}).get("contents") or []
        if not filenames:
            return [], None, {"display": "none"}, "", ""

        options = []
        enabled_values = []
        for i, (name, content) in enumerate(zip(filenames, contents)):
            kind = _peek_kind_safe(name, content)
            value = str(i)
            if kind in _IMPORTABLE_KINDS:
                options.append({"label": f"{name}  ({kind})", "value": value})
                enabled_values.append(value)
            elif kind is None:
                options.append(
                    {
                        "label": f"{name}  (unreadable)",
                        "value": value,
                        "disabled": True,
                    }
                )
            else:
                # Reachable only for a kind pycsamt.format itself does
                # not recognise (a hand-crafted or future-format file --
                # peek_kind reads whatever string is there without
                # validating it against GEOMETRY_KINDS).
                options.append(
                    {
                        "label": f"{name}  ({kind} — not a recognised PCSF kind)",
                        "value": value,
                        "disabled": True,
                    }
                )

        value = enabled_values[0] if len(enabled_values) == 1 else None
        n = len(filenames)
        count_text = f"{n} file{'s' if n != 1 else ''} found"
        if not enabled_values:
            status = "⚠ None of these are importable here yet — see labels below."
        elif len(enabled_values) == 1:
            status = "✓ Ready to import."
        else:
            status = f"{len(enabled_values)} importable — pick one below."
        return options, value, {"display": "block"}, count_text, status


# ── resolve: chosen candidate -> the single-file store confirm() reads ──


def _register_resolve_pick(app) -> None:
    @app.callback(
        Output(IDs.INV_RESOLVED_STORE, "data"),
        Input(IDs.INV_CANDIDATE_PICKER, "value"),
        State(IDs.INV_CANDIDATES_STORE, "data"),
        prevent_initial_call=True,
    )
    def resolve(value, candidates):
        if value is None:
            return {}
        filenames = (candidates or {}).get("filenames") or []
        contents = (candidates or {}).get("contents") or []
        idx = int(value)
        if idx < 0 or idx >= len(filenames):
            return {}
        return {"filenames": [filenames[idx]], "contents": [contents[idx]]}


def _register_confirm(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Output(IDs.INV_STATUS, "children"),
        Output(IDs.MODAL_LOAD, "is_open", allow_duplicate=True),
        Output(IDs.DATA_BADGE_TEXT, "children", allow_duplicate=True),
        Output(IDs.DATA_BADGE, "className", allow_duplicate=True),
        Input(IDs.BTN_INV_CONFIRM, "n_clicks"),
        State(IDs.INV_RESOLVED_STORE, "data"),
        State(IDs.CK_INV_KNOWN_STA, "value"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.LOAD_MODE_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def confirm(_n, staged, use_known, session_id, mode, theme):
        filenames = (staged or {}).get("filenames") or []
        contents = (staged or {}).get("contents") or []
        if not filenames:
            return (
                no_update,
                "⚠ Browse to a folder or drop a .pcsf/.pcsm/.pcsm.gz "
                "file, then pick one to import, first.",
                no_update,
                no_update,
                no_update,
            )
        if not session_id:
            return (
                no_update,
                "⚠ Session not initialised — refresh.",
                no_update,
                no_update,
                no_update,
            )

        theme = theme or "light"
        old = get_view(session_id)
        known_stations = None
        if use_known and old is not None:
            known_stations = old.data.stations

        tmpdir = None
        try:
            from pycsamt.map import MapView

            tmp_path = _decode_to_tempfile(filenames[0], contents[0])
            tmpdir = os.path.dirname(tmp_path)
            view = MapView.from_pcsf(
                tmp_path,
                known_stations=known_stations,
                theme=theme,
            )
            if view.n_stations == 0:
                return (
                    no_update,
                    "⚠ No stations could be parsed from that file.",
                    no_update,
                    no_update,
                    no_update,
                )

            if mode == _MODE_APPEND and old is not None:
                view = merge_views(old, view)
            set_view(session_id, view)
            store = store_from_view(view)
            n_s, n_l = store["n_stations"], store["n_lines"]
            feedback = (
                f"✓ Imported {n_s} station(s) from {n_l} line(s) — "
                f"{filenames[0]}."
            )
            badge = f"{n_s} stations · {n_l} line(s)"
            return store, feedback, False, badge, "mv-data-badge visible"
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            return (
                no_update,
                f"✗ Error: {exc}",
                no_update,
                no_update,
                no_update,
            )
        finally:
            if tmpdir:
                shutil.rmtree(tmpdir, ignore_errors=True)
