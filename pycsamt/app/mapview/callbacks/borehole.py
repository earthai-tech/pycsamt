# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Borehole Studio — import, edit, and apply PCBH boreholes in Map View."""

from __future__ import annotations

import base64
import binascii
import io
import json
from datetime import datetime, timezone
from typing import Any

from dash import Input, Output, State, ctx, dcc, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app._borehole import (
    builder_draft_from_document,
    decode_pcbh_upload,
    document_from_builder_draft,
    document_from_store,
    strip_log_figure,
)

from .._ids import IDs

_COLLAR_COLS = (
    "id", "name", "kind", "status", "x", "y", "z", "total_depth_md",
)
_LAYER_COLS = (
    "borehole_id", "from_md", "to_md", "code", "label",
    "resistivity_ohm_m", "description",
)
_SURVEY_COLS = ("borehole_id", "md", "azimuth_deg", "inclination_deg")


def register_borehole(app) -> None:
    _register_accordion_upload(app)
    _register_modal_toggle(app)
    _register_imports(app)
    _register_xlsx(app)
    _register_points(app)
    _register_tables(app)
    _register_apply(app)
    _register_export(app)


def _register_points(app) -> None:
    @app.callback(
        Output(IDs.PCPT_STORE, "data"),
        Output(IDs.PCPT_IMPORT_INFO, "children"),
        Input(IDs.PCPT_IMPORT, "contents"),
        State(IDs.PCPT_IMPORT, "filename"),
        prevent_initial_call=True,
    )
    def load_points(contents, filename):
        if not contents:
            raise PreventUpdate
        from pycsamt.app._points import decode_points_upload

        try:
            store = decode_points_upload(contents, filename)
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error))
        return store, _ok(
            f"{store['filename']} — {store['n_points']} point(s)"
        )


# ---------------------------------------------------------------------------
# the slim upload in the 3-D "Boreholes (PCBH)" accordion
# ---------------------------------------------------------------------------


def _register_accordion_upload(app) -> None:
    @app.callback(
        Output(IDs.PCBH_STORE, "data", allow_duplicate=True),
        Output(IDs.PCBH_UPLOAD_INFO, "children"),
        Input(IDs.PCBH_UPLOAD, "contents"),
        State(IDs.PCBH_UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def load_pcbh(contents, filename):
        if not contents:
            raise PreventUpdate
        try:
            store = decode_pcbh_upload(contents, filename)
        except (TypeError, ValueError) as error:
            return no_update, html.Span(str(error), className="text-danger")
        return store, html.Span(
            f"{store['filename']} — {store['n_boreholes']} borehole(s)",
            className="text-success",
        )


# ---------------------------------------------------------------------------
# modal open / close
# ---------------------------------------------------------------------------


def _register_modal_toggle(app) -> None:
    @app.callback(
        Output(IDs.BH_STUDIO_MODAL, "is_open"),
        Input(IDs.BTN_BH_STUDIO, "n_clicks"),
        Input(IDs.BTN_BH_STUDIO_TB, "n_clicks"),
        State(IDs.BH_STUDIO_MODAL, "is_open"),
        prevent_initial_call=True,
    )
    def toggle(n_accordion, n_toolbar, is_open):
        if not (n_accordion or n_toolbar):
            raise PreventUpdate
        return not is_open


# ---------------------------------------------------------------------------
# imports (PCBH json / combined CSV)
# ---------------------------------------------------------------------------


def _register_imports(app) -> None:
    @app.callback(
        Output(IDs.PCBH_STORE, "data", allow_duplicate=True),
        Output(IDs.PCBH_DRAFT_STORE, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_COLLARS, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_LAYERS, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_SURVEY, "data", allow_duplicate=True),
        Output(IDs.BH_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.BH_IMPORT_PCBH, "contents"),
        Input(IDs.BH_IMPORT_CSV, "contents"),
        State(IDs.BH_IMPORT_PCBH, "filename"),
        State(IDs.BH_IMPORT_CSV, "filename"),
        prevent_initial_call=True,
    )
    def do_import(pcbh_contents, csv_contents, pcbh_name, csv_name):
        trigger = ctx.triggered_id
        try:
            if trigger == IDs.BH_IMPORT_PCBH and pcbh_contents:
                store = decode_pcbh_upload(pcbh_contents, pcbh_name)
                document = document_from_store(store)
            elif trigger == IDs.BH_IMPORT_CSV and csv_contents:
                from pycsamt.format.borehole import (
                    boreholes_from_csv,
                    pcbh_to_dict,
                )

                raw = _decode(csv_contents)
                path = _tmp_write(raw, csv_name or "boreholes.csv")
                document, _ = boreholes_from_csv(path, strict=False)
                store = {
                    "filename": csv_name,
                    "document": pcbh_to_dict(document),
                    "n_boreholes": len(document.boreholes),
                }
            else:
                raise PreventUpdate
        except PreventUpdate:
            raise
        except Exception as error:  # noqa: BLE001
            return (
                no_update, no_update, no_update, no_update, no_update,
                _err(str(error)),
            )
        draft = builder_draft_from_document(document)
        collars, layers, survey = _tables_from_draft(draft)
        return (
            store, draft, collars, layers, survey,
            _ok(f"Imported {len(document.boreholes)} borehole(s)"),
        )


# ---------------------------------------------------------------------------
# spreadsheet import
# ---------------------------------------------------------------------------


def _register_xlsx(app) -> None:
    @app.callback(
        Output(IDs.BH_XLSX_RAW, "data"),
        Output(IDs.BH_XLSX_SHEET, "options"),
        Output(IDs.BH_XLSX_SHEET, "value"),
        Output(IDs.BH_XLSX_PREVIEW, "children"),
        Input(IDs.BH_IMPORT_XLSX, "contents"),
        prevent_initial_call=True,
    )
    def inspect(contents):
        if not contents:
            raise PreventUpdate
        from pycsamt.format.borehole import inspect_workbook

        try:
            raw = _decode(contents)
            outline = inspect_workbook(raw)
        except Exception as error:  # noqa: BLE001
            return None, [], None, _err(str(error))
        options = [
            {"label": s.name, "value": s.name} for s in outline.sheets
        ]
        first = outline.sheets[0]
        preview = _preview_table(first.preview)
        payload = {"b64": base64.b64encode(raw).decode()}
        return payload, options, first.name, preview

    @app.callback(
        Output(IDs.BH_XLSX_MAP_ID, "options"),
        Output(IDs.BH_XLSX_MAP_LITH, "options"),
        Output(IDs.BH_XLSX_MAP_FROM, "options"),
        Output(IDs.BH_XLSX_MAP_TO, "options"),
        Output(IDs.BH_XLSX_MAP_ID, "value"),
        Output(IDs.BH_XLSX_MAP_LITH, "value"),
        Output(IDs.BH_XLSX_MAP_FROM, "value"),
        Output(IDs.BH_XLSX_MAP_TO, "value"),
        Output(IDs.BH_XLSX_PREVIEW, "children", allow_duplicate=True),
        Input(IDs.BH_XLSX_SHEET, "value"),
        Input(IDs.BH_XLSX_HEADER, "value"),
        State(IDs.BH_XLSX_RAW, "data"),
        prevent_initial_call=True,
    )
    def populate_columns(sheet, header, payload):
        if not payload or not sheet:
            raise PreventUpdate
        from pycsamt.format.borehole.xlsxio import sheet_rows

        try:
            raw = base64.b64decode(payload["b64"], validate=True)
            header_row = (
                int(header) - 1 if header not in (None, "") else None
            )
            _name, headers, rows = sheet_rows(
                raw, sheet=sheet, header_row=header_row
            )
        except Exception:  # noqa: BLE001
            raise PreventUpdate
        blank = {"label": "— auto —", "value": ""}
        options = [blank] + [
            {"label": h, "value": h} for h in headers
        ]
        guess = _guess_columns(headers)
        preview = _preview_table([headers, *rows[:10]])
        return (
            options, options, options, options,
            guess.get("id", ""),
            guess.get("lithology", ""),
            guess.get("from_md", ""),
            guess.get("to_md", ""),
            preview,
        )

    @app.callback(
        Output(IDs.PCBH_STORE, "data", allow_duplicate=True),
        Output(IDs.PCBH_DRAFT_STORE, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_COLLARS, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_LAYERS, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_SURVEY, "data", allow_duplicate=True),
        Output(IDs.BH_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.BH_XLSX_APPLY, "n_clicks"),
        State(IDs.BH_XLSX_RAW, "data"),
        State(IDs.BH_XLSX_SHEET, "value"),
        State(IDs.BH_XLSX_HEADER, "value"),
        State(IDs.BH_XLSX_COLLAR, "value"),
        State(IDs.BH_XLSX_MAP_ID, "value"),
        State(IDs.BH_XLSX_MAP_LITH, "value"),
        State(IDs.BH_XLSX_MAP_FROM, "value"),
        State(IDs.BH_XLSX_MAP_TO, "value"),
        prevent_initial_call=True,
    )
    def apply_sheet(
        n_clicks, payload, sheet, header, collar_text,
        map_id, map_lith, map_from, map_to,
    ):
        if not n_clicks or not payload:
            raise PreventUpdate
        from pycsamt.format.borehole import (
            boreholes_from_xlsx,
            pcbh_to_dict,
        )

        try:
            raw = base64.b64decode(payload["b64"], validate=True)
            collars = _parse_collar(collar_text)
            header_row = (
                int(header) - 1 if header not in (None, "") else None
            )
            constants: dict[str, Any] = {}
            if collars is None or not collars.get("crs"):
                constants["crs.horizontal"] = "LOCAL:studio-grid"
            columns = {
                canonical: source
                for canonical, source in (
                    ("borehole.id", map_id),
                    ("interval.lithology", map_lith),
                    ("interval.from_md", map_from),
                    ("interval.to_md", map_to),
                )
                if source
            } or None
            document, report = boreholes_from_xlsx(
                raw,
                sheet=sheet,
                header_row=header_row,
                columns=columns,
                constants=constants,
                collars=collars,
                default_borehole_id=_slug(sheet) or "HOLE-1",
                strict=False,
            )
        except Exception as error:  # noqa: BLE001
            return (
                no_update, no_update, no_update, no_update, no_update,
                _err(str(error)),
            )
        store = {
            "filename": f"{sheet}.xlsx",
            "document": pcbh_to_dict(document),
            "n_boreholes": len(document.boreholes),
        }
        draft = builder_draft_from_document(document)
        collars_t, layers_t, survey_t = _tables_from_draft(draft)
        msg = (
            f"Imported {len(document.boreholes)} hole(s), "
            f"{report.rows_accepted} interval(s)"
        )
        if report.rows_rejected:
            msg += f" — {report.rows_rejected} row(s) rejected"
        return store, draft, collars_t, layers_t, survey_t, _ok(msg)


# ---------------------------------------------------------------------------
# table editing -> draft -> validation + preview
# ---------------------------------------------------------------------------


def _register_tables(app) -> None:
    @app.callback(
        Output(IDs.PCBH_DRAFT_STORE, "data", allow_duplicate=True),
        Output(IDs.BH_STUDIO_VALIDATION, "children"),
        Output(IDs.BH_STUDIO_PREVIEW, "figure"),
        Input(IDs.BH_TABLE_COLLARS, "data"),
        Input(IDs.BH_TABLE_LAYERS, "data"),
        Input(IDs.BH_TABLE_SURVEY, "data"),
        State(IDs.PCBH_DRAFT_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def rebuild(collars, layers, survey, draft, theme):
        import plotly.graph_objects as go

        draft = _draft_from_tables(collars, layers, survey, draft)
        try:
            document = document_from_builder_draft(draft)
        except Exception as error:  # noqa: BLE001
            return draft, _err(str(error)), go.Figure()
        figure = strip_log_figure(document, theme=theme or "light")
        return (
            draft,
            _ok(f"Valid — {len(document.boreholes)} borehole(s)"),
            figure,
        )

    @app.callback(
        Output(IDs.BH_TABLE_COLLARS, "data", allow_duplicate=True),
        Output(IDs.BH_TABLE_LAYERS, "data", allow_duplicate=True),
        Input(IDs.BH_BTN_ADD_ROW, "n_clicks"),
        State(IDs.BH_STUDIO_TABS, "active_tab"),
        State(IDs.BH_TABLE_COLLARS, "data"),
        State(IDs.BH_TABLE_LAYERS, "data"),
        prevent_initial_call=True,
    )
    def add_row(n_clicks, active, collars, layers):
        if not n_clicks:
            raise PreventUpdate
        collars = list(collars or [])
        layers = list(layers or [])
        if active == "bh-layers":
            layers.append({key: None for key in _LAYER_COLS})
        else:
            collars.append({key: None for key in _COLLAR_COLS})
        return collars, layers


# ---------------------------------------------------------------------------
# apply to the map / 3-D view
# ---------------------------------------------------------------------------


def _register_apply(app) -> None:
    @app.callback(
        Output(IDs.PCBH_STORE, "data", allow_duplicate=True),
        Output(IDs.BH_STUDIO_MODAL, "is_open", allow_duplicate=True),
        Output(IDs.BH_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.BH_BTN_APPLY, "n_clicks"),
        State(IDs.PCBH_DRAFT_STORE, "data"),
        prevent_initial_call=True,
    )
    def apply(n_clicks, draft):
        if not n_clicks:
            raise PreventUpdate
        from pycsamt.format.borehole import pcbh_to_dict

        try:
            document = document_from_builder_draft(draft)
        except Exception as error:  # noqa: BLE001
            return no_update, no_update, _err(str(error))
        if document is None:
            return no_update, no_update, _err("nothing to apply")
        store = {
            "filename": "studio.pcbh.json",
            "document": pcbh_to_dict(document),
            "n_boreholes": len(document.boreholes),
        }
        return store, False, _ok("Applied to view")


def _register_export(app) -> None:
    @app.callback(
        Output(IDs.BH_EXPORT_DL, "data"),
        Input(IDs.BH_BTN_EXPORT, "n_clicks"),
        State(IDs.PCBH_DRAFT_STORE, "data"),
        State(IDs.PCBH_STORE, "data"),
        prevent_initial_call=True,
    )
    def export(n_clicks, draft, store):
        if not n_clicks:
            raise PreventUpdate
        from pycsamt.format.borehole import pcbh_to_dict

        document = None
        try:
            document = document_from_builder_draft(draft)
        except Exception:  # noqa: BLE001
            document = document_from_store(store)
        if document is None:
            raise PreventUpdate
        payload = json.dumps(
            pcbh_to_dict(document), ensure_ascii=False, indent=2
        )
        name = str(document.document_id).replace(":", "-") or "boreholes"
        return dcc.send_string(payload + "\n", f"{name}.pcbh.json")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _decode(contents: str) -> bytes:
    try:
        return base64.b64decode(contents.split(",", 1)[1], validate=True)
    except (binascii.Error, IndexError, ValueError) as error:
        raise ValueError("could not decode the uploaded file") from error


def _tmp_write(raw: bytes, name: str) -> str:
    import tempfile
    from pathlib import Path

    folder = tempfile.mkdtemp(prefix="pcbh-")
    path = Path(folder) / name
    path.write_bytes(raw)
    return str(path)


def _slug(text: str | None) -> str:
    import re

    return re.sub(r"[^A-Za-z0-9_-]+", "-", str(text or "")).strip("-")[:40]


def _guess_columns(headers: list[str]) -> dict[str, str]:
    """Alias-match id / lithology / from / to against sheet headers."""
    from pycsamt.format.borehole import ImportReport
    from pycsamt.format.borehole.mapping import resolve_csv_columns

    report = ImportReport(
        source="", source_sha256="", delimiter="xlsx", strict=False
    )
    resolved = resolve_csv_columns(
        headers, explicit=None, constants={}, report=report
    )
    out: dict[str, str] = {}
    for canonical, key in (
        ("borehole.id", "id"),
        ("interval.lithology", "lithology"),
        ("interval.from_md", "from_md"),
        ("interval.to_md", "to_md"),
    ):
        if canonical in resolved:
            out[key] = resolved[canonical]
    return out


def _parse_collar(text: str | None):
    if not text or not str(text).strip():
        return None
    parts = [p for p in str(text).replace(";", ",").split(",") if p.strip()]
    try:
        values = [float(p) for p in parts]
    except ValueError:
        return None
    if len(values) < 3:
        return None
    return {"x": values[0], "y": values[1], "z": values[2]}


def _tables_from_draft(draft: dict[str, Any] | None):
    draft = draft or {}
    collars = [
        {key: row.get(key) for key in _COLLAR_COLS}
        for row in draft.get("boreholes", [])
    ]
    layers = [
        {key: row.get(key) for key in _LAYER_COLS}
        for row in draft.get("intervals", [])
        if row.get("family", "lithology") == "lithology"
    ]
    survey = [
        {key: row.get(key) for key in _SURVEY_COLS}
        for row in draft.get("surveys", [])
    ]
    return collars, layers, survey


def _draft_from_tables(collars, layers, survey, previous):
    draft = dict(previous or {})
    project = dict(draft.get("project") or {})
    project.setdefault("document_id", "pcbh:studio")
    project.setdefault("created_by", "pyCSAMT Borehole Studio")
    project.setdefault("crs_horizontal", "LOCAL:studio-grid")
    project.setdefault("crs_vertical", "unknown")
    project.setdefault("coordinate_unit", "m")
    project.setdefault("depth_unit", "m")
    project.setdefault("diameter_unit", "m")
    draft["project"] = project
    draft["boreholes"] = [
        {**{k: r.get(k) for k in _COLLAR_COLS}, "north_reference": "true"}
        for r in (collars or [])
        if r.get("id")
    ]
    draft["intervals"] = [
        {**{k: r.get(k) for k in _LAYER_COLS}, "family": "lithology"}
        for r in (layers or [])
        if r.get("borehole_id") and r.get("from_md") is not None
    ]
    draft["surveys"] = [
        {k: r.get(k) for k in _SURVEY_COLS}
        for r in (survey or [])
        if r.get("borehole_id")
    ]
    for key in (
        "structures", "water", "construction", "samples", "assays",
        "mapping_profiles",
    ):
        draft.setdefault(key, [])
    draft["created_at"] = (
        datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    )
    return draft


def _preview_table(rows: list[list[str]]):
    from dash import dash_table

    if not rows:
        return html.Span("empty sheet", className="text-muted")
    width = max(len(r) for r in rows)
    columns = [{"name": str(i), "id": str(i)} for i in range(width)]
    data = [
        {str(i): (row[i] if i < len(row) else "") for i in range(width)}
        for row in rows[:12]
    ]
    return dash_table.DataTable(
        columns=columns,
        data=data,
        page_action="none",
        style_table={"overflowX": "auto", "maxHeight": "220px"},
        style_cell={"fontSize": "10px", "padding": "2px", "maxWidth": "120px"},
    )


def _ok(text: str):
    return html.Span(text, className="text-success")


def _err(text: str):
    return html.Span(text, className="text-danger")


# retained for import compatibility with older call sites
def _bytes_io(raw: bytes) -> io.BytesIO:
    return io.BytesIO(raw)
