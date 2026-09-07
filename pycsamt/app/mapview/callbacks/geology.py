# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Interpretation Studio — import, auto-suggest, edit, and apply a PCGL
geology legend to Map View's 3-D block / fence / depth / iso views."""

from __future__ import annotations

import json

from dash import Input, Output, State, dcc, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app._borehole import document_from_store
from pycsamt.app._geology import (
    TABLE_COLUMNS,
    auto_suggest_legend,
    borehole_lithology_colors,
    decode_geology_upload,
    legend_from_store,
    legend_from_table_rows,
    legend_preview_figure,
    store_from_legend,
    sync_legend_colors_with_boreholes,
    table_rows_from_legend,
)
from pycsamt.app._structure import (
    FAULT_COLUMNS,
    LINEAR_COLUMNS,
    PLANAR_COLUMNS,
    decode_structure_csv_upload,
    decode_structure_json_upload,
    model_from_table_rows,
    store_from_structure,
    structure_from_store,
    structure_section_figure,
    table_rows_from_model,
)
from pycsamt.format.geology import GeologyLegend, legend_to_dict
from pycsamt.format.structure import StructModel, structure_to_dict

from .._ids import IDs
from ..cache import get_view


def register_geology(app) -> None:
    _register_modal_toggle(app)
    _register_imports(app)
    _register_quick_upload(app)
    _register_default(app)
    _register_auto_suggest(app)
    _register_table(app)
    _register_add_row(app)
    _register_sync_colors(app)
    _register_apply(app)
    _register_export(app)
    _register_strip(app)
    _register_structure_imports(app)
    _register_structure_table(app)
    _register_structure_add_rows(app)
    _register_structure_apply(app)
    _register_structure_export(app)
    _register_structure_strip(app)


# ---------------------------------------------------------------------------
# modal open / close
# ---------------------------------------------------------------------------


def _register_modal_toggle(app) -> None:
    @app.callback(
        Output(IDs.GEO_STUDIO_MODAL, "is_open"),
        Input(IDs.BTN_GEO_STUDIO, "n_clicks"),
        State(IDs.GEO_STUDIO_MODAL, "is_open"),
        prevent_initial_call=True,
    )
    def toggle(n_clicks, is_open):
        if not n_clicks:
            raise PreventUpdate
        return not is_open


# ---------------------------------------------------------------------------
# imports (PCGL json / legend CSV)
# ---------------------------------------------------------------------------


def _register_imports(app) -> None:
    @app.callback(
        Output(IDs.GEO_TABLE_LEGEND, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_IMPORT_PCGL, "contents"),
        Input(IDs.GEO_IMPORT_CSV, "contents"),
        State(IDs.GEO_IMPORT_PCGL, "filename"),
        State(IDs.GEO_IMPORT_CSV, "filename"),
        prevent_initial_call=True,
    )
    def do_import(pcgl_contents, csv_contents, pcgl_name, csv_name):
        from dash import ctx

        trigger = ctx.triggered_id
        contents = pcgl_contents if trigger == IDs.GEO_IMPORT_PCGL \
            else csv_contents
        name = pcgl_name if trigger == IDs.GEO_IMPORT_PCGL else csv_name
        if not contents:
            raise PreventUpdate
        try:
            store = decode_geology_upload(contents, name)
            legend = legend_from_store(store)
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error))
        return (
            table_rows_from_legend(legend),
            _ok(f"Imported {len(legend.entries)} entrie(s)"),
        )


def _register_quick_upload(app) -> None:
    """The slim direct-to-view upload in the Geology rail panel (GRP_GEO)
    -- mirrors PCBH_UPLOAD's role for boreholes: skips the Studio table
    entirely and writes GEO_STORE straight away."""

    @app.callback(
        Output(IDs.GEO_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_QUICK_UPLOAD_INFO, "children"),
        Input(IDs.GEO_QUICK_UPLOAD, "contents"),
        State(IDs.GEO_QUICK_UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def quick_upload(contents, filename):
        if not contents:
            raise PreventUpdate
        try:
            store = decode_geology_upload(contents, filename)
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error))
        return store, _ok(
            f"{store['filename']} — {store['n_entries']} entrie(s)"
        )


def _register_default(app) -> None:
    @app.callback(
        Output(IDs.GEO_TABLE_LEGEND, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_LOAD_DEFAULT, "n_clicks"),
        prevent_initial_call=True,
    )
    def load_default(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        legend = GeologyLegend.from_rock_database(title="Default rock database")
        return (
            table_rows_from_legend(legend),
            _ok(f"Loaded {len(legend.entries)} default entrie(s)"),
        )


def _register_auto_suggest(app) -> None:
    @app.callback(
        Output(IDs.GEO_TABLE_LEGEND, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_AUTO_SUGGEST, "n_clicks"),
        State(IDs.GEO_AUTO_NBINS, "value"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def auto_suggest(n_clicks, n_bins, session_id):
        if not n_clicks:
            raise PreventUpdate
        from pycsamt.map import all_resistivity_values

        view = get_view(session_id) if session_id else None
        if view is None or getattr(view, "data", None) is None:
            return no_update, _err(
                "load survey lines or an inversion result first"
            )
        rho = all_resistivity_values(view.data)
        try:
            legend = auto_suggest_legend(
                rho, n_bins=int(n_bins) if n_bins else 8
            )
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error))
        return (
            table_rows_from_legend(legend),
            _ok(f"Suggested {len(legend.entries)} entrie(s)"),
        )


# ---------------------------------------------------------------------------
# table editing -> validation + preview
# ---------------------------------------------------------------------------


def _register_table(app) -> None:
    @app.callback(
        Output(IDs.GEO_DRAFT_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_VALIDATION, "children"),
        Output(IDs.GEO_STUDIO_PREVIEW, "figure"),
        Input(IDs.GEO_TABLE_LEGEND, "data"),
        prevent_initial_call=True,
    )
    def rebuild(rows):
        import plotly.graph_objects as go

        try:
            legend = legend_from_table_rows(rows)
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error)), go.Figure()
        return (
            legend_to_dict(legend),
            _ok(f"Valid — {len(legend.entries)} entrie(s)"),
            legend_preview_figure(legend),
        )


def _register_add_row(app) -> None:
    @app.callback(
        Output(IDs.GEO_TABLE_LEGEND, "data", allow_duplicate=True),
        Input(IDs.GEO_BTN_ADD_ROW, "n_clicks"),
        State(IDs.GEO_TABLE_LEGEND, "data"),
        prevent_initial_call=True,
    )
    def add_row(n_clicks, rows):
        if not n_clicks:
            raise PreventUpdate
        rows = list(rows or [])
        rows.append({key: None for key in TABLE_COLUMNS})
        return rows


def _register_sync_colors(app) -> None:
    """"Sync colours with loaded boreholes" -- overwrite each Legend
    row's colour with the matching (exact-name) borehole lithology
    colour, so the same unit reads identically in both overlays. See
    ``pycsamt.app._geology.sync_legend_colors_with_boreholes``."""

    @app.callback(
        Output(IDs.GEO_TABLE_LEGEND, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_SYNC_COLORS, "n_clicks"),
        State(IDs.GEO_TABLE_LEGEND, "data"),
        State(IDs.PCBH_STORE, "data"),
        prevent_initial_call=True,
    )
    def sync_colors(n_clicks, rows, pcbh_store):
        if not n_clicks:
            raise PreventUpdate
        document = document_from_store(pcbh_store) if pcbh_store else None
        if document is None:
            return no_update, _err(
                "no borehole loaded — import one in the Boreholes "
                "section first"
            )
        colors = borehole_lithology_colors(document)
        if not colors:
            return no_update, _err(
                "the loaded borehole(s) have no lithology log to match"
            )
        updated, n_matched = sync_legend_colors_with_boreholes(rows, colors)
        if n_matched == 0:
            return updated, _err(
                "no legend row name matched a borehole lithology exactly"
            )
        return updated, _ok(
            f"Synced {n_matched} row colour(s) to the loaded borehole(s)"
        )


# ---------------------------------------------------------------------------
# apply to the 3-D view
# ---------------------------------------------------------------------------


def _register_apply(app) -> None:
    @app.callback(
        Output(IDs.GEO_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_MODAL, "is_open", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_APPLY, "n_clicks"),
        State(IDs.GEO_TABLE_LEGEND, "data"),
        prevent_initial_call=True,
    )
    def apply(n_clicks, rows):
        if not n_clicks:
            raise PreventUpdate
        try:
            legend = legend_from_table_rows(rows)
        except Exception as error:  # noqa: BLE001
            return no_update, no_update, _err(str(error))
        store = store_from_legend(legend, filename="studio.pcgl.json")
        return store, False, _ok(
            "Applied — turn on “Apply geology legend to this view” in "
            "the Interpretation accordion"
        )


def _register_export(app) -> None:
    @app.callback(
        Output(IDs.GEO_EXPORT_DL, "data"),
        Input(IDs.GEO_BTN_EXPORT, "n_clicks"),
        State(IDs.GEO_TABLE_LEGEND, "data"),
        prevent_initial_call=True,
    )
    def export(n_clicks, rows):
        if not n_clicks:
            raise PreventUpdate
        try:
            legend = legend_from_table_rows(rows)
        except Exception:  # noqa: BLE001
            raise PreventUpdate
        payload = json.dumps(legend_to_dict(legend), ensure_ascii=False, indent=2)
        name = str(legend.document_id).replace(":", "-") or "legend"
        return dcc.send_string(payload + "\n", f"{name}.pcgl.json")


# ---------------------------------------------------------------------------
# accordion mini legend-strip preview
# ---------------------------------------------------------------------------


def _register_strip(app) -> None:
    @app.callback(
        Output(IDs.GEO_STRIP, "children"),
        Input(IDs.GEO_STORE, "data"),
        prevent_initial_call=False,
    )
    def strip(store):
        legend = legend_from_store(store) if store else None
        if legend is None or not legend.entries:
            return html.Span("No legend applied", className="text-muted")
        chips = [
            html.Span(
                entry.name,
                title=f"{entry.rho_min:g}–{entry.rho_max:g} Ω·m",
                style={
                    "backgroundColor": entry.color,
                    "display": "inline-block",
                    "padding": "1px 6px",
                    "borderRadius": "3px",
                    "fontSize": "10px",
                    "marginRight": "3px",
                    "marginBottom": "3px",
                    "color": "#111",
                },
            )
            for entry in legend.entries
        ]
        return html.Div(chips)


# ---------------------------------------------------------------------------
# structural geology (PCGS) — import
# ---------------------------------------------------------------------------


def _register_structure_imports(app) -> None:
    @app.callback(
        Output(IDs.GEO_TABLE_PLANAR, "data", allow_duplicate=True),
        Output(IDs.GEO_TABLE_LINEAR, "data", allow_duplicate=True),
        Output(IDs.GEO_TABLE_FAULTS, "data", allow_duplicate=True),
        Output(IDs.GEO_STRUCT_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_STRUCT_IMPORT_JSON, "contents"),
        Input(IDs.GEO_STRUCT_IMPORT_PLANAR, "contents"),
        Input(IDs.GEO_STRUCT_IMPORT_LINEAR, "contents"),
        Input(IDs.GEO_STRUCT_IMPORT_FAULTS, "contents"),
        State(IDs.GEO_STRUCT_IMPORT_JSON, "filename"),
        State(IDs.GEO_STRUCT_IMPORT_PLANAR, "filename"),
        State(IDs.GEO_STRUCT_IMPORT_LINEAR, "filename"),
        State(IDs.GEO_STRUCT_IMPORT_FAULTS, "filename"),
        State(IDs.GEO_TABLE_PLANAR, "data"),
        State(IDs.GEO_TABLE_LINEAR, "data"),
        State(IDs.GEO_TABLE_FAULTS, "data"),
        prevent_initial_call=True,
    )
    def do_import_structure(
        json_c, planar_c, linear_c, faults_c,
        json_n, planar_n, linear_n, faults_n,
        planar_rows, linear_rows, faults_rows,
    ):
        from dash import ctx

        trigger = ctx.triggered_id
        try:
            if trigger == IDs.GEO_STRUCT_IMPORT_JSON and json_c:
                store = decode_structure_json_upload(json_c, json_n)
                model = structure_from_store(store)
                rows = table_rows_from_model(model.model if model else None)
                n = (model and len(model)) or 0
                return (
                    rows["planar"], rows["linear"], rows["faults"],
                    _ok(f"Imported {n} item(s)"),
                )
            if trigger == IDs.GEO_STRUCT_IMPORT_PLANAR and planar_c:
                rows = decode_structure_csv_upload(
                    planar_c, planar_n, kind="planar"
                )
                return (
                    rows, no_update, no_update,
                    _ok(f"Imported {len(rows)} planar measurement(s)"),
                )
            if trigger == IDs.GEO_STRUCT_IMPORT_LINEAR and linear_c:
                rows = decode_structure_csv_upload(
                    linear_c, linear_n, kind="linear"
                )
                return (
                    no_update, rows, no_update,
                    _ok(f"Imported {len(rows)} linear measurement(s)"),
                )
            if trigger == IDs.GEO_STRUCT_IMPORT_FAULTS and faults_c:
                rows = decode_structure_csv_upload(
                    faults_c, faults_n, kind="faults"
                )
                return (
                    no_update, no_update, rows,
                    _ok(f"Imported {len(rows)} fault trace(s)"),
                )
            raise PreventUpdate
        except PreventUpdate:
            raise
        except Exception as error:  # noqa: BLE001
            return no_update, no_update, no_update, _err(str(error))


# ---------------------------------------------------------------------------
# structural geology — table editing -> validation + preview
# ---------------------------------------------------------------------------


def _register_structure_table(app) -> None:
    @app.callback(
        Output(IDs.STRUCT_DRAFT_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_STRUCT_VALIDATION, "children"),
        Output(IDs.GEO_STRUCT_PREVIEW, "figure"),
        Input(IDs.GEO_TABLE_PLANAR, "data"),
        Input(IDs.GEO_TABLE_LINEAR, "data"),
        Input(IDs.GEO_TABLE_FAULTS, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def rebuild_structure(planar_rows, linear_rows, fault_rows, theme):
        import plotly.graph_objects as go

        try:
            model = model_from_table_rows(planar_rows, linear_rows, fault_rows)
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error)), go.Figure()
        doc = StructModel.from_structural_model(model)
        n = len(doc)
        return (
            structure_to_dict(doc),
            _ok(f"Valid — {n} item(s)"),
            structure_section_figure(model, theme=theme or "light"),
        )


def _register_structure_add_rows(app) -> None:
    def _add(columns, name):
        def _cb(n_clicks, rows):
            if not n_clicks:
                raise PreventUpdate
            rows = list(rows or [])
            rows.append({key: None for key in columns})
            return rows

        _cb.__name__ = name
        return _cb

    app.callback(
        Output(IDs.GEO_TABLE_PLANAR, "data", allow_duplicate=True),
        Input(IDs.GEO_BTN_ADD_PLANAR, "n_clicks"),
        State(IDs.GEO_TABLE_PLANAR, "data"),
        prevent_initial_call=True,
    )(_add(PLANAR_COLUMNS, "add_planar_row"))
    app.callback(
        Output(IDs.GEO_TABLE_LINEAR, "data", allow_duplicate=True),
        Input(IDs.GEO_BTN_ADD_LINEAR, "n_clicks"),
        State(IDs.GEO_TABLE_LINEAR, "data"),
        prevent_initial_call=True,
    )(_add(LINEAR_COLUMNS, "add_linear_row"))
    app.callback(
        Output(IDs.GEO_TABLE_FAULTS, "data", allow_duplicate=True),
        Input(IDs.GEO_BTN_ADD_FAULT, "n_clicks"),
        State(IDs.GEO_TABLE_FAULTS, "data"),
        prevent_initial_call=True,
    )(_add(FAULT_COLUMNS, "add_fault_row"))


# ---------------------------------------------------------------------------
# structural geology — apply / export
# ---------------------------------------------------------------------------


def _register_structure_apply(app) -> None:
    @app.callback(
        Output(IDs.STRUCT_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_STUDIO_MODAL, "is_open", allow_duplicate=True),
        Output(IDs.GEO_STRUCT_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_APPLY_STRUCT, "n_clicks"),
        State(IDs.GEO_TABLE_PLANAR, "data"),
        State(IDs.GEO_TABLE_LINEAR, "data"),
        State(IDs.GEO_TABLE_FAULTS, "data"),
        prevent_initial_call=True,
    )
    def apply_structure(n_clicks, planar_rows, linear_rows, fault_rows):
        if not n_clicks:
            raise PreventUpdate
        try:
            model = model_from_table_rows(planar_rows, linear_rows, fault_rows)
        except Exception as error:  # noqa: BLE001
            return no_update, no_update, _err(str(error))
        doc = StructModel.from_structural_model(
            model, document_id="pcgs:studio"
        )
        store = store_from_structure(doc, filename="studio.pcgs.json")
        return store, False, _ok(
            "Applied — turn on “Apply structure” in the Interpretation "
            "accordion"
        )


def _register_structure_export(app) -> None:
    @app.callback(
        Output(IDs.GEO_STRUCT_EXPORT_DL, "data"),
        Input(IDs.GEO_BTN_STRUCT_EXPORT, "n_clicks"),
        State(IDs.GEO_TABLE_PLANAR, "data"),
        State(IDs.GEO_TABLE_LINEAR, "data"),
        State(IDs.GEO_TABLE_FAULTS, "data"),
        prevent_initial_call=True,
    )
    def export_structure(n_clicks, planar_rows, linear_rows, fault_rows):
        if not n_clicks:
            raise PreventUpdate
        try:
            model = model_from_table_rows(planar_rows, linear_rows, fault_rows)
        except Exception:  # noqa: BLE001
            raise PreventUpdate
        doc = StructModel.from_structural_model(
            model, document_id="pcgs:studio"
        )
        payload = json.dumps(
            structure_to_dict(doc), ensure_ascii=False, indent=2
        )
        return dcc.send_string(payload + "\n", "structure.pcgs.json")


# ---------------------------------------------------------------------------
# accordion / panel mini structure-status strip
# ---------------------------------------------------------------------------


def _register_structure_strip(app) -> None:
    @app.callback(
        Output(IDs.STRUCT_STRIP, "children"),
        Input(IDs.STRUCT_STORE, "data"),
        prevent_initial_call=False,
    )
    def strip_structure(store):
        model = structure_from_store(store) if store else None
        if model is None or len(model) == 0:
            return html.Span("No structure applied", className="text-muted")
        return html.Span(
            f"{len(model.faults)} fault(s), {len(model.planar)} planar, "
            f"{len(model.linear)} linear"
        )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _ok(text: str):
    return html.Span(text, className="text-success")


def _err(text: str):
    return html.Span(text, className="text-danger")
