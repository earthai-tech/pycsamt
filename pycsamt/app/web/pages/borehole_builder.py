# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Accessible PCBH builder/import page."""

from __future__ import annotations

import base64
import json
import tempfile
from pathlib import Path

import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash import (
    ALL,
    Input,
    Output,
    State,
    ctx,
    dash_table,
    dcc,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app._borehole import pcbh_plotly_traces
from pycsamt.format import read_pcsf, write_pcsf
from pycsamt.format.borehole import (
    CSVMappingProfile,
    boreholes_from_csv,
    document_to_builder,
    embed_pcbh,
    new_builder_draft,
    pcbh_to_dict,
    preview_csv_mapping,
    validate_builder,
)

PAGE_ID = "borehole-builder"

_DRAFT = "pcbh-builder-draft"
_DOCUMENT = "pcbh-builder-document"
_PROFILES = "pcbh-builder-profiles"
_CSV_PREVIEW = "pcbh-builder-csv-preview"
_VALIDATION = "pcbh-builder-validation"
_LOG_GRAPH = "pcbh-builder-log-graph"
_TRAJECTORY_GRAPH = "pcbh-builder-trajectory-graph"

_TABLES = {
    "boreholes": [
        "id",
        "name",
        "kind",
        "status",
        "x",
        "y",
        "z",
        "total_depth_md",
        "diameter",
        "north_reference",
    ],
    "surveys": ["borehole_id", "md", "azimuth_deg", "inclination_deg"],
    "intervals": [
        "borehole_id",
        "family",
        "from_md",
        "to_md",
        "code",
        "label",
        "color",
        "data_nature",
        "description",
    ],
    "structures": [
        "borehole_id",
        "kind",
        "at_md",
        "from_md",
        "to_md",
        "orientation_representation",
        "dip_deg",
        "dip_direction_deg",
    ],
    "water": ["borehole_id", "kind", "at_md", "from_md", "to_md", "value"],
    "construction": [
        "borehole_id",
        "kind",
        "from_md",
        "to_md",
        "diameter",
        "material",
    ],
    "samples": ["borehole_id", "sample_id", "from_md", "to_md", "type"],
    "assays": ["borehole_id", "sample_id", "analyte", "value", "unit"],
}


def _table(name: str):
    return dash_table.DataTable(
        id=f"pcbh-builder-{name}",
        columns=[
            {"name": item.replace("_", " ").title(), "id": item}
            for item in _TABLES[name]
        ],
        data=[],
        editable=True,
        row_deletable=True,
        virtualization=True,
        page_action="none",
        fixed_rows={"headers": True},
        style_table={"height": "330px", "overflowY": "auto"},
        style_cell={"minWidth": "110px", "fontSize": "12px"},
        tooltip_header={item: item for item in _TABLES[name]},
    )


def _project_panel():
    fields = [
        ("Document ID", "document_id", "pcbh:new-project"),
        ("Created by", "created_by", "pyCSAMT Borehole Builder"),
        ("Horizontal CRS", "crs_horizontal", "EPSG:4326"),
        ("Vertical CRS", "crs_vertical", "unknown"),
        ("Coordinate unit", "coordinate_unit", "m"),
        ("Depth unit", "depth_unit", "m"),
        ("Diameter unit", "diameter_unit", "m"),
    ]
    return dbc.Row(
        [
            dbc.Col(
                [
                    dbc.Label(label, html_for=f"pcbh-builder-project-{key}"),
                    dbc.Input(
                        id=f"pcbh-builder-project-{key}",
                        value=value,
                        debounce=True,
                    ),
                ],
                md=4,
                className="mb-2",
            )
            for label, key, value in fields
        ]
    )


def _csv_panel():
    return html.Div(
        [
            dcc.Upload(
                id="pcbh-builder-csv-upload",
                children=html.Div("Drop or browse a combined borehole CSV"),
                accept=".csv,.tsv,.txt",
                multiple=False,
                style={
                    "border": "1px dashed #7f849c",
                    "padding": "18px",
                    "textAlign": "center",
                    "borderRadius": "6px",
                },
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Dropdown(
                            id="pcbh-builder-profile-select",
                            placeholder="Use a saved mapping profile",
                            clearable=True,
                        )
                    ),
                    dbc.Col(
                        dbc.Input(
                            id="pcbh-builder-profile-name",
                            placeholder="Mapping profile name",
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            "Save mapping profile",
                            id="pcbh-builder-save-profile",
                            title="Save inferred CSV mapping in this browser",
                        ),
                        width="auto",
                    ),
                ],
                className="my-2",
            ),
            html.Div(id="pcbh-builder-csv-status", role="status"),
            dash_table.DataTable(
                id="pcbh-builder-csv-table",
                page_size=10,
                style_table={"overflowX": "auto"},
            ),
        ]
    )


def layout() -> html.Div:
    tabs = [
        dbc.Tab(_project_panel(), label="Project / CRS", tab_id="project"),
        dbc.Tab(_table("boreholes"), label="Collars", tab_id="boreholes"),
        dbc.Tab(_table("surveys"), label="Survey", tab_id="surveys"),
        dbc.Tab(_table("intervals"), label="Intervals", tab_id="intervals"),
        dbc.Tab(_table("structures"), label="Structures", tab_id="structures"),
        dbc.Tab(_table("water"), label="Water", tab_id="water"),
        dbc.Tab(
            _table("construction"), label="Construction", tab_id="construction"
        ),
        dbc.Tab(_table("samples"), label="Samples", tab_id="samples"),
        dbc.Tab(_table("assays"), label="Assays", tab_id="assays"),
        dbc.Tab(_csv_panel(), label="CSV import", tab_id="csv"),
    ]
    return html.Div(
        [
            dcc.Store(
                id=_DRAFT, storage_type="local", data=new_builder_draft()
            ),
            dcc.Store(id=_DOCUMENT, storage_type="memory"),
            dcc.Store(id=_PROFILES, storage_type="local", data={}),
            dcc.Store(id=_CSV_PREVIEW, storage_type="memory"),
            dcc.Download(id="pcbh-builder-download"),
            dcc.Download(id="pcbh-builder-pcsf-download"),
            html.Div(
                [
                    html.Div(
                        [
                            html.H4("Borehole Builder"),
                            html.Small("Canonical PCBH · autosaved locally"),
                        ]
                    ),
                    html.Div(
                        [
                            dbc.Button(
                                "Add row",
                                id="pcbh-builder-add-row",
                                title="Add a row to the active editor",
                            ),
                            dbc.Button(
                                "Restore autosave",
                                id="pcbh-builder-restore",
                                title=(
                                    "Restore the locally saved builder draft"
                                ),
                            ),
                            dbc.Button(
                                "Export PCBH",
                                id="pcbh-builder-export",
                                color="primary",
                            ),
                            dcc.Upload(
                                id="pcbh-builder-pcsf-upload",
                                children=dbc.Button(
                                    "Attach to PCSF", color="secondary"
                                ),
                                accept=".pcsf",
                                multiple=False,
                            ),
                        ],
                        className="d-flex gap-2",
                    ),
                ],
                className=(
                    "d-flex justify-content-between align-items-center p-3"
                ),
            ),
            dbc.Tabs(tabs, id="pcbh-builder-tabs", active_tab="project"),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H6("Validation"),
                            html.Div(id=_VALIDATION, role="alert"),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        dcc.Graph(id=_LOG_GRAPH, figure=go.Figure()), md=4
                    ),
                    dbc.Col(
                        dcc.Graph(id=_TRAJECTORY_GRAPH, figure=go.Figure()),
                        md=4,
                    ),
                ],
                className="p-2 flex-grow-1",
            ),
        ],
        style={
            "display": "flex",
            "flexDirection": "column",
            "height": "100%",
            "overflow": "auto",
        },
    )


def register_callbacks(app) -> None:
    _register_add_row(app)
    _register_restore(app)
    _register_builder(app)
    _register_issue_navigation(app)
    _register_csv(app)
    _register_profile(app)
    _register_export(app)
    _register_pcsf(app)


def _register_add_row(app):
    outputs = [Output(f"pcbh-builder-{name}", "data") for name in _TABLES]

    @app.callback(
        outputs,
        Input("pcbh-builder-add-row", "n_clicks"),
        State("pcbh-builder-tabs", "active_tab"),
        *[State(f"pcbh-builder-{name}", "data") for name in _TABLES],
        prevent_initial_call=True,
    )
    def add_row(n_clicks, active, *tables):
        if not n_clicks or active not in _TABLES:
            raise PreventUpdate
        result = [list(value or []) for value in tables]
        index = list(_TABLES).index(active)
        result[index].append({key: None for key in _TABLES[active]})
        return result


def _register_restore(app):
    project_keys = [
        "document_id",
        "created_by",
        "crs_horizontal",
        "crs_vertical",
        "coordinate_unit",
        "depth_unit",
        "diameter_unit",
    ]

    @app.callback(
        *[
            Output(
                f"pcbh-builder-project-{key}",
                "value",
                allow_duplicate=True,
            )
            for key in project_keys
        ],
        *[
            Output(f"pcbh-builder-{name}", "data", allow_duplicate=True)
            for name in _TABLES
        ],
        Input("pcbh-builder-restore", "n_clicks"),
        State(_DRAFT, "data"),
        prevent_initial_call=True,
    )
    def restore(n_clicks, draft):
        if not n_clicks or not draft:
            raise PreventUpdate
        project = draft.get("project", {})
        return (
            *[project.get(key) for key in project_keys],
            *[draft.get(name, []) for name in _TABLES],
        )


def _register_builder(app):
    project_keys = [
        "document_id",
        "created_by",
        "crs_horizontal",
        "crs_vertical",
        "coordinate_unit",
        "depth_unit",
        "diameter_unit",
    ]

    @app.callback(
        Output(_DRAFT, "data"),
        Output(_DOCUMENT, "data"),
        Output(_VALIDATION, "children"),
        Output(_LOG_GRAPH, "figure"),
        Output(_TRAJECTORY_GRAPH, "figure"),
        *[
            Input(f"pcbh-builder-project-{key}", "value")
            for key in project_keys
        ],
        *[Input(f"pcbh-builder-{name}", "data") for name in _TABLES],
        State(_DRAFT, "data"),
        prevent_initial_call=True,
    )
    def build(*values):
        project_values = values[: len(project_keys)]
        table_values = values[
            len(project_keys) : len(project_keys) + len(_TABLES)
        ]
        previous = values[-1] or new_builder_draft()
        draft = dict(previous)
        draft["project"] = dict(zip(project_keys, project_values))
        for name, rows in zip(_TABLES, table_values):
            draft[name] = rows or []
        validation = validate_builder(draft)
        if not validation.ok:
            messages = [
                html.Button(
                    _diagnostic_label(item),
                    id={
                        "type": "pcbh-builder-issue",
                        "editor": item.editor,
                        "row": item.row or 0,
                    },
                    className=(
                        "btn btn-link text-danger text-start p-0 d-block"
                    ),
                    title=item.path,
                )
                for item in validation.diagnostics[:100]
            ]
            return draft, None, messages, go.Figure(), go.Figure()
        document = validation.document
        canonical = pcbh_to_dict(document)
        log_figure = _log_figure(document)
        trajectory = go.Figure()
        for trace in pcbh_plotly_traces(document):
            trajectory.add_trace(trace)
        trajectory.update_layout(title="3-D trajectory preview", height=350)
        message = dbc.Alert(
            f"Valid PCBH · {len(document.boreholes)} borehole(s)",
            color="success",
        )
        return draft, canonical, message, log_figure, trajectory


def _register_issue_navigation(app):
    @app.callback(
        Output("pcbh-builder-tabs", "active_tab", allow_duplicate=True),
        *[
            Output(
                f"pcbh-builder-{name}",
                "active_cell",
                allow_duplicate=True,
            )
            for name in _TABLES
        ],
        Input(
            {"type": "pcbh-builder-issue", "editor": ALL, "row": ALL},
            "n_clicks",
        ),
        prevent_initial_call=True,
    )
    def navigate(_clicks):
        target = ctx.triggered_id
        if not isinstance(target, dict):
            raise PreventUpdate
        editor = target.get("editor", "project")
        row = int(target.get("row", 0))
        cells = []
        for name in _TABLES:
            cells.append(
                {"row": row, "column": 0} if name == editor else no_update
            )
        return editor, *cells


def _register_csv(app):
    @app.callback(
        Output(_CSV_PREVIEW, "data"),
        Output("pcbh-builder-csv-table", "data"),
        Output("pcbh-builder-csv-table", "columns"),
        Output("pcbh-builder-csv-status", "children"),
        Output(_DRAFT, "data", allow_duplicate=True),
        Input("pcbh-builder-csv-upload", "contents"),
        State("pcbh-builder-csv-upload", "filename"),
        State("pcbh-builder-project-crs_horizontal", "value"),
        State("pcbh-builder-profile-select", "value"),
        State(_PROFILES, "data"),
        prevent_initial_call=True,
    )
    def import_csv(contents, filename, crs, profile_name, profiles):
        if not contents:
            raise PreventUpdate
        try:
            raw = base64.b64decode(contents.split(",", 1)[1], validate=True)
            profile = None
            if profile_name and profile_name in (profiles or {}):
                profile = CSVMappingProfile.from_dict(profiles[profile_name])
            preview = preview_csv_mapping(raw, profile=profile)
            with tempfile.TemporaryDirectory() as folder:
                source = Path(folder) / (filename or "boreholes.csv")
                source.write_bytes(raw)
                constants = (
                    dict(profile.constants if profile else {})
                    if "crs.horizontal" in preview["mapping"]
                    else {
                        **(profile.constants if profile else {}),
                        "crs.horizontal": crs,
                    }
                )
                document, report = boreholes_from_csv(
                    source,
                    columns=profile.columns if profile else None,
                    constants=constants,
                    strict=False,
                )
            draft = document_to_builder(document)
            status = dbc.Alert(
                f"Imported {len(document.boreholes)} borehole(s); "
                f"{report.rows_rejected} rejected row(s).",
                color="success" if not report.errors else "warning",
            )
        except Exception as error:
            return (
                None,
                [],
                [],
                dbc.Alert(str(error), color="danger"),
                no_update,
            )
        columns = [{"name": name, "id": name} for name in preview["headers"]]
        return preview, preview["rows"], columns, status, draft


def _register_profile(app):
    @app.callback(
        Output(_PROFILES, "data"),
        Input("pcbh-builder-save-profile", "n_clicks"),
        State("pcbh-builder-profile-name", "value"),
        State(_CSV_PREVIEW, "data"),
        State(_PROFILES, "data"),
        prevent_initial_call=True,
    )
    def save_profile(n_clicks, name, preview, profiles):
        if not n_clicks or not name or not preview:
            raise PreventUpdate
        profile = CSVMappingProfile(name, preview.get("mapping", {}))
        result = dict(profiles or {})
        result[name] = profile.to_dict()
        return result

    @app.callback(
        Output("pcbh-builder-profile-select", "options"),
        Input(_PROFILES, "data"),
    )
    def profile_options(profiles):
        return [
            {"label": name, "value": name}
            for name in sorted((profiles or {}).keys())
        ]


def _register_export(app):
    @app.callback(
        Output("pcbh-builder-download", "data"),
        Input("pcbh-builder-export", "n_clicks"),
        State(_DOCUMENT, "data"),
        prevent_initial_call=True,
    )
    def export(n_clicks, document):
        if not n_clicks or not document:
            raise PreventUpdate
        name = str(document.get("document_id", "boreholes")).replace(":", "-")
        return dcc.send_string(
            json.dumps(document, ensure_ascii=False, indent=2) + "\n",
            f"{name}.pcbh.json",
        )


def _register_pcsf(app):
    @app.callback(
        Output("pcbh-builder-pcsf-download", "data"),
        Input("pcbh-builder-pcsf-upload", "contents"),
        State("pcbh-builder-pcsf-upload", "filename"),
        State(_DOCUMENT, "data"),
        prevent_initial_call=True,
    )
    def attach(contents, filename, document_data):
        if not contents or not document_data:
            raise PreventUpdate
        raw = base64.b64decode(contents.split(",", 1)[1], validate=True)
        from pycsamt.format.borehole import pcbh_from_dict

        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "source.pcsf"
            target = Path(folder) / "with-boreholes.pcsf"
            source.write_bytes(raw)
            model = embed_pcbh(
                read_pcsf(source), pcbh_from_dict(document_data)
            )
            write_pcsf(model, target)
            result = target.read_bytes()
        output_name = Path(filename or "model.pcsf").stem + "-pcbh.pcsf"
        return dcc.send_bytes(result, output_name)


def _log_figure(document):
    figure = go.Figure()
    for index, hole in enumerate(document.boreholes):
        for interval in hole.interval_logs.get("lithology", []):
            figure.add_trace(
                go.Bar(
                    x=[1],
                    y=[interval.to_md - interval.from_md],
                    base=interval.from_md,
                    name=interval.label or interval.code or "Unknown",
                    marker_color="#808080",
                    offsetgroup=str(index),
                    hovertemplate=(
                        f"{hole.id}<br>{interval.from_md:g}–"
                        f"{interval.to_md:g} MD<extra></extra>"
                    ),
                )
            )
    figure.update_layout(
        title="2-D lithology log preview",
        barmode="stack",
        yaxis={"autorange": "reversed", "title": "Measured depth"},
        height=350,
        showlegend=False,
    )
    return figure


def _diagnostic_label(item):
    location = item.editor
    if item.row is not None:
        location += f" row {item.row + 1}"
    return f"{location}: {item.message}"
