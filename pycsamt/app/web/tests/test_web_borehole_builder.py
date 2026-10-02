"""Borehole Builder page and integration contracts."""

from __future__ import annotations

import base64
import binascii
import json
from types import SimpleNamespace

import numpy as np
import plotly.graph_objects as go
import pytest
from dash import no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.web.pages import borehole_builder
from pycsamt.format import Grid3DGeometry, PCSFModel, read_pcsf, write_pcsf
from pycsamt.format.borehole import (
    BuilderDiagnostic,
    BuilderValidation,
    CSVMappingProfile,
    Collar,
    CoordinateReferenceSystem,
    PCBHBorehole,
    PCBHDocument,
    document_from_builder,
    document_to_builder,
    new_builder_draft,
    pcbh_to_dict,
    validate_builder,
)


class _CaptureApp:
    """Records callback bodies so they can be invoked without a browser."""

    def __init__(self):
        self.functions = []

    def callback(self, *_args, **_kwargs):
        def decorate(fn):
            self.functions.append(fn)
            return fn

        return decorate

    def get(self, name, index=0):
        return [fn for fn in self.functions if fn.__name__ == name][index]


def _capture(register_name):
    app = _CaptureApp()
    getattr(borehole_builder, register_name)(app)
    return app


_PROJECT_KEYS = [
    "document_id",
    "created_by",
    "crs_horizontal",
    "crs_vertical",
    "coordinate_unit",
    "depth_unit",
    "diameter_unit",
]


def _full_draft():
    draft = new_builder_draft()
    draft["project"].update(
        {
            "document_id": "builder:test",
            "crs_horizontal": "EPSG:32629",
        }
    )
    draft["boreholes"] = [
        {
            "id": "BH-1",
            "name": "Builder hole",
            "kind": "water",
            "status": "completed",
            "x": 500000,
            "y": 600000,
            "z": 100,
            "total_depth_md": 50,
            "diameter": 0.2,
            "north_reference": "grid",
        }
    ]
    draft["surveys"] = [
        {
            "borehole_id": "BH-1",
            "md": 0,
            "azimuth_deg": 0,
            "inclination_deg": 0,
        },
        {
            "borehole_id": "BH-1",
            "md": 50,
            "azimuth_deg": 10,
            "inclination_deg": 15,
        },
    ]
    draft["intervals"] = [
        {
            "borehole_id": "BH-1",
            "family": "lithology",
            "from_md": 0,
            "to_md": 50,
            "code": "SAND",
            "label": "Sand",
            "color": "#D4B483",
        }
    ]
    draft["structures"] = [
        {"borehole_id": "BH-1", "kind": "fracture", "at_md": 25}
    ]
    draft["water"] = [
        {"borehole_id": "BH-1", "at_md": 12, "kind": "water_strike"}
    ]
    draft["construction"] = [
        {"borehole_id": "BH-1", "from_md": 0, "to_md": 20, "kind": "casing"}
    ]
    draft["samples"] = [
        {
            "borehole_id": "BH-1",
            "sample_id": "S-1",
            "from_md": 10,
            "to_md": 11,
        }
    ]
    draft["assays"] = [
        {
            "borehole_id": "BH-1",
            "sample_id": "S-1",
            "analyte": "Au",
            "value": 1.2,
        }
    ]
    return draft


def _project_values(draft):
    return [draft["project"].get(key) for key in _PROJECT_KEYS]


def _table_values(draft):
    return [draft.get(name, []) for name in borehole_builder._TABLES]


def _pcsf_model():
    geometry = Grid3DGeometry(
        x=np.array([25.0, 75.0]),
        y=np.array([25.0, 75.0]),
        z=np.array([20.0, 60.0, 100.0]),
        x_nodes=np.array([0.0, 50.0, 100.0]),
        y_nodes=np.array([0.0, 50.0, 100.0]),
        z_nodes=np.array([0.0, 40.0, 80.0, 120.0]),
        origin=np.array([500000.0, 4500000.0, 0.0]),
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=np.full(geometry.resistivity_shape, 100.0),
        crs="EPSG:32629",
    )


def _pcsf_document():
    return PCBHDocument(
        document_id="test:pcsf-association",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:32629", vertical="EPSG:5714"),
        boreholes=[
            PCBHBorehole(
                id="BH-1",
                name="Borehole 1",
                kind="mining_exploration",
                status="completed",
                collar=Collar(500050.0, 4500050.0, 0.0),
                total_depth_md=80.0,
            )
        ],
    )


def _b64_contents(raw: bytes, mime: str = "text/csv") -> str:
    return f"data:{mime};base64," + base64.b64encode(raw).decode("ascii")


def test_builder_page_exposes_all_editors_and_accessible_feedback():
    rendered = str(borehole_builder.layout())

    for name in (
        "boreholes",
        "surveys",
        "intervals",
        "structures",
        "water",
        "construction",
        "samples",
        "assays",
    ):
        assert f"pcbh-builder-{name}" in rendered
    assert "pcbh-builder-validation" in rendered
    assert "role='alert'" in rendered or '"role":"alert"' in rendered
    assert "pcbh-builder-draft" in rendered
    assert "storage_type='local'" in rendered or '"storage_type":"local"' in rendered


def test_builder_callbacks_are_registered(web_app):
    callback_text = str(web_app.callback_map)
    assert "pcbh-builder-document" in callback_text
    assert "pcbh-builder-download" in callback_text
    assert "pcbh-builder-pcsf-download" in callback_text
    assert "pcbh-builder-csv-preview" in callback_text


def test_navigation_contains_builder_page(web_app):
    layout_text = str(web_app.layout)
    assert "nav-btn-borehole-builder" in layout_text
    assert "page-borehole-builder" in layout_text


def test_large_builder_draft_constructs_many_holes():
    draft = new_builder_draft()
    draft["project"]["crs_horizontal"] = "EPSG:32629"
    draft["boreholes"] = [
        {
            "id": f"BH-{index:04d}",
            "name": f"Hole {index}",
            "kind": "mining_exploration",
            "status": "planned",
            "x": 500000 + index,
            "y": 600000,
            "z": 100,
            "total_depth_md": 10,
        }
        for index in range(1_000)
    ]

    document = document_from_builder(draft)
    assert len(document.boreholes) == 1_000


# ── _register_add_row ───────────────────────────────────────────────────────


class TestAddRow:
    def test_adds_blank_row_to_active_tab_only(self):
        add_row = _capture("_register_add_row").get("add_row")
        tables = [[] for _ in borehole_builder._TABLES]

        result = add_row(1, "surveys", *tables)

        names = list(borehole_builder._TABLES)
        for name, rows in zip(names, result):
            if name == "surveys":
                assert rows == [
                    {key: None for key in borehole_builder._TABLES["surveys"]}
                ]
            else:
                assert rows == []

    def test_preserves_existing_rows_in_active_tab(self):
        add_row = _capture("_register_add_row").get("add_row")
        tables = [[] for _ in borehole_builder._TABLES]
        boreholes_index = list(borehole_builder._TABLES).index("boreholes")
        tables[boreholes_index] = [{"id": "BH-1"}]

        result = add_row(1, "boreholes", *tables)

        assert result[boreholes_index][0] == {"id": "BH-1"}
        assert len(result[boreholes_index]) == 2

    def test_inactive_tab_raises_prevent_update(self):
        add_row = _capture("_register_add_row").get("add_row")
        tables = [[] for _ in borehole_builder._TABLES]

        with pytest.raises(PreventUpdate):
            add_row(1, "project", *tables)

    def test_no_clicks_raises_prevent_update(self):
        add_row = _capture("_register_add_row").get("add_row")
        tables = [[] for _ in borehole_builder._TABLES]

        with pytest.raises(PreventUpdate):
            add_row(None, "boreholes", *tables)


# ── _register_restore ───────────────────────────────────────────────────────


class TestRestore:
    def test_restores_project_fields_and_tables(self):
        restore = _capture("_register_restore").get("restore")
        draft = _full_draft()

        result = restore(1, draft)

        project_out = result[: len(_PROJECT_KEYS)]
        table_out = result[len(_PROJECT_KEYS) :]
        assert list(project_out) == _project_values(draft)
        assert list(table_out) == _table_values(draft)

    def test_no_draft_raises_prevent_update(self):
        restore = _capture("_register_restore").get("restore")

        with pytest.raises(PreventUpdate):
            restore(1, None)

    def test_no_clicks_raises_prevent_update(self):
        restore = _capture("_register_restore").get("restore")

        with pytest.raises(PreventUpdate):
            restore(None, _full_draft())


# ── _register_builder ───────────────────────────────────────────────────────


class TestBuilder:
    def test_valid_draft_produces_canonical_document_and_figures(self):
        build = _capture("_register_builder").get("build")
        draft = _full_draft()

        draft_out, canonical, message, log_fig, traj_fig = build(
            *_project_values(draft), *_table_values(draft), None
        )

        assert draft_out["boreholes"] == draft["boreholes"]
        assert canonical["document_id"] == "builder:test"
        assert "Valid PCBH" in message.children
        assert "1 borehole(s)" in message.children
        assert isinstance(log_fig, go.Figure)
        assert isinstance(traj_fig, go.Figure)
        assert len(log_fig.data) == 1
        assert len(traj_fig.data) >= 1

    def test_invalid_draft_without_borehole_or_crs_yields_unique_issue_ids(
        self,
    ):
        build = _capture("_register_builder").get("build")
        draft = new_builder_draft()
        draft["project"]["crs_horizontal"] = ""
        draft["boreholes"] = []

        draft_out, document, messages, log_fig, traj_fig = build(
            *_project_values(draft), *_table_values(draft), None
        )

        assert document is None
        assert isinstance(log_fig, go.Figure) and not log_fig.data
        assert isinstance(traj_fig, go.Figure) and not traj_fig.data
        assert len(messages) >= 2

        ids = [button.id for button in messages]
        # All issues here are project-level (no specific row) - historically
        # `item.row or 0` collapsed every such issue onto the same
        # {"editor": "project", "row": 0} id as a genuine row-0 issue would
        # use, so Dash rendered duplicate component ids. Confirm they are
        # now distinguishable.
        assert len(ids) == len({tuple(sorted(i.items())) for i in ids})
        for issue_id in ids:
            assert issue_id["editor"] == "project"
            assert issue_id["row"] == -1

    def test_row_zero_issue_is_distinguishable_from_no_row_issue(self):
        # Directly exercises the id-collision fix: a diagnostic tied to the
        # first row of a table (row=0) must not collide with a diagnostic
        # that has no row at all (row=None) under the same editor.
        build = _capture("_register_builder").get("build")

        crafted = BuilderValidation(
            None,
            (
                BuilderDiagnostic(
                    "error", "x.missing", "no row here", "$", "boreholes"
                ),
                BuilderDiagnostic(
                    "error", "y.missing", "row zero", "$.boreholes[0]",
                    "boreholes", 0,
                ),
            ),
        )

        import pycsamt.app.web.pages.borehole_builder as module

        original = module.validate_builder
        module.validate_builder = lambda draft: crafted
        try:
            draft = new_builder_draft()
            _, _, messages, _, _ = build(
                *_project_values(draft), *_table_values(draft), None
            )
        finally:
            module.validate_builder = original

        assert len(messages) == 2
        rows = [button.id["row"] for button in messages]
        assert -1 in rows
        assert 0 in rows
        ids = [tuple(sorted(button.id.items())) for button in messages]
        assert len(ids) == len(set(ids))

    def test_empty_tables_edge_case_still_builds(self):
        build = _capture("_register_builder").get("build")
        draft = _full_draft()
        for name in (
            "surveys",
            "structures",
            "water",
            "construction",
            "samples",
            "assays",
        ):
            draft[name] = []

        draft_out, canonical, message, log_fig, traj_fig = build(
            *_project_values(draft), *_table_values(draft), None
        )

        assert canonical is not None
        assert "Valid PCBH" in message.children

    def test_previous_draft_state_used_when_falsy(self):
        build = _capture("_register_builder").get("build")
        draft = _full_draft()

        draft_out, *_ = build(
            *_project_values(draft), *_table_values(draft), None
        )

        assert draft_out["project"]["document_id"] == "builder:test"


# ── _register_issue_navigation ──────────────────────────────────────────────


class TestIssueNavigation:
    def test_non_dict_triggered_id_raises_prevent_update(self, monkeypatch):
        navigate = _capture("_register_issue_navigation").get("navigate")
        monkeypatch.setattr(
            borehole_builder, "ctx", SimpleNamespace(triggered_id=None)
        )

        with pytest.raises(PreventUpdate):
            navigate([1])

    def test_valid_issue_switches_tab_and_sets_active_cell(self, monkeypatch):
        navigate = _capture("_register_issue_navigation").get("navigate")
        monkeypatch.setattr(
            borehole_builder,
            "ctx",
            SimpleNamespace(
                triggered_id={
                    "type": "pcbh-builder-issue",
                    "editor": "intervals",
                    "row": 3,
                    "index": 1,
                }
            ),
        )

        result = navigate([1])

        active_tab, *cells = result
        assert active_tab == "intervals"
        names = list(borehole_builder._TABLES)
        for name, cell in zip(names, cells):
            if name == "intervals":
                assert cell == {"row": 3, "column": 0}
            else:
                assert cell is no_update

    def test_no_row_issue_falls_back_to_first_row_cell(self, monkeypatch):
        navigate = _capture("_register_issue_navigation").get("navigate")
        monkeypatch.setattr(
            borehole_builder,
            "ctx",
            SimpleNamespace(
                triggered_id={
                    "type": "pcbh-builder-issue",
                    "editor": "boreholes",
                    "row": -1,
                    "index": 0,
                }
            ),
        )

        active_tab, *cells = navigate([1])

        assert active_tab == "boreholes"
        names = list(borehole_builder._TABLES)
        boreholes_cell = cells[names.index("boreholes")]
        assert boreholes_cell == {"row": 0, "column": 0}


# ── _register_csv ────────────────────────────────────────────────────────────


class TestCsvImport:
    def test_csv_column_crs_takes_precedence_over_project_field(self):
        import_csv = _capture("_register_csv").get("import_csv")
        raw = (
            "borehole_id,x,y,z,from_md,to_md,lithology,crs\n"
            "BH-1,500000,600000,100,0,10,SAND,EPSG:32629\n"
            "BH-1,500000,600000,100,10,20,CLAY,EPSG:32629\n"
        ).encode("utf-8")

        preview, rows, columns, status, draft = import_csv(
            _b64_contents(raw), "holes.csv", "EPSG:9999", None, {}
        )

        assert "crs.horizontal" in preview["mapping"]
        assert draft["project"]["crs_horizontal"] == "EPSG:32629"
        assert "Imported 1 borehole(s)" in status.children
        assert len(rows) == 2

    def test_csv_without_crs_column_uses_project_field_constant(self):
        import_csv = _capture("_register_csv").get("import_csv")
        raw = (
            "borehole_id,x,y,z,from_md,to_md,lithology\n"
            "BH-2,500100,600100,110,0,5,SAND\n"
            "BH-2,500100,600100,110,5,15,CLAY\n"
        ).encode("utf-8")

        preview, rows, columns, status, draft = import_csv(
            _b64_contents(raw), "holes.csv", "EPSG:4326", None, {}
        )

        assert "crs.horizontal" not in preview["mapping"]
        assert draft["project"]["crs_horizontal"] == "EPSG:4326"

    def test_csv_with_saved_profile_reuses_mapping(self):
        import_csv = _capture("_register_csv").get("import_csv")
        raw = (
            "HOLE,E,N,RL,START,END,ROCK\n"
            "W1,10,20,30,0,15,Clay\n"
            "W1,10,20,30,15,25,Sand\n"
        ).encode("utf-8")
        profile = CSVMappingProfile(
            "wells",
            {
                "borehole.id": "HOLE",
                "collar.x": "E",
                "collar.y": "N",
                "collar.z": "RL",
                "interval.from_md": "START",
                "interval.to_md": "END",
                "interval.lithology": "ROCK",
            },
        )
        profiles = {"wells": profile.to_dict()}

        preview, rows, columns, status, draft = import_csv(
            _b64_contents(raw), "wells.csv", "EPSG:4326", "wells", profiles
        )

        assert preview["mapping"]["borehole.id"] == "HOLE"
        assert draft["boreholes"][0]["id"] == "W1"

    def test_no_contents_raises_prevent_update(self):
        import_csv = _capture("_register_csv").get("import_csv")

        with pytest.raises(PreventUpdate):
            import_csv(None, "holes.csv", "EPSG:4326", None, {})

    def test_malformed_contents_returns_error_alert(self):
        import_csv = _capture("_register_csv").get("import_csv")

        preview, rows, columns, status, draft = import_csv(
            "data:text/csv;base64,not-valid-base64!!",
            "holes.csv",
            "EPSG:4326",
            None,
            {},
        )

        assert preview is None
        assert rows == []
        assert columns == []
        assert status.color == "danger"
        assert draft is no_update


# ── _register_profile ───────────────────────────────────────────────────────


class TestProfile:
    def test_save_profile_stores_mapping_under_name(self):
        save_profile = _capture("_register_profile").get("save_profile")
        preview = {"mapping": {"borehole.id": "HOLE"}}

        result = save_profile(1, "wells", preview, {})

        assert result["wells"] == CSVMappingProfile(
            "wells", {"borehole.id": "HOLE"}
        ).to_dict()

    def test_save_profile_merges_into_existing_profiles(self):
        save_profile = _capture("_register_profile").get("save_profile")
        preview = {"mapping": {"borehole.id": "HOLE"}}
        existing = {"other": {"name": "other", "columns": {}, "constants": {}}}

        result = save_profile(1, "wells", preview, existing)

        assert set(result) == {"other", "wells"}

    def test_save_profile_missing_name_raises_prevent_update(self):
        save_profile = _capture("_register_profile").get("save_profile")

        with pytest.raises(PreventUpdate):
            save_profile(1, "", {"mapping": {}}, {})

    def test_save_profile_missing_preview_raises_prevent_update(self):
        save_profile = _capture("_register_profile").get("save_profile")

        with pytest.raises(PreventUpdate):
            save_profile(1, "wells", None, {})

    def test_save_profile_no_clicks_raises_prevent_update(self):
        save_profile = _capture("_register_profile").get("save_profile")

        with pytest.raises(PreventUpdate):
            save_profile(None, "wells", {"mapping": {}}, {})

    def test_profile_options_sorted_by_name(self):
        profile_options = _capture("_register_profile").get(
            "profile_options"
        )

        result = profile_options({"zeta": {}, "alpha": {}, "mid": {}})

        assert result == [
            {"label": "alpha", "value": "alpha"},
            {"label": "mid", "value": "mid"},
            {"label": "zeta", "value": "zeta"},
        ]

    def test_profile_options_empty_when_no_profiles(self):
        profile_options = _capture("_register_profile").get(
            "profile_options"
        )

        assert profile_options(None) == []


# ── _register_export ─────────────────────────────────────────────────────────


class TestExport:
    def test_export_produces_pcbh_json_download(self):
        export = _capture("_register_export").get("export")
        document = {"document_id": "pcbh:my-project", "boreholes": []}

        result = export(1, document)

        assert result["filename"] == "pcbh-my-project.pcbh.json"
        assert json.loads(result["content"]) == document

    def test_export_no_document_raises_prevent_update(self):
        export = _capture("_register_export").get("export")

        with pytest.raises(PreventUpdate):
            export(1, None)

    def test_export_no_clicks_raises_prevent_update(self):
        export = _capture("_register_export").get("export")

        with pytest.raises(PreventUpdate):
            export(None, {"document_id": "pcbh:x"})


# ── _register_pcsf ───────────────────────────────────────────────────────────


class TestPcsfAttach:
    def test_attach_embeds_document_and_round_trips(self, tmp_path):
        attach = _capture("_register_pcsf").get("attach")
        source_path = write_pcsf(_pcsf_model(), tmp_path / "source.pcsf")
        raw = source_path.read_bytes()
        document = _pcsf_document()
        document_data = pcbh_to_dict(document)

        result = attach(
            _b64_contents(raw, "application/octet-stream"),
            "source.pcsf",
            document_data,
        )

        assert result["filename"] == "source-pcbh.pcsf"
        content = base64.b64decode(result["content"])
        out_path = tmp_path / "roundtrip.pcsf"
        out_path.write_bytes(content)
        restored = read_pcsf(out_path)
        assert restored.boreholes.embedded.document_id == document.document_id

    def test_attach_default_filename_when_missing(self, tmp_path):
        attach = _capture("_register_pcsf").get("attach")
        source_path = write_pcsf(_pcsf_model(), tmp_path / "source.pcsf")
        raw = source_path.read_bytes()
        document_data = pcbh_to_dict(_pcsf_document())

        result = attach(
            _b64_contents(raw, "application/octet-stream"),
            None,
            document_data,
        )

        assert result["filename"] == "model-pcbh.pcsf"

    def test_attach_no_contents_raises_prevent_update(self):
        attach = _capture("_register_pcsf").get("attach")

        with pytest.raises(PreventUpdate):
            attach(None, "source.pcsf", {"document_id": "pcbh:x"})

    def test_attach_no_document_raises_prevent_update(self, tmp_path):
        attach = _capture("_register_pcsf").get("attach")
        source_path = write_pcsf(_pcsf_model(), tmp_path / "source.pcsf")
        raw = source_path.read_bytes()

        with pytest.raises(PreventUpdate):
            attach(
                _b64_contents(raw, "application/octet-stream"),
                "source.pcsf",
                None,
            )


# ── _log_figure ──────────────────────────────────────────────────────────────


class TestLogFigure:
    def test_multiple_boreholes_and_intervals_produce_one_trace_each(self):
        draft = _full_draft()
        draft["boreholes"].append(
            {
                "id": "BH-2",
                "name": "Second hole",
                "kind": "water",
                "status": "completed",
                "x": 500010,
                "y": 600010,
                "z": 100,
                "total_depth_md": 30,
            }
        )
        draft["intervals"].append(
            {
                "borehole_id": "BH-2",
                "family": "lithology",
                "from_md": 0,
                "to_md": 15,
                "code": "CLAY",
                "label": "Clay",
            }
        )
        draft["intervals"].append(
            {
                "borehole_id": "BH-2",
                "family": "lithology",
                "from_md": 15,
                "to_md": 30,
                "code": "SILT",
                "label": None,
            }
        )
        document = document_from_builder(draft)

        figure = borehole_builder._log_figure(document)

        assert len(figure.data) == 3
        names = [trace.name for trace in figure.data]
        assert names == ["Sand", "Clay", "SILT"]
        assert figure.data[0].offsetgroup == "0"
        assert figure.data[1].offsetgroup == "1"

    def test_borehole_without_lithology_yields_no_traces(self):
        draft = new_builder_draft()
        draft["project"]["crs_horizontal"] = "EPSG:4326"
        draft["boreholes"] = [
            {
                "id": "BH-DRY",
                "name": "Dry hole",
                "kind": "water",
                "status": "abandoned",
                "x": 1,
                "y": 2,
                "z": 3,
                "total_depth_md": 5,
            }
        ]
        document = document_from_builder(draft)

        figure = borehole_builder._log_figure(document)

        assert len(figure.data) == 0
        assert figure.layout.title.text == "2-D lithology log preview"


# ── _diagnostic_label ────────────────────────────────────────────────────────


class TestDiagnosticLabel:
    def test_label_with_row_is_one_indexed(self):
        item = BuilderDiagnostic(
            "error", "x.missing", "value is required", "$.boreholes[2]",
            "boreholes", 2,
        )

        assert (
            borehole_builder._diagnostic_label(item)
            == "boreholes row 3: value is required"
        )

    def test_label_without_row_omits_row_suffix(self):
        item = BuilderDiagnostic(
            "error", "crs.horizontal", "must be set", "$.crs.horizontal",
            "project", None,
        )

        assert (
            borehole_builder._diagnostic_label(item)
            == "project: must be set"
        )
