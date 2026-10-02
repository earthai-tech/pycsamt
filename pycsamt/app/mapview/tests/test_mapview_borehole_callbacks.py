# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Behavioral coverage for the Map View Borehole Studio callback module."""

from __future__ import annotations

import base64
import io
import json
import types

import pytest
from dash.exceptions import PreventUpdate

from pycsamt.app.mapview.callbacks import borehole as bh
from pycsamt.app.mapview._ids import IDs
from pycsamt.format.borehole import pcbh_to_dict
from pycsamt.format.borehole.builder import document_from_builder, new_builder_draft


def _capture():
    captured: dict = {}

    class _App:
        def callback(self, *a, **k):
            def deco(fn):
                captured[fn.__name__] = fn
                return fn

            return deco

    bh.register_borehole(_App())
    return captured


def _set_trigger(monkeypatch, triggered_id):
    """``ctx.triggered_id`` is only populated by a real Dash dispatch;
    unit-testing the unwrapped callback needs the module's bound ``ctx``
    name replaced directly (see memory: monkeypatch the consuming
    module's import, not Dash's own global)."""
    monkeypatch.setattr(
        bh, "ctx", types.SimpleNamespace(triggered_id=triggered_id)
    )


def _b64(raw: bytes, mime: str = "application/json") -> str:
    return f"data:{mime};base64," + base64.b64encode(raw).decode()


def _valid_draft() -> dict:
    draft = new_builder_draft()
    draft["boreholes"] = [
        {
            "id": "BH-1",
            "name": "Hole 1",
            "kind": "water",
            "status": "completed",
            "x": 10.0,
            "y": 20.0,
            "z": 100.0,
            "total_depth_md": 30.0,
        }
    ]
    draft["intervals"] = [
        {
            "borehole_id": "BH-1",
            "from_md": 0.0,
            "to_md": 15.0,
            "code": "SAP",
            "label": "Saprolite",
            "family": "lithology",
        }
    ]
    return draft


def _valid_document():
    return document_from_builder(_valid_draft())


def _pcbh_upload_url() -> str:
    payload = pcbh_to_dict(_valid_document())
    raw = json.dumps(payload).encode("utf-8")
    return _b64(raw)


def _combined_csv_url() -> str:
    csv = (
        "borehole_id,x,y,z,crs,kind,status,total_depth_md,from_md,to_md,"
        "lithology,resistivity_ohm_m\n"
        "BH-9,1.0,2.0,3.0,EPSG:32629,water,completed,20,0,10,Clay,50\n"
    )
    return _b64(csv.encode("utf-8"), mime="text/csv")


# ---------------------------------------------------------------------------
# registration
# ---------------------------------------------------------------------------


def test_register_borehole_wires_all_eight_groups():
    captured = _capture()
    assert {
        "load_points", "load_pcbh", "toggle", "do_import", "inspect",
        "populate_columns", "apply_sheet", "rebuild", "add_row", "apply",
        "export",
    } <= set(captured)


# ---------------------------------------------------------------------------
# _register_points
# ---------------------------------------------------------------------------


def test_load_points_no_contents_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["load_points"](None, "f.csv")


def test_load_points_success_reports_point_count():
    captured = _capture()
    csv = "Target,Longitude,Latitude\nP1,103.1,22.2\n"
    url = _b64(csv.encode(), mime="text/csv")
    store, info = captured["load_points"](url, "targets.csv")
    assert store["n_points"] == 1
    assert "point(s)" in info.children


def test_load_points_bad_payload_returns_error_span():
    captured = _capture()
    store, info = captured["load_points"]("not-a-data-url", "f.csv")
    assert store is bh.no_update
    assert info.__class__.__name__ == "Span"
    assert "text-danger" in info.className


# ---------------------------------------------------------------------------
# _register_accordion_upload
# ---------------------------------------------------------------------------


def test_load_pcbh_no_contents_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["load_pcbh"](None, "f.pcbh.json")


def test_load_pcbh_success_reports_borehole_count():
    captured = _capture()
    store, info = captured["load_pcbh"](
        _pcbh_upload_url(), "boreholes.pcbh.json"
    )
    assert store["n_boreholes"] == 1
    assert "1 borehole(s)" in info.children
    assert "text-success" in info.className


def test_load_pcbh_bad_extension_returns_error_span():
    captured = _capture()
    store, info = captured["load_pcbh"](_pcbh_upload_url(), "boreholes.txt")
    assert store is bh.no_update
    assert "text-danger" in info.className


# ---------------------------------------------------------------------------
# _register_modal_toggle
# ---------------------------------------------------------------------------


def test_toggle_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["toggle"](None, None, False)


@pytest.mark.parametrize("n_accordion, n_toolbar", [(1, None), (None, 1), (1, 1)])
def test_toggle_flips_open_state(n_accordion, n_toolbar):
    captured = _capture()
    assert captured["toggle"](n_accordion, n_toolbar, False) is True
    assert captured["toggle"](n_accordion, n_toolbar, True) is False


# ---------------------------------------------------------------------------
# _register_imports (do_import) -- ctx.triggered_id branching
# ---------------------------------------------------------------------------


def test_do_import_no_matching_trigger_prevents_update(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, "some-other-id")
    with pytest.raises(PreventUpdate):
        captured["do_import"](None, None, None, None)


def test_do_import_pcbh_branch_populates_tables(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.BH_IMPORT_PCBH)
    store, draft, collars, layers, survey, status = captured["do_import"](
        _pcbh_upload_url(), None, "boreholes.pcbh.json", None
    )
    assert store["n_boreholes"] == 1
    assert len(collars) == 1 and collars[0]["id"] == "BH-1"
    assert len(layers) == 1
    assert "Imported 1 borehole(s)" in status.children


def test_do_import_csv_branch_populates_tables(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.BH_IMPORT_CSV)
    store, draft, collars, layers, survey, status = captured["do_import"](
        None, _combined_csv_url(), None, "boreholes.csv"
    )
    assert store["n_boreholes"] == 1
    assert store["filename"] == "boreholes.csv"
    assert len(layers) == 1
    assert "Imported" in status.children


def test_do_import_pcbh_trigger_but_no_contents_prevents_update(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.BH_IMPORT_PCBH)
    with pytest.raises(PreventUpdate):
        captured["do_import"](None, None, "boreholes.pcbh.json", None)


def test_do_import_bad_pcbh_returns_error_status_not_crash(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.BH_IMPORT_PCBH)
    result = captured["do_import"](
        "not-a-data-url", None, "boreholes.pcbh.json", None
    )
    *stores, status = result
    assert all(item is bh.no_update for item in stores)
    assert "text-danger" in status.className


def test_do_import_bad_csv_returns_error_status(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.BH_IMPORT_CSV)
    bad_url = "data:text/csv;base64,not-valid-base64!!"
    result = captured["do_import"](None, bad_url, None, "boreholes.csv")
    *stores, status = result
    assert all(item is bh.no_update for item in stores)
    assert "text-danger" in status.className


# ---------------------------------------------------------------------------
# _register_xlsx
# ---------------------------------------------------------------------------


def _xlsx_url() -> str:
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.append(["Rock name", "From_m", "To_m", "Resistivity"])
    ws.append(["granodiorite", 0, 12.5, 240])
    ws.append(["hornblende", 12.5, 30, 85])
    buffer = io.BytesIO()
    wb.save(buffer)
    return "data:x;base64," + base64.b64encode(buffer.getvalue()).decode()


def test_inspect_no_contents_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["inspect"](None)


def test_inspect_bad_workbook_returns_error_preview():
    captured = _capture()
    payload, options, value, preview = captured["inspect"]("data:x;base64,AAAA")
    assert payload is None
    assert options == []
    assert value is None
    assert "text-danger" in preview.className


def test_inspect_success_returns_sheet_options_and_preview():
    captured = _capture()
    payload, options, value, preview = captured["inspect"](_xlsx_url())
    assert payload["b64"]
    assert options == [{"label": "Sheet1", "value": "Sheet1"}]
    assert value == "Sheet1"
    assert preview.__class__.__name__ == "DataTable"


def test_populate_columns_no_payload_or_sheet_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["populate_columns"]("Sheet1", 1, None)
    with pytest.raises(PreventUpdate):
        captured["populate_columns"](None, 1, {"b64": "AAAA"})


def test_populate_columns_bad_payload_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["populate_columns"]("Sheet1", 1, {"b64": "not-base64!"})


def test_populate_columns_guesses_headers():
    captured = _capture()
    payload, *_ = captured["inspect"](_xlsx_url())
    result = captured["populate_columns"]("Sheet1", 1, payload)
    (
        id_opts, lith_opts, from_opts, to_opts,
        id_val, lith_val, from_val, to_val, preview,
    ) = result
    assert {"", "Rock name", "From_m", "To_m"} <= {o["value"] for o in lith_opts}
    assert lith_val == "Rock name"
    assert from_val == "From_m"
    assert to_val == "To_m"
    assert preview.__class__.__name__ == "DataTable"


def test_populate_columns_no_header_row_uses_auto():
    captured = _capture()
    payload, *_ = captured["inspect"](_xlsx_url())
    result = captured["populate_columns"]("Sheet1", None, payload)
    assert result[8].__class__.__name__ == "DataTable"


def test_apply_sheet_no_clicks_or_payload_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["apply_sheet"](
            None, {"b64": "AAAA"}, "Sheet1", 1, "", "", "", "", ""
        )
    with pytest.raises(PreventUpdate):
        captured["apply_sheet"](1, None, "Sheet1", 1, "", "", "", "", "")


def test_apply_sheet_success_builds_store_and_tables():
    captured = _capture()
    payload, *_ = captured["inspect"](_xlsx_url())
    store, draft, collars, layers, survey, status = captured["apply_sheet"](
        1, payload, "Sheet1", 1, "350000, 3000000, 1200", "", "", "", "",
    )
    assert store["n_boreholes"] == 1
    assert len(layers) == 2
    assert collars[0]["z"] == 1200.0
    assert "Imported" in status.children


class _FakeReport:
    rows_accepted = 1
    rows_rejected = 1


def test_apply_sheet_reports_rejected_rows_when_present(monkeypatch):
    captured = _capture()
    payload, *_ = captured["inspect"](_xlsx_url())

    def _fake_boreholes_from_xlsx(*a, **k):
        return _valid_document(), _FakeReport()

    import pycsamt.format.borehole as fb

    monkeypatch.setattr(fb, "boreholes_from_xlsx", _fake_boreholes_from_xlsx)
    store, draft, collars, layers, survey, status = captured["apply_sheet"](
        1, payload, "Sheet1", 1, "0,0,0", "", "", "", "",
    )
    assert "rejected" in status.children


def test_apply_sheet_skips_crs_override_when_collar_has_crs(monkeypatch):
    captured = _capture()
    payload, *_ = captured["inspect"](_xlsx_url())
    monkeypatch.setattr(
        bh, "_parse_collar",
        lambda text: {"x": 0.0, "y": 0.0, "z": 0.0, "crs": "EPSG:4326"},
    )
    seen: dict = {}

    def _fake_boreholes_from_xlsx(raw, **kwargs):
        seen.update(kwargs)
        return _valid_document(), _FakeReport()

    import pycsamt.format.borehole as fb

    monkeypatch.setattr(fb, "boreholes_from_xlsx", _fake_boreholes_from_xlsx)
    captured["apply_sheet"](
        1, payload, "Sheet1", 1, "irrelevant", "", "", "", "",
    )
    assert "crs.horizontal" not in seen["constants"]


def test_apply_sheet_error_returns_error_status():
    captured = _capture()
    store, draft, collars, layers, survey, status = captured["apply_sheet"](
        1, {"b64": "not-base64!"}, "Sheet1", 1, "", "", "", "", "",
    )
    assert store is bh.no_update
    assert "text-danger" in status.className


# ---------------------------------------------------------------------------
# _register_tables
# ---------------------------------------------------------------------------


def test_rebuild_success_returns_valid_status_and_figure():
    captured = _capture()
    collars, layers, survey = bh._tables_from_draft(_valid_draft())
    new_draft, validation, figure = captured["rebuild"](
        collars, layers, survey, _valid_draft(), "dark"
    )
    assert "Valid" in validation.children
    assert len(figure.data) >= 1


def test_rebuild_empty_tables_reports_validation_error():
    captured = _capture()
    new_draft, validation, figure = captured["rebuild"]([], [], [], None, "light")
    assert "text-danger" in validation.className
    assert figure.data == ()


def test_add_row_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["add_row"](None, "bh-collars", [], [])


def test_add_row_layers_tab_appends_layer_row():
    captured = _capture()
    collars, layers = captured["add_row"](1, "bh-layers", None, None)
    assert collars == []
    assert len(layers) == 1
    assert set(layers[0]) == set(bh._LAYER_COLS)


def test_add_row_default_tab_appends_collar_row():
    captured = _capture()
    collars, layers = captured["add_row"](1, "bh-collars", [], [])
    assert len(collars) == 1
    assert layers == []


# ---------------------------------------------------------------------------
# _register_apply
# ---------------------------------------------------------------------------


def test_apply_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["apply"](None, _valid_draft())


def test_apply_empty_draft_reports_nothing_to_apply():
    captured = _capture()
    store, is_open, status = captured["apply"](1, None)
    assert store is bh.no_update and is_open is bh.no_update
    assert "nothing to apply" in status.children


def test_apply_invalid_nonempty_draft_reports_error():
    captured = _capture()
    store, is_open, status = captured["apply"](1, {"boreholes": []})
    assert store is bh.no_update
    assert "text-danger" in status.className


def test_apply_success_writes_store_and_closes_modal():
    captured = _capture()
    store, is_open, status = captured["apply"](1, _valid_draft())
    assert store["n_boreholes"] == 1
    assert is_open is False
    assert "Applied" in status.children


# ---------------------------------------------------------------------------
# _register_export
# ---------------------------------------------------------------------------


def test_export_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["export"](None, _valid_draft(), None)


def test_export_no_draft_and_no_store_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["export"](1, None, None)


def test_export_from_draft_returns_download():
    captured = _capture()
    result = captured["export"](1, _valid_draft(), None)
    assert result["filename"].endswith(".pcbh.json")
    body = json.loads(result["content"])
    assert body["boreholes"][0]["id"] == "BH-1"


def test_export_falls_back_to_store_when_draft_is_empty():
    """Real bug fix: an empty/None draft used to short-circuit straight
    to PreventUpdate even when the PCBH store held a valid document
    (e.g. after an accordion-only upload with no studio draft built
    yet); it must fall back to the store instead."""
    captured = _capture()
    store = {
        "filename": "x.pcbh.json",
        "document": pcbh_to_dict(_valid_document()),
        "n_boreholes": 1,
    }
    result = captured["export"](1, None, store)
    body = json.loads(result["content"])
    assert body["boreholes"][0]["id"] == "BH-1"


def test_export_falls_back_to_store_when_draft_is_invalid():
    captured = _capture()
    store = {
        "filename": "x.pcbh.json",
        "document": pcbh_to_dict(_valid_document()),
        "n_boreholes": 1,
    }
    result = captured["export"](1, {"boreholes": [{"id": None}]}, store)
    body = json.loads(result["content"])
    assert body["boreholes"][0]["id"] == "BH-1"


# ---------------------------------------------------------------------------
# module helpers
# ---------------------------------------------------------------------------


def test_decode_rejects_missing_comma():
    with pytest.raises(ValueError, match="could not decode"):
        bh._decode("no-comma-here")


def test_decode_rejects_bad_base64():
    with pytest.raises(ValueError, match="could not decode"):
        bh._decode("data:text/csv;base64,not-valid-base64!!")


def test_decode_roundtrips_good_payload():
    assert bh._decode(_b64(b"hello")) == b"hello"


def test_tmp_write_writes_bytes_to_disk():
    from pathlib import Path

    path = bh._tmp_write(b"payload", "sample.csv")
    assert Path(path).read_bytes() == b"payload"
    assert Path(path).name == "sample.csv"


@pytest.mark.parametrize(
    "text, expected",
    [
        ("Sheet 1!!", "Sheet-1"),
        ("  ..weird//name..  ", "weird-name"),
        ("", ""),
        (None, ""),
    ],
)
def test_slug_normalizes_text(text, expected):
    assert bh._slug(text) == expected


def test_slug_truncates_to_40_chars():
    assert len(bh._slug("a" * 100)) == 40


@pytest.mark.parametrize(
    "text, expected",
    [
        (None, None),
        ("", None),
        ("   ", None),
        ("1,2", None),
        ("a,2,3", None),
        ("1,2,3", {"x": 1.0, "y": 2.0, "z": 3.0}),
        ("1;2;3", {"x": 1.0, "y": 2.0, "z": 3.0}),
        ("1, 2, 3, 999", {"x": 1.0, "y": 2.0, "z": 3.0}),
    ],
)
def test_parse_collar(text, expected):
    assert bh._parse_collar(text) == expected


def test_tables_from_draft_none_returns_empty_tables():
    collars, layers, survey = bh._tables_from_draft(None)
    assert collars == [] and layers == [] and survey == []


def test_tables_from_draft_filters_non_lithology_intervals():
    draft = {
        "boreholes": [{"id": "BH-1", "x": 1, "y": 2, "z": 3}],
        "intervals": [
            {"borehole_id": "BH-1", "from_md": 0, "to_md": 1, "family": "lithology"},
            {"borehole_id": "BH-1", "from_md": 1, "to_md": 2, "family": "structure"},
        ],
        "surveys": [{"borehole_id": "BH-1", "md": 0, "azimuth_deg": 10, "inclination_deg": -80}],
    }
    collars, layers, survey = bh._tables_from_draft(draft)
    assert len(collars) == 1
    assert len(layers) == 1
    assert len(survey) == 1


def test_draft_from_tables_sets_project_defaults_when_no_previous():
    draft = bh._draft_from_tables([], [], [], None)
    assert draft["project"]["document_id"] == "pcbh:studio"
    assert draft["project"]["crs_horizontal"] == "LOCAL:studio-grid"
    for key in (
        "structures", "water", "construction", "samples", "assays",
        "mapping_profiles",
    ):
        assert draft[key] == []
    assert draft["created_at"].endswith("Z")


def test_draft_from_tables_keeps_previous_project_values():
    previous = {"project": {"document_id": "pcbh:keep-me"}}
    draft = bh._draft_from_tables([], [], [], previous)
    assert draft["project"]["document_id"] == "pcbh:keep-me"


def test_draft_from_tables_filters_rows_missing_required_keys():
    collars = [{"id": "BH-1", "x": 1}, {"id": None, "x": 2}, {"id": "", "x": 3}]
    layers = [
        {"borehole_id": "BH-1", "from_md": 0.0},
        {"borehole_id": "BH-1", "from_md": None},
        {"borehole_id": None, "from_md": 1.0},
    ]
    survey = [{"borehole_id": "BH-1"}, {"borehole_id": None}]
    draft = bh._draft_from_tables(collars, layers, survey, None)
    assert len(draft["boreholes"]) == 1
    assert len(draft["intervals"]) == 1
    assert len(draft["surveys"]) == 1


def test_preview_table_empty_shows_empty_message():
    result = bh._preview_table([])
    assert result.__class__.__name__ == "Span"
    assert "empty sheet" in result.children


def test_preview_table_pads_ragged_rows():
    table = bh._preview_table([["a", "b", "c"], ["1"]])
    assert len(table.columns) == 3
    assert table.data[1] == {"0": "1", "1": "", "2": ""}


def test_preview_table_truncates_to_twelve_rows():
    rows = [["h"]] + [[str(i)] for i in range(20)]
    table = bh._preview_table(rows)
    assert len(table.data) == 12


def test_ok_and_err_spans():
    ok = bh._ok("good")
    err = bh._err("bad")
    assert ok.children == "good" and "text-success" in ok.className
    assert err.children == "bad" and "text-danger" in err.className


def test_bytes_io_wraps_raw_bytes():
    buf = bh._bytes_io(b"xyz")
    assert isinstance(buf, io.BytesIO)
    assert buf.getvalue() == b"xyz"


def test_guess_columns_matches_known_aliases():
    guess = bh._guess_columns(["ID", "Lithology", "From_m", "To_m"])
    assert guess["id"] == "ID"
    assert guess["lithology"] == "Lithology"
    assert guess["from_md"] == "From_m"
    assert guess["to_md"] == "To_m"


def test_guess_columns_no_match_returns_empty_dict():
    assert bh._guess_columns(["Foo", "Bar"]) == {}
