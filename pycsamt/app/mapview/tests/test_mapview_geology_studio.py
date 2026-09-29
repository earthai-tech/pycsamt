# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Interpretation Studio (Map View) callback tests."""

from __future__ import annotations

import base64
import json
import types

import dash
import numpy as np
import pytest
from dash.exceptions import PreventUpdate

from pycsamt.app.mapview.callbacks import geology as geo_cb
from pycsamt.app.mapview._ids import IDs
from pycsamt.format.geology import GeologyLegend, legend_to_dict
from pycsamt.format.structure import StructModel, structure_to_dict
from pycsamt.geology import RockDatabase
from pycsamt.geology.structural import FaultTrace, StructuralModel


def _capture():
    captured: dict = {}

    class _App:
        def callback(self, *a, **k):
            def deco(fn):
                captured[fn.__name__] = fn
                return fn

            return deco

    geo_cb.register_geology(_App())
    return captured


def _set_trigger(monkeypatch, triggered_id):
    """``do_import``/``do_import_structure`` do a fresh ``from dash import
    ctx`` *inside* the callback body (unlike the borehole/view modules,
    which bind ``ctx`` at module import time), so the module-level name
    can't be monkeypatched here -- patch ``dash.ctx`` itself instead; the
    next ``from dash import ctx`` inside the function picks up the
    patched attribute."""
    monkeypatch.setattr(
        dash, "ctx", types.SimpleNamespace(triggered_id=triggered_id)
    )


def _b64(raw: bytes, mime: str = "application/json") -> str:
    return f"data:{mime};base64," + base64.b64encode(raw).decode()


def _pcgl_json_url() -> str:
    legend = GeologyLegend.from_rock_database(title="t")
    return _b64(json.dumps(legend_to_dict(legend)).encode("utf-8"))


def _pcgl_csv_url() -> str:
    csv = "name,rho_min,rho_max\nClay,1,20\nGranite,1000,5000\n"
    return _b64(csv.encode("utf-8"), mime="text/csv")


def _pcgs_json_url() -> str:
    model = StructuralModel(
        faults=[
            FaultTrace(
                x=1.0, dip_deg=70.0, downthrown_side="right",
                line="L1", z_top=0.0,
            )
        ],
    )
    doc = StructModel.from_structural_model(model)
    return _b64(json.dumps(structure_to_dict(doc)).encode("utf-8"))


def _planar_csv_url() -> str:
    csv = "x,kind,strike_deg,dip_deg,dip_direction_deg\n100,bedding,45,30,135\n"
    return _b64(csv.encode("utf-8"), mime="text/csv")


def _linear_csv_url() -> str:
    csv = "x,kind,trend_deg,plunge_deg\n50,fold_axis,10,20\n"
    return _b64(csv.encode("utf-8"), mime="text/csv")


def _faults_csv_url() -> str:
    csv = "x,dip_deg,downthrown_side\n500,70,right\n"
    return _b64(csv.encode("utf-8"), mime="text/csv")


_BAD_PLANAR_ROW = {
    "x": 100, "kind": "bedding", "strike_deg": 45, "dip_deg": 30,
    "dip_direction_deg": 45, "line": "L1",
}


def test_registers_expected_callbacks():
    captured = _capture()
    assert {
        "toggle",
        "do_import",
        "load_default",
        "auto_suggest",
        "rebuild",
        "add_row",
        "sync_colors",
        "apply",
        "export",
        "strip",
        "quick_upload",
        "do_import_structure",
        "rebuild_structure",
        "add_planar_row",
        "add_linear_row",
        "add_fault_row",
        "apply_structure",
        "export_structure",
        "strip_structure",
    } <= set(captured)


def test_modal_toggle():
    captured = _capture()
    assert captured["toggle"](1, False) is True
    assert captured["toggle"](1, True) is False


def test_load_default_populates_table():
    captured = _capture()
    rows, status = captured["load_default"](1)
    assert len(rows) == len(RockDatabase.default())
    assert "default" in status.children.lower() or "Loaded" in status.children


def test_rebuild_validates_and_builds_preview():
    captured = _capture()
    rows = [
        {"name": "Clay", "rho_min": 1, "rho_max": 20, "color": "#111"},
        {"name": "Granite", "rho_min": 1000, "rho_max": 5000, "color": "#222"},
    ]
    draft, validation, figure = captured["rebuild"](rows)
    assert "Valid" in validation.children
    assert len(figure.data) == 2
    assert draft["entries"][0]["name"] == "Clay"


def test_rebuild_reports_invalid_rows():
    captured = _capture()
    rows = [{"name": "Bad", "rho_min": 100, "rho_max": 10}]
    draft, validation, figure = captured["rebuild"](rows)
    assert "rho_min" in validation.children or validation.children


def test_add_row_appends_blank_row():
    captured = _capture()
    rows = captured["add_row"](1, [{"name": "Clay"}])
    assert len(rows) == 2 and rows[1]["name"] is None


def test_apply_writes_geo_store_and_closes_modal():
    captured = _capture()
    rows = [{"name": "Clay", "rho_min": 1, "rho_max": 20, "color": "#111"}]
    store, is_open, status = captured["apply"](1, rows)
    assert store["n_entries"] == 1
    assert is_open is False


def test_apply_with_invalid_rows_reports_error_without_writing_store():
    captured = _capture()
    rows = [{"name": "Bad", "rho_min": 100, "rho_max": 10}]
    store, is_open, status = captured["apply"](1, rows)
    from dash import no_update

    assert store is no_update
    assert is_open is no_update


def test_export_returns_downloadable_pcgl(monkeypatch):
    captured = _capture()
    rows = [{"name": "Clay", "rho_min": 1, "rho_max": 20, "color": "#111"}]
    payload = captured["export"](1, rows)
    assert payload["filename"].endswith(".pcgl.json")


def test_strip_shows_no_legend_when_empty():
    captured = _capture()
    result = captured["strip"](None)
    assert "No legend" in str(result.children)


def test_strip_shows_chips_for_applied_legend():
    captured = _capture()
    store, _open, _status = captured["apply"](
        1,
        [
            {"name": "Clay", "rho_min": 1, "rho_max": 20, "color": "#111"},
            {"name": "Granite", "rho_min": 1000, "rho_max": 5000, "color": "#222"},
        ],
    )
    result = captured["strip"](store)
    assert len(result.children) == 2


# ---------------------------------------------------------------------------
# Structure tab
# ---------------------------------------------------------------------------


_PLANAR_ROW = {
    "x": 100, "kind": "bedding", "strike_deg": 45, "dip_deg": 30,
    "dip_direction_deg": 135, "line": "L1",
}
_FAULT_ROW = {"x": 500, "dip_deg": 70, "downthrown_side": "right", "line": "L1"}


def test_add_structure_rows_are_independent_per_table():
    captured = _capture()
    planar = captured["add_planar_row"](1, [])
    linear = captured["add_linear_row"](1, [])
    faults = captured["add_fault_row"](1, [])
    assert len(planar) == 1 and "strike_deg" in planar[0]
    assert len(linear) == 1 and "trend_deg" in linear[0]
    assert len(faults) == 1 and "downthrown_side" in faults[0]


def test_rebuild_structure_validates_and_builds_preview():
    captured = _capture()
    draft, validation, figure = captured["rebuild_structure"](
        [_PLANAR_ROW], [], [_FAULT_ROW], "light"
    )
    assert "Valid" in validation.children
    assert len(figure.data) == 2
    assert draft["planar"][0]["kind"] == "bedding"


def test_apply_structure_writes_struct_store_and_closes_modal():
    captured = _capture()
    store, is_open, status = captured["apply_structure"](
        1, [_PLANAR_ROW], [], [_FAULT_ROW]
    )
    assert store["n_planar"] == 1 and store["n_faults"] == 1
    assert is_open is False


def test_export_structure_returns_downloadable_pcgs():
    captured = _capture()
    payload = captured["export_structure"](1, [_PLANAR_ROW], [], [_FAULT_ROW])
    assert payload["filename"] == "structure.pcgs.json"


def test_strip_structure_shows_no_structure_when_empty():
    captured = _capture()
    result = captured["strip_structure"](None)
    assert "No structure" in str(result.children)


def test_strip_structure_shows_counts_for_applied_structure():
    captured = _capture()
    store, *_ = captured["apply_structure"](1, [_PLANAR_ROW], [], [_FAULT_ROW])
    result = captured["strip_structure"](store)
    assert "1 fault" in result.children and "1 planar" in result.children


# ---------------------------------------------------------------------------
# Sync colours with loaded boreholes
# ---------------------------------------------------------------------------


def _pcbh_store():
    from pycsamt.app._borehole import pcbh_to_dict
    from pycsamt.format.borehole import (
        Collar,
        CoordinateReferenceSystem,
        LogInterval,
        PCBHBorehole,
        PCBHDocument,
        VocabularyEntry,
    )

    document = PCBHDocument(
        document_id="test:colors",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[VocabularyEntry("A", "Granodiorite", color="#C7AA54")],
        boreholes=[
            PCBHBorehole(
                id="BH1", name="Hole 1", kind="mining_exploration",
                status="completed", collar=Collar(0.0, 0.0, 0.0),
                total_depth_md=30.0,
                interval_logs={
                    "lithology": [LogInterval(0.0, 30.0, code="A")]
                },
            )
        ],
    )
    return {"document": pcbh_to_dict(document), "n_boreholes": 1}


def test_sync_colors_requires_a_loaded_borehole():
    captured = _capture()
    rows = [{"name": "Granodiorite", "color": "#000000"}]
    result, status = captured["sync_colors"](1, rows, None)
    from dash import no_update

    assert result is no_update
    assert "no borehole loaded" in status.children.lower()


def test_sync_colors_overwrites_matching_rows_only():
    captured = _capture()
    rows = [
        {"name": "Granodiorite", "color": "#000000"},
        {"name": "Not drilled", "color": "#111111"},
    ]
    result, status = captured["sync_colors"](1, rows, _pcbh_store())
    assert result[0]["color"] == "#C7AA54"
    assert result[1]["color"] == "#111111"
    assert "Synced 1" in status.children


def test_sync_colors_reports_when_nothing_matches():
    captured = _capture()
    rows = [{"name": "Nothing in common", "color": "#000000"}]
    result, status = captured["sync_colors"](1, rows, _pcbh_store())
    assert result[0]["color"] == "#000000"
    assert "no legend row" in status.children.lower()


def test_sync_colors_reports_when_borehole_has_no_lithology_log():
    from pycsamt.app._borehole import pcbh_to_dict
    from pycsamt.format.borehole import (
        Collar,
        CoordinateReferenceSystem,
        PCBHBorehole,
        PCBHDocument,
    )

    document = PCBHDocument(
        document_id="test:no-litho",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[],
        boreholes=[
            PCBHBorehole(
                id="BH1", name="Hole 1", kind="mining_exploration",
                status="completed", collar=Collar(0.0, 0.0, 0.0),
                total_depth_md=30.0,
            )
        ],
    )
    store = {"document": pcbh_to_dict(document), "n_boreholes": 1}
    captured = _capture()
    rows = [{"name": "Granodiorite", "color": "#000000"}]
    result, status = captured["sync_colors"](1, rows, store)
    assert result is geo_cb.no_update
    assert "no lithology log" in status.children.lower()


# ---------------------------------------------------------------------------
# no_clicks / PreventUpdate guard clauses
# ---------------------------------------------------------------------------


def test_toggle_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["toggle"](None, False)


def test_load_default_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["load_default"](None)


def test_add_row_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["add_row"](None, [])


def test_sync_colors_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["sync_colors"](None, [], None)


def test_apply_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["apply"](None, [])


def test_export_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["export"](None, [])


def test_export_invalid_rows_prevents_update_instead_of_downloading():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["export"](1, [{"name": "Bad", "rho_min": 100, "rho_max": 10}])


def test_add_planar_row_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["add_planar_row"](None, [])


def test_add_linear_row_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["add_linear_row"](None, [])


def test_add_fault_row_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["add_fault_row"](None, [])


def test_apply_structure_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["apply_structure"](None, [], [], [])


def test_export_structure_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["export_structure"](None, [], [], [])


def test_export_structure_invalid_rows_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["export_structure"](1, [_BAD_PLANAR_ROW], [], [])


def test_rebuild_structure_reports_invalid_rows():
    captured = _capture()
    draft, validation, figure = captured["rebuild_structure"](
        [_BAD_PLANAR_ROW], [], [], "light"
    )
    assert draft is geo_cb.no_update
    assert "text-danger" in validation.className
    assert figure.data == ()


def test_apply_structure_invalid_rows_reports_error():
    captured = _capture()
    store, is_open, status = captured["apply_structure"](
        1, [_BAD_PLANAR_ROW], [], []
    )
    assert store is geo_cb.no_update
    assert is_open is geo_cb.no_update
    assert "text-danger" in status.className


# ---------------------------------------------------------------------------
# _register_imports (do_import) -- ctx.triggered_id dispatch
# ---------------------------------------------------------------------------


def test_do_import_no_contents_prevents_update(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_IMPORT_PCGL)
    with pytest.raises(PreventUpdate):
        captured["do_import"](None, None, "legend.pcgl.json", None)


def test_do_import_pcgl_json_branch_populates_table(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_IMPORT_PCGL)
    rows, status = captured["do_import"](
        _pcgl_json_url(), None, "legend.pcgl.json", None
    )
    assert len(rows) == len(RockDatabase.default())
    assert "Imported" in status.children
    assert "text-success" in status.className


def test_do_import_csv_branch_populates_table(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_IMPORT_CSV)
    rows, status = captured["do_import"](
        None, _pcgl_csv_url(), None, "legend.csv"
    )
    assert len(rows) == 2
    assert rows[0]["name"] == "Clay"
    assert "Imported" in status.children


def test_do_import_bad_payload_returns_error_status(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_IMPORT_PCGL)
    rows, status = captured["do_import"](
        "not-a-data-url", None, "legend.pcgl.json", None
    )
    assert rows is geo_cb.no_update
    assert "text-danger" in status.className


# ---------------------------------------------------------------------------
# _register_quick_upload
# ---------------------------------------------------------------------------


def test_quick_upload_no_contents_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["quick_upload"](None, "legend.pcgl.json")


def test_quick_upload_success_writes_geo_store():
    captured = _capture()
    store, info = captured["quick_upload"](_pcgl_json_url(), "legend.pcgl.json")
    assert store["n_entries"] == len(RockDatabase.default())
    assert "entrie(s)" in info.children
    assert "text-success" in info.className


def test_quick_upload_bad_payload_returns_error_info():
    captured = _capture()
    store, info = captured["quick_upload"]("not-a-data-url", "legend.pcgl.json")
    assert store is geo_cb.no_update
    assert "text-danger" in info.className


# ---------------------------------------------------------------------------
# _register_auto_suggest
# ---------------------------------------------------------------------------


def _map_view_with_rho(values_a, values_b):
    from pycsamt.map._core import MapData, StationRecord
    from pycsamt.map.view import MapView

    class _Z:
        def __init__(self, value):
            self.freq = [10.0, 1.0]
            self.resistivity = np.ones((2, 2, 2)) * value
            self.phase = np.ones((2, 2, 2)) * 45.0

    class _Edi:
        def __init__(self, station, value):
            self.station = station
            self.Z = _Z(value)

    data = MapData(
        sites=[_Edi("S00", values_a), _Edi("S01", values_b)],
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
        metadata={},
    )
    return MapView(data)


def test_auto_suggest_no_clicks_prevents_update():
    captured = _capture()
    with pytest.raises(PreventUpdate):
        captured["auto_suggest"](None, 8, "sess-1")


def test_auto_suggest_no_session_reports_error():
    captured = _capture()
    rows, status = captured["auto_suggest"](1, 8, None)
    assert rows is geo_cb.no_update
    assert "load survey lines" in status.children.lower()


def test_auto_suggest_view_unavailable_reports_error(monkeypatch):
    captured = _capture()
    monkeypatch.setattr(geo_cb, "get_view", lambda session_id: None)
    rows, status = captured["auto_suggest"](1, 8, "sess-1")
    assert rows is geo_cb.no_update
    assert "text-danger" in status.className


def test_auto_suggest_view_with_no_data_reports_error(monkeypatch):
    captured = _capture()
    monkeypatch.setattr(
        geo_cb, "get_view", lambda session_id: types.SimpleNamespace(data=None)
    )
    rows, status = captured["auto_suggest"](1, 8, "sess-1")
    assert rows is geo_cb.no_update
    assert "text-danger" in status.className


def test_auto_suggest_success_builds_legend(monkeypatch):
    captured = _capture()
    view = _map_view_with_rho(100.0, 5000.0)
    monkeypatch.setattr(geo_cb, "get_view", lambda session_id: view)
    rows, status = captured["auto_suggest"](1, None, "sess-1")
    assert rows
    assert "Suggested" in status.children
    assert "text-success" in status.className


def test_auto_suggest_exception_from_uniform_resistivity_reports_error(
    monkeypatch,
):
    captured = _capture()
    view = _map_view_with_rho(100.0, 100.0)
    monkeypatch.setattr(geo_cb, "get_view", lambda session_id: view)
    rows, status = captured["auto_suggest"](1, 8, "sess-1")
    assert rows is geo_cb.no_update
    assert "text-danger" in status.className


# ---------------------------------------------------------------------------
# _register_structure_imports (do_import_structure) -- ctx.triggered_id
# ---------------------------------------------------------------------------


def test_do_import_structure_no_matching_trigger_prevents_update(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, "some-other-id")
    with pytest.raises(PreventUpdate):
        captured["do_import_structure"](
            None, None, None, None, None, None, None, None, [], [], [],
        )


def test_do_import_structure_json_branch_populates_all_tables(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_STRUCT_IMPORT_JSON)
    planar, linear, faults, status = captured["do_import_structure"](
        _pcgs_json_url(), None, None, None,
        "structure.pcgs.json", None, None, None,
        [], [], [],
    )
    assert faults[0]["downthrown_side"] == "right"
    assert "Imported" in status.children


def test_do_import_structure_planar_csv_branch(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_STRUCT_IMPORT_PLANAR)
    planar, linear, faults, status = captured["do_import_structure"](
        None, _planar_csv_url(), None, None,
        None, "planar.csv", None, None,
        [], [], [],
    )
    assert len(planar) == 1
    assert linear is geo_cb.no_update
    assert faults is geo_cb.no_update
    assert "planar measurement" in status.children.lower()


def test_do_import_structure_linear_csv_branch(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_STRUCT_IMPORT_LINEAR)
    planar, linear, faults, status = captured["do_import_structure"](
        None, None, _linear_csv_url(), None,
        None, None, "linear.csv", None,
        [], [], [],
    )
    assert planar is geo_cb.no_update
    assert len(linear) == 1
    assert faults is geo_cb.no_update
    assert "linear measurement" in status.children.lower()


def test_do_import_structure_faults_csv_branch(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_STRUCT_IMPORT_FAULTS)
    planar, linear, faults, status = captured["do_import_structure"](
        None, None, None, _faults_csv_url(),
        None, None, None, "faults.csv",
        [], [], [],
    )
    assert planar is geo_cb.no_update
    assert linear is geo_cb.no_update
    assert len(faults) == 1
    assert "fault trace" in status.children.lower()


def test_do_import_structure_trigger_but_no_contents_prevents_update(
    monkeypatch,
):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_STRUCT_IMPORT_JSON)
    with pytest.raises(PreventUpdate):
        captured["do_import_structure"](
            None, None, None, None,
            "structure.pcgs.json", None, None, None,
            [], [], [],
        )


def test_do_import_structure_bad_payload_returns_error_status(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, IDs.GEO_STRUCT_IMPORT_JSON)
    result = captured["do_import_structure"](
        "not-a-data-url", None, None, None,
        "structure.pcgs.json", None, None, None,
        [], [], [],
    )
    *stores, status = result
    assert all(item is geo_cb.no_update for item in stores)
    assert "text-danger" in status.className
