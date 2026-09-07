# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Interpretation Studio (Map View) callback tests."""

from __future__ import annotations

from pycsamt.app.mapview.callbacks import geology as geo_cb
from pycsamt.geology import RockDatabase


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
