# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Borehole Studio (Map View) callback tests."""

from __future__ import annotations

import base64
import io

import pytest

openpyxl = pytest.importorskip("openpyxl")

from pycsamt.app.mapview.callbacks import borehole as bh  # noqa: E402


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


def _clean_xlsx() -> str:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.append(["Rock name", "From_m", "To_m", "Resistivity"])
    ws.append(["granodiorite", 0, 12.5, 240])
    ws.append(["hornblende", 12.5, 30, 85])
    buffer = io.BytesIO()
    wb.save(buffer)
    return "data:x;base64," + base64.b64encode(buffer.getvalue()).decode()


def test_registers_expected_callbacks():
    captured = _capture()
    assert {
        "load_pcbh",
        "toggle",
        "do_import",
        "inspect",
        "populate_columns",
        "apply_sheet",
        "rebuild",
        "add_row",
        "apply",
        "export",
        "load_points",
    } <= set(captured)


def test_populate_columns_offers_headers_and_auto_guesses():
    captured = _capture()
    payload, *_ = captured["inspect"](_clean_xlsx())
    result = captured["populate_columns"]("Sheet1", 1, payload)
    id_opts, lith_opts = result[0], result[1]
    lith_value, from_value, to_value = result[5], result[6], result[7]
    assert {opt["value"] for opt in lith_opts} >= {"", "Rock name", "From_m"}
    assert lith_value == "Rock name"
    assert from_value == "From_m" and to_value == "To_m"


def test_load_points_populates_pcpt_store():
    import base64

    captured = _capture()
    csv = "Target,Longitude,Latitude\nP1,103.1,22.2\n"
    url = "data:text/csv;base64," + base64.b64encode(csv.encode()).decode()
    store, info = captured["load_points"](url, "targets.csv")
    assert store["n_points"] == 1


def test_xlsx_inspect_then_apply_populates_tables_and_store():
    captured = _capture()
    payload, options, _value, _preview = captured["inspect"](_clean_xlsx())
    assert [opt["value"] for opt in options] == ["Sheet1"]

    (
        store,
        draft,
        collars,
        layers,
        survey,
        status,
    ) = captured["apply_sheet"](
        1, payload, "Sheet1", 1, "350000, 3000000, 1200", "", "", "", ""
    )
    assert store["n_boreholes"] == 1
    assert len(layers) == 2
    assert collars[0]["z"] == 1200.0
    assert "Imported" in status.children


def test_tables_rebuild_validates_and_previews():
    captured = _capture()
    payload, *_ = captured["inspect"](_clean_xlsx())
    _s, draft, collars, layers, survey, _st = captured["apply_sheet"](
        1, payload, "Sheet1", 1, "0,0,0", "", "", "", ""
    )
    new_draft, validation, figure = captured["rebuild"](
        collars, layers, survey, draft, "light"
    )
    assert "Valid" in validation.children
    assert len(figure.data) == 2


def test_apply_writes_pcbh_store_and_closes_modal():
    captured = _capture()
    payload, *_ = captured["inspect"](_clean_xlsx())
    _s, draft, collars, layers, survey, _st = captured["apply_sheet"](
        1, payload, "Sheet1", 1, "0,0,0", "", "", "", ""
    )
    draft2, *_ = captured["rebuild"](collars, layers, survey, draft, "light")
    store, is_open, status = captured["apply"](1, draft2)
    assert store["n_boreholes"] == 1
    assert is_open is False


def test_add_row_targets_the_active_tab():
    captured = _capture()
    collars, layers = captured["add_row"](1, "bh-layers", [], [])
    assert len(layers) == 1 and len(collars) == 0
    collars, layers = captured["add_row"](1, "bh-collars", [], [])
    assert len(collars) == 1 and len(layers) == 0
