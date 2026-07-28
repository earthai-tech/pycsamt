"""Direct callback tests for MapView.

Dash normally hides callback bodies behind its dispatcher.  A tiny registrar
captures the original functions so their state/validation contracts can be
tested without a browser.
"""

from __future__ import annotations

import base64
import inspect
import json
from types import SimpleNamespace

import pytest
from dash import no_update


class CaptureApp:
    def __init__(self):
        self.functions = []
        self.clientside = []

    def callback(self, *_args, **_kwargs):
        def decorate(fn):
            self.functions.append(fn)
            return fn

        return decorate

    def clientside_callback(self, *args, **kwargs):
        self.clientside.append((args, kwargs))

    def get(self, name, index=0):
        return [fn for fn in self.functions if fn.__name__ == name][index]


def _capture(module, register_name):
    app = CaptureApp()
    getattr(module, register_name)(app)
    return app


def test_controls_callbacks(monkeypatch):
    from pycsamt.app.mapview._ids import IDs
    from pycsamt.app.mapview.callbacks import controls

    app = _capture(controls, "register_controls")
    monkeypatch.setattr(controls, "ctx", SimpleNamespace(triggered_id=IDs.BTN_DEPTH_1K))
    assert app.get("apply_preset")() == (0, 1000)
    monkeypatch.setattr(controls, "ctx", SimpleNamespace(triggered_id=IDs.BTN_RHO_COND))
    assert app.get("apply_rho_preset")() == (1, 100)
    assert app.get("populate")({}) == (0, 0, None, 0)
    assert app.get("populate")({"frequencies": [0.1, 1, 1000]})[1] == 2

    gather = app.get("gather")
    values = [None] * len(inspect.signature(gather).parameters)
    values[-1] = {"frequencies": [0.1, 10.0]}
    values[3] = 99
    result, label = gather(*values)
    assert result["frequency"] == 10.0
    assert result["overlay"] == "index"
    assert "Hz" in label


def test_export_and_view_callbacks(monkeypatch):
    from pycsamt.app.mapview.callbacks import export, view

    app = _capture(export, "register_export")
    cb = app.get("export")
    monkeypatch.setattr(export, "get_view", lambda _sid: None)
    assert cb(1, "sid", "map", {}, "light", {}) is no_update
    fig = SimpleNamespace(to_html=lambda **_: "<html/>")
    monkeypatch.setattr(export, "get_view", lambda _sid: object())
    monkeypatch.setattr(export, "figure_for", lambda *a, **k: fig)
    payload = cb(1, "sid", "map3d", {}, "dark", {"active": ["L1"]})
    assert payload["filename"] == "mapview_map3d.html"

    app = _capture(view, "register_view")
    monkeypatch.setattr(view, "ctx", SimpleNamespace(triggered_id="unknown"))
    assert app.get("switch")()[0] == "map"
    render = app.get("render")
    monkeypatch.setattr(view, "empty_figure", lambda *a: ("empty", a))
    assert render({}, "map", {}, "light", {}, 0, [], "sid")[1]["display"] == "flex"
    monkeypatch.setattr(view, "get_view", lambda _sid: object())
    monkeypatch.setattr(view, "figure_for", lambda *a, **k: "figure")
    assert render({"n_stations": 1}, "map", {}, None, {}, 2, [], "sid") == (
        "figure",
        {"display": "none"},
    )


def test_inversion_import_decode_and_confirm(monkeypatch, tmp_path):
    import pycsamt.map as map_pkg
    from pycsamt.app.mapview.callbacks import inversion_import as inv

    monkeypatch.setattr(inv.tempfile, "mkdtemp", lambda **_: str(tmp_path))
    content = "data:application/octet-stream;base64," + base64.b64encode(b"x").decode()
    inv._decode_to_tempdir(["folder/model.dat", "bad"], [content, "?"])
    assert (tmp_path / "model.dat").read_bytes() == b"x"

    app = _capture(inv, "register_inversion_import")
    confirm = app.get("confirm")
    assert confirm(1, {}, False, "sid", "replace", "light")[0] is no_update
    assert (
        "Session" in confirm(1, {"filenames": ["x"]}, False, "", "replace", "light")[1]
    )

    fake_view = SimpleNamespace(n_stations=2, data=SimpleNamespace(stations=["known"]))
    monkeypatch.setattr(inv, "_decode_to_tempdir", lambda *a: str(tmp_path))
    monkeypatch.setattr(
        map_pkg,
        "MapView",
        SimpleNamespace(from_inversion_results=lambda *a, **k: fake_view),
    )
    monkeypatch.setattr(inv, "get_view", lambda _sid: fake_view)
    monkeypatch.setattr(inv, "merge_views", lambda old, new: new)
    monkeypatch.setattr(inv, "set_view", lambda *a: None)
    monkeypatch.setattr(
        inv,
        "store_from_view",
        lambda _v: {"n_stations": 2, "n_lines": 1},
    )
    out = confirm(
        1,
        {"filenames": ["x"], "contents": [content]},
        True,
        "sid",
        "append",
        None,
    )
    assert out[0]["n_stations"] == 2 and out[2] is False


def test_settings_callbacks(monkeypatch):
    from pycsamt.app.mapview.callbacks import settings

    app = _capture(settings, "register_settings")
    assert app.get("toggle")(1, False) is True
    assert app.get("clear")(1) == []
    store = {
        "n_stations": 2,
        "lines": ["L1", "L2"],
        "station_records": [
            {"ID": "A", "Line": "L1"},
            {"ID": "B", "Line": "L2"},
        ],
    }
    assert app.get("mask_hidden")(1, store, {"active": ["L1"]}, []) == ["B"]
    opts, active = app.get("populate")(store, {"all": ["L1", "L2"], "active": ["L1"]})
    assert len(opts) == 2 and active == ["L1"]
    assert app.get("set_active")(["L2"], {"all": ["L1", "L2"], "active": ["L1"]})[
        "active"
    ] == ["L2"]
    assert "No data" in str(app.get("summary")({}, {}, [])[0])
    summary, chips = app.get("summary")(store, {"active": ["L1"]}, ["A"])
    assert "masked" in str(summary) and "A" in str(chips)

    monkeypatch.setattr(
        settings,
        "ctx",
        SimpleNamespace(triggered_id={"index": "A"}, triggered=[{"value": 1}]),
    )
    assert app.get("toggle_mask")([], []) == ["A"]
    assert app.get("toggle_mask")([], ["A"]) == []
    assert app.get("unmask")([], ["A", "B"]) == ["B"]


def test_session_callbacks():
    from pycsamt.app.mapview.callbacks import session

    with pytest.raises(ValueError):
        session._validate_snapshot([])
    with pytest.raises(ValueError):
        session._validate_snapshot({"version": "2"})
    snap = session._build_snapshot(
        "map",
        {"a": 1},
        {"active": ["L"]},
        None,
        ["S"],
        "dark",
        {"n_stations": 1},
        "note",
    )
    session._validate_snapshot(snap)
    app = _capture(session, "register_session")
    assert app.get("_toggle")(1, False) is True
    assert app.get("_auto_snapshot")("map", {}, {}, None, [], "light", {}, "") == (
        no_update,
        no_update,
    )
    saved, chip = app.get("_auto_snapshot")(
        "map", {}, {}, None, [], "light", {"n_stations": 1}, ""
    )
    assert saved["app"] == "mapview" and "Auto-saved" in str(chip)
    assert app.get("_download")(0, *([None] * 8)) == (no_update, no_update)
    download = app.get("_download")(
        1, "map", {}, {}, None, [], "light", {"n_stations": 1}, "n"
    )
    assert download[0]["filename"].endswith(".json")

    encoded = base64.b64encode(json.dumps(snap).encode()).decode()
    restored = app.get("_upload_restore")("data:json;base64," + encoded, "s.json")
    assert restored[0] == "map"
    assert app.get("_upload_restore")("bad", "bad.json")[0] is no_update
    assert app.get("_browser_restore")(1, snap)[0] == "map"
    assert app.get("_browser_restore")(1, None)[0] is no_update
    assert app.get("_clear")(1)[0] is None


def test_toolbar_callbacks(monkeypatch):
    from pycsamt.app.mapview._ids import IDs
    from pycsamt.app.mapview.callbacks import toolbar

    app = _capture(toolbar, "register_toolbar")
    assert app.get("info")({}) == ""
    assert "2 stations" in app.get("info")({"n_stations": 2, "n_lines": 1})
    monkeypatch.setattr(toolbar, "ctx", SimpleNamespace(triggered_id=IDs.TB_BM_DARK))
    assert app.get("quickswitch", 0)() == "carto-darkmatter"
    monkeypatch.setattr(
        toolbar, "ctx", SimpleNamespace(triggered_id=IDs.TB3D_MODE_DEPTH)
    )
    assert app.get("quickswitch", 1)() == "depth"
    assert "EPSG" in str(app.get("crs_info")("geo", None, None, None))


def test_topography_callbacks(monkeypatch, tmp_path):
    import tempfile

    import pycsamt.map as map_pkg
    from pycsamt.app.mapview.callbacks import topo

    app = _capture(topo, "register_topo")
    parse = app.get("parse")
    assert parse(None, None) == (no_update, no_update)
    monkeypatch.setattr(map_pkg, "parse_elevation_file", lambda *a: {"S1": 10})
    assert parse("encoded", "e.csv")[0] == {"S1": 10}
    apply = app.get("apply")
    monkeypatch.setattr(topo, "get_view", lambda _sid: None)
    assert "Load lines" in str(apply(1, "upload", {}, None, "sid")[1])

    station = SimpleNamespace(id="S1", elevation=10)
    new_view = SimpleNamespace(data=SimpleNamespace(stations=[station]))
    old_view = SimpleNamespace(
        data=SimpleNamespace(stations=[station]),
        with_elevations=lambda _e: new_view,
        fetch_elevations=lambda **_: {"S1": 20},
    )
    monkeypatch.setattr(topo, "get_view", lambda _sid: old_view)
    monkeypatch.setattr(topo, "set_view", lambda *a: None)
    monkeypatch.setattr(topo, "store_from_view", lambda _v: {"n_stations": 1})
    assert apply(1, "upload", {"S1": 20}, None, "sid")[2] is True
    assert apply(1, "fetch", {}, "api", "sid")[2] is True
    assert apply(1, "stations", {}, None, "sid")[2] is True

    export = app.get("export")
    monkeypatch.setattr(tempfile, "mkdtemp", lambda **_: str(tmp_path))

    def write(path, fmt):
        mode = "wb" if fmt != "csv" else "w"
        with open(path, mode) as fh:
            fh.write(b"x" if "b" in mode else "x")

    old_view.export_topography = write
    assert export(1, "csv", "sid")[0]["filename"] == "topography.csv"


def test_lines_callbacks(monkeypatch):
    import numpy as np

    import pycsamt.app.mapview._render as render_mod
    from pycsamt.app.mapview._ids import IDs
    from pycsamt.app.mapview.callbacks import lines

    app = _capture(lines, "register_lines")
    records = [
        {
            "ID": "A",
            "Line": "L1",
            "Latitude": 1,
            "Longitude": 2,
            "Elevation": 100,
            "Index": 0,
        },
        {
            "ID": "B",
            "Line": "L2",
            "Latitude": 3,
            "Longitude": 4,
            "Elevation": 200,
            "Index": 1,
        },
    ]
    store = {"lines": ["L1", "L2"], "station_records": records}
    state = app.get("init_lines")(store, {"active": ["L2"]})
    assert state["active"] == ["L2"]
    pills, rows = app.get("render")(store, state, {"station_id": "B"}, ["B"])
    assert len(pills) == 3 and "mv-sta-masked" in str(rows)
    assert "No lines" in str(app.get("render")({}, {}, {}, [])[1])
    assert "100 m" in lines._row_meta(records[0])
    assert lines._fmt_coord("bad", "°").endswith("—")

    monkeypatch.setattr(
        lines,
        "ctx",
        SimpleNamespace(
            triggered_id={"type": "mv-pill", "index": "L2"},
            triggered=[{"value": 1}],
        ),
    )
    assert app.get("toggle")([], state)["active"] == ["L1", "L2"]
    monkeypatch.setattr(
        lines,
        "ctx",
        SimpleNamespace(triggered_id=IDs.CANVAS_GRAPH, triggered=[{"value": 1}]),
    )
    assert app.get("select")([], {"points": [{"customdata": ["A"]}]}) == {
        "station_id": "A"
    }
    inspected = app.get("inspect")({"station_id": "A"}, {}, store)
    assert "Elevation" in str(inspected)
    assert "Station not found" in str(
        app.get("inspect")({"station_id": "X"}, {}, store)
    )
    data, cols = app.get("fill_table")(
        store, {"active": ["L1"]}, {"crs_mode": "geo"}, []
    )
    assert [r["ID"] for r in data] == ["A"]
    monkeypatch.setattr(
        render_mod,
        "project_to_crs",
        lambda *a: (np.array([10.0]), np.array([20.0]), 32650),
    )
    projected, cols = lines._add_projected_columns(
        [dict(records[0])], {"crs_mode": "utm"}, cols
    )
    assert projected[0]["E (32650)"] == 10.0


def test_load_helpers_and_callbacks(monkeypatch):
    from pycsamt.app.mapview._ids import IDs
    from pycsamt.app.mapview.callbacks import load

    assert load._sanitize("../bad:name.edi", "x") == "bad_name.edi"
    entries = load._entries(
        ["a", "b"], ["root/L1/a.edi", "root/L2/readme.txt"], source="folder"
    )
    assert len(entries) == 1
    assert load._line_from_path("root/L1/a.edi") == "L1"
    assert load._infer_lines(["root/L1/a.edi", "root/L1/b.edi"]) == {"L1": 2}
    assert load._filtered_entries(entries, ["L1"]) == entries
    assert load._filtered_entries(entries, ["L2"]) == []
    assert "Add Lines" in str(load._modal_title("append"))
    assert "merged" in load._mode_hint("append")

    app = _capture(load, "register_load")
    monkeypatch.setattr(load, "ctx", SimpleNamespace(triggered_id=IDs.BTN_ADD_LINE))
    assert app.get("open_modal")(1, 0, 1)[1] == "append"
    monkeypatch.setattr(load, "ctx", SimpleNamespace(triggered_id=IDs.FOLDER_STORE))
    staged = {"contents": ["x"], "filenames": ["root/L1/a.edi"]}
    assert app.get("capture")(None, staged, None)[0]["source"] == "folder"
    assert app.get("populate")([])[2] == {"display": "none"}
    multi = load._entries(
        ["a", "b"], ["root/L1/a.edi", "root/L2/b.edi"], source="folder"
    )
    opts, selected, style = app.get("populate")(multi)
    assert len(opts) == 2 and selected == ["L1", "L2"]
    assert "No files" in str(app.get("preflight")([], [], "replace"))
    assert "2 file" in str(app.get("preflight")(multi, selected, "append"))

    confirm = app.get("confirm")
    assert confirm(1, [], [], "sid", "replace", None, None)[0] is no_update
    assert confirm(1, multi, selected, "", "replace", None, None)[0] is no_update
    fake_view = SimpleNamespace(n_stations=2)
    old_view = SimpleNamespace(n_stations=1)
    monkeypatch.setattr(load, "_decode", lambda *a: ({"L1": ["a"]}, None))
    monkeypatch.setattr(load, "_build_view", lambda *a, **k: fake_view)
    monkeypatch.setattr(load, "get_view", lambda _sid: old_view)
    monkeypatch.setattr(load, "merge_views", lambda old, new: new)
    monkeypatch.setattr(load, "set_view", lambda *a: None)
    monkeypatch.setattr(
        load,
        "store_from_view",
        lambda _view: {"n_stations": 2, "n_lines": 1},
    )
    result = confirm(1, multi, selected, "sid", "append", None, None)
    assert result[0]["n_stations"] == 2 and result[2] is False
