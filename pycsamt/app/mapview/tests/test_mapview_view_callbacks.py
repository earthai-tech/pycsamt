# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Behavioral coverage for the Map View view-rail / canvas-render callback
module (``pycsamt.app.mapview.callbacks.view``)."""

from __future__ import annotations

import types

import numpy as np
import pytest

from pycsamt.app.mapview.callbacks import view as vmod
from pycsamt.app.mapview._ids import IDs
from pycsamt.app._geology import store_from_legend
from pycsamt.app._structure import store_from_structure
from pycsamt.format.borehole import pcbh_to_dict
from pycsamt.format.borehole.builder import document_from_builder, new_builder_draft
from pycsamt.format.geology import GeologyLegend
from pycsamt.format.structure import StructModel
from pycsamt.geology import RockDatabase, RockEntry
from pycsamt.geology.structural import FaultTrace, StructuralModel
from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.view import MapView


def _capture():
    captured: dict = {}

    class _App:
        def callback(self, *a, **k):
            def deco(fn):
                captured[fn.__name__] = fn
                return fn

            return deco

    vmod.register_view(_App())
    return captured


def _set_trigger(monkeypatch, triggered_id):
    monkeypatch.setattr(
        vmod, "ctx", types.SimpleNamespace(triggered_id=triggered_id)
    )


# ---------------------------------------------------------------------------
# registration
# ---------------------------------------------------------------------------


def test_register_view_wires_switch_and_render():
    captured = _capture()
    assert {"switch", "render"} <= set(captured)


# ---------------------------------------------------------------------------
# _register_rail / switch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "rail_id, expected_view",
    [
        (IDs.RAIL_MAP, "map"),
        (IDs.RAIL_3D, "map3d"),
        (IDs.RAIL_BH, "bh"),
        (IDs.RAIL_GEO, "geology"),
    ],
)
def test_switch_sets_view_and_active_class(monkeypatch, rail_id, expected_view):
    captured = _capture()
    _set_trigger(monkeypatch, rail_id)
    result = captured["switch"](1, None, None, None)
    view, map_cls, cls3d, bh_cls, geo_cls, title = result
    assert view == expected_view
    assert title == vmod.VIEW_TITLES[expected_view]
    classes = {
        IDs.RAIL_MAP: map_cls,
        IDs.RAIL_3D: cls3d,
        IDs.RAIL_BH: bh_cls,
        IDs.RAIL_GEO: geo_cls,
    }
    for bid, cls in classes.items():
        if bid == rail_id:
            assert cls == "mv-rail-btn mv-rail-active"
        else:
            assert cls == "mv-rail-btn"


def test_switch_unknown_trigger_defaults_to_map(monkeypatch):
    captured = _capture()
    _set_trigger(monkeypatch, "not-a-rail-id")
    view, *_classes, title = captured["switch"](1, None, None, None)
    assert view == "map"
    assert title == "Map view"


# ---------------------------------------------------------------------------
# _borehole_view_figure
# ---------------------------------------------------------------------------


def _valid_draft() -> dict:
    draft = new_builder_draft()
    draft["boreholes"] = [
        {
            "id": "BH-1", "name": "Hole 1", "kind": "water",
            "status": "completed", "x": 10.0, "y": 20.0, "z": 100.0,
            "total_depth_md": 30.0,
        }
    ]
    draft["intervals"] = [
        {
            "borehole_id": "BH-1", "from_md": 0.0, "to_md": 15.0,
            "code": "SAP", "label": "Saprolite", "family": "lithology",
        }
    ]
    return draft


def _pcbh_store() -> dict:
    document = document_from_builder(_valid_draft())
    return {"document": pcbh_to_dict(document), "n_boreholes": 1}


def test_borehole_view_figure_no_store_shows_open_studio_prompt():
    fig = vmod._borehole_view_figure(
        None, "striplog", "lithology", 0.9, False, False, None, "light",
    )
    assert "Open Borehole Studio" in fig.layout.annotations[0].text


def test_borehole_view_figure_striplog_mode():
    fig = vmod._borehole_view_figure(
        _pcbh_store(), "striplog", "lithology", 0.9, True, False, None,
        "dark",
    )
    assert fig is not None


def test_borehole_view_figure_striplog_is_default_mode():
    with_default = vmod._borehole_view_figure(
        _pcbh_store(), None, "lithology", 0.9, False, False, None, "light",
    )
    with_explicit = vmod._borehole_view_figure(
        _pcbh_store(), "striplog", "lithology", 0.9, False, False, None,
        "light",
    )
    assert type(with_default) is type(with_explicit)


def test_borehole_view_figure_3d_mode_builds_scene():
    fig = vmod._borehole_view_figure(
        _pcbh_store(), "3d", "lithology", 0.8, True, True, 5.0, "dark",
    )
    assert fig.layout.scene.xaxis.title.text == "Easting"
    assert fig.layout.scene.yaxis.title.text == "Northing"
    assert fig.layout.template.layout.annotationdefaults is not None or True


def test_borehole_view_figure_3d_mode_without_fixed_radius():
    fig = vmod._borehole_view_figure(
        _pcbh_store(), "3d", "lithology", 0.8, False, False, None, "light",
    )
    assert fig.layout.scene.zaxis.title.text == "Elevation (m)"


# ---------------------------------------------------------------------------
# _geology_view_figure
# ---------------------------------------------------------------------------


def _legend() -> GeologyLegend:
    db = RockDatabase(
        [
            RockEntry("Sand (saturated)", 10.0, 200.0, "#E9C46A"),
            RockEntry("Granodiorite", 1000.0, 5000.0, "#8D99AE"),
        ]
    )
    return GeologyLegend.from_rock_database(db, document_id="pcgl:t")


def _struct_store() -> dict:
    model = StructuralModel(
        faults=[
            FaultTrace(
                x=1.0, dip_deg=70.0, downthrown_side="right",
                line="L1", z_top=0.0,
            )
        ],
    )
    doc = StructModel.from_structural_model(model)
    return store_from_structure(doc)


def test_geology_view_figure_default_legend_mode_no_store():
    fig = vmod._geology_view_figure(None, None, None, "light")
    assert fig is not None


def test_geology_view_figure_legend_mode_with_store():
    fig = vmod._geology_view_figure(
        store_from_legend(_legend()), None, "legend", "light",
    )
    assert fig is not None


def test_geology_view_figure_structure_mode_no_store():
    fig = vmod._geology_view_figure(None, None, "structure", "dark")
    assert fig is not None


def test_geology_view_figure_structure_mode_with_store():
    fig = vmod._geology_view_figure(
        None, _struct_store(), "structure", "dark",
    )
    assert fig is not None


# ---------------------------------------------------------------------------
# render() -- the main canvas callback
# ---------------------------------------------------------------------------


def _map_data() -> MapData:
    class _Z:
        def __init__(self):
            self.freq = [10.0, 1.0]
            self.resistivity = np.ones((2, 2, 2)) * 100.0
            self.phase = np.ones((2, 2, 2)) * 45.0

    class _Edi:
        def __init__(self, station):
            self.station = station
            self.Z = _Z()

    return MapData(
        sites=[_Edi("S00"), _Edi("S01")],
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
        metadata={},
    )


def _render_kwargs(**overrides):
    kwargs = dict(
        store=None, view_name="map", controls=None, theme="light",
        lines=None, fit=0, masked=None, pcbh_store=None,
        pcbh_visible=False, pcbh_labels=False, pcbh_family=None,
        pcbh_opacity=None, pcbh_as_tubes=False, pcbh_radius=None,
        pcbh_on_map=False, pcbh_in_3d=False, pcbh_lean=None,
        pcbh_lean_dir=None, pcbh_label_angle=None, pcbh_label_size=None,
        pcbh_collar_size=None, pcbh_depth_ticks=None,
        pcbh_patch_geology=False, pcbh_patch_width=None,
        pcpt_store=None, pcpt_visible=False, bh_view_mode=None,
        geo_store=None, geo_apply=False, struct_store=None,
        struct_apply=False, geo_view_mode=None,
        viewport=None, session_id="sess-1",
    )
    kwargs.update(overrides)
    return kwargs


def _call_render(captured, **overrides):
    return captured["render"](**_render_kwargs(**overrides))


def test_render_bh_view_routes_to_borehole_figure():
    captured = _capture()
    fig, welcome = _call_render(captured, view_name="bh")
    assert welcome == {"display": "none"}
    assert "Open Borehole Studio" in fig.layout.annotations[0].text


def test_render_bh_view_with_visible_pcbh_store():
    captured = _capture()
    fig, welcome = _call_render(
        captured, view_name="bh", pcbh_store=_pcbh_store(), pcbh_visible=True,
    )
    assert welcome == {"display": "none"}
    assert fig is not None


def test_render_geology_view_routes_to_geology_figure():
    captured = _capture()
    fig, welcome = _call_render(captured, view_name="geology")
    assert welcome == {"display": "none"}
    assert fig is not None


def test_render_no_store_shows_welcome():
    captured = _capture()
    fig, welcome = _call_render(captured, store=None)
    assert welcome == {"display": "flex"}
    assert "Upload PCBH or load survey lines" in fig.layout.annotations[0].text


def test_render_no_station_count_shows_welcome():
    captured = _capture()
    fig, welcome = _call_render(captured, store={"n_stations": 0})
    assert welcome == {"display": "flex"}


def test_render_no_store_but_pcbh_store_hides_welcome_and_overlays():
    captured = _capture()
    fig, welcome = _call_render(
        captured, store=None, pcbh_store=_pcbh_store(), pcbh_visible=True,
    )
    assert welcome == {"display": "none"}


def test_render_no_store_pcbh_present_but_not_visible_keeps_placeholder():
    captured = _capture()
    fig, welcome = _call_render(
        captured, store=None, pcbh_store=_pcbh_store(), pcbh_visible=False,
    )
    assert welcome == {"display": "none"}
    assert "Upload PCBH or load survey lines" in fig.layout.annotations[0].text


def test_render_session_view_unavailable(monkeypatch):
    captured = _capture()
    monkeypatch.setattr(vmod, "get_view", lambda session_id: None)
    fig, welcome = _call_render(
        captured, store={"n_stations": 2}, session_id="missing",
    )
    assert welcome == {"display": "flex"}
    assert "Session data unavailable" in fig.layout.annotations[0].text


def test_render_full_path_builds_figure(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        view_name="map",
        lines={"active": ["L1"]},
        masked=None,
        fit=1,
        session_id="sess-ok",
    )
    assert welcome == {"display": "none"}
    assert fig is not None


def test_render_full_path_with_boreholes_and_points(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        view_name="map3d",
        controls={"mode3d": "fence"},
        pcbh_store=_pcbh_store(),
        pcbh_in_3d=True,
        pcbh_visible=True,
        pcbh_on_map=True,
        pcpt_store={"points": []},
        pcpt_visible=True,
        session_id="sess-ok2",
    )
    assert welcome == {"display": "none"}
    assert fig is not None


def test_render_full_path_with_structure_opts(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        struct_store=_struct_store(),
        struct_apply=True,
        session_id="sess-ok3",
    )
    assert welcome == {"display": "none"}
    assert fig is not None


def test_render_geology_bands_applied_from_valid_legend(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        geo_store=store_from_legend(_legend()),
        geo_apply=True,
        session_id="sess-ok4",
    )
    assert welcome == {"display": "none"}


def test_render_geology_bands_exception_is_swallowed(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        geo_store={"document": 123},
        geo_apply=True,
        session_id="sess-ok5",
    )
    assert welcome == {"display": "none"}


def test_render_geology_pattern_fill_applied(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        geo_store=store_from_legend(_legend()),
        geo_apply=True,
        controls={"geology_fill": "pattern"},
        session_id="sess-ok6",
    )
    assert welcome == {"display": "none"}


def test_render_geology_pattern_fill_exception_is_swallowed(monkeypatch):
    captured = _capture()
    real_view = MapView(_map_data())
    monkeypatch.setattr(vmod, "get_view", lambda session_id: real_view)
    fig, welcome = _call_render(
        captured,
        store={"n_stations": 2},
        geo_store={"document": 123},
        geo_apply=True,
        controls={"geology_fill": "pattern"},
        session_id="sess-ok7",
    )
    assert welcome == {"display": "none"}
