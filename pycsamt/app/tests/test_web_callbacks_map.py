# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.map (station map, overlays, contours).

Strategy
--------
``build_station_map``/``build_contour_overlay`` (pycsamt.app.web.utils)
are real, fast Plotly/matplotlib figure builders driven by plain
DataFrames -- exercised for real with WILLY station coordinates rather
than mocked. The one genuinely external/heavy dependency is a full
Occam2D/ModEM inversion run, so the "inversion overlay" map-type is
exercised against a small duck-typed fake result object (same shape the
real ``InversionResult.to_resistivity_model()`` produces) instead of
running a real inversion.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from dash import no_update

import pycsamt.app.web.callbacks.map as map_mod
from pycsamt.app.web.cache import cache_set, cache_set_inversion_result
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture(scope="session")
def willy_records(willy_sites):
    recs = []
    for i, ed in enumerate(willy_sites.as_list()):
        summary = willy_sites.get(ed.station).summary()
        recs.append(
            {
                "ID": ed.station,
                "Latitude": summary["lat"],
                "Longitude": summary["lon"],
                "Elevation": summary["elev"],
                "Line": "L1" if i % 2 == 0 else "L2",
            }
        )
    return recs


@pytest.fixture
def store_data_willy(willy_records):
    return {
        "station_records": willy_records,
        "n_stations": len(willy_records),
        "n_lines": 2,
    }


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-map-session"
    cache_set(session_id, willy_sites)
    return session_id


class _FakeResModel:
    def __init__(self, rho_2d, z_centers, station_names):
        self.rho_2d = rho_2d
        self.z_centers = z_centers
        self.station_names = station_names


class _FakeInvResult:
    def __init__(self, model, method="occam2d", dimension="2d"):
        self._model = model
        self.method = method
        self.dimension = dimension

    def to_resistivity_model(self):
        return self._model


@pytest.fixture
def cached_inv_result(willy_records):
    n_x = len(willy_records)
    n_z = 5
    rho_2d = np.random.default_rng(0).uniform(1.0, 500.0, size=(n_z, n_x))
    z_centers = np.linspace(10, 1000, n_z)
    station_names = [r["ID"] for r in willy_records]
    model = _FakeResModel(rho_2d, z_centers, station_names)
    session_id = "test-map-inv-session"
    cache_set_inversion_result(session_id, _FakeInvResult(model))
    return session_id


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(k for k in web_app.callback_map if all(s in k for s in substrings))
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(AttributeDict(triggered_inputs=[{"prop_id": prop_id}]))


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


# ── Pure helpers ─────────────────────────────────────────────────────────────


class TestStationDataframe:
    def test_no_records_returns_empty(self):
        assert map_mod._station_dataframe(None).empty
        assert map_mod._station_dataframe({}).empty

    def test_alias_mapping(self):
        store = {
            "station_records": [
                {"station": "S1", "lat": 1.0, "long": 2.0, "profile": "P1"}
            ]
        }
        df = map_mod._station_dataframe(store)
        assert {"ID", "Latitude", "Longitude", "Line"}.issubset(df.columns)
        assert df["ID"].iloc[0] == "S1"

    def test_line_default_dash_for_missing(self):
        store = {"station_records": [{"ID": "S1", "Line": ""}]}
        df = map_mod._station_dataframe(store)
        assert df["Line"].iloc[0] == "—"

    def test_numeric_coercion(self):
        store = {
            "station_records": [{"ID": "S1", "Latitude": "bad", "Longitude": "1.5"}]
        }
        df = map_mod._station_dataframe(store)
        assert pd.isna(df["Latitude"].iloc[0])
        assert df["Longitude"].iloc[0] == 1.5


class TestUsableLineFilter:
    def test_no_filter_returns_none(self):
        df = pd.DataFrame({"Line": ["L1", "L2"]})
        assert map_mod._usable_line_filter(None, df) is None
        assert map_mod._usable_line_filter([], df) is None

    def test_no_line_column_returns_none(self):
        df = pd.DataFrame({"ID": ["S1"]})
        assert map_mod._usable_line_filter(["L1"], df) is None

    def test_keeps_only_matching_lines(self):
        df = pd.DataFrame({"Line": ["L1", "L2", "L1"]})
        assert map_mod._usable_line_filter(["L1", "L9"], df) == ["L1"]

    def test_no_match_returns_none(self):
        df = pd.DataFrame({"Line": ["L1"]})
        assert map_mod._usable_line_filter(["ZZ"], df) is None


class TestLayoutHasId:
    def test_none_layout(self):
        assert map_mod._layout_has_id(None, "x") is False

    def test_finds_nested_id(self):
        class _Node:
            def __init__(self, id_=None, children=None):
                self.id = id_
                self.children = children

        tree = _Node(children=[_Node(id_="a"), _Node(children=[_Node(id_="target")])])
        assert map_mod._layout_has_id(tree, "target") is True
        assert map_mod._layout_has_id(tree, "missing") is False

    def test_list_layout(self):
        class _Node:
            def __init__(self, id_=None, children=None):
                self.id = id_
                self.children = children

        assert map_mod._layout_has_id([_Node(id_="a"), _Node(id_="b")], "b") is True


class TestResolveCrsInfo:
    def test_geo(self):
        assert "4326" in map_mod._resolve_crs_info("geo")

    def test_utm_north(self):
        out = map_mod._resolve_crs_info("utm", zone=50, hem="N")
        assert "32650" in out

    def test_utm_south(self):
        out = map_mod._resolve_crs_info("utm", zone=50, hem="S")
        assert "32750" in out

    def test_custom_valid_epsg(self):
        out = map_mod._resolve_crs_info("custom", epsg="4326")
        assert "4326" in out

    def test_custom_invalid_epsg_falls_back(self):
        out = map_mod._resolve_crs_info("custom", epsg="not-a-code")
        assert "cannot validate" in out


class TestExtractRhoAtFreq:
    def test_real_extraction(self, willy_sites):
        ed = willy_sites.as_list()[0]
        f0 = float(np.asarray(ed.Z.freq, float)[0])
        result = map_mod._extract_rho_at_freq(willy_sites, f0, "xy")
        assert isinstance(result, dict)
        assert len(result) > 0
        assert all(v > 0 for v in result.values())

    def test_det_component(self, willy_sites):
        ed = willy_sites.as_list()[0]
        f0 = float(np.asarray(ed.Z.freq, float)[0])
        result = map_mod._extract_rho_at_freq(willy_sites, f0, "det")
        assert isinstance(result, dict)


class TestExtractDepthAtFreq:
    def test_skin_depth_mode(self, willy_sites):
        ed = willy_sites.as_list()[0]
        f0 = float(np.asarray(ed.Z.freq, float)[0])
        result = map_mod._extract_depth_at_freq(willy_sites, f0, "xy")
        assert isinstance(result, dict)
        assert len(result) > 0

    def test_target_depth_mode_returns_period(self, willy_sites):
        ed = willy_sites.as_list()[0]
        f0 = float(np.asarray(ed.Z.freq, float)[0])
        result = map_mod._extract_depth_at_freq(
            willy_sites, f0, "xy", target_depth_m=100.0
        )
        assert isinstance(result, dict)
        assert len(result) > 0


# ── Callbacks ─────────────────────────────────────────────────────────────


class TestHomeMap:
    """update_home_map is only registered when the layout actually contains
    both HOME_MAP_GRAPH and HOME_MAP_OVERLAY (see _layout_has_id gating in
    register_map). The real app.layout doesn't currently render the
    home-page mini-map, so a standalone minimal app is built here to force
    registration and exercise the callback body directly."""

    @pytest.fixture
    def home_map_app(self):
        import dash
        from dash import dcc, html

        app = dash.Dash(__name__)
        app.layout = html.Div(
            [
                dcc.Graph(id=IDs.HOME_MAP_GRAPH),
                dcc.Dropdown(id=IDs.HOME_MAP_OVERLAY),
                dcc.Store(id=IDs.STORE_DATA),
                dcc.Store(id=IDs.STORE_SELECTION),
                dcc.Store(id=IDs.STORE_THEME),
                dcc.Graph(id=IDs.MAP_GRAPH),
                dcc.Store(id=IDs.SESSION_ID),
                dcc.Store(id=IDs.STORE_ACTIVE_LINES),
            ]
        )
        map_mod.register_map(app)
        return app

    def _fn(self, app):
        return _cb(app, f"{IDs.HOME_MAP_GRAPH}.figure")

    def test_no_store_data(self, home_map_app):
        fig = self._fn(home_map_app)(None, "Index", None, "dark")
        assert fig is not None

    def test_empty_records(self, home_map_app):
        fig = self._fn(home_map_app)({"station_records": []}, "Index", None, "dark")
        assert fig is not None

    def test_real_data(self, home_map_app, store_data_willy):
        fig = self._fn(home_map_app)(
            store_data_willy, "Index", {"station_id": "18-001A"}, "light"
        )
        assert fig is not None


class TestUpdateMap:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP_GRAPH}.figure")

    def _base_args(self, **overrides):
        args = dict(
            nav_section="map",
            store_data=None,
            overlay="Index",
            map_type="station",
            freq_sel=None,
            comp="xy",
            depth_target=None,
            selection=None,
            theme="dark",
            line_filter=None,
            marker_size=10,
            basemap_style=None,
            cmap="plasma",
            opacity_pct=90,
            show_opts=[],
            crs_mode="geo",
            utm_zone=50,
            utm_hem="N",
            epsg="4326",
            contour_enable=[],
            contour_mode="filled+lines",
            contour_cmap="jet",
            contour_levels=12,
            contour_opacity_pct=60,
            contour_interp="cubic",
            contour_extra=0.12,
            contour_smooth=1.0,
            contour_lwidth=1.0,
            contour_opts=[],
            contour_res=150,
            inv_depth_m=500,
            inv_line_sel=None,
            bearing=0,
            _refresh=0,
            _fit=0,
            session_id=None,
            active_lines_store=None,
        )
        args.update(overrides)
        return args

    def _call(self, web_app, trigger=f"{IDs.STORE_DATA}.data", **overrides):
        fn = self._fn(web_app)
        _set_triggered(trigger)
        try:
            return fn(**self._base_args(**overrides))
        finally:
            _clear_triggered()

    def test_nav_skip_when_leaving_map(self, web_app):
        out = self._call(web_app, trigger=f"{IDs.NAV_SECTION}.data", nav_section="home")
        assert out is no_update

    def test_no_store_data(self, web_app):
        fig = self._call(web_app)
        assert fig is not None

    def test_empty_df(self, web_app):
        fig = self._call(web_app, store_data={"station_records": []})
        assert fig is not None

    def test_station_type_basic(self, web_app, store_data_willy):
        fig = self._call(web_app, store_data=store_data_willy)
        assert fig is not None

    def test_elevation_type(self, web_app, store_data_willy):
        fig = self._call(web_app, store_data=store_data_willy, map_type="elevation")
        assert fig is not None

    def test_resistivity_type(
        self, web_app, store_data_willy, cached_session, willy_sites
    ):
        ed = willy_sites.as_list()[0]
        f0 = str(float(np.asarray(ed.Z.freq, float)[0]))
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            map_type="resistivity",
            freq_sel=f0,
            comp="xy",
            session_id=cached_session,
        )
        assert fig is not None

    def test_depth_type(self, web_app, store_data_willy, cached_session, willy_sites):
        ed = willy_sites.as_list()[0]
        f0 = str(float(np.asarray(ed.Z.freq, float)[0]))
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            map_type="depth",
            freq_sel=f0,
            comp="xy",
            depth_target=100.0,
            session_id=cached_session,
        )
        assert fig is not None

    def test_resistivity_type_no_session_falls_back(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            map_type="resistivity",
            freq_sel="1000.0",
            session_id=None,
        )
        assert fig is not None

    def test_inv_profile_type(self, web_app, store_data_willy, cached_inv_result):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            map_type="inv_profile",
            session_id=cached_inv_result,
        )
        assert fig is not None

    def test_inv_depth_slice_type_forces_contour(
        self, web_app, store_data_willy, cached_inv_result
    ):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            map_type="inv_depth_slice",
            inv_depth_m=200,
            session_id=cached_inv_result,
        )
        assert fig is not None

    def test_inv_type_no_result_falls_back(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            map_type="inv_profile",
            session_id="no-such-session",
        )
        assert fig is not None

    def test_crs_utm_reprojection(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            crs_mode="utm",
            utm_zone=50,
            utm_hem="N",
        )
        assert fig is not None

    def test_crs_custom_reprojection(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            crs_mode="custom",
            epsg="32650",
        )
        assert fig is not None

    def test_contour_overlay_enabled(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            contour_enable=["on"],
        )
        assert fig is not None

    def test_line_filter_and_muted_lines(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            line_filter=["L1"],
            active_lines_store={"active": ["L1"], "all": ["L1", "L2"]},
        )
        assert fig is not None

    def test_show_labels_and_profiles(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            store_data=store_data_willy,
            show_opts=["labels", "profiles"],
        )
        assert fig is not None

    def test_fit_button_triggers_uirevision_reset(self, web_app, store_data_willy):
        fig = self._call(
            web_app,
            trigger="map-tb-fit.n_clicks",
            store_data=store_data_willy,
            _fit=3,
        )
        assert fig is not None


class TestPopulateFreqOptions:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP_PAGE_FREQ_SEL}.options")

    def test_no_store_data(self, web_app):
        opts, val = self._fn(web_app)(None, "hz", "s")
        assert opts == []
        assert val is None

    def test_no_session(self, web_app, store_data_willy):
        opts, val = self._fn(web_app)(store_data_willy, "hz", None)
        assert opts == []

    def test_hz_units(self, web_app, store_data_willy, cached_session):
        opts, val = self._fn(web_app)(store_data_willy, "hz", cached_session)
        assert opts
        assert "Hz" in opts[0]["label"]

    def test_period_units(self, web_app, store_data_willy, cached_session):
        opts, val = self._fn(web_app)(store_data_willy, "period", cached_session)
        assert opts
        assert "s" in opts[0]["label"]


class TestUpdateDepthInfo:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.MAP_PAGE_DEPTH_INFO}.children")

    def test_missing_inputs_returns_empty(self, web_app):
        assert self._fn(web_app)(None, None, "s") == ""
        assert self._fn(web_app)("1000", None, "s") == ""
        assert self._fn(web_app)("1000", "100", None) == ""

    def test_valid_inputs(self, web_app):
        out = self._fn(web_app)("1000", "100", "s")
        assert "Hz" in out
        assert "δ" in out


class TestPopulateInvLines:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP_PAGE_INV_LINE_SEL}.options")

    def test_no_session(self, web_app):
        opts, val = self._fn(web_app)(None, None)
        assert opts == []
        assert val is None

    def test_no_result_cached(self, web_app):
        opts, val = self._fn(web_app)(None, "no-such-session")
        assert opts[0]["disabled"] is True

    def test_result_cached(self, web_app, cached_inv_result):
        opts, val = self._fn(web_app)(None, cached_inv_result)
        assert val == "default"
        assert "occam2d" in opts[0]["label"].lower()


class TestUpdateInvDepthInfo:
    def _fn(self, web_app):
        return _cb(web_app, "map-inv-depth-info.children")

    def test_missing_inputs(self, web_app):
        assert self._fn(web_app)(None, "s") == ""
        assert self._fn(web_app)("100", None) == ""

    def test_no_result_cached(self, web_app):
        out = self._fn(web_app)("100", "no-such-session")
        assert "No inversion result cached." in out

    def test_result_cached(self, web_app, cached_inv_result):
        out = self._fn(web_app)("200", cached_inv_result)
        assert "Layer" in out
        assert "stations" in out


class TestUpdateCrsInfo:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.MAP_PAGE_CRS_INFO}.children")

    def test_geo_default(self, web_app):
        out = self._fn(web_app)(None, None, None, None)
        assert "4326" in out

    def test_utm(self, web_app):
        out = self._fn(web_app)("utm", 50, "N", "4326")
        assert "UTM" in out


class TestSelectionSync:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.STORE_SELECTION}.data", IDs.MAP_GRAPH)

    # NOTE: the real app registers ``sync_selection`` with
    # ``include_home_map=False`` (the home mini-map isn't in the current
    # layout -- see TestHomeMap), so the registered callback takes exactly
    # 4 positional args: (click_data, selected_rows, _row_clicks, store_data).

    def test_map_click(self, web_app):
        _set_triggered(f"{IDs.MAP_GRAPH}.clickData")
        try:
            cd = {"points": [{"customdata": "18-001A"}]}
            out = self._fn(web_app)(cd, [], [], None)
            assert out == {"station_id": "18-001A"}
        finally:
            _clear_triggered()

    def test_table_row_selected(self, web_app, store_data_willy):
        _set_triggered(f"{IDs.STATION_TABLE}.selected_rows")
        try:
            out = self._fn(web_app)(None, [0], [], store_data_willy)
            assert out == {"station_id": store_data_willy["station_records"][0]["ID"]}
        finally:
            _clear_triggered()

    def test_row_click_pattern_id(self, web_app):
        _set_triggered('{"type":"sta-row","index":"18-003A"}.n_clicks')
        try:
            out = self._fn(web_app)(None, [], [1], None)
            assert out == {"station_id": "18-003A"}
        finally:
            _clear_triggered()

    def test_no_matching_trigger_no_update(self, web_app):
        _set_triggered("")
        try:
            out = self._fn(web_app)(None, [], [], None)
            assert out is no_update
        finally:
            _clear_triggered()


class TestUpdateStatusStation:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.STATUS_STATION}.children")

    def test_no_selection(self, web_app):
        assert self._fn(web_app)(None, None) == "No station selected"

    def test_no_station_id(self, web_app):
        assert self._fn(web_app)({}, None) == "No station selected"

    def test_match_with_coords(self, web_app, store_data_willy):
        sid = store_data_willy["station_records"][0]["ID"]
        out = self._fn(web_app)({"station_id": sid}, store_data_willy)
        assert sid in out
        assert "°N" in out

    def test_no_match_falls_back_to_bare_label(self, web_app, store_data_willy):
        out = self._fn(web_app)({"station_id": "UNKNOWN"}, store_data_willy)
        assert out == "Station: UNKNOWN"


class TestToolbarInfo:
    def _fn(self, web_app):
        return _cb(web_app, "map-tb-info.children")

    def test_no_data(self, web_app):
        out = self._fn(web_app)(None)
        assert "No data loaded" in out.children

    def test_with_data(self, web_app, store_data_willy):
        out = self._fn(web_app)(store_data_willy)
        texts = [getattr(c, "children", c) for c in out]
        assert any(str(store_data_willy["n_stations"]) in str(t) for t in texts)


# ── _inject_contour direct test ──────────────────────────────────────────────


class TestInjectContour:
    def test_too_few_points_returns_fig_unchanged(self, willy_records):
        import plotly.graph_objects as go

        fig = go.Figure()
        df = pd.DataFrame(willy_records[:2])
        out = map_mod._inject_contour(
            fig,
            df,
            "Index",
            None,
            contour_mode="filled+lines",
            contour_cmap="jet",
            n_levels=12,
            c_opacity=0.6,
            interp_method="cubic",
            extra_factor=0.12,
            smooth_sigma=1.0,
            line_width=1.0,
            contour_opts=[],
            grid_res=50,
            dark=True,
        )
        assert out is fig

    def test_real_contour_injects_layer(self, willy_records):
        import plotly.graph_objects as go

        fig = go.Figure()
        fig.update_layout(map={"style": "carto-darkmatter"})
        df = pd.DataFrame(willy_records)
        out = map_mod._inject_contour(
            fig,
            df,
            "Index",
            None,
            contour_mode="filled+lines",
            contour_cmap="jet",
            n_levels=8,
            c_opacity=0.6,
            interp_method="linear",
            extra_factor=0.12,
            smooth_sigma=1.0,
            line_width=1.0,
            contour_opts=[],
            grid_res=40,
            dark=True,
        )
        assert out.layout.map.layers
