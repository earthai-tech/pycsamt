# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.tools -- strike / validator / converter /
coords tool panels, plus the shared helpers near the top of the module that
all tool panels depend on.

Strategy
--------
Real domain objects are exercised wherever practical instead of mocks:
``estimate_strike_*`` (``pycsamt.emtools.strike``) run against real WILLY
L18PLT EDI stations (``data/AMT/WILLY_DATA/L18PLT``, 28 stations), the format
converter runs a real AVG -> EDI conversion via ``ConversionController``
against ``data/avg/K1.AVG`` + ``K1.stn``, and the coordinate transformer uses
real ``pyproj``/GDAL-backed projections from ``pycsamt.gis.utils``. Only
Dash's own dispatch context (``callback_context`` / ``set_props``) is faked,
since it genuinely requires a live request to exist.
"""

from __future__ import annotations

import base64
from pathlib import Path

import pandas as pd
import pytest
from dash import html, no_update
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.tools as tools_mod
from pycsamt.app.web.cache import cache_get, cache_set

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))
_K1_AVG = _ROOT / "data" / "avg" / "K1.AVG"
_K1_STN = _ROOT / "data" / "avg" / "K1.stn"
_HAS_K1 = _K1_AVG.exists() and _K1_STN.exists()


# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


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


def _raises_prevent_update():
    return pytest.raises(PreventUpdate)


# ── EDI fixtures ─────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture(scope="session")
def willy_ids(willy_sites):
    from pycsamt.emtools._core import _iter_items, _name

    return [_name(ed, i) for i, ed in enumerate(_iter_items(willy_sites))]


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-tools-strike-session"
    cache_set(session_id, willy_sites)
    return session_id


@pytest.fixture
def store_data_willy(willy_ids):
    half = len(willy_ids) // 2
    records = [
        {"ID": sid, "Line": "L1" if i < half else "L2"}
        for i, sid in enumerate(willy_ids)
    ]
    return {
        "station_records": records,
        "n_stations": len(records),
        "n_lines": 2,
    }


ACTIVE_ALL = {"active": ["L1", "L2"], "all": ["L1", "L2"]}


def _records(pairs):
    return [{"ID": sid, "Line": ln} for sid, ln in pairs]


# ═══════════════════════════════════════════════════════════════════════════
# 1. Shared helpers (lines ~49-561)
# ═══════════════════════════════════════════════════════════════════════════


class TestStyledTable:
    def test_builds_datatable_from_dataframe(self):
        df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
        tbl = tools_mod._styled_table(df)
        assert tbl.columns == [
            {"name": "a", "id": "a"},
            {"name": "b", "id": "b"},
        ]
        assert tbl.data == [{"a": 1, "b": "x"}, {"a": 2, "b": "y"}]


class TestNoDataMsg:
    def test_contains_hint_text(self):
        msg = tools_mod._no_data_msg()
        assert "Load survey data first." in msg.children


class TestRunBtn:
    def test_defaults(self):
        btn = tools_mod._run_btn("my-btn", "Run It")
        assert btn.id == "my-btn"
        assert btn.n_clicks == 0
        assert "btn-primary" in btn.className
        assert "Run It" in btn.children

    def test_custom_icon(self):
        btn = tools_mod._run_btn("id2", "Go", icon="rocket")
        assert "bi-rocket" in btn.children[0].className


class TestActionBar:
    def test_wraps_children(self):
        bar = tools_mod._action_bar(html.Span("a"), html.Span("b"))
        assert bar.className == "tool-action-bar"
        assert len(bar.children) == 2


class TestPlotlyConfig:
    def test_default_filename(self):
        cfg = tools_mod._plotly_config()
        assert cfg["toImageButtonOptions"]["filename"] == "pycsamt_plot"
        assert cfg["displaylogo"] is False

    def test_custom_filename(self):
        cfg = tools_mod._plotly_config("custom_name")
        assert cfg["toImageButtonOptions"]["filename"] == "custom_name"


class TestOutArea:
    def test_default_no_children(self):
        area = tools_mod._out_area("out-1")
        assert area.id == "out-1"
        assert area.children is None
        assert area.style == {"minHeight": "160px"}

    def test_explicit_children_used(self):
        area = tools_mod._out_area("out-2", children=html.Div("hi"))
        assert area.children.children == "hi"

    def test_last_output_restores(self):
        stored = {"type": "text", "content": "hello", "cls": "small"}
        area = tools_mod._out_area("out-3", last_output=stored)
        assert area.children.children == "hello"

    def test_children_takes_priority_over_last_output(self):
        explicit = html.Div("explicit")
        area = tools_mod._out_area(
            "out-4", children=explicit, last_output={"type": "text"}
        )
        assert area.children is explicit


class TestRestoreFromStore:
    def test_none_returns_none(self):
        assert tools_mod._restore_from_store(None) is None

    def test_empty_dict_returns_none(self):
        assert tools_mod._restore_from_store({}) is None

    def test_unknown_type_returns_none(self):
        assert tools_mod._restore_from_store({"type": "bogus"}) is None

    def test_fig_json(self):
        import plotly.graph_objects as go

        fig = go.Figure(go.Scatter(x=[1, 2], y=[3, 4]))
        stored = {"type": "fig_json", "fig": fig.to_json()}
        out = tools_mod._restore_from_store(stored)
        assert out.__class__.__name__ == "Graph"

    def test_fig_json_malformed_returns_none(self):
        # missing "fig" key -> KeyError swallowed by the broad except
        out = tools_mod._restore_from_store({"type": "fig_json"})
        assert out is None

    def test_strike_payload_delegates_to_renderer(self):
        payload = {
            "scope_label": "All",
            "n_stations": 1,
            "median_strike": 10.0,
            "median_iqr": 2.0,
            "records": [{"station": "S1", "line": "L1", "ang_axial": 10.0}],
            "table": [{"Station": "S1"}],
            "columns": ["Station"],
            "page_size": 10,
            "method_label": "Consensus strike",
        }
        out = tools_mod._restore_from_store(
            {"type": "strike_payload", "payload": payload, "theme": "dark"}
        )
        assert isinstance(out, html.Div)

    def test_validator_rows_empty_returns_none(self):
        out = tools_mod._restore_from_store(
            {"type": "validator_rows", "rows": [], "cols": []}
        )
        assert out is None

    def test_validator_rows_builds_chips_and_table(self):
        stored = {
            "type": "validator_rows",
            "rows": [{"Station": "S1", "Status": "PASS"}],
            "cols": ["Station", "Status"],
            "summary": {
                "scope": "All",
                "total": 1,
                "pass": 1,
                "warn": 0,
                "fail": 0,
            },
        }
        out = tools_mod._restore_from_store(stored)
        assert isinstance(out, html.Div)
        assert "(restored)" in str(out)

    def test_text_type(self):
        out = tools_mod._restore_from_store(
            {
                "type": "text",
                "content": "hi there",
                "cls": "small text-warning",
            }
        )
        assert out.children == "hi there"
        assert out.className == "small text-warning"

    def test_html_type_no_imgs_shows_placeholder_alert(self):
        out = tools_mod._restore_from_store({"type": "html", "imgs": []})
        assert "click Run to regenerate" in str(out)

    def test_html_type_with_imgs(self):
        out = tools_mod._restore_from_store({"type": "html", "imgs": ["AAAA", "BBBB"]})
        imgs = [c for c in out.children if getattr(c, "src", None)]
        assert len(imgs) == 2
        assert imgs[0].src.startswith("data:image/png;base64,AAAA")

    def test_lm_type_with_name_and_layers(self):
        stored = {
            "type": "lm",
            "model": {"name": "test-model", "resistivity": [10.0, 100.0]},
            "imgs": ["CCCC"],
            "solver": "mt1d",
        }
        out = tools_mod._restore_from_store(stored)
        assert isinstance(out, html.Div)
        assert "2-layer model" in str(out)

    def test_lm_type_without_name_or_imgs(self):
        stored = {"type": "lm", "model": {}, "imgs": [], "solver": ""}
        out = tools_mod._restore_from_store(stored)
        assert "0-layer model" in str(out)


class TestPathRow:
    def test_input_only(self):
        row = tools_mod._path_row("in-1", "placeholder")
        assert len(row.children) == 1

    def test_file_and_folder_buttons(self):
        row = tools_mod._path_row(
            "in-2",
            "ph",
            file_btn_id={"type": "x", "target": "in-2", "mode": "file"},
            folder_btn_id={"type": "x", "target": "in-2", "mode": "folder"},
        )
        assert len(row.children) == 3


class TestActiveLineNames:
    def test_active_present_returned_verbatim(self):
        out = tools_mod._active_line_names(
            {"station_records": _records([("A", "L1")])},
            {"active": ["L2"], "all": ["L1", "L2"]},
        )
        assert out == ["L2"]

    def test_no_active_falls_back_to_sorted_unique_lines(self):
        out = tools_mod._active_line_names(
            {"station_records": _records([("A", "L2"), ("B", "L1")])},
            None,
        )
        assert out == ["L1", "L2"]

    def test_no_data_returns_empty(self):
        assert tools_mod._active_line_names(None, None) == []


class TestStationOptionsForLines:
    STORE = {"station_records": _records([("A1", "L1"), ("A2", "L1"), ("B1", "L2")])}

    def test_all_active_no_selection(self):
        opts = tools_mod._station_options_for_lines(
            self.STORE, {"active": ["L1", "L2"], "all": ["L1", "L2"]}
        )
        assert {o["value"] for o in opts} == {"A1", "A2", "B1"}
        assert opts[0]["label"] == "A1 · L1"

    def test_active_line_filters(self):
        opts = tools_mod._station_options_for_lines(
            self.STORE, {"active": ["L1"], "all": ["L1", "L2"]}
        )
        assert {o["value"] for o in opts} == {"A1", "A2"}

    def test_selected_lines_further_narrows(self):
        opts = tools_mod._station_options_for_lines(
            self.STORE,
            {"active": ["L1", "L2"], "all": ["L1", "L2"]},
            selected_lines=["L2"],
        )
        assert {o["value"] for o in opts} == {"B1"}

    def test_falls_back_to_station_or_name_key(self):
        store = {"station_records": [{"station": "S1", "Line": "L1"}]}
        opts = tools_mod._station_options_for_lines(store, None)
        assert opts == [{"label": "S1 · L1", "value": "S1"}]


class TestStationLineMap:
    def test_maps_id_to_line(self):
        out = tools_mod._station_line_map({"station_records": _records([("A1", "L1")])})
        assert out == {"A1": "L1"}

    def test_missing_line_defaults_unassigned(self):
        out = tools_mod._station_line_map(
            {"station_records": [{"ID": "A1", "Line": ""}]}
        )
        assert out == {"A1": "Unassigned"}

    def test_no_data_returns_empty_dict(self):
        assert tools_mod._station_line_map(None) == {}


class TestFilterSitesByStationIds:
    def test_no_ids_returns_sites_unchanged(self):
        sentinel = object()
        assert tools_mod._filter_sites_by_station_ids(sentinel, None) is sentinel

    def test_exception_path_returns_sites_unchanged(self, monkeypatch):
        # A matching item is found (so the "no items" early-return isn't
        # taken); force ensure_sites() itself to raise so the broad except
        # falls back to *sites* unchanged.
        import pycsamt.emtools._core as core_mod

        class FakeEd:
            station = "match"

        def _boom(*a, **k):
            raise RuntimeError("not a real site")

        monkeypatch.setattr(core_mod, "ensure_sites", _boom)
        sites = [FakeEd()]
        out = tools_mod._filter_sites_by_station_ids(sites, ["match"])
        assert out is sites

    def test_real_subset_returns_filtered_sites(self, willy_sites, willy_ids):
        wanted = willy_ids[:3]
        out = tools_mod._filter_sites_by_station_ids(willy_sites, wanted)
        assert out is not None
        from pycsamt.emtools._core import _iter_items, _name

        got = {_name(ed, i) for i, ed in enumerate(_iter_items(out))}
        assert got == set(wanted)

    def test_real_no_match_returns_none(self, willy_sites):
        out = tools_mod._filter_sites_by_station_ids(willy_sites, ["does-not-exist"])
        assert out is None


class TestScopedSites:
    def test_session_expired(self, store_data_willy):
        sites, msg = tools_mod._scoped_sites(
            "no-such-session", store_data_willy, ACTIVE_ALL
        )
        assert sites is None
        assert "Session expired" in msg

    def test_no_active_lines_left(self, cached_session, store_data_willy):
        sites, msg = tools_mod._scoped_sites(
            cached_session,
            store_data_willy,
            {"active": ["NOPE"], "all": ["L1"]},
        )
        assert sites is None
        assert "No active survey lines" in msg

    def test_selected_lines_empty_after_filter(self, cached_session, store_data_willy):
        sites, msg = tools_mod._scoped_sites(
            cached_session,
            store_data_willy,
            ACTIVE_ALL,
            selected_lines=["NOPE"],
        )
        assert sites is None
        assert "selected line scope" in msg

    def test_selected_stations_empty_after_filter(
        self, cached_session, store_data_willy
    ):
        sites, msg = tools_mod._scoped_sites(
            cached_session,
            store_data_willy,
            ACTIVE_ALL,
            selected_stations=["does-not-exist"],
        )
        assert sites is None
        assert "selected station scope" in msg

    def test_success_returns_sites(self, cached_session, store_data_willy):
        sites, msg = tools_mod._scoped_sites(
            cached_session, store_data_willy, ACTIVE_ALL
        )
        assert sites is not None
        assert msg is None


class TestBuildToolBody:
    def test_strike_delegates(self):
        out = tools_mod._build_tool_body("strike", None)
        assert "Load survey data first." in out.children

    def test_validator_delegates(self):
        out = tools_mod._build_tool_body("validator", None)
        assert "Load survey data first." in out.children

    def test_converter_delegates(self):
        out = tools_mod._build_tool_body("converter", None)
        assert isinstance(out, html.Div)

    def test_coords_delegates(self):
        out = tools_mod._build_tool_body("coords", None)
        assert isinstance(out, html.Div)

    def test_unknown_tool_falls_back_to_placeholder(self):
        out = tools_mod._build_tool_body("not-a-real-tool", None)
        assert "tool-placeholder" in out.className

    def test_strike_with_data_builds_full_body(self, store_data_willy):
        out = tools_mod._build_tool_body(
            "strike", store_data_willy, ACTIVE_ALL, None, "dark", None
        )
        assert isinstance(out, html.Div)
        assert len(out.children) == 6  # P, row, row, row, action_bar, out_area


class TestMetricChip:
    def test_structure(self):
        chip = tools_mod._metric_chip("Label", "Value")
        assert chip.className == "tool-metric-chip"
        assert chip.children[0].children == "Label"
        assert chip.children[1].children == "Value"


# ═══════════════════════════════════════════════════════════════════════════
# 2. Strike panel
# ═══════════════════════════════════════════════════════════════════════════


class TestStrikeBody:
    def test_no_data_shows_hint(self):
        out = tools_mod._strike_body(0, tools_mod._no_data_msg(), None, None)
        assert "Load survey data first." in out.children

    def test_with_data_builds_controls(self, store_data_willy):
        out = tools_mod._strike_body(
            store_data_willy["n_stations"],
            tools_mod._no_data_msg(),
            store_data_willy,
            ACTIVE_ALL,
        )
        line_dd = out.children[1].children[0].children[1]
        assert line_dd.id == "tool-strike-lines"
        assert {o["value"] for o in line_dd.options} == {"L1", "L2"}


class TestSyncStrikeStationScope:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "tool-strike-stations.options", "tool-strike-lines"
        )

    def test_returns_options_and_resets_value(self, web_app, store_data_willy):
        opts, val = self._fn(web_app)(None, store_data_willy, ACTIVE_ALL)
        assert val is None
        assert len(opts) == store_data_willy["n_stations"]

    def test_narrows_to_selected_line(self, web_app, store_data_willy):
        opts, _ = self._fn(web_app)(["L1"], store_data_willy, ACTIVE_ALL)
        assert len(opts) == len(store_data_willy["station_records"]) // 2


class TestRunStrike:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-strike-out.children", "tool-strike-run")

    def test_no_clicks_returns_triple_no_update(self, web_app):
        out = self._fn(web_app)(
            None,
            None,
            None,
            "consensus",
            "all",
            1.0,
            40,
            None,
            None,
            "sid",
            "dark",
            None,
            None,
        )
        assert out == (no_update, no_update, no_update)

    def test_session_expired_returns_two_tuple_not_three(
        self, web_app, store_data_willy
    ):
        """
        Documents a real bug: the "scope resolution failed" early-return
        inside ``_run_strike`` (``if msg: return _warn(msg), no_update``)
        only returns 2 values, but the callback declares 3 Outputs
        (tool-strike-out.children, tools-run-store.data,
        tool-outputs-store.data). Under a real Dash dispatch this raises;
        here (calling the unwrapped function directly) it simply surfaces
        as a too-short tuple.
        """
        out = self._fn(web_app)(
            1,
            None,
            None,
            "consensus",
            "all",
            1.0,
            40,
            store_data_willy,
            ACTIVE_ALL,
            "no-such-session",
            "dark",
            None,
            None,
        )
        assert len(out) == 2
        assert "Session expired" in str(out[0])

    def test_real_run_consensus_default(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:4],
            "consensus",
            "all",
            10.0,
            40,
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            "dark",
            None,
            None,
        )
        assert len(out) == 3
        children, run_store, saved = out
        assert isinstance(children, html.Div)
        assert run_store["strike"]["method_label"] == "Consensus strike"
        assert saved["strike"]["type"] == "strike_payload"
        # 4 metric chips, rose graph, box graph, table
        assert len(children.children) == 4

    def test_real_run_sweep_method(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:3],
            "sweep",
            "high",
            10.0,
            10,
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            "light",
            None,
            None,
        )
        _, run_store, _ = out
        assert run_store["strike"]["method_label"] == "Z tensor rotation sweep"

    def test_real_run_phase_tensor_method(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:3],
            "pt",
            "all",
            10.0,
            10,
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            "dark",
            None,
            None,
        )
        _, run_store, _ = out
        assert run_store["strike"]["method_label"] == "Phase tensor strike"

    def test_empty_dataframe_result_warns(
        self, monkeypatch, web_app, cached_session, store_data_willy, willy_ids
    ):
        import pycsamt.emtools.strike as strike_mod

        monkeypatch.setattr(
            strike_mod,
            "estimate_strike_consensus",
            lambda *a, **k: pd.DataFrame(),
        )
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:2],
            "consensus",
            "all",
            5.0,
            10,
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            "dark",
            None,
            None,
        )
        # same missing-3rd-output bug as the session-expired branch
        assert len(out) == 2
        assert "No strike estimates were produced" in str(out[0])

    def test_all_nonfinite_result_warns(
        self, monkeypatch, web_app, cached_session, store_data_willy, willy_ids
    ):
        import pycsamt.emtools.strike as strike_mod

        fake_df = pd.DataFrame(
            {
                "station": ["S1"],
                "ang": [float("nan")],
                "iqr": [0.1],
                "lo": [1e-3],
                "hi": [10.0],
                "n": [5],
            }
        )
        monkeypatch.setattr(
            strike_mod, "estimate_strike_consensus", lambda *a, **k: fake_df
        )
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:2],
            "consensus",
            "all",
            5.0,
            10,
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            "dark",
            None,
            None,
        )
        assert len(out) == 2
        assert "non-finite" in str(out[0])

    def test_exception_path_returns_error_and_no_update(
        self, monkeypatch, web_app, cached_session, store_data_willy, willy_ids
    ):
        import pycsamt.emtools.strike as strike_mod

        def _boom(*a, **k):
            raise RuntimeError("kaboom")

        monkeypatch.setattr(strike_mod, "estimate_strike_sweep", _boom)
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:2],
            "sweep",
            "all",
            5.0,
            10,
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            "dark",
            None,
            None,
        )
        assert len(out) == 3
        assert "kaboom" in str(out[0])
        assert out[1] is no_update
        assert out[2] is no_update


class TestRenderStoredToolResult:
    def test_no_payload_returns_none(self):
        assert tools_mod._render_stored_tool_result("strike", None, "dark") is None

    def test_unknown_tool_returns_none(self):
        assert (
            tools_mod._render_stored_tool_result("other", {"other": {}}, "dark") is None
        )

    def test_strike_payload_renders(self):
        payload = {
            "scope_label": "All active lines (2)",
            "n_stations": 2,
            "median_strike": 45.0,
            "median_iqr": float("nan"),
            "records": [
                {"station": "S1", "line": "L1", "ang_axial": 30.0},
                {"station": "S2", "line": "L2", "ang_axial": 60.0},
            ],
            "table": [{"Station": "S1", "Strike (°)": 30.0}],
            "columns": ["Station", "Strike (°)"],
            "page_size": 10,
            "method_label": "Consensus strike",
        }
        out = tools_mod._render_stored_tool_result(
            "strike", {"strike": payload}, "dark"
        )
        assert isinstance(out, html.Div)


class TestStrikeBandFromChoice:
    def test_high(self):
        assert tools_mod._strike_band_from_choice("high") == (0.0, 0.01)

    def test_low(self):
        assert tools_mod._strike_band_from_choice("low") == (1.0, 1.0e12)

    def test_all_and_default(self):
        assert tools_mod._strike_band_from_choice("all") is None
        assert tools_mod._strike_band_from_choice(None) is None
        assert tools_mod._strike_band_from_choice("unknown") is None


class TestStrikeScopeLabel:
    def test_selected_stations_wins(self):
        out = tools_mod._strike_scope_label(None, None, ["L1"], ["A1", "A2"])
        assert out == "2 selected station(s)"

    def test_selected_lines_short_list(self):
        out = tools_mod._strike_scope_label(None, None, ["L1", "L2"], None)
        assert out == "L1, L2"

    def test_selected_lines_truncated_with_count(self):
        out = tools_mod._strike_scope_label(None, None, ["L1", "L2", "L3", "L4"], None)
        assert out == "L1, L2, L3 +1"

    def test_active_lines_fallback(self, store_data_willy):
        out = tools_mod._strike_scope_label(store_data_willy, ACTIVE_ALL, None, None)
        assert out == "All active lines (2)"

    def test_no_scope_at_all(self):
        assert tools_mod._strike_scope_label(None, None, None, None) == (
            "All loaded stations"
        )


class TestStrikeFigures:
    DF = pd.DataFrame(
        {
            "line": ["L1", "L1", "L2"],
            "ang_axial": [10.0, 20.0, 170.0],
        }
    )

    def test_rose_figure_dark(self):
        fig = tools_mod._strike_rose_figure(self.DF, "Title", "dark")
        assert fig.layout.title.text == "Title"
        assert len(fig.data) == 2  # one trace per line

    def test_rose_figure_light(self):
        fig = tools_mod._strike_rose_figure(self.DF, "Title", "light")
        assert fig.layout.paper_bgcolor == "#ffffff"

    def test_box_figure_dark(self):
        fig = tools_mod._strike_box_figure(self.DF, "dark")
        assert len(fig.data) == 2
        assert fig.layout.yaxis.title.text == "Axial strike (0-180°)"

    def test_box_figure_light(self):
        fig = tools_mod._strike_box_figure(self.DF, "light")
        assert fig.layout.paper_bgcolor == "#ffffff"


# ═══════════════════════════════════════════════════════════════════════════
# 3. Validator panel
# ═══════════════════════════════════════════════════════════════════════════


class TestValidatorBody:
    def test_no_data_shows_hint(self):
        out = tools_mod._validator_body(0, tools_mod._no_data_msg(), None, None)
        assert "Load survey data first." in out.children

    def test_with_data_builds_controls(self, store_data_willy):
        out = tools_mod._validator_body(
            store_data_willy["n_stations"],
            tools_mod._no_data_msg(),
            store_data_willy,
            ACTIVE_ALL,
        )
        assert isinstance(out, html.Div)


class TestSyncValidatorStationScope:
    def test_returns_options(self, web_app, store_data_willy):
        fn = _cb_by_input(web_app, "tool-valid-stations.options", "tool-valid-lines")
        opts, val = fn(["L2"], store_data_willy, ACTIVE_ALL)
        assert val is None
        assert len(opts) == len(store_data_willy["station_records"]) // 2


class TestRunValidator:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-valid-out.children", "tool-valid-run")

    def test_no_clicks_returns_pair_no_update(self, web_app):
        out = self._fn(web_app)(
            None, None, None, "all", 3, ["coords"], None, None, "sid", None
        )
        assert out == (no_update, no_update)

    def test_session_expired(self, web_app, store_data_willy):
        out = self._fn(web_app)(
            1,
            None,
            None,
            "all",
            3,
            ["coords"],
            store_data_willy,
            ACTIVE_ALL,
            "no-such-session",
            None,
        )
        assert len(out) == 2
        assert "Session expired" in str(out[0])

    def test_real_run_all_checks_default_severity(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:6],
            "all",
            3,
            ["coords", "z", "errors", "freq", "nan", "duplicates", "tipper"],
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            None,
        )
        children, saved = out
        assert isinstance(children, html.Div)
        summary = saved["validator"]["summary"]
        assert summary["total"] == 6
        assert summary["pass"] + summary["warn"] + summary["fail"] == 6

    def test_severity_filters_issues_only(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:6],
            "issues",
            3,
            ["coords", "z"],
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            None,
        )
        _, saved = out
        rows = saved["validator"]["rows"]
        assert all(r["Status"] in ("WARN", "FAIL") for r in rows)

    def test_severity_pass_only(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:6],
            "pass",
            1,
            ["z"],
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            None,
        )
        _, saved = out
        rows = saved["validator"]["rows"]
        assert all(r["Status"] == "PASS" for r in rows)

    def test_severity_fail_only_with_high_minfreq(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        # An absurdly high min-freq threshold cannot itself cause FAIL (only
        # WARN), so with only the "freq" check enabled this should yield an
        # empty (but valid) fail-only table -- exercises the "fail" branch.
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:4],
            "fail",
            999999,
            ["freq"],
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            None,
        )
        _, saved = out
        assert saved["validator"]["rows"] == []

    def test_exception_path(
        self, monkeypatch, web_app, cached_session, store_data_willy, willy_ids
    ):
        monkeypatch.setattr(
            tools_mod,
            "_validate_sites_rows",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        out = self._fn(web_app)(
            1,
            None,
            willy_ids[:2],
            "all",
            3,
            ["coords"],
            store_data_willy,
            ACTIVE_ALL,
            cached_session,
            None,
        )
        assert "boom" in str(out[0])
        assert out[1] is no_update


class TestValidatorStatusBar:
    def test_all_zero_still_renders_no_segments(self):
        out = tools_mod._validator_status_bar(0, 0, 0)
        assert out.children == []

    def test_mixed_counts_render_segments_with_widths(self):
        out = tools_mod._validator_status_bar(3, 1, 0)
        assert len(out.children) == 2
        widths = {seg.style["width"] for seg in out.children}
        assert widths == {"75.00%", "25.00%"}


class TestValidateSitesRows:
    def test_real_data_all_checks_pass_shape(self, willy_sites, willy_ids):
        rows = tools_mod._validate_sites_rows(
            willy_sites,
            store_data={"station_records": _records([(willy_ids[0], "L1")])},
            checks={"coords", "z", "errors", "freq", "nan", "tipper"},
            min_freq=3,
        )
        assert len(rows) == len(willy_ids)
        assert rows[0]["Station"] == willy_ids[0]
        assert rows[0]["Line"] == "L1"

    def test_empty_object_flags_missing_everything(self):
        class _Empty:
            pass

        rows = tools_mod._validate_sites_rows(
            [_Empty()],
            store_data={},
            checks={"coords", "z", "errors", "freq", "nan", "tipper"},
            min_freq=3,
        )
        row = rows[0]
        assert row["Status"] == "FAIL"
        assert row["Z ok"] == "NO"
        assert row["Errors"] == "NO"
        assert row["Tipper"] == "NO"
        assert row["N freq"] == 0
        assert "Missing impedance tensor Z." in row["Issues"]
        assert "Missing positive frequency values." in row["Issues"]

    def test_duplicates_check_flags_repeated_station_names(self, willy_sites):
        from pycsamt.emtools._core import _iter_items

        one = next(iter(_iter_items(willy_sites)))
        rows = tools_mod._validate_sites_rows(
            [one, one],
            store_data={},
            checks={"duplicates"},
            min_freq=1,
        )
        assert len(rows) == 2
        assert all("Duplicate station identifier." in r["Issues"] for r in rows)
        assert all(r["Status"] != "FAIL" for r in rows)

    def test_checks_disabled_still_reports_freq_and_coords_when_present(
        self, willy_sites, willy_ids
    ):
        rows = tools_mod._validate_sites_rows(
            willy_sites, store_data={}, checks=set(), min_freq=3
        )
        # "freq"/"coords" not in checks -> best-effort population without
        # contributing issues/warnings.
        assert rows[0]["Status"] == "PASS"


# ═══════════════════════════════════════════════════════════════════════════
# 4. Converter panel
# ═══════════════════════════════════════════════════════════════════════════


class TestConverterBody:
    def test_builds_controls_and_default_option_visibility(self):
        out = tools_mod._converter_body()
        assert isinstance(out, html.Div)
        dropdown = out.children[2]
        assert dropdown.id == "tool-conv-type"
        assert dropdown.value == "AVG -> EDI"


class TestSyncConverterOptions:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            "tool-conv-avg-options.style",
            "tool-conv-spectra-options.style",
        )

    def test_avg_shows_avg_hides_spectra(self, web_app):
        avg_style, spec_style = self._fn(web_app)("AVG -> EDI")
        assert avg_style == {"display": "block"}
        assert spec_style == {"display": "none"}

    def test_spectra_shows_spectra_hides_avg(self, web_app):
        avg_style, spec_style = self._fn(web_app)("Spectra -> EDI")
        assert avg_style == {"display": "none"}
        assert spec_style == {"display": "block"}

    def test_j_hides_both(self, web_app):
        avg_style, spec_style = self._fn(web_app)("J -> EDI")
        assert avg_style == {"display": "none"}
        assert spec_style == {"display": "none"}

    def test_none_hides_both(self, web_app):
        avg_style, spec_style = self._fn(web_app)(None)
        assert avg_style == {"display": "none"}
        assert spec_style == {"display": "none"}


@pytest.mark.skipif(not _HAS_K1, reason="K1 AVG/STN fixtures are not available")
class TestRunConverter:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-conv-out.children", "tool-conv-run")

    def test_no_clicks_returns_pair_no_update(self, web_app):
        out = self._fn(web_app)(
            None,
            "AVG -> EDI",
            "",
            "",
            "",
            "",
            "",
            True,
            "EX,EY",
            "HX,HY",
            "",
            [],
            "descending",
            None,
            ["compute_z"],
            "",
            None,
            "sid",
        )
        assert out == (no_update, no_update)

    def test_empty_source_path_warns(self, web_app):
        out = self._fn(web_app)(
            1,
            "AVG -> EDI",
            "  ",
            "",
            "",
            "",
            "",
            True,
            "EX,EY",
            "HX,HY",
            "",
            [],
            "descending",
            None,
            ["compute_z"],
            "",
            None,
            "sid",
        )
        assert "Enter a source AVG/J/Spectra path first." in str(out[0])
        assert out[1] is no_update

    def test_nonexistent_source_warns(self, web_app):
        out = self._fn(web_app)(
            1,
            "AVG -> EDI",
            "/no/such/path.AVG",
            "",
            "",
            "",
            "",
            True,
            "EX,EY",
            "HX,HY",
            "",
            [],
            "descending",
            None,
            ["compute_z"],
            "",
            None,
            "sid",
        )
        assert "Source path does not exist" in str(out[0])
        assert out[1] is no_update

    def test_real_avg_conversion_no_session_load(self, web_app, tmp_path):
        out_dir = tmp_path / "converted"
        out = self._fn(web_app)(
            1,
            "AVG -> EDI",
            str(_K1_AVG),
            str(out_dir),
            str(_K1_STN),
            "49N",
            "",
            True,
            "EX,EY",
            "HX,HY",
            "",
            [],
            "descending",
            None,
            ["compute_z", "compute_rho_phi"],
            "",
            None,
            "some-session-id",
        )
        children, new_store = out
        assert isinstance(children, html.Div)
        assert new_store is no_update
        # NB: str(<dash component>) reprs its string children (so a
        # Windows path's backslashes get doubled) -- compare against the
        # folder name only rather than the full backslash-laden path.
        assert "Written to" in str(children)
        assert out_dir.name in str(children)

    def test_real_avg_conversion_with_load_session(self, web_app, tmp_path):
        session_id = "test-converter-load-session"
        out = self._fn(web_app)(
            1,
            "AVG -> EDI",
            str(_K1_AVG),
            "",
            str(_K1_STN),
            "49N",
            "",
            True,
            "EX,EY",
            "HX,HY",
            "",
            [],
            "descending",
            None,
            ["compute_z", "compute_rho_phi", "load_session"],
            "",
            None,
            session_id,
        )
        children, new_store = out
        assert new_store is not no_update
        assert new_store["converted"] is True
        assert new_store["n_stations"] >= 1
        assert cache_get(session_id) is not None

    def test_exception_path_bad_avg_file(self, web_app, tmp_path):
        bad = tmp_path / "bad.AVG"
        bad.write_text("this is not a valid AVG file\n", encoding="utf-8")
        out = self._fn(web_app)(
            1,
            "AVG -> EDI",
            str(bad),
            "",
            "",
            "",
            "",
            True,
            "EX,EY",
            "HX,HY",
            "",
            [],
            "descending",
            None,
            ["compute_z", "compute_rho_phi"],
            "",
            None,
            "sid",
        )
        children, new_store = out
        assert "✗" in str(children)
        assert new_store is no_update


class TestConverterOptionsFromValues:
    def test_defaults(self):
        opts = tools_mod._converter_options_from_values(
            output_dir=None,
            stn_path=None,
            utm_zone=None,
            epsg=None,
            convert_coords=False,
            e_labels=None,
            h_labels=None,
            suffix=None,
            spectra_flags=None,
            freq_order=None,
            freq_tol=None,
            core_flags=None,
            station_name=None,
        )
        assert opts["e_labels"] == "EX,EY"
        assert opts["h_labels"] == "HX,HY"
        assert opts["freq_order"] == "descending"
        assert opts["compute_z"] is False
        assert "freq_tol" not in opts

    def test_flags_and_freq_tol(self):
        opts = tools_mod._converter_options_from_values(
            output_dir="/out",
            stn_path="/stn",
            utm_zone="49N",
            epsg="32649",
            convert_coords=True,
            e_labels="A,B",
            h_labels="C,D",
            suffix="_X",
            spectra_flags=["estimate_error", "use_remote"],
            freq_order="ascending",
            freq_tol="1e-6",
            core_flags=["compute_z", "compute_rho_phi"],
            station_name="K1",
        )
        assert opts["estimate_error"] is True
        assert opts["use_remote"] is True
        assert opts["skip_errors"] is False
        assert opts["compute_z"] is True
        assert opts["freq_tol"] == 1e-6
        assert opts["station_name"] == "K1"


class TestConvertedSitesToStore:
    def test_builds_store_without_old_store(self, willy_sites):
        store = tools_mod._converted_sites_to_store(
            willy_sites, old_store=None, data_dir="/tmp/out"
        )
        assert store["converted"] is True
        assert store["data_dir"] == "/tmp/out"
        assert store["n_stations"] == len(store["station_records"])

    def test_falls_back_when_no_data_dir(self, willy_sites):
        store = tools_mod._converted_sites_to_store(
            willy_sites, old_store=None, data_dir=""
        )
        assert store["data_dir"] == "[converted]"

    def test_reuses_line_assignments_from_old_store(self, willy_sites, willy_ids):
        old_store = {"station_records": _records([(willy_ids[0], "West Line")])}
        store = tools_mod._converted_sites_to_store(
            willy_sites, old_store=old_store, data_dir="/tmp"
        )
        lines = {r["ID"]: r.get("Line") for r in store["station_records"]}
        assert lines.get(willy_ids[0]) == "West Line"


class TestConverterResultView:
    def test_no_rows_shows_warning(self):
        out = tools_mod._converter_result_view(
            stats={}, failures=[], output_dir="", loaded=False
        )
        assert "No EDI objects were produced." in str(out)

    def test_rows_and_failures_render_tables(self):
        stats = {
            "rows": [
                {
                    "station": "S1",
                    "n_freqs": 10,
                    "f_min": 0.1,
                    "f_max": 100.0,
                    "lat": 26.0,
                    "lon": 113.0,
                    "elev": 500.0,
                    "has_Z": True,
                    "has_tipper": False,
                }
            ],
            "n_total": 1,
            "n_failures": 1,
        }
        failures = [{"source": "bad.avg", "error": "parse error"}]
        out = tools_mod._converter_result_view(
            stats=stats, failures=failures, output_dir="/out", loaded=True
        )
        text = str(out)
        assert "Written to /out" in text
        assert "Loaded into session." in text
        assert "bad.avg" in text

    def test_no_output_dir_message(self):
        out = tools_mod._converter_result_view(
            stats={"rows": [], "n_total": 0, "n_failures": 0},
            failures=[],
            output_dir="",
            loaded=False,
        )
        assert "No output folder provided" in str(out)


class TestFmtFloatAndFiniteNumber:
    def test_fmt_float_finite(self):
        assert tools_mod._fmt_float(3.14159, nd=4) == "3.142"

    def test_fmt_float_nan(self):
        assert tools_mod._fmt_float(float("nan")) == "—"

    def test_fmt_float_non_numeric(self):
        assert tools_mod._fmt_float("abc") == "—"

    def test_fmt_float_high_precision_uses_f_format(self):
        assert tools_mod._fmt_float(1.23456789, nd=5) == "1.23457"

    def test_is_finite_number_true(self):
        assert tools_mod._is_finite_number(1.5) is True

    def test_is_finite_number_false_for_nan_and_bad(self):
        assert tools_mod._is_finite_number(float("nan")) is False
        assert tools_mod._is_finite_number("abc") is False


class TestZBlockWithErrors:
    def test_real_edi_returns_four_tuple(self, willy_sites):
        from pycsamt.emtools._core import _iter_items

        ed = next(iter(_iter_items(willy_sites)))
        out = tools_mod._z_block_with_errors(ed)
        assert len(out) == 4


class TestExtractLatLon:
    def test_direct_lat_lon_attrs(self):
        class Obj:
            lat = 26.05
            lon = 113.48

        assert tools_mod._extract_lat_lon(Obj()) == (26.05, 113.48)

    def test_latitude_longitude_attrs(self):
        class Obj:
            latitude = 1.0
            longitude = 2.0

        assert tools_mod._extract_lat_lon(Obj()) == (1.0, 2.0)

    def test_coords_attribute_fallback(self):
        class Obj:
            coords = (5.0, 6.0)

        assert tools_mod._extract_lat_lon(Obj()) == (5.0, 6.0)

    def test_wrapped_used_when_raw_has_nothing(self):
        class Raw:
            pass

        class Wrapped:
            lat = 9.0
            lon = 10.0

        assert tools_mod._extract_lat_lon(Raw(), Wrapped()) == (9.0, 10.0)

    def test_head_attribute_fallback(self):
        class Head:
            lat = 11.0
            lon = 12.0

        class Raw:
            Header = Head()

        assert tools_mod._extract_lat_lon(Raw()) == (11.0, 12.0)

    def test_dict_head_fallback(self):
        class Raw:
            Head = {"lat": 13.0, "lon": 14.0}

        assert tools_mod._extract_lat_lon(Raw()) == (13.0, 14.0)

    def test_nothing_found_returns_none_none(self):
        class Raw:
            pass

        assert tools_mod._extract_lat_lon(Raw()) == (None, None)


class TestErrWarn:
    def test_err_span(self):
        out = tools_mod._err(ValueError("bad"))
        assert "✗ bad" in out.children

    def test_warn_span(self):
        out = tools_mod._warn("careful")
        assert "⚠ careful" in out.children


class TestDd2Dms:
    def test_lat_positive_is_north(self):
        assert tools_mod._dd2dms(40.7128, "lat").endswith("N")

    def test_lat_negative_is_south(self):
        assert tools_mod._dd2dms(-33.87, "lat").endswith("S")

    def test_lon_positive_is_east(self):
        assert tools_mod._dd2dms(151.2, "lon").endswith("E")

    def test_lon_negative_is_west(self):
        assert tools_mod._dd2dms(-74.006, "lon").endswith("W")

    def test_known_value_format(self):
        # 40.7128 deg = 40 deg 42' 46.08"
        s = tools_mod._dd2dms(40.7128, "lat")
        assert s.startswith("40°42′")


# ═══════════════════════════════════════════════════════════════════════════
# 5. Coordinate transformer panel
# ═══════════════════════════════════════════════════════════════════════════


class TestCoordsBody:
    def test_builds_all_mode_panels(self):
        out = tools_mod._coords_body()
        assert isinstance(out, html.Div)
        text = str(out)
        for cid in (
            "tool-coord-mode-single",
            "tool-coord-mode-batch",
            "tool-coord-mode-survey",
            "tool-coord-direction",
            "tool-coord-easting",
            "tool-coord-lat",
            "tool-coord-epsg-src",
            "tool-coord-upload",
            "tool-coord-zone-survey",
            "tool-coord-run",
            "tool-coord-export-btn",
        ):
            assert cid in text

    def test_export_button_disabled_initially(self):
        out = tools_mod._coords_body()
        action_bar = out.children[8]
        export_btn = action_bar.children[2]
        assert export_btn.disabled is True


class TestCoordParseCsv:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "tool-coord-upload-store.data", "tool-coord-upload"
        )

    def test_no_contents_returns_all_no_update(self, web_app):
        out = self._fn(web_app)(None, "f.csv")
        assert out == (no_update,) * 5

    def test_parses_valid_csv(self, web_app):
        raw = "Easting,Northing,ID\n500000,4500000,ST1\n510000,4510000,ST2\n"
        contents = "data:text/csv;base64," + base64.b64encode(raw.encode()).decode()
        data, info, e_opts, n_opts, id_opts = self._fn(web_app)(contents, "points.csv")
        assert len(data) == 2
        assert "points.csv" in str(info)
        assert {o["value"] for o in e_opts} == {"Easting", "Northing", "ID"}
        assert id_opts[0] == {"label": "— none —", "value": ""}

    def test_malformed_contents_returns_error(self, web_app):
        data, info, e_opts, n_opts, id_opts = self._fn(web_app)(
            "no-comma-here", "f.csv"
        )
        assert data is None
        assert "✗" in str(info)
        assert e_opts == [] and n_opts == [] and id_opts == []


class TestRunCoords:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-coord-out.children", "tool-coord-run")

    def _call(self, web_app, **over):
        base = dict(
            n=1,
            mode="single",
            direction="utm2ll",
            datum="WGS84",
            easting=None,
            northing=None,
            zone=None,
            lat=None,
            lon=None,
            target_zone=None,
            epsg_src=None,
            epsg_dst=None,
            csv_data=None,
            col_e=None,
            col_n=None,
            col_id=None,
            zone_batch=None,
            zone_survey=None,
            survey_store=None,
        )
        base.update(over)
        fn = self._fn(web_app)
        return fn(*base.values())

    def test_no_clicks_returns_six_no_update(self, web_app):
        out = self._call(web_app, n=None)
        assert out == (no_update,) * 6

    def test_single_utm2ll_missing_fields_warns(self, web_app):
        out = self._call(web_app, direction="utm2ll")
        assert "Fill in Easting, Northing and UTM Zone." in str(out[0])
        assert out[5] is True

    def test_single_utm2ll_success(self, web_app):
        out = self._call(
            web_app,
            direction="utm2ll",
            easting=500000.0,
            northing=2881263.0,
            zone="49N",
        )
        text, result, pct, label, style, disabled = out
        assert pct == 100
        assert label == "Done"
        assert disabled is False
        assert result[0]["Zone"] == "49N"
        assert -90 <= result[0]["Latitude"] <= 90

    def test_single_ll2utm_missing_fields_warns(self, web_app):
        out = self._call(web_app, direction="ll2utm")
        assert "Fill in Latitude and Longitude." in str(out[0])

    def test_single_ll2utm_success(self, web_app):
        out = self._call(web_app, direction="ll2utm", lat=26.05, lon=113.48)
        text, result, pct, label, style, disabled = out
        assert pct == 100
        assert result[0]["Latitude"] == 26.05
        assert result[0]["Easting(m)"] > 0

    def test_single_epsg2epsg_missing_fields_warns(self, web_app):
        out = self._call(web_app, direction="epsg2epsg")
        assert "Fill in X, Y and both EPSG codes." in str(out[0])

    def test_single_epsg2epsg_success(self, web_app):
        out = self._call(
            web_app,
            direction="epsg2epsg",
            easting=500000.0,
            northing=2881263.0,
            epsg_src=32649,
            epsg_dst=4326,
        )
        text, result, pct, label, style, disabled = out
        assert pct == 100
        assert (
            "EPSG:32649" in result[0]["EPSG_src"]
            if isinstance(result[0]["EPSG_src"], str)
            else result[0]["EPSG_src"] == 32649
        )

    def test_batch_no_csv_warns(self, web_app):
        out = self._call(web_app, mode="batch")
        assert "Upload a CSV file first." in str(out[0])

    def test_batch_missing_columns_selection_warns(self, web_app):
        out = self._call(web_app, mode="batch", csv_data=[{"E": 1, "N": 2}])
        assert "Select the Easting and Northing columns." in str(out[0])

    def test_batch_columns_not_found_warns(self, web_app):
        out = self._call(
            web_app,
            mode="batch",
            csv_data=[{"E": 1, "N": 2}],
            col_e="Easting",
            col_n="Northing",
        )
        assert "not found in data" in str(out[0])

    def test_batch_utm2ll_success_with_zone_and_id(self, web_app):
        out = self._call(
            web_app,
            mode="batch",
            direction="utm2ll",
            csv_data=[
                {"Easting": 500000, "Northing": 2881263, "ID": "P1"},
                {"Easting": 501000, "Northing": 2882263, "ID": "P2"},
            ],
            col_e="Easting",
            col_n="Northing",
            col_id="ID",
            zone_batch="49N",
        )
        text, result, pct, label, style, disabled = out
        assert len(result) == 2
        assert result[0]["ID"] == "P1"
        assert disabled is False

    def test_batch_ll2utm_success(self, web_app):
        out = self._call(
            web_app,
            mode="batch",
            direction="ll2utm",
            csv_data=[{"Lat": 26.05, "Lon": 113.48}],
            col_e="Lat",
            col_n="Lon",
        )
        _, result, _, _, _, _ = out
        assert result[0]["Easting(m)"] > 0

    def test_batch_epsg2epsg_success(self, web_app):
        out = self._call(
            web_app,
            mode="batch",
            direction="epsg2epsg",
            csv_data=[{"X": 500000, "Y": 2881263}],
            col_e="X",
            col_n="Y",
            epsg_src=32649,
            epsg_dst=4326,
        )
        _, result, _, _, _, _ = out
        assert "X_out" in result[0]

    def test_batch_row_conversion_failure_recorded(self, web_app):
        out = self._call(
            web_app,
            mode="batch",
            direction="utm2ll",
            csv_data=[{"Easting": "not-a-number", "Northing": 2881263}],
            col_e="Easting",
            col_n="Northing",
            zone_batch="49N",
        )
        _, result, _, _, _, _ = out
        assert "Error" in result[0]

    def test_survey_mode_with_loaded_data_still_warns_no_survey_data(
        self, web_app, store_data_willy
    ):
        """
        Documents a real bug: "Survey" mode reads
        ``survey_store.get("stations")``, but the shared app-wide
        STORE_DATA store (populated by callbacks/data.py, see
        ``_active_line_names``/``_station_line_map`` elsewhere in this
        module) never has a "stations" key -- only "station_records". So
        Survey mode always falls into the "No survey data loaded." branch
        even when a survey is actively loaded, making it permanently
        unreachable in the running app.
        """
        out = self._call(web_app, mode="survey", survey_store=store_data_willy)
        assert "No survey data loaded." in str(out[0])
        assert out[5] is True

    def test_survey_mode_no_store_warns(self, web_app):
        out = self._call(web_app, mode="survey", survey_store=None)
        assert "No survey data loaded." in str(out[0])

    def test_survey_mode_success_with_synthetic_stations_key(self, web_app):
        # Exercises the (currently unreachable in the live app -- see bug
        # documented above) success path directly by supplying the
        # "stations" key the code actually looks for.
        out = self._call(
            web_app,
            mode="survey",
            survey_store={
                "stations": [
                    {"name": "ST1", "easting": 500000, "northing": 2881263},
                    {"name": "ST2", "x": 501000, "y": 2882263, "zone": "49N"},
                    {"name": "ST3"},
                ]
            },
        )
        _, result, pct, label, _, disabled = out
        assert len(result) == 3
        assert result[2]["Error"] == "no coords"
        assert disabled is False

    def test_epsg_conversion_without_pyproj_available(self, monkeypatch, web_app):
        import builtins

        real_import = builtins.__import__

        def _fake_import(name, *a, **k):
            if name == "pyproj":
                raise ImportError("no pyproj")
            return real_import(name, *a, **k)

        monkeypatch.setattr(builtins, "__import__", _fake_import)
        out = self._call(
            web_app,
            direction="epsg2epsg",
            easting=1.0,
            northing=2.0,
            epsg_src=32649,
            epsg_dst=4326,
        )
        assert "pyproj is required" in str(out[0])

    def test_exception_path_returns_error(self, web_app):
        # zone=None with direction utm2ll normally warns before reaching
        # project_point_utm2ll; force the generic except by passing a zone
        # string that project_point_utm2ll cannot parse at all.
        out = self._call(
            web_app,
            direction="utm2ll",
            easting=1.0,
            northing=2.0,
            zone="!!!not-a-zone!!!",
        )
        assert "✗" in str(out[0]) or "Fill in" in str(out[0])


class TestCoordMap:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-coord-out", "tool-coord-map-btn")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, [{"Latitude": 1, "Longitude": 2}]) is no_update

    def test_no_results_no_update(self, web_app):
        assert self._fn(web_app)(1, None) is no_update

    def test_no_latlon_in_results_shows_hint(self, web_app):
        out = self._fn(web_app)(1, [{"Easting(m)": 1}])
        assert "run a UTM→LL conversion first" in str(out)

    def test_builds_map_figure(self, web_app):
        out = self._fn(web_app)(
            1, [{"Latitude": 26.0, "Longitude": 113.0, "Station": "S1"}]
        )
        assert out.__class__.__name__ == "Graph"


class TestCoordExport:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "tool-coord-download.data", "tool-coord-export-btn"
        )

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, [{"a": 1}]) is no_update

    def test_no_results_no_update(self, web_app):
        assert self._fn(web_app)(1, None) is no_update

    def test_builds_csv_download(self, web_app):
        out = self._fn(web_app)(1, [{"Latitude": 1.0, "Longitude": 2.0}])
        assert out["filename"] == "coordinates_converted.csv"
        assert out["type"] == "text/csv"
        assert "Latitude" in out["content"]


# ═══════════════════════════════════════════════════════════════════════════
# 6. Server-side path browser (used by the converter panel's File/Folder
#    buttons) -- exercised via monkeypatched ``tools_mod.ctx``/``set_props``
#    since it genuinely needs Dash's live callback context.
# ═══════════════════════════════════════════════════════════════════════════


class _FakeCtx:
    def __init__(self, triggered_id):
        self.triggered_id = triggered_id


class TestToolsPathBrowse:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tools-path-modal.is_open", "tools-path-browse")

    def test_no_clicks_raises_prevent_update(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)([None, 0])

    def test_non_dict_trigger_raises_prevent_update(self, monkeypatch, web_app):
        monkeypatch.setattr(tools_mod, "ctx", _FakeCtx("some-string-id"))
        with _raises_prevent_update():
            self._fn(web_app)([1])

    def test_opens_with_target_metadata(self, monkeypatch, web_app):
        monkeypatch.setattr(
            tools_mod,
            "ctx",
            _FakeCtx(
                {
                    "type": "tools-path-browse",
                    "target": "tool-conv-source",
                    "mode": "file",
                    "title": "Pick a file",
                }
            ),
        )
        is_open, current, target, title = self._fn(web_app)([1])
        assert is_open is True
        assert target["target"] == "tool-conv-source"
        assert title == "Pick a file"


class TestRefreshToolsPathListing:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tools-path-list.children")

    def test_lists_folder_contents(self, web_app, tmp_path):
        (tmp_path / "sub").mkdir()
        (tmp_path / "notes.txt").write_text("x", encoding="utf-8")
        (tmp_path / "ignored.bin").write_bytes(b"\x00")
        rows, cwd = self._fn(web_app)(str(tmp_path), {"mode": "file"})
        assert cwd == str(tmp_path)
        assert "sub" in str(rows) and "notes.txt" in str(rows)

    def test_folder_mode_excludes_files(self, web_app, tmp_path):
        (tmp_path / "notes.txt").write_text("x", encoding="utf-8")
        rows, _ = self._fn(web_app)(str(tmp_path), {"mode": "folder"})
        assert "No sub-folders here." in str(rows) or "notes.txt" not in str(rows)

    def test_empty_folder_shows_hint(self, web_app, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        rows, _ = self._fn(web_app)(str(empty), {"mode": "file"})
        assert "No folders or supported files here." in str(rows)

    def test_unreadable_path_reports_error(self, web_app):
        rows, _ = self._fn(web_app)(
            str(Path("Z:/definitely/not/a/real/path")), {"mode": "file"}
        )
        assert "Cannot read this folder" in str(rows)


class TestToolsPathEnterDir:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tools-path-current.data", "tools-path-dir")

    def test_no_clicks_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)([None])

    def test_enters_clicked_dir(self, monkeypatch, web_app):
        monkeypatch.setattr(
            tools_mod,
            "ctx",
            _FakeCtx({"type": "tools-path-dir", "path": "/x/y"}),
        )
        assert self._fn(web_app)([1]) == "/x/y"

    def test_non_dict_trigger_raises(self, monkeypatch, web_app):
        monkeypatch.setattr(tools_mod, "ctx", _FakeCtx(None))
        with _raises_prevent_update():
            self._fn(web_app)([1])


class TestToolsPathParent:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tools-path-current.data", "tools-path-up")

    def test_no_clicks_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)(None, "/a/b")

    def test_goes_to_parent(self, web_app, tmp_path):
        child = tmp_path / "child"
        child.mkdir()
        parent = self._fn(web_app)(1, str(child))
        assert Path(parent) == tmp_path

    def test_root_raises_prevent_update(self, web_app):
        root = Path(_fs_root())
        with _raises_prevent_update():
            self._fn(web_app)(1, str(root))


def _fs_root() -> str:
    return Path(__file__).drive + "\\" if Path(__file__).drive else "/"


class TestToolsPathMkdir:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tools-path-current.data", "tools-path-mkdir")

    def test_no_clicks_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)(None, "newdir", "/tmp")

    def test_no_name_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)(1, "", "/tmp")

    def test_creates_directory(self, web_app, tmp_path):
        new_path, cleared = self._fn(web_app)(1, "sub folder/../weird", str(tmp_path))
        assert cleared == ""
        assert Path(new_path).exists()

    def test_blank_after_sanitizing_raises(self, web_app, tmp_path):
        with _raises_prevent_update():
            self._fn(web_app)(1, "   ", str(tmp_path))


class TestToolsPathCancel:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tools-path-modal.is_open", "tools-path-cancel")

    def test_no_clicks_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)(None)

    def test_click_closes(self, web_app):
        assert self._fn(web_app)(1) is False


class TestToolsPathSelectFolder:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "tools-path-modal.is_open", "tools-path-select-folder"
        )

    def test_no_clicks_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)(None, "/x", {"target": "tool-conv-source"})

    def test_no_target_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)(1, "/x", None)

    def test_selects_and_calls_set_props(self, monkeypatch, web_app, tmp_path):
        calls = []
        monkeypatch.setattr(
            tools_mod,
            "set_props",
            lambda cid, props: calls.append((cid, props)),
        )
        is_open = self._fn(web_app)(1, str(tmp_path), {"target": "tool-conv-source"})
        assert is_open is False
        assert calls[0][0] == "tool-conv-source"
        assert calls[0][1]["value"] == str(tmp_path)


class TestToolsPathSelectFile:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tools-path-modal.is_open", "tools-path-file")

    def test_no_clicks_raises(self, web_app):
        with _raises_prevent_update():
            self._fn(web_app)([None], {"target": "tool-conv-source"})

    def test_non_dict_trigger_raises(self, monkeypatch, web_app):
        monkeypatch.setattr(tools_mod, "ctx", _FakeCtx("str-id"))
        with _raises_prevent_update():
            self._fn(web_app)([1], {"target": "tool-conv-source"})

    def test_selects_file_and_calls_set_props(self, monkeypatch, web_app):
        calls = []
        monkeypatch.setattr(
            tools_mod,
            "set_props",
            lambda cid, props: calls.append((cid, props)),
        )
        monkeypatch.setattr(
            tools_mod,
            "ctx",
            _FakeCtx({"type": "tools-path-file", "path": "/a/b/K1.AVG"}),
        )
        is_open = self._fn(web_app)([1], {"target": "tool-conv-source"})
        assert is_open is False
        assert calls[0][0] == "tool-conv-source"
        assert "K1.AVG" in calls[0][1]["value"]
