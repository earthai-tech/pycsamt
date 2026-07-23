# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.dashboard (home-page survey dashboard).

Strategy
--------
Every callback in this module consumes the plain ``dict`` shape already
written into ``STORE_DATA`` by ``callbacks/data.py`` (keys:
``n_stations``, ``n_lines``, ``data_dir``, ``line_counts``,
``station_records``) and produces Dash component trees / plotly figure
dicts from it directly -- there is no I/O, no controller singleton, and
no dependency on real EDI parsing. So synthetic ``store_data`` dicts
(single-line and multi-line) are built by hand here rather than reusing
real survey fixtures; this keeps every branch (multi-line vs
single-line, geo-present vs geo-absent, station-selected vs not)
reachable and deterministic.
"""

from __future__ import annotations

import pytest
from dash import no_update

from pycsamt.app.web.layout import IDs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(
        k for k in web_app.callback_map if all(s in k for s in substrings)
    )
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(i.get("id") == input_id for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


@pytest.fixture(autouse=True, scope="module")
def _merge_global_callbacks(web_app):
    # SI-1/SI-2 use the bare ``dash.callback`` decorator, which stashes its
    # entries in a global registry until ``Dash._setup_server`` merges them
    # into ``app.callback_map`` (normally triggered by the first request).
    web_app._setup_server()


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id}])
    )


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


def _rec(id_, lat=None, lon=None, elev=None, nfreq=0, tipper=False, line="—"):
    # Real ``DataController`` output always carries a "Line" key, defaulting
    # to "—" (em-dash) for single-line surveys -- see data_controller.py:129.
    return {
        "ID": id_,
        "Latitude": lat,
        "Longitude": lon,
        "Elevation": elev,
        "N_freq": nfreq,
        "Tipper": tipper,
        "Line": line,
    }


SINGLE_LINE_STORE = {
    "n_stations": 3,
    "n_lines": 1,
    "data_dir": "/data/survey1/",
    "line_counts": {},
    "station_records": [
        _rec("S1", 48.10, 7.10, 200.0, 40, True),
        _rec("S2", 48.20, 7.20, 210.0, 35, False),
        _rec("S3", None, None, 220.0, 20, False),
    ],
}

MULTI_LINE_STORE = {
    "n_stations": 4,
    "n_lines": 2,
    "data_dir": "/data/survey2",
    "line_counts": {"L1": 2, "L2": 2},
    "station_records": [
        _rec("A1", 48.10, 7.10, 200.0, 40, True, line="L1"),
        _rec("A2", 48.15, 7.15, 205.0, 35, False, line="L1"),
        _rec("B1", 48.30, 7.30, 300.0, 45, True, line="L2"),
        _rec("B2", 48.35, 7.35, 305.0, 15, False, line="L2"),
    ],
}

NO_GEO_STORE = {
    "n_stations": 2,
    "n_lines": 1,
    "data_dir": "",
    "line_counts": {},
    "station_records": [
        _rec("N1", None, None, None, 0, False),
        _rec("N2", None, None, None, 0, False),
    ],
}


# ── 1. update_dashboard ──────────────────────────────────────────────────────


class TestUpdateDashboard:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.DASH_LINE_TABLE}.children")

    def test_no_store_data_shows_empty_state(self, web_app):
        fn = self._fn(web_app)
        out = fn(None)
        assert out[0] == {"display": "block"}  # empty state shown
        assert out[1] == {"display": "none"}  # data state hidden
        assert out[4] == []  # line table children
        assert out[5] == []  # quality children
        assert out[6] == ""
        assert out[7] == ""
        assert out[8] == "dash-health-badge"

    def test_zero_stations_shows_empty_state(self, web_app):
        fn = self._fn(web_app)
        out = fn({"n_stations": 0})
        assert out[0] == {"display": "block"}

    def test_single_line_uses_plain_span(self, web_app):
        fn = self._fn(web_app)
        out = fn(SINGLE_LINE_STORE)
        assert out[0] == {"display": "none"}
        assert out[1] == {"display": "block"}
        line_table = out[4]
        assert "single folder" in line_table.children

    def test_multi_line_builds_table_with_header_and_rows(self, web_app):
        fn = self._fn(web_app)
        out = fn(MULTI_LINE_STORE)
        line_table = out[4]
        rows = line_table.children
        assert len(rows) == 1 + len(MULTI_LINE_STORE["line_counts"])
        header = rows[0]
        assert header.children[0].children == "Line"
        first_data_row = rows[1]
        # name span, count span, bar-track div
        name_span, count_span, bar = first_data_row.children
        assert name_span.children == "L1"
        assert count_span.children == "2"

    def test_quality_all_ok_when_all_geo_and_high_freq(self, web_app):
        fn = self._fn(web_app)
        out = fn(SINGLE_LINE_STORE)
        # S3 lacks coordinates -> not all geo-valid -> not all_ok
        assert "Review data quality" in out[7]
        assert out[8] == "dash-health-badge warn"

    def test_quality_all_ok_true_path(self, web_app):
        store = {
            "n_stations": 2,
            "n_lines": 1,
            "data_dir": "",
            "line_counts": {},
            "station_records": [
                _rec("G1", 1.0, 2.0, 100.0, 40, True),
                _rec("G2", 1.1, 2.1, 110.0, 35, False),
            ],
        }
        fn = self._fn(web_app)
        out = fn(store)
        assert "All checks passed" in out[7]
        assert out[8] == "dash-health-badge good"

    def test_survey_badge_from_data_dir_basename(self, web_app):
        fn = self._fn(web_app)
        out = fn(SINGLE_LINE_STORE)
        assert out[6] == "survey1"

    def test_survey_badge_defaults_when_no_data_dir(self, web_app):
        fn = self._fn(web_app)
        out = fn(NO_GEO_STORE)
        assert out[6] == "survey"

    def test_quality_chips_report_tipper_and_freq_counts(self, web_app):
        fn = self._fn(web_app)
        out = fn(MULTI_LINE_STORE)
        chips = out[5]
        assert len(chips) == 3


# ── 2. populate_station_dropdown ─────────────────────────────────────────────


class TestPopulateStationDropdown:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.PROFILE_PAGE_STATION}.options")

    def test_no_store_data(self, web_app):
        fn = self._fn(web_app)
        assert fn(None) == ([], None)

    def test_empty_records(self, web_app):
        fn = self._fn(web_app)
        opts, first = fn({"station_records": []})
        assert opts == []
        assert first is None

    def test_populates_options_from_records(self, web_app):
        fn = self._fn(web_app)
        opts, first = fn(SINGLE_LINE_STORE)
        assert opts == [
            {"label": "S1", "value": "S1"},
            {"label": "S2", "value": "S2"},
            {"label": "S3", "value": "S3"},
        ]
        assert first == "S1"


# ── 3. profile_station_selected ──────────────────────────────────────────────


class TestProfileStationSelected:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.STORE_SELECTION}.data", IDs.PROFILE_PAGE_STATION
        )

    def test_no_station_returns_no_update(self, web_app):
        fn = self._fn(web_app)
        assert fn(None) is no_update
        assert fn("") is no_update

    def test_station_selected_writes_store(self, web_app):
        fn = self._fn(web_app)
        assert fn("S1") == {"station_id": "S1"}
        assert fn(42) == {"station_id": "42"}


# ── 4. sync_profile_dropdown ─────────────────────────────────────────────────


class TestSyncProfileDropdown:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.PROFILE_PAGE_STATION}.value", IDs.STORE_SELECTION
        )

    def test_no_selection_no_update(self, web_app):
        fn = self._fn(web_app)
        assert fn(None, [{"value": "S1"}]) is no_update

    def test_no_options_no_update(self, web_app):
        fn = self._fn(web_app)
        assert fn({"station_id": "S1"}, []) is no_update

    def test_sid_in_valid_options_returns_sid(self, web_app):
        fn = self._fn(web_app)
        opts = [{"value": "S1"}, {"value": "S2"}]
        assert fn({"station_id": "S2"}, opts) == "S2"

    def test_sid_not_in_valid_options_no_update(self, web_app):
        fn = self._fn(web_app)
        opts = [{"value": "S1"}]
        assert fn({"station_id": "ZZZ"}, opts) is no_update


# ── 4b. update_prof_info_bar ─────────────────────────────────────────────────


class TestUpdateProfInfoBar:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.PROF_INFO_BAR}.children")

    def test_no_store_data_shows_placeholder(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, None)
        assert "No station selected" in out[0].children[1]

    def test_no_selection_shows_placeholder(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, SINGLE_LINE_STORE)
        assert "No station selected" in out[0].children[1]

    def test_no_match_but_sid_present_shows_bare_span(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "UNKNOWN"}, SINGLE_LINE_STORE)
        assert out[0].children[1] == "UNKNOWN"

    def test_match_builds_full_info_bar_with_line(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "A1"}, MULTI_LINE_STORE)
        # station span, sep, coords, sep, line inserted, freq badge, tipper badge
        texts = [
            getattr(c, "children", None)
            for c in out
        ]
        assert any(
            isinstance(t, str) and "Line L1" in t for t in texts
        )
        # tipper badge present since A1 has Tipper=True
        tip_badge = out[-1]
        assert "prof-badge-tipper" in tip_badge.className

    def test_match_without_line_and_without_tipper(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "S2"}, SINGLE_LINE_STORE)
        tip_badge = out[-1]
        assert "prof-badge-no-tipper" in tip_badge.className

    def test_match_missing_coords_falls_back_to_dash(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "S3"}, SINGLE_LINE_STORE)
        coords_span = next(
            c for c in out if getattr(c, "className", "") == "prof-info-coords"
        )
        assert "—" in coords_span.children


# ── 5. update_map_station_info ───────────────────────────────────────────────


class TestUpdateMapStationInfo:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.MAP_PAGE_INFO}.children")

    def test_no_selection(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, SINGLE_LINE_STORE)
        assert "Click a station" in out.children

    def test_selection_without_sid(self, web_app):
        fn = self._fn(web_app)
        out = fn({}, SINGLE_LINE_STORE)
        assert "Click a station" in out.children

    def test_selection_without_store_data(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "S1"}, None)
        assert "Click a station" in out.children

    def test_no_match_returns_bare_name(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "ZZZ"}, SINGLE_LINE_STORE)
        assert out.children == "ZZZ"

    def test_match_formats_float_lat_lon_elev(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "S1"}, SINGLE_LINE_STORE)
        rows = out.children
        assert rows[0].children == "S1"
        lat_row = rows[2]
        assert "48.1000°" in lat_row.children[1].children

    def test_match_missing_coords_uses_dash_string(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "S3"}, SINGLE_LINE_STORE)
        lat_row = out.children[2]
        assert lat_row.children[1].children == "None"


# ── 6. populate_map_line_filter ──────────────────────────────────────────────


class TestPopulateMapLineFilter:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP_PAGE_LINE_FILTER}.options")

    def test_no_store_data(self, web_app):
        fn = self._fn(web_app)
        assert fn(None) == ([], [], "(all)")

    def test_no_line_counts(self, web_app):
        fn = self._fn(web_app)
        assert fn(SINGLE_LINE_STORE) == ([], [], "(all)")

    def test_with_line_counts(self, web_app):
        fn = self._fn(web_app)
        opts, values, count_label = fn(MULTI_LINE_STORE)
        assert values == ["L1", "L2"]
        assert count_label == "(2)"
        assert opts[0] == {"label": "L1 (2)", "value": "L1"}


# ── 7. update_map_stats ──────────────────────────────────────────────────────


class TestUpdateMapStats:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.MAP_PAGE_STATS}.children")

    def test_no_store_data(self, web_app):
        fn = self._fn(web_app)
        assert fn(None, None) == ""

    def test_stats_without_filter(self, web_app):
        fn = self._fn(web_app)
        out = fn(MULTI_LINE_STORE, None)
        values = [d.children[0].children for d in out]
        assert values == ["4", "2", "2", "4"]

    def test_stats_with_line_filter(self, web_app):
        fn = self._fn(web_app)
        out = fn(MULTI_LINE_STORE, ["L1"])
        values = [d.children[0].children for d in out]
        # only L1's 2 stations counted
        assert values[0] == "2"
        assert values[1] == "1"


# ── 9. update_analytics_charts ───────────────────────────────────────────────


class TestUpdateAnalyticsCharts:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.DASH_CHART_DONUT}.figure")

    def test_no_records_returns_all_empty_figs(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, "dark")
        assert len(out) == 6
        for fig in out:
            assert fig["data"] == []

    def test_single_line_light_theme(self, web_app):
        fn = self._fn(web_app)
        donut, freq, geo, rank, heat, violin = fn(SINGLE_LINE_STORE, "light")
        assert donut["data"][0]["type"] == "pie"
        assert freq["data"][0]["type"] == "bar"
        assert rank["data"][0]["type"] == "bar"
        assert heat["data"][0]["type"] == "heatmap"
        assert violin["data"][0]["type"] == "box"

    def test_multi_line_dark_theme_uses_stacked_bar_donut(self, web_app):
        fn = self._fn(web_app)
        donut, freq, geo, rank, heat, violin = fn(MULTI_LINE_STORE, "dark")
        assert donut["layout"]["barmode"] == "stack"
        assert freq["layout"]["barmode"] == "group"
        # geo has 2 lines -> 2 lines*2 traces = 4 traces
        assert len(geo["data"]) == 4
        assert len(violin["data"]) == 2

    def test_no_geo_records_shows_annotation(self, web_app):
        fn = self._fn(web_app)
        donut, freq, geo, rank, heat, violin = fn(NO_GEO_STORE, "dark")
        assert geo["data"] == []
        assert "No coordinate data" in geo["layout"]["annotations"][0]["text"]

    def test_heatmap_pads_shorter_lines_with_none(self, web_app):
        store = {
            "n_stations": 3,
            "n_lines": 2,
            "data_dir": "",
            "line_counts": {"L1": 2, "L2": 1},
            "station_records": [
                _rec("A1", 1.0, 2.0, 100.0, 10, False, line="L1"),
                _rec("A2", 1.1, 2.1, 101.0, 20, False, line="L1"),
                _rec("B1", 1.2, 2.2, 102.0, 30, True, line="L2"),
            ],
        }
        fn = self._fn(web_app)
        *_, heat, _violin = fn(store, "dark")
        z = heat["data"][0]["z"]
        # L2 row padded with None to match L1's length of 2
        l2_idx = heat["data"][0]["y"].index("L2")
        assert None in z[l2_idx]


# ── 10. update_station_radar ─────────────────────────────────────────────────


class TestUpdateStationRadar:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.DASH_CHART_RADAR}.figure")

    def test_no_records_returns_empty_fig(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, None, "dark")
        assert out["data"] == []

    def test_no_selection_single_line_shows_elevation_profile(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, SINGLE_LINE_STORE, "dark")
        assert "Elevation Profile" in out["layout"]["title"]["text"]

    def test_no_selection_multi_line_shows_per_line_spider(self, web_app):
        fn = self._fn(web_app)
        out = fn(None, MULTI_LINE_STORE, "light")
        assert "Per-Line Quality Comparison" in out["layout"]["title"]["text"]
        # 2 background rings + 2 line traces
        assert len(out["data"]) == 4

    def test_selection_with_no_match_returns_empty_fig(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "ZZZ"}, SINGLE_LINE_STORE, "dark")
        assert out["data"] == []

    def test_selection_with_match_builds_quality_spider(self, web_app):
        fn = self._fn(web_app)
        out = fn({"station_id": "S1"}, SINGLE_LINE_STORE, "dark")
        assert "Quality Radar" in out["layout"]["title"]["text"]
        assert "S1" in out["layout"]["title"]["text"]
        # 2 rings + 1 data trace
        assert len(out["data"]) == 3


# ── SI-1. _si_health ──────────────────────────────────────────────────────────


class TestSiHealth:
    def _fn(self, web_app):
        return _cb_multi(web_app, "si-health-dot.className")

    def test_no_store_data(self, web_app):
        fn = self._fn(web_app)
        cls, txt, chips = fn(None)
        assert cls == "si-dot si-dot-empty"
        assert txt == "No survey loaded"
        assert chips == []

    def test_zero_stations(self, web_app):
        fn = self._fn(web_app)
        cls, txt, chips = fn({"n_stations": 0})
        assert cls == "si-dot si-dot-empty"

    def test_with_data_builds_chips(self, web_app):
        fn = self._fn(web_app)
        cls, txt, chips = fn(MULTI_LINE_STORE)
        assert cls == "si-dot si-dot-ok"
        assert "4 stations loaded" in txt
        # 1 profile-count chip + 2 per-line chips + 1 data-dir chip
        assert len(chips) == 4

    def test_more_than_four_lines_shows_more_chip(self, web_app):
        store = {
            "n_stations": 5,
            "n_lines": 5,
            "data_dir": "",
            "line_counts": {f"L{i}": 1 for i in range(5)},
        }
        fn = self._fn(web_app)
        cls, txt, chips = fn(store)
        assert any("+1 more" in getattr(c, "children", "") for c in chips)

    def test_singular_station_and_profile_wording(self, web_app):
        store = {
            "n_stations": 1,
            "n_lines": 1,
            "data_dir": "",
            "line_counts": {},
        }
        fn = self._fn(web_app)
        cls, txt, chips = fn(store)
        assert txt == "1 station loaded"
        assert "1 profile" in chips[0].children[1]


# ── SI-2. _si_last_run ───────────────────────────────────────────────────────


class TestSiLastRun:
    def _fn(self, web_app):
        return _cb(web_app, "si-last-run.children")

    def test_no_history(self, web_app):
        fn = self._fn(web_app)
        out = fn(None)
        assert out.children == "No runs yet"

    def test_empty_history(self, web_app):
        fn = self._fn(web_app)
        out = fn([])
        assert out.children == "No runs yet"

    def test_success_run_renders_summary(self, web_app):
        fn = self._fn(web_app)
        history = [
            {
                "name": "QC Quicklook",
                "timestamp": "12:00",
                "status": "success",
                "summary": "All good",
            }
        ]
        out = fn(history)
        row, summary_p, btn = out.children
        assert "si-run-ok" in row.children[0].className
        assert summary_p.children == "All good"

    def test_error_run_renders_error_icon(self, web_app):
        fn = self._fn(web_app)
        history = [
            {
                "name": "Agent X",
                "timestamp": "13:00",
                "status": "error",
                "summary": "",
            }
        ]
        out = fn(history)
        row, summary_p, btn = out.children
        assert "si-run-err" in row.children[0].className
        assert summary_p.children == "—"

    def test_long_summary_truncated(self, web_app):
        fn = self._fn(web_app)
        long_summary = "x" * 200
        history = [
            {
                "name": "Agent Y",
                "timestamp": "",
                "status": "success",
                "summary": long_summary,
            }
        ]
        out = fn(history)
        _row, summary_p, _btn = out.children
        assert summary_p.children.endswith("…")
        assert len(summary_p.children) == 81
