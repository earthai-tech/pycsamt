# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for the pure/light helpers in
pycsamt.app.web.callbacks.tools (Tools off-canvas panel).

Complements test_web_tools_validator.py (``_validate_sites_rows``) and
test_web_tools_converter.py (converter/batch-export helpers), which
already cover the heavier pipeline functions.
"""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pandas as pd


class _FakeZ:
    def __init__(self, freq, z, z_err=None):
        self.freq = np.asarray(freq, dtype=float)
        self.z = np.asarray(z, dtype=complex)
        self.z_err = z_err


class _FakeEDI:
    def __init__(self, station, lat=None, lon=None, z_obj=None):
        self.station = station
        self.lat = lat
        self.lon = lon
        self.Z = z_obj


def _z(freq):
    arr = np.zeros((len(freq), 2, 2), dtype=complex)
    arr[:, 0, 1] = 1.0 + 1.0j
    arr[:, 1, 0] = -1.0 - 1.0j
    return arr


# ── trivial Dash component builders ─────────────────────────────────────────


def test_no_data_msg_is_a_div():
    from dash import html

    from pycsamt.app.web.callbacks.tools import _no_data_msg

    div = _no_data_msg()
    assert isinstance(div, html.Div)
    assert div.className == "tool-no-data"


def test_run_btn_sets_id_and_label():
    from pycsamt.app.web.callbacks.tools import _run_btn

    btn = _run_btn("run-strike", "Run strike", icon="magic")
    assert btn.id == "run-strike"
    assert btn.n_clicks == 0
    assert btn.children[1] == "Run strike"


def test_action_bar_wraps_children():
    from pycsamt.app.web.callbacks.tools import _action_bar

    bar = _action_bar("a", "b")
    assert bar.children == ["a", "b"]
    assert bar.className == "tool-action-bar"


def test_plotly_config_sets_filename():
    from pycsamt.app.web.callbacks.tools import _plotly_config

    cfg = _plotly_config("my_plot")
    assert cfg["toImageButtonOptions"]["filename"] == "my_plot"
    assert cfg["displaylogo"] is False


def test_path_row_includes_file_and_folder_buttons_when_requested():
    from pycsamt.app.web.callbacks.tools import _path_row

    row = _path_row(
        "input-1",
        "Choose a path",
        file_btn_id="file-btn",
        folder_btn_id="folder-btn",
    )
    child_ids = [getattr(c, "id", None) for c in row.children if hasattr(c, "id")]
    assert "input-1" in child_ids
    assert "file-btn" in child_ids
    assert "folder-btn" in child_ids


def test_path_row_omits_buttons_when_not_requested():
    from pycsamt.app.web.callbacks.tools import _path_row

    row = _path_row("input-1", "Choose a path")
    assert len(row.children) == 1


def test_placeholder_body_uses_tool_label():
    from pycsamt.app.web.callbacks.tools import _placeholder_body

    div = _placeholder_body("strike")
    assert div.className == "tool-placeholder"


def test_metric_chip_shows_label_and_value():
    from pycsamt.app.web.callbacks.tools import _metric_chip

    chip = _metric_chip("Stations", "12")
    assert chip.children[0].children == "Stations"
    assert chip.children[1].children == "12"


def test_err_and_warn_spans():
    from pycsamt.app.web.callbacks.tools import _err, _warn

    e = _err(ValueError("boom"))
    w = _warn("careful")
    assert "boom" in e.children
    assert e.className == "text-danger small"
    assert "careful" in w.children
    assert w.className == "text-warning small"


def test_coming_soon_mentions_tool_name():
    from pycsamt.app.web.callbacks.tools import _coming_soon

    div = _coming_soon("3-D Inversion")
    assert "3-D Inversion" in div.children[1].children


# ── active line / station scoping (pure dict logic) ─────────────────────────


def test_active_line_names_prefers_active_store():
    from pycsamt.app.web.callbacks.tools import _active_line_names

    store = {
        "station_records": [
            {"Line": "L1"},
            {"Line": "L2"},
        ]
    }
    assert _active_line_names(store, {"active": ["L2"]}) == ["L2"]


def test_active_line_names_falls_back_to_all_records():
    from pycsamt.app.web.callbacks.tools import _active_line_names

    store = {
        "station_records": [
            {"Line": "L2"},
            {"Line": "L1"},
            {"Line": "L1"},
        ]
    }
    assert _active_line_names(store, None) == ["L1", "L2"]


def test_station_options_for_lines_filters_by_active_and_selection():
    from pycsamt.app.web.callbacks.tools import _station_options_for_lines

    store = {
        "station_records": [
            {"ID": "S1", "Line": "L1"},
            {"ID": "S2", "Line": "L2"},
            {"ID": "S3", "Line": "L1"},
        ]
    }
    opts = _station_options_for_lines(store, {"active": ["L1"]})
    assert [o["value"] for o in opts] == ["S1", "S3"]

    opts2 = _station_options_for_lines(
        store, {"active": ["L1", "L2"]}, selected_lines=["L2"]
    )
    assert [o["value"] for o in opts2] == ["S2"]


def test_station_line_map_defaults_to_unassigned():
    from pycsamt.app.web.callbacks.tools import _station_line_map

    store = {
        "station_records": [
            {"ID": "S1", "Line": "L1"},
            {"ID": "S2"},
        ]
    }
    mapping = _station_line_map(store)
    assert mapping == {"S1": "L1", "S2": "Unassigned"}


def test_filter_sites_by_station_ids_passthrough_when_empty():
    from pycsamt.app.web.callbacks.tools import _filter_sites_by_station_ids

    sentinel = object()
    assert _filter_sites_by_station_ids(sentinel, None) is sentinel
    assert _filter_sites_by_station_ids(sentinel, []) is sentinel


def test_filter_sites_by_station_ids_returns_none_when_nothing_matches():
    from pycsamt.app.web.callbacks.tools import _filter_sites_by_station_ids

    sites = [_FakeEDI("S1"), _FakeEDI("S2")]
    assert _filter_sites_by_station_ids(sites, ["ZZZ"]) is None


def test_filter_sites_by_station_ids_falls_back_on_error():
    from pycsamt.app.web.callbacks.tools import _filter_sites_by_station_ids

    # A lookup failure (e.g. _iter_items itself raising) should degrade to
    # returning the original sites, not propagate the exception.
    sentinel = object()
    with patch(
        "pycsamt.emtools._core._iter_items",
        side_effect=RuntimeError("boom"),
    ):
        result = _filter_sites_by_station_ids(sentinel, ["S1"])
    assert result is sentinel


def test_scoped_sites_reports_expired_session():
    from pycsamt.app.web.callbacks.tools import _scoped_sites

    with patch("pycsamt.app.web.cache.cache_get", return_value=None):
        sites, err = _scoped_sites("sess-1", {}, {})
    assert sites is None
    assert "expired" in err.lower()


def test_scoped_sites_returns_sites_when_no_scoping_requested():
    from pycsamt.app.web.callbacks.tools import _scoped_sites

    sentinel = object()
    with patch("pycsamt.app.web.cache.cache_get", return_value=sentinel):
        sites, err = _scoped_sites("sess-1", {}, {})
    assert sites is sentinel
    assert err is None


# ── theming / formatting ─────────────────────────────────────────────────────


def test_tool_theme_light_vs_dark():
    from pycsamt.app.web.callbacks.tools import _tool_theme

    assert _tool_theme("light")["bg"] == "#ffffff"
    assert _tool_theme("dark")["bg"] == "#1e1e2e"
    assert _tool_theme(None)["bg"] == "#1e1e2e"


def test_tool_table_styles_returns_three_parts():
    from pycsamt.app.web.callbacks.tools import _tool_table_styles

    cell, header, conditional = _tool_table_styles("light")
    assert cell["backgroundColor"] == "#ffffff"
    assert header["backgroundColor"] == "#e6e9ef"
    assert isinstance(conditional, list) and len(conditional) == 2


def test_style_fig_sets_layout_height():
    import plotly.graph_objects as go

    from pycsamt.app.web.callbacks.tools import _style_fig

    fig = go.Figure()
    _style_fig(fig, "Value", height=300)
    assert fig.layout.height == 300
    assert fig.layout.yaxis.title.text == "Value"


def test_safe_stem_strips_unsafe_characters():
    from pycsamt.app.web.callbacks.tools import _safe_stem

    assert _safe_stem("My Figure #1 / 2024!") == "My_Figure_1_2024"
    assert _safe_stem("") == "figure"
    assert _safe_stem(None) == "figure"


def test_fmt_float_handles_non_finite_and_precision():
    from pycsamt.app.web.callbacks.tools import _fmt_float

    assert _fmt_float(float("nan")) == "—"
    assert _fmt_float("not a number") == "—"
    assert _fmt_float(3.14159, nd=2) == "3.1"
    assert _fmt_float(3.14159, nd=6) == "3.141590"


def test_is_finite_number():
    from pycsamt.app.web.callbacks.tools import _is_finite_number

    assert _is_finite_number(1.5) is True
    assert _is_finite_number(float("nan")) is False
    assert _is_finite_number("oops") is False


def test_dd2dms_north_and_south():
    from pycsamt.app.web.callbacks.tools import _dd2dms

    north = _dd2dms(40.7127, axis="lat")
    south = _dd2dms(-33.8688, axis="lat")
    east = _dd2dms(7.75, axis="lon")

    assert north.endswith("N")
    assert south.endswith("S")
    assert east.endswith("E")
    assert north.startswith("40°")


# ── station count / lat-lon extraction (fake EDI objects) ───────────────────


def test_station_count_uses_iter_items():
    from pycsamt.app.web.callbacks.tools import _station_count

    sites = [_FakeEDI("S1"), _FakeEDI("S2"), _FakeEDI("S3")]
    assert _station_count(sites) == 3


def test_station_count_degrades_to_zero_on_error():
    from pycsamt.app.web.callbacks.tools import _station_count

    assert _station_count(object()) == 0


def test_extract_lat_lon_reads_direct_attributes():
    from pycsamt.app.web.callbacks.tools import _extract_lat_lon

    edi = _FakeEDI("S1", lat=10.5, lon=20.5)
    lat, lon = _extract_lat_lon(edi)
    assert lat == 10.5
    assert lon == 20.5


def test_extract_lat_lon_returns_none_when_missing():
    from pycsamt.app.web.callbacks.tools import _extract_lat_lon

    lat, lon = _extract_lat_lon(object())
    assert lat is None
    assert lon is None


def test_first_existing_attr_returns_first_present():
    from pycsamt.app.web.callbacks.tools import _first_existing_attr

    obj = _FakeEDI("S1", lat=None, lon=5.0)
    assert _first_existing_attr(obj, ("lat", "lon")) == 5.0
    assert _first_existing_attr(obj, ("missing",)) is None


def test_z_block_with_errors_reads_z_err():
    from pycsamt.app.web.callbacks.tools import _z_block_with_errors

    freq = np.asarray([1.0, 10.0])
    z_err = np.ones((2, 2, 2))
    edi = _FakeEDI("S1", z_obj=_FakeZ(freq, _z(freq), z_err=z_err))

    Z, z, fr, ze = _z_block_with_errors(edi)
    assert Z is edi.Z
    np.testing.assert_allclose(fr, freq)
    np.testing.assert_allclose(ze, z_err)


# ── strike scope / plotting helpers ──────────────────────────────────────────


def test_strike_band_from_choice():
    from pycsamt.app.web.callbacks.tools import _strike_band_from_choice

    assert _strike_band_from_choice("high") == (0.0, 0.01)
    assert _strike_band_from_choice("low") == (1.0, 1.0e12)
    assert _strike_band_from_choice("all") is None
    assert _strike_band_from_choice(None) is None


def test_strike_scope_label_prioritises_stations_then_lines_then_active():
    from pycsamt.app.web.callbacks.tools import _strike_scope_label

    assert (
        _strike_scope_label(None, None, None, ["S1", "S2"]) == "2 selected station(s)"
    )
    assert _strike_scope_label(None, None, ["L1", "L2"], None) == "L1, L2"
    assert (
        _strike_scope_label(None, None, ["L1", "L2", "L3", "L4"], None)
        == "L1, L2, L3 +1"
    )

    store = {"station_records": [{"Line": "L1"}, {"Line": "L2"}]}
    assert (
        _strike_scope_label(store, {"active": ["L1", "L2"]}, None, None)
        == "All active lines (2)"
    )
    assert _strike_scope_label(None, None, None, None) == "All loaded stations"


def test_validator_status_bar_only_renders_nonzero_segments():
    from pycsamt.app.web.callbacks.tools import _validator_status_bar

    bar = _validator_status_bar(8, 2, 0)
    assert len(bar.children) == 2
    widths = {seg.className.split()[-1]: seg.style["width"] for seg in bar.children}
    assert widths["tool-valid-pass"] == "80.00%"
    assert widths["tool-valid-warn"] == "20.00%"


def _strike_df():
    return pd.DataFrame(
        {
            "line": ["L1", "L1", "L2"],
            "ang_axial": [10.0, 20.0, 170.0],
        }
    )


def test_strike_rose_figure_has_one_trace_per_line():
    from pycsamt.app.web.callbacks.tools import _strike_rose_figure

    fig = _strike_rose_figure(_strike_df(), "Strike", theme="dark")
    assert len(fig.data) == 2
    assert {tr.name for tr in fig.data} == {"L1", "L2"}


def test_strike_box_figure_has_one_trace_per_line():
    from pycsamt.app.web.callbacks.tools import _strike_box_figure

    fig = _strike_box_figure(_strike_df(), theme="light")
    assert len(fig.data) == 2


def test_render_strike_result_builds_full_layout():
    from pycsamt.app.web.callbacks.tools import _render_strike_result

    payload = {
        "records": [
            {"line": "L1", "ang_axial": 10.0},
            {"line": "L1", "ang_axial": 20.0},
        ],
        "method_label": "Consensus",
        "scope_label": "L1",
        "n_stations": 2,
        "median_strike": 15.0,
        "median_iqr": 5.0,
        "columns": ["Station", "Strike"],
        "table": [{"Station": "S1", "Strike": 10.0}],
        "page_size": 5,
    }
    div = _render_strike_result(payload, "dark")
    # metric row + 2 graphs + table
    assert len(div.children) == 4


def test_render_stored_tool_result_returns_none_for_unknown_tool():
    from pycsamt.app.web.callbacks.tools import _render_stored_tool_result

    assert _render_stored_tool_result("strike", None, "dark") is None
    assert _render_stored_tool_result("elevation", {}, "dark") is None


# ── restore-from-store dispatch ──────────────────────────────────────────────


def test_restore_from_store_returns_none_when_empty():
    from pycsamt.app.web.callbacks.tools import _restore_from_store

    assert _restore_from_store(None) is None
    assert _restore_from_store({}) is None


def test_restore_from_store_text_type():
    from pycsamt.app.web.callbacks.tools import _restore_from_store

    span = _restore_from_store(
        {"type": "text", "content": "hello", "cls": "small text-muted"}
    )
    assert span.children == "hello"
    assert span.className == "small text-muted"


def test_restore_from_store_unknown_type_returns_none():
    from pycsamt.app.web.callbacks.tools import _restore_from_store

    assert _restore_from_store({"type": "nonsense"}) is None


def test_restore_from_store_html_type_with_images():
    from pycsamt.app.web.callbacks.tools import _restore_from_store

    div = _restore_from_store({"type": "html", "imgs": ["Zm9v"]})
    assert div is not None


def test_restore_from_store_fig_json_type():
    import plotly.graph_objects as go

    from pycsamt.app.web.callbacks.tools import _restore_from_store

    fig = go.Figure(data=[go.Scatter(x=[1, 2], y=[3, 4])])
    graph = _restore_from_store({"type": "fig_json", "fig": fig.to_json()})
    assert list(graph.figure["data"][0]["y"]) == [3, 4]


# ── elevation summary table ──────────────────────────────────────────────────


def test_elev_summary_table_counts_ok_stations():
    from pycsamt.app.web.callbacks.tools import _elev_summary_table

    df = pd.DataFrame(
        {
            "station": ["S1", "S2"],
            "lat": [10.0, 11.0],
            "lon": [20.0, 21.0],
            "elev": [500.0, float("nan")],
        }
    )
    div = _elev_summary_table(df, title="Elevation fetch")
    header_text = div.children[0].children[1]
    assert "1/2 station(s)" in header_text
    assert "Elevation fetch" in header_text
