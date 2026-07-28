# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the dimensionality / frequency-editor / batch / elevation tool
panels in ``pycsamt.app.web.callbacks.tools``.

Strategy
--------
Real WILLY EDI data (``data/AMT/WILLY_DATA/L18PLT``) is loaded once via
``pycsamt.emtools.ensure_sites`` and cached with
``pycsamt.app.web.cache.cache_set`` so every callback exercised here runs
against real ``Sites``/EDI objects and real ``pycsamt.emtools`` analysis
routines (``classify_dimensionality``, ``edit_frequencies_by_confidence``,
``select_band``, ``regrid_logspace``, ``align_grid``, ``smooth_mavg``,
``decimate_step``, ``drop_duplicates``, the various pseudosection/heatmap
plot functions...) rather than mocks. File-writing tool paths (batch
figure export, edited-EDI export) write into ``tmp_path`` rather than
mocking the filesystem.

Body-builder functions (``_dimensionality_body``, ``_freq_editor_body``,
``_batch_body``, ``_elevation_body``) are pure/static Dash-layout
constructors -- they are called directly (no Dash server needed) since
they are not themselves registered as callbacks. The actual wired
callbacks (``_run_dimensionality``, ``_run_freq_editor``, ``_run_batch``,
the elevation-panel callbacks, ...) are looked up from
``web_app.callback_map`` using the shared helper pattern established in
``test_web_callbacks_lines.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from dash import no_update

import pycsamt.app.web.callbacks.tools as tools_mod
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


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


# ── Real WILLY EDI fixtures ──────────────────────────────────────────────────


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
def active_lines_both():
    return {"active": ["L1", "L2"], "all": ["L1", "L2"]}


@pytest.fixture
def cached_session(willy_sites):
    from pycsamt.app.web.cache import cache_set

    session_id = "test-tools-dim-session"
    cache_set(session_id, willy_sites)
    return session_id


NO_DATA = tools_mod._no_data_msg()


# ══════════════════════════════════════════════════════════════════════════
# 1. Dimensionality panel
# ══════════════════════════════════════════════════════════════════════════


class TestDimensionalityBody:
    def test_zero_clicks_returns_no_data(self):
        assert tools_mod._dimensionality_body(0, NO_DATA) is NO_DATA

    def test_none_clicks_returns_no_data(self):
        assert tools_mod._dimensionality_body(None, NO_DATA) is NO_DATA

    def test_builds_body_with_line_and_station_options(
        self, store_data_willy, active_lines_both
    ):
        body = tools_mod._dimensionality_body(
            1, NO_DATA, store_data_willy, active_lines_both, None
        )
        # dropdown line options reflect the two synthetic lines
        line_dd = body.children[1].children[0].children[1]
        values = [o["value"] for o in line_dd.options]
        assert set(values) == {"L1", "L2"}

    def test_builds_body_without_store_data(self):
        body = tools_mod._dimensionality_body(1, NO_DATA, None, None, None)
        line_dd = body.children[1].children[0].children[1]
        assert line_dd.options == []


class TestSyncDimStationScope:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-dim-stations.options", "tool-dim-lines")

    def test_filters_stations_to_selected_line(
        self, web_app, store_data_willy, active_lines_both
    ):
        opts, value = self._fn(web_app)(["L1"], store_data_willy, active_lines_both)
        assert value is None
        assert opts  # non-empty
        assert all("L1" in o["label"] for o in opts)

    def test_no_selection_returns_all_active_stations(
        self, web_app, store_data_willy, active_lines_both
    ):
        opts, _ = self._fn(web_app)(None, store_data_willy, active_lines_both)
        assert len(opts) == len(store_data_willy["station_records"])


class TestRunDimensionality:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-dim-out.children", "tool-dim-run")

    def test_no_clicks_no_update(self, web_app):
        out, saved = self._fn(web_app)(
            None,
            None,
            None,
            "psection",
            3.0,
            0.2,
            None,
            None,
            10.0,
            24,
            None,
            None,
            "sess",
            "dark",
            None,
        )
        assert out is no_update
        assert saved is no_update

    def test_session_expired_shows_warning(self, web_app, store_data_willy):
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            "psection",
            3.0,
            0.2,
            None,
            None,
            10.0,
            24,
            store_data_willy,
            None,
            "not-a-real-session",
            "dark",
            None,
        )
        assert "Session expired" in out.children
        assert saved is no_update

    def test_psection_view_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            "psection",
            3.0,
            0.2,
            None,
            None,
            10.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        assert saved["dimensionality"]["type"] == "html"
        assert saved["dimensionality"]["imgs"]  # at least one PNG produced
        # summary bar always appended alongside the requested view
        assert len(out.children) >= 2

    def test_occupancy_view_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            "occupancy",
            3.0,
            0.2,
            None,
            None,
            10.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "light",
            {"other-tool": {"foo": "bar"}},
        )
        assert saved["dimensionality"]["imgs"]
        assert "other-tool" in saved  # merges into existing saved-outputs

    def test_map_view_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            "map",
            3.0,
            0.2,
            None,
            None,
            15.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        assert saved["dimensionality"]["imgs"]

    def test_scatter_view_with_period_filter_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            "scatter",
            3.0,
            0.2,
            1e-3,
            1e3,
            10.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        # scatter view produces a dcc.Graph + summary bar -> 2 children
        assert len(out.children) == 2

    def test_scatter_view_period_filter_excludes_everything(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        """Period range far outside the survey's actual periods empties
        df_filt -> hits the 'No data in selected period range' warning
        branch, and the always-on summary bar is skipped too since it
        also filters on df_filt."""
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            "scatter",
            3.0,
            0.2,
            1e6,
            1e7,
            10.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        texts = [getattr(c, "children", "") for c in out.children]
        assert any("No data in selected period range" in str(t) for t in texts)

    def test_all_views_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved = self._fn(web_app)(
            1,
            ["L1"],
            None,
            "all",
            3.0,
            0.2,
            None,
            None,
            10.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        # psection + occupancy + map + scatter + summary bar == 5 sections
        assert len(out.children) == 5
        assert len(saved["dimensionality"]["imgs"]) == 3  # 3 mpl figures

    def test_defaults_used_when_thresholds_none(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        """skew_th/ellipt_th/map_period/n_bands fall back to hardcoded
        defaults when the Input values are None (e.g. cleared by the
        user)."""
        out, saved = self._fn(web_app)(
            1,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        assert saved["dimensionality"]["type"] == "html"

    def test_station_scope_selection(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        first_id = store_data_willy["station_records"][0]["ID"]
        out, saved = self._fn(web_app)(
            1,
            None,
            [first_id],
            "psection",
            3.0,
            0.2,
            None,
            None,
            10.0,
            24,
            store_data_willy,
            active_lines_both,
            cached_session,
            "dark",
            None,
        )
        assert saved["dimensionality"]["imgs"]


# ══════════════════════════════════════════════════════════════════════════
# 2. Frequency editor panel
# ══════════════════════════════════════════════════════════════════════════


class TestFreqEditorBody:
    def test_zero_clicks_returns_no_data(self):
        assert tools_mod._freq_editor_body(0, NO_DATA) is NO_DATA

    def test_builds_body_with_scope_options(self, store_data_willy, active_lines_both):
        body = tools_mod._freq_editor_body(
            1, NO_DATA, store_data_willy, active_lines_both, None
        )
        assert body.children  # non-trivial layout produced

    def test_last_output_restores_from_store(self, store_data_willy, active_lines_both):
        stored = {"type": "text", "content": "hello", "cls": "small"}
        body = tools_mod._freq_editor_body(
            1, NO_DATA, store_data_willy, active_lines_both, stored
        )
        # last output area is the very last child of the panel
        out_area = body.children[-1]
        assert out_area.children.children == "hello"


class TestFreqBeforeAfterHeatmaps:
    def test_produces_alert_and_two_heatmap_sections(self, willy_sites, cached_session):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dash import html as dhtml

        from pycsamt.emtools.frequency import plot_coverage_quality_heatmap

        def _mpl_wrap(plot_fn, *args, **kw):
            fig, ax = plt.subplots()
            plot_fn(*args, ax=ax, **kw)
            plt.close(fig)
            return dhtml.Img(src="data:image/png;base64,")

        def _style_ax(ax, fig):
            pass

        children = tools_mod._freq_before_after_heatmaps(
            willy_sites,
            willy_sites,
            _mpl_wrap,
            _style_ax,
            "#000",
            "#fff",
            "#333",
            "#111",
            plot_coverage_quality_heatmap,
        )
        assert len(children) == 3  # alert + before + after

    def test_counts_helper_swallows_exceptions(self):
        # A non-Sites object makes the internal ``_counts`` closure raise
        # and fall back to (0, 0) -- still produces a valid alert message.
        children = tools_mod._freq_before_after_heatmaps(
            object(),
            object(),
            lambda *a, **k: None,
            lambda *a, **k: None,
            "#000",
            "#fff",
            "#333",
            "#111",
            lambda *a, **k: None,
        )
        assert "0 stations" in children[0].children[1]


class TestSyncFreqStationScope:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-freq-stations.options", "tool-freq-lines")

    def test_filters_by_selected_lines(
        self, web_app, store_data_willy, active_lines_both
    ):
        opts, value = self._fn(web_app)(["L2"], store_data_willy, active_lines_both)
        assert value is None
        assert all("L2" in o["label"] for o in opts)


class TestSyncFreqSource:
    def _fn(self, web_app):
        return _cb(web_app, "tool-freq-folder-row.style")

    def test_folder_source_shows_row(self, web_app):
        assert self._fn(web_app)("folder") == {"display": "block"}

    def test_session_source_hides_row(self, web_app):
        assert self._fn(web_app)("session") == {"display": "none"}

    def test_none_source_hides_row(self, web_app):
        assert self._fn(web_app)(None) == {"display": "none"}


class TestSyncFreqActionParams:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-freq-conf-params.style", "tool-freq-action")

    def test_confidence_shows_conf_params_only(self, web_app):
        styles = self._fn(web_app)("confidence")
        assert styles[0] == {"display": "block"}
        assert all(s == {"display": "none"} for s in styles[1:])

    def test_band_shows_band_params_only(self, web_app):
        styles = self._fn(web_app)("band")
        assert styles[1] == {"display": "block"}
        assert styles[0] == {"display": "none"}

    def test_regrid_shows_regrid_params_only(self, web_app):
        styles = self._fn(web_app)("regrid")
        assert styles[2] == {"display": "block"}

    def test_align_shows_align_params_only(self, web_app):
        styles = self._fn(web_app)("align")
        assert styles[3] == {"display": "block"}

    def test_smooth_shows_smooth_params_only(self, web_app):
        styles = self._fn(web_app)("smooth")
        assert styles[4] == {"display": "block"}

    def test_diagnostics_hides_all_params(self, web_app):
        styles = self._fn(web_app)("diagnostics")
        assert all(s == {"display": "none"} for s in styles)


class TestRunFreqEditor:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-freq-out.children", "tool-freq-run")

    def _call(self, fn, action, session_id, store_data, active_lines, **kw):
        defaults = dict(
            n=1,
            source="session",
            folder=None,
            sel_lines=None,
            sel_stations=None,
            action=action,
            mode="recover",
            method="composite",
            thresh=0.5,
            ci_hi=0.9,
            ci_lo=0.5,
            also="both",
            interp="linear",
            reject="drop",
            fmin=None,
            fmax=None,
            keep_str=None,
            rg_fmin=None,
            rg_fmax=None,
            rg_ppd=6,
            rg_method="nearest",
            align_mode="union",
            align_ref=None,
            smooth_op="smooth",
            smooth_k=3,
            smooth_step=2,
            store_data=store_data,
            active_lines_store=active_lines,
            session_id=session_id,
            theme="dark",
            saved_outputs=None,
        )
        defaults.update(kw)
        return fn(
            defaults["n"],
            defaults["source"],
            defaults["folder"],
            defaults["sel_lines"],
            defaults["sel_stations"],
            defaults["action"],
            defaults["mode"],
            defaults["method"],
            defaults["thresh"],
            defaults["ci_hi"],
            defaults["ci_lo"],
            defaults["also"],
            defaults["interp"],
            defaults["reject"],
            defaults["fmin"],
            defaults["fmax"],
            defaults["keep_str"],
            defaults["rg_fmin"],
            defaults["rg_fmax"],
            defaults["rg_ppd"],
            defaults["rg_method"],
            defaults["align_mode"],
            defaults["align_ref"],
            defaults["smooth_op"],
            defaults["smooth_k"],
            defaults["smooth_step"],
            defaults["store_data"],
            defaults["active_lines_store"],
            defaults["session_id"],
            defaults["theme"],
            defaults["saved_outputs"],
        )

    def test_no_clicks_no_update(self, web_app):
        out, saved, disabled = self._fn(web_app)(
            None,
            "session",
            None,
            None,
            None,
            "diagnostics",
            "recover",
            "composite",
            0.5,
            0.9,
            0.5,
            "both",
            "linear",
            "drop",
            None,
            None,
            None,
            None,
            None,
            6,
            "nearest",
            "union",
            None,
            "smooth",
            3,
            2,
            None,
            None,
            "sess",
            "dark",
            None,
        )
        assert out is no_update and saved is no_update and disabled is no_update

    def test_folder_source_missing_dir(self, web_app):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "diagnostics",
            "sess",
            None,
            None,
            source="folder",
            folder="Z:/definitely/not/a/real/path",
        )
        assert "Folder not found" in out.children
        assert disabled is True

    def test_session_expired(self, web_app, store_data_willy):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "diagnostics",
            "no-such-session",
            store_data_willy,
            None,
        )
        assert "Session expired" in out.children
        assert disabled is True

    def test_diagnostics_action_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "diagnostics",
            cached_session,
            store_data_willy,
            active_lines_both,
        )
        assert saved["freq-editor"]["imgs"]
        assert disabled is True  # diagnostics never sets has_edit

    def test_folder_source_real_edi_dir(self, web_app, active_lines_both):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "diagnostics",
            "unused-session",
            None,
            active_lines_both,
            source="folder",
            folder=str(_WILLY_L18),
        )
        assert saved["freq-editor"]["imgs"]

    def test_confidence_action_recover_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "confidence",
            cached_session,
            store_data_willy,
            active_lines_both,
            mode="recover",
        )
        assert saved["freq-editor"]["imgs"]
        assert disabled is False  # an edit was applied -> export enabled
        assert "Mode:" in out.children[0].children[1]

    def test_confidence_action_drop_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "confidence",
            cached_session,
            store_data_willy,
            active_lines_both,
            mode="drop",
            reject="mask",
        )
        assert disabled is False

    def test_band_action_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "band",
            cached_session,
            store_data_willy,
            active_lines_both,
            fmin=1.0,
            fmax=1000.0,
        )
        assert disabled is False
        assert saved["freq-editor"]["imgs"]

    def test_band_action_unparseable_keep_list(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "band",
            cached_session,
            store_data_willy,
            active_lines_both,
            keep_str="not,a,number",
        )
        texts = [str(getattr(c, "children", c)) for c in out.children]
        assert any("Could not parse" in t for t in texts)

    def test_regrid_action_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "regrid",
            cached_session,
            store_data_willy,
            active_lines_both,
            rg_ppd=4,
            rg_method="nearest",
        )
        assert disabled is False
        assert saved["freq-editor"]["imgs"]

    def test_align_action_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "align",
            cached_session,
            store_data_willy,
            active_lines_both,
            align_mode="union",
        )
        assert disabled is False
        assert saved["freq-editor"]["imgs"]

    def test_smooth_op_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "smooth",
            cached_session,
            store_data_willy,
            active_lines_both,
            smooth_op="smooth",
            smooth_k=3,
        )
        assert disabled is False

    def test_decimate_op_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "smooth",
            cached_session,
            store_data_willy,
            active_lines_both,
            smooth_op="decimate",
            smooth_step=2,
        )
        assert disabled is False

    def test_drop_duplicates_op_real_data(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "smooth",
            cached_session,
            store_data_willy,
            active_lines_both,
            smooth_op="dedupe",
        )
        assert disabled is False

    def test_saved_outputs_merges_with_existing(
        self, web_app, store_data_willy, active_lines_both, cached_session
    ):
        out, saved, disabled = self._call(
            self._fn(web_app),
            "diagnostics",
            cached_session,
            store_data_willy,
            active_lines_both,
            saved_outputs={"dimensionality": {"type": "html", "imgs": []}},
        )
        assert "dimensionality" in saved
        assert "freq-editor" in saved


class TestExportFreqEdis:
    def _fn(self, web_app):
        return _cb(web_app, "tool-freq-export-out.children")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, "x", "{station}.edi", "sess") is no_update

    def test_no_out_path_warns(self, web_app):
        out = self._fn(web_app)(1, None, "{station}.edi", "sess")
        assert "No output folder" in out.children

    def test_no_edited_data_in_cache_warns(self, web_app, tmp_path):
        out = self._fn(web_app)(
            1, str(tmp_path), "{station}.edi", "no-freq-edit-cached"
        )
        assert "No edited data in cache" in out.children

    def test_exports_edited_edis_to_tmp_path(
        self, web_app, willy_sites, cached_session, tmp_path
    ):
        from pycsamt.app.web.cache import cache_set_freq_edit

        cache_set_freq_edit(cached_session, willy_sites)
        out_dir = tmp_path / "edited_out"
        out = self._fn(web_app)(1, str(out_dir), "{station}.edi", cached_session)
        assert "Exported" in out.children
        written = list(out_dir.glob("*.edi"))
        assert len(written) == len(willy_sites.as_list())


# ══════════════════════════════════════════════════════════════════════════
# 3. Batch export panel
# ══════════════════════════════════════════════════════════════════════════


class TestBatchBody:
    def test_zero_clicks_returns_no_data(self):
        assert tools_mod._batch_body(0, NO_DATA) is NO_DATA

    def test_builds_full_body(self):
        body = tools_mod._batch_body(1, NO_DATA)
        assert body.children


class TestRunBatch:
    def _fn(self, web_app):
        return _cb(web_app, "tool-batch-out.children")

    def test_no_clicks_no_update(self, web_app):
        assert (
            self._fn(web_app)(
                None, "/tmp/x", "png", 150, ["map"], [], None, None, "sess"
            )
            is no_update
        )

    def test_no_out_path_warns(self, web_app):
        out = self._fn(web_app)(1, "  ", "png", 150, ["map"], [], None, None, "sess")
        assert "output folder" in out.children

    def test_no_items_warns(self, web_app, tmp_path):
        out = self._fn(web_app)(
            1, str(tmp_path), "png", 150, [], [], None, None, "sess"
        )
        assert "at least one figure type" in out.children

    def test_session_expired_warns(self, web_app, tmp_path):
        out = self._fn(web_app)(
            1,
            str(tmp_path),
            "png",
            150,
            ["map"],
            [],
            None,
            None,
            "no-such-session",
        )
        assert "Session expired" in out.children

    def test_exports_map_only_to_tmp_path(
        self,
        web_app,
        store_data_willy,
        active_lines_both,
        cached_session,
        tmp_path,
    ):
        self._fn(web_app)(
            1,
            str(tmp_path),
            "png",
            96,
            ["map"],
            ["manifest"],
            store_data_willy,
            active_lines_both,
            cached_session,
        )
        # _batch_result_view renders a metric row + DataTable
        assert (tmp_path / "station_map.png").exists()
        assert (tmp_path / "pycsamt_batch_export_manifest.csv").exists()

    def test_exports_curves_and_pseudosections_to_tmp_path(
        self,
        web_app,
        store_data_willy,
        active_lines_both,
        cached_session,
        tmp_path,
    ):
        self._fn(web_app)(
            1,
            str(tmp_path),
            "png",
            96,
            ["curves", "rho_pseudo", "phi_pseudo"],
            [],
            store_data_willy,
            active_lines_both,
            cached_session,
        )
        assert (tmp_path / "mt_curves_rho_phase.png").exists()
        assert (tmp_path / "rho_xy_pseudosection.png").exists()
        assert (tmp_path / "phase_xy_pseudosection.png").exists()

    def test_exports_strike_dimensionality_to_tmp_path(
        self,
        web_app,
        store_data_willy,
        active_lines_both,
        cached_session,
        tmp_path,
    ):
        self._fn(web_app)(
            1,
            str(tmp_path),
            "svg",
            72,
            ["strike"],
            [],
            store_data_willy,
            active_lines_both,
            cached_session,
        )
        assert (tmp_path / "skew_traffic_pseudosection.svg").exists()
        assert (tmp_path / "dimensionality_pseudosection.svg").exists()

    def test_skip_existing_when_not_overwriting(
        self,
        web_app,
        store_data_willy,
        active_lines_both,
        cached_session,
        tmp_path,
    ):
        (tmp_path / "station_map.png").write_bytes(b"placeholder")
        out = self._fn(web_app)(
            1,
            str(tmp_path),
            "png",
            96,
            ["map"],
            [],
            store_data_willy,
            active_lines_both,
            cached_session,
        )
        table = out.children[-1]
        rows = table.data
        assert rows[0]["Status"] == "SKIPPED"

    def test_overwrite_flag_replaces_existing(
        self,
        web_app,
        store_data_willy,
        active_lines_both,
        cached_session,
        tmp_path,
    ):
        (tmp_path / "station_map.png").write_bytes(b"placeholder")
        out = self._fn(web_app)(
            1,
            str(tmp_path),
            "png",
            96,
            ["map"],
            ["overwrite"],
            store_data_willy,
            active_lines_both,
            cached_session,
        )
        table = out.children[-1]
        assert table.data[0]["Status"] == "SAVED"


# ══════════════════════════════════════════════════════════════════════════
# 4. Elevation panel
# ══════════════════════════════════════════════════════════════════════════


class TestElevationBody:
    def test_zero_clicks_returns_no_data(self):
        assert tools_mod._elevation_body(0, NO_DATA) is NO_DATA

    def test_builds_full_body_with_station_options(self, store_data_willy):
        body = tools_mod._elevation_body(1, NO_DATA, store_data_willy)
        assert body.children

    def test_builds_body_without_store_data(self):
        body = tools_mod._elevation_body(1, NO_DATA, None)
        assert body.children


class TestElevUpdateStations:
    def _fn(self, web_app):
        return _cb(web_app, "tool-elev-stations.options")

    def test_no_filter_returns_all(self, web_app, store_data_willy):
        opts = self._fn(web_app)(None, store_data_willy)
        assert len(opts) == len(store_data_willy["station_records"])

    def test_filters_by_selected_lines(self, web_app, store_data_willy):
        opts = self._fn(web_app)(["L1"], store_data_willy)
        assert all("L1" in o["label"] for o in opts)


class TestElevParseUpload:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "tool-elev-upload-store.data", "tool-elev-upload-input"
        )

    def test_no_contents_no_update(self, web_app):
        result = self._fn(web_app)(None, "f.csv")
        assert result == (no_update,) * 10

    def test_parses_csv_with_autodetected_columns(self, web_app):
        import base64

        csv_text = "Station,Lat,Lon,Elevation\nS1,10.0,20.0,300.0\nS2,11.0,21.0,310.0\n"
        b64 = base64.b64encode(csv_text.encode()).decode()
        contents = f"data:text/csv;base64,{b64}"
        (
            data,
            info,
            opt_id,
            opt_lat,
            opt_lon,
            opt_elev,
            v_id,
            v_lat,
            v_lon,
            v_elev,
        ) = self._fn(web_app)(contents, "stations.csv")
        assert v_id == "Station"
        assert v_lat == "Lat"
        assert v_lon == "Lon"
        assert v_elev == "Elevation"
        assert len(data["records"]) == 2

    def test_parse_error_returns_error_span(self, web_app):
        result = self._fn(web_app)("not-a-valid-data-uri", "bad.csv")
        data, info = result[0], result[1]
        assert data is None
        assert "Parse error" in info.children


class TestElevLoadFile:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-elev-raw-store.data", "tool-elev-load-btn")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(
            None, {"records": [], "columns": []}, None, None, None, None
        ) == (no_update, no_update)

    def test_no_upload_data_no_update(self, web_app):
        assert self._fn(web_app)(1, None, None, None, None, None) == (
            no_update,
            no_update,
        )

    def test_loads_rows_with_explicit_columns(self, web_app):
        upload_data = {
            "records": [
                {
                    "Station": "S1",
                    "Lat": 10.0,
                    "Lon": 20.0,
                    "Elevation": 300.0,
                },
                {
                    "Station": "S2",
                    "Lat": 11.0,
                    "Lon": 21.0,
                    "Elevation": 310.0,
                },
            ],
            "columns": ["Station", "Lat", "Lon", "Elevation"],
        }
        rows, out = self._fn(web_app)(
            1, upload_data, "Station", "Lat", "Lon", "Elevation"
        )
        assert len(rows) == 2
        assert rows[0] == {
            "station": "S1",
            "lat": 10.0,
            "lon": 20.0,
            "elev": 300.0,
        }
        assert "loaded" in out.children[0].children[1]

    def test_missing_elevation_column_warns(self, web_app):
        upload_data = {
            "records": [{"Station": "S1"}],
            "columns": ["Station"],
        }
        rows, out = self._fn(web_app)(1, upload_data, "Station", None, None, None)
        assert rows is no_update
        assert "Could not find an elevation column" in out.children


class TestElevFetch:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-elev-out.children", "tool-elev-run")

    def test_no_clicks_no_update(self, web_app):
        result = self._fn(web_app)(
            None, "session", "WGS84", "all", None, None, None, "sess"
        )
        assert result == (no_update,) * 5

    def test_session_source_uses_edi_header_elevations(
        self, web_app, store_data_willy, cached_session
    ):
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "session",
            "WGS84",
            "all",
            None,
            None,
            store_data_willy,
            cached_session,
        )
        assert style == {}
        assert prog == 100
        assert len(rows) == len(store_data_willy["station_records"])
        assert "stations have elevation" in label

    def test_scope_lines_filters_records(
        self, web_app, store_data_willy, cached_session
    ):
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "session",
            "WGS84",
            "lines",
            ["L1"],
            None,
            store_data_willy,
            cached_session,
        )
        assert len(rows) == sum(
            1 for r in store_data_willy["station_records"] if r["Line"] == "L1"
        )

    def test_scope_selected_stations_filters_records(
        self, web_app, store_data_willy, cached_session
    ):
        first_id = store_data_willy["station_records"][0]["ID"]
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "session",
            "WGS84",
            "sel",
            None,
            [first_id],
            store_data_willy,
            cached_session,
        )
        assert len(rows) == 1

    def test_scope_selected_falls_back_to_lines_when_no_stations_picked(
        self, web_app, store_data_willy, cached_session
    ):
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "session",
            "WGS84",
            "sel",
            ["L2"],
            None,
            store_data_willy,
            cached_session,
        )
        assert len(rows) == sum(
            1 for r in store_data_willy["station_records"] if r["Line"] == "L2"
        )

    def test_no_records_after_scope_filter_warns(
        self, web_app, store_data_willy, cached_session
    ):
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "session",
            "WGS84",
            "lines",
            ["NoSuchLine"],
            None,
            store_data_willy,
            cached_session,
        )
        assert "No stations matched the scope" in out.children
        assert style == {"display": "none"}

    def test_online_source_calls_elevation_api(
        self, web_app, store_data_willy, cached_session, monkeypatch
    ):
        """Network APIs are stubbed out -- this exercises the online-fetch
        branch (lat/lon extraction, api_name mapping) without a real HTTP
        call."""
        import numpy as np

        import pycsamt.gis.utils as gis_utils

        def _fake_get_elevation(lats, lons, api_name=None):
            assert api_name == "open_elevation"
            return np.full(len(lats), 123.0)

        monkeypatch.setattr(gis_utils, "get_elevation_from_api", _fake_get_elevation)
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "open-elevation",
            "WGS84",
            "all",
            None,
            None,
            store_data_willy,
            cached_session,
        )
        assert "fetched" in label
        assert any(r.get("elev") == 123.0 for r in rows)

    def test_session_expired_warns(self, web_app, store_data_willy):

        # session_id that was never cache_set -> cache_get returns None
        out, rows, prog, label, style = self._fn(web_app)(
            1,
            "session",
            "WGS84",
            "all",
            None,
            None,
            store_data_willy,
            "definitely-not-cached-session-id",
        )
        assert "Session expired" in out.children


class TestElevCorrect:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-elev-out.children", "tool-elev-correct-btn")

    RAW = [
        {"station": "S1", "lat": 1.0, "lon": 1.0, "elev": 100.0},
        {"station": "S2", "lat": 1.1, "lon": 1.1, "elev": None},
        {"station": "S3", "lat": 1.2, "lon": 1.2, "elev": 300.0},
        {"station": "S4", "lat": 1.3, "lon": 1.3, "elev": 9999.0},
        {"station": "S5", "lat": 1.4, "lon": 1.4, "elev": 110.0},
    ]

    def test_no_clicks_no_update(self, web_app):
        result = self._fn(web_app)(None, self.RAW, ["smooth"], "loess", 5, 0.0, [], 3.0)
        assert result == (no_update,) * 5

    def test_no_raw_data_no_update(self, web_app):
        result = self._fn(web_app)(1, None, ["smooth"], "loess", 5, 0.0, [], 3.0)
        assert result == (no_update,) * 5

    def test_missing_elev_column_warns(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, [{"station": "S1"}], ["smooth"], "loess", 5, 0.0, [], 3.0
        )
        assert "No elevation column" in out.children

    def test_fill_op_interpolates_nan(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["fill"], "loess", 5, 0.0, [], 3.0
        )
        rows = {r["station"]: r["elev_corrected"] for r in corrected}
        assert rows["S2"] == pytest.approx(200.0)  # midpoint of neighbours

    def test_outlier_op_flags_and_interpolates(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["outlier"], "loess", 5, 0.0, [], 1.0
        )
        rows = {r["station"]: r["elev_corrected"] for r in corrected}
        assert rows["S4"] != 9999.0  # the extreme outlier was replaced

    def test_smooth_loess(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["smooth"], "loess", 3, 0.0, [], 3.0
        )
        assert "loess" in label

    def test_smooth_mean(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["smooth"], "mean", 3, 0.0, [], 3.0
        )
        assert "mean" in label

    def test_smooth_savgol(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["smooth"], "savgol", 3, 0.0, [], 3.0
        )
        assert "savgol" in label

    def test_smooth_gaussian(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["smooth"], "gaussian", 3, 0.0, [], 3.0
        )
        assert "gaussian" in label

    def test_smooth_median(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["smooth"], "median", 3, 0.0, [], 3.0
        )
        assert "median" in label

    def test_shift_op_adds_delta(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, ["shift"], "loess", 5, 12.5, [], 3.0
        )
        rows = {r["station"]: r["elev_corrected"] for r in corrected}
        assert rows["S1"] == pytest.approx(112.5)
        assert "shifted" in label

    def test_no_ops_applied(self, web_app):
        out, corrected, prog, label, style = self._fn(web_app)(
            1, self.RAW, [], "loess", 5, 0.0, [], 3.0
        )
        assert "No ops applied" in out.children[0].children[1]

    def test_fill_with_zeros_toggle(self, web_app):
        raw = [
            {"station": "S1", "lat": 1.0, "lon": 1.0, "elev": 0.0},
            {"station": "S2", "lat": 1.1, "lon": 1.1, "elev": 200.0},
        ]
        out, corrected, prog, label, style = self._fn(web_app)(
            1, raw, ["fill"], "loess", 5, 0.0, ["zeros"], 3.0
        )
        rows = {r["station"]: r["elev_corrected"] for r in corrected}
        # zero treated as missing -> filled by (limited) interpolation
        assert rows["S1"] != 0.0 or rows["S1"] == pytest.approx(200.0)


class TestElevExportCsv:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-elev-dl-csv.data", "tool-elev-export-csv")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, None, None) == (no_update, no_update)

    def test_no_data_warns(self, web_app):
        dl, status = self._fn(web_app)(1, None, None)
        assert dl is no_update
        assert "fetch or load first" in status.children

    def test_exports_corrected_csv(self, web_app):
        corrected = [
            {
                "station": "S1",
                "lat": 1.0,
                "lon": 1.0,
                "elev": 100.0,
                "elev_corrected": 101.0,
            },
        ]
        dl, status = self._fn(web_app)(1, corrected, None)
        assert dl["filename"] == "pycsamt_elevation.csv"
        assert "S1" in dl["content"]
        assert "1 rows" in status.children[1]

    def test_exports_raw_when_no_corrected(self, web_app):
        raw = [{"station": "S1", "lat": 1.0, "lon": 1.0, "elev": 100.0}]
        dl, status = self._fn(web_app)(1, None, raw)
        assert dl["filename"] == "pycsamt_elevation.csv"


class TestElevExportH5:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-elev-dl-h5.data", "tool-elev-export-h5")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, None, None) == (no_update, no_update)

    def test_no_data_warns(self, web_app):
        dl, status = self._fn(web_app)(1, None, None)
        assert dl is no_update
        assert "fetch or load first" in status.children

    def test_exports_h5_or_npz(self, web_app):
        corrected = [
            {
                "station": "S1",
                "lat": 1.0,
                "lon": 1.0,
                "elev": 100.0,
                "elev_corrected": 101.0,
            },
        ]
        dl, status = self._fn(web_app)(1, corrected, None)
        assert dl["base64"] is True
        assert dl["filename"] in (
            "pycsamt_elevation.h5",
            "pycsamt_elevation.npz",
        )
        assert "exported" in status.children[1]


class TestElevUpdateSession:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app,
            "tool-elev-export-status.children",
            "tool-elev-update-session",
        )

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, [{"station": "S1"}], "sess") is no_update

    def test_no_corrected_no_update(self, web_app):
        assert self._fn(web_app)(1, None, "sess") is no_update

    def test_session_expired_warns(self, web_app):
        out = self._fn(web_app)(
            1,
            [{"station": "S1", "elev_corrected": 100.0}],
            "no-such-session-xyz",
        )
        assert "Session expired" in out.children

    def test_updates_edi_headers_in_cached_session(
        self, web_app, willy_sites, cached_session
    ):
        first = willy_sites.as_list()[0]
        corrected = [
            {"station": first.station, "elev_corrected": 4242.0},
            {"station": "not-a-real-station", "elev_corrected": 1.0},
        ]
        out = self._fn(web_app)(1, corrected, cached_session)
        assert "station(s) updated in session" in out.children[1]
        assert "1 " in out.children[1]  # only the real station matched


class TestRestoreElevOutput:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-elev-out.children", IDs.ACTIVE_TOOL)

    def test_wrong_tool_no_update(self, web_app):
        result = self._fn(web_app)("strike", None, None)
        assert result == (no_update,) * 4

    def test_no_stored_data_no_update(self, web_app):
        result = self._fn(web_app)("elevation", None, None)
        assert result == (no_update,) * 4

    def test_restores_corrected_data(self, web_app):
        corrected = [
            {
                "station": "S1",
                "lat": 1.0,
                "lon": 1.0,
                "elev": 100.0,
                "elev_corrected": 101.0,
            },
        ]
        out, style, prog, label = self._fn(web_app)("elevation", None, corrected)
        assert style == {}
        assert prog == 100
        assert "1/1 stations" in label

    def test_restores_raw_data_when_no_corrected(self, web_app):
        raw = [{"station": "S1", "lat": 1.0, "lon": 1.0, "elev": 100.0}]
        out, style, prog, label = self._fn(web_app)("elevation", raw, None)
        assert prog == 100


# ══════════════════════════════════════════════════════════════════════════
# 5. Shared tail helpers
# ══════════════════════════════════════════════════════════════════════════


class TestStyleFig:
    def test_updates_layout_in_place(self):
        import plotly.graph_objects as go

        fig = go.Figure()
        tools_mod._style_fig(fig, "Y label", height=333)
        assert fig.layout.yaxis.title.text == "Y label"
        assert fig.layout.height == 333


class TestToolTheme:
    def test_dark_default(self):
        p = tools_mod._tool_theme(None)
        assert p["bg"] == "#1e1e2e"

    def test_dark_explicit(self):
        p = tools_mod._tool_theme("dark")
        assert p["bg"] == "#1e1e2e"

    def test_light(self):
        p = tools_mod._tool_theme("light")
        assert p["bg"] == "#ffffff"

    def test_light_case_insensitive(self):
        p = tools_mod._tool_theme("LIGHT")
        assert p["bg"] == "#ffffff"


class TestToolTableStyles:
    def test_dark(self):
        cell, header, conditional = tools_mod._tool_table_styles("dark")
        assert cell["backgroundColor"] == "#1e1e2e"
        assert conditional[0]["backgroundColor"] == "#242438"

    def test_light(self):
        cell, header, conditional = tools_mod._tool_table_styles("light")
        assert cell["backgroundColor"] == "#ffffff"
        assert conditional[0]["backgroundColor"] == "#f6f7fb"


class TestSafeStem:
    def test_strips_unsafe_characters(self):
        assert tools_mod._safe_stem("a/b:c*d?.png") == "a_b_c_d_.png"

    def test_empty_label_falls_back_to_figure(self):
        assert tools_mod._safe_stem("") == "figure"
        assert tools_mod._safe_stem(None) == "figure"

    def test_truncates_long_labels(self):
        long_label = "x" * 300
        assert len(tools_mod._safe_stem(long_label)) == 96


class TestStationCount:
    def test_real_sites_count(self, willy_sites):
        n = tools_mod._station_count(willy_sites)
        assert n == len(willy_sites.as_list())

    def test_bad_object_returns_zero(self):
        assert tools_mod._station_count(object()) == 0


class TestSaveMplFigure:
    def _fig(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        return fig

    def test_saves_new_file(self, tmp_path):
        path = tmp_path / "fig.png"
        row = tools_mod._save_mpl_figure(
            self._fig(), path, dpi=72, transparent=False, overwrite=False
        )
        assert row["Status"] == "SAVED"
        assert path.exists()

    def test_skips_existing_without_overwrite(self, tmp_path):
        path = tmp_path / "fig.png"
        path.write_bytes(b"x")
        row = tools_mod._save_mpl_figure(
            self._fig(), path, dpi=72, transparent=False, overwrite=False
        )
        assert row["Status"] == "SKIPPED"

    def test_overwrites_existing(self, tmp_path):
        path = tmp_path / "fig.png"
        path.write_bytes(b"x")
        row = tools_mod._save_mpl_figure(
            self._fig(), path, dpi=72, transparent=True, overwrite=True
        )
        assert row["Status"] == "SAVED"


class TestPlotStationMapMpl:
    def test_empty_store_shows_placeholder_text(self):
        fig = tools_mod._plot_station_map_mpl(None, None)
        ax = fig.axes[0]
        assert ax.texts[0].get_text() == "No station records"

    def test_plots_real_lat_lon_records(self, store_data_willy):
        fig = tools_mod._plot_station_map_mpl(store_data_willy, None)
        ax = fig.axes[0]
        assert ax.get_xlabel() == "Longitude"

    def test_active_lines_filter_can_empty_result(self, store_data_willy):
        fig = tools_mod._plot_station_map_mpl(
            store_data_willy, {"active": ["NoSuchLine"]}
        )
        ax = fig.axes[0]
        assert "No active station records" in ax.texts[0].get_text()

    def test_falls_back_to_easting_northing(self):
        store = {
            "station_records": [
                {"ID": "S1", "Easting": 500000.0, "Northing": 4000000.0},
                {"ID": "S2", "Easting": 500100.0, "Northing": 4000100.0},
            ]
        }
        fig = tools_mod._plot_station_map_mpl(store, None)
        ax = fig.axes[0]
        assert ax.get_xlabel() == "Easting (m)"

    def test_falls_back_to_station_index(self):
        store = {"station_records": [{"ID": "S1"}, {"ID": "S2"}]}
        fig = tools_mod._plot_station_map_mpl(store, None)
        ax = fig.axes[0]
        assert ax.get_xlabel() == "Station index"


class TestBatchExportToolFigures:
    def test_unsupported_format_raises(self, willy_sites, tmp_path):
        with pytest.raises(ValueError, match="Unsupported export format"):
            tools_mod._batch_export_tool_figures(
                willy_sites,
                dest=tmp_path,
                fmt="bmp",
                dpi=100,
                items=["map"],
                flags=set(),
                store_data={},
                active_lines_store={},
            )

    def test_pseudo_alias_expands_to_both_pseudosections(self, willy_sites, tmp_path):
        result = tools_mod._batch_export_tool_figures(
            willy_sites,
            dest=tmp_path,
            fmt="png",
            dpi=72,
            items=["pseudo"],
            flags=set(),
            store_data={},
            active_lines_store={},
        )
        assert (tmp_path / "rho_xy_pseudosection.png").exists()
        assert (tmp_path / "phase_xy_pseudosection.png").exists()
        assert result["saved"] == 2

    def test_no_manifest_flag_skips_csv(self, willy_sites, tmp_path):
        tools_mod._batch_export_tool_figures(
            willy_sites,
            dest=tmp_path,
            fmt="png",
            dpi=72,
            items=["map"],
            flags=set(),
            store_data={},
            active_lines_store={},
        )
        assert not (tmp_path / "pycsamt_batch_export_manifest.csv").exists()


class TestBatchResultView:
    def test_no_rows_shows_warning(self):
        out = tools_mod._batch_result_view(
            {
                "rows": [],
                "saved": 0,
                "skipped": 0,
                "failed": 0,
                "n_stations": 0,
                "fmt": "png",
                "dpi": 150,
                "dest": "x",
            }
        )
        assert "No figures were generated" in out.children[-1].children

    def test_with_rows_shows_table_and_manifest(self):
        result = {
            "rows": [
                {
                    "Figure": "a",
                    "Status": "SAVED",
                    "File": "a.png",
                    "Message": "",
                },
            ],
            "saved": 1,
            "skipped": 0,
            "failed": 0,
            "n_stations": 5,
            "fmt": "png",
            "dpi": 150,
            "dest": "x",
            "manifest": "x/manifest.csv",
        }
        out = tools_mod._batch_result_view(result)
        table = out.children[-1]
        assert table.data[0]["Status"] == "SAVED"


class TestRenderStoredToolResult:
    def test_no_payload_returns_none(self):
        assert (
            tools_mod._render_stored_tool_result("dimensionality", None, "dark") is None
        )
        assert (
            tools_mod._render_stored_tool_result(
                "dimensionality", {"strike": {}}, "dark"
            )
            is None
        )

    def test_non_strike_tool_returns_none(self):
        assert (
            tools_mod._render_stored_tool_result(
                "freq-editor", {"freq-editor": {"foo": "bar"}}, "dark"
            )
            is None
        )

    def test_strike_tool_delegates_to_render_strike_result(self):
        payload = {
            "records": [
                {"line": "L1", "ang_axial": 10.0},
                {"line": "L1", "ang_axial": 15.0},
            ],
            "method_label": "MB04",
            "median_iqr": 1.0,
        }
        out = tools_mod._render_stored_tool_result(
            "strike", {"strike": payload}, "dark"
        )
        assert out is not None


class TestComingSoon:
    def test_returns_placeholder_div(self):
        div = tools_mod._coming_soon("Some Tool")
        assert "not yet available in web mode" in div.children[1].children


class TestElevSummaryTable:
    def test_builds_table_and_header_counts(self):
        import pandas as pd

        df = pd.DataFrame(
            [
                {"station": "S1", "lat": 1.0, "lon": 2.0, "elev": 100.0},
                {"station": "S2", "lat": 1.1, "lon": 2.1, "elev": None},
            ]
        )
        out = tools_mod._elev_summary_table(df, "My Title")
        header, tbl = out.children
        assert "My Title: 1/2 station(s)" in header.children[1]
        assert len(tbl.data) == 2
        assert tbl.data[1]["elev"] == "—"
