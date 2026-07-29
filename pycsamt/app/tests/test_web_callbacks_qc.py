# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.qc (QC page v2 callbacks).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ -- real EDI-derived Sites, loaded once via
``pycsamt.emtools.ensure_sites`` and cached with
``pycsamt.app.web.cache.cache_set`` so ``run_qc``/``run_qc_overview``
exercise the real ``QCController.draw()`` -> real ``pycsamt.emtools``
plot functions -> real ``fig_to_src`` pipeline, rather than mocks.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from dash import no_update

import pycsamt.app.web.callbacks.qc as qc_mod
from pycsamt.app.desktop.controllers.qc_controller import ALL_GROUPS
from pycsamt.app.web.cache import cache_set
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]  # pycsamt/
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
    """station_records using the *real* STORE_DATA schema (capitalized
    "ID"/"Line" -- see pycsamt.app.web.callbacks.data)."""
    recs = []
    for i, ed in enumerate(willy_sites.as_list()):
        recs.append(
            {
                "ID": ed.station,
                "Line": "L1" if i % 2 == 0 else "L2",
            }
        )
    return recs


@pytest.fixture
def store_data_willy(willy_records):
    return {"station_records": willy_records}


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-qc-session"
    cache_set(session_id, willy_sites)
    return session_id


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


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(AttributeDict(triggered_inputs=[{"prop_id": prop_id}]))


def _set_no_trigger():
    """Valid (in-callback) context with nothing triggered -- as opposed
    to no context at all, which dash 4.x raises on ``ctx.triggered_id``."""
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(AttributeDict(triggered_inputs=[]))


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


# ── 1. _plot_options ─────────────────────────────────────────────────────────


class TestPlotOptions:
    def test_known_group_returns_catalogue_options(self):
        opts = qc_mod._plot_options("Coverage")
        labels = {o["label"] for o in opts}
        values = {o["value"] for o in opts}
        assert "Coverage (per-site presence)" in labels
        assert "plot_coverage" in values
        assert len(opts) == 3  # COVERAGE_PLOTS has 3 entries

    def test_unknown_group_returns_empty(self):
        assert qc_mod._plot_options("Nope") == []

    def test_all_groups_produce_matching_option_counts(self):
        for gname, plots in ALL_GROUPS:
            opts = qc_mod._plot_options(gname)
            assert len(opts) == len(plots)


# ── 2. _has_ax ────────────────────────────────────────────────────────────────


class TestHasAx:
    def test_single_axes_fn_is_true(self):
        assert qc_mod._has_ax("plot_coverage") is True

    def test_multi_axes_fn_is_false(self):
        assert qc_mod._has_ax("plot_qc_quicklook") is False

    def test_unknown_fn_defaults_true(self):
        assert qc_mod._has_ax("does_not_exist") is True


# ── 3. set_active_group ──────────────────────────────────────────────────────


class TestSetActiveGroup:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.QC_GROUP_SEL}.data")

    def test_no_trigger_no_update(self, web_app):
        _set_no_trigger()
        try:
            out = self._fn(web_app)([1], "Overview")
        finally:
            _clear_triggered()
        assert out == (no_update, no_update)

    def test_non_dict_trigger_no_update(self, web_app):
        _set_triggered(f"{IDs.BTN_QC_RUN}.n_clicks")
        try:
            out = self._fn(web_app)([1], "Overview")
        finally:
            _clear_triggered()
        assert out == (no_update, no_update)

    def test_pill_click_sets_group_and_classes(self, web_app):
        _set_triggered('{"index":"Coverage","type":"qc-grp-btn"}.n_clicks')
        try:
            new_group, classes = self._fn(web_app)([1, 0], "Overview")
        finally:
            _clear_triggered()
        assert new_group == "Coverage"
        names = [g for g, _ in ALL_GROUPS]
        assert classes[names.index("Coverage")] == "qc-grp-btn active"
        assert classes[names.index("Overview")] == "qc-grp-btn"


# ── 4. sync_qc_plots ─────────────────────────────────────────────────────────


class TestSyncQcPlots:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.QC_PLOT_SELECT}.options")

    def test_none_defaults_to_overview(self, web_app):
        opts, val = self._fn(web_app)(None)
        assert val == "plot_qc_quicklook"
        assert len(opts) == 5

    def test_known_group(self, web_app):
        opts, val = self._fn(web_app)("Static Shift")
        assert val == "ss_qc_psection"
        assert len(opts) == 4

    def test_unknown_group_empty_value_none(self, web_app):
        opts, val = self._fn(web_app)("Bogus")
        assert opts == []
        assert val is None


# ── 5. populate_qc_lines ─────────────────────────────────────────────────────


class TestPopulateQcLines:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.QC_LINE_SEL}.options")

    def test_no_store_data(self, web_app):
        assert self._fn(web_app)(None) == []

    def test_lowercase_line_key_populates(self, web_app):
        """``populate_qc_lines`` reads ``r.get("line")`` (lowercase)."""
        store = {
            "station_records": [
                {"line": "L2"},
                {"line": "L1"},
                {"line": ""},
                {},
            ]
        }
        opts = self._fn(web_app)(store)
        assert opts == [
            {"label": "L1", "value": "L1"},
            {"label": "L2", "value": "L2"},
        ]

    def test_real_store_data_uses_capitalized_line_key_bug(
        self, web_app, store_data_willy
    ):
        """
        Real bug: ``populate_qc_lines`` looks up ``r.get("line", "")``
        (lowercase), but the real STORE_DATA schema built by
        ``pycsamt.app.web.callbacks.data`` always uses a capitalized
        ``"Line"`` key (see e.g. ``data.py`` -- ``r.get("Line", "")``,
        ``by_id[r["ID"]]``). Because of this key-casing mismatch, the
        QC page's "Lines" multi-dropdown is *always empty* when fed
        real survey data, even though the survey clearly has lines
        (``store_data_willy`` here has real "L1"/"L2" line assignments).
        Not fixed here, per instructions -- documented via this test.
        """
        assert "Line" in store_data_willy["station_records"][0]
        opts = self._fn(web_app)(store_data_willy)
        assert opts == []  # should list L1/L2, but the key mismatch empties it


# ── 6. populate_qc_stns ──────────────────────────────────────────────────────


class TestPopulateQcStns:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.QC_STN_SEL}.options")

    def test_no_store_data(self, web_app):
        assert self._fn(web_app)(None, None) == []

    def test_lowercase_station_key_all_lines(self, web_app):
        store = {
            "station_records": [
                {"line": "L1", "station": "B1"},
                {"line": "L1", "station": "A1"},
                {"line": "L2", "station": "C1"},
            ]
        }
        opts = self._fn(web_app)(None, store)
        assert opts == [
            {"label": "A1", "value": "A1"},
            {"label": "B1", "value": "B1"},
            {"label": "C1", "value": "C1"},
        ]

    def test_filters_by_selected_lines(self, web_app):
        store = {
            "station_records": [
                {"line": "L1", "station": "A1"},
                {"line": "L2", "station": "C1"},
            ]
        }
        opts = self._fn(web_app)(["L2"], store)
        assert opts == [{"label": "C1", "value": "C1"}]

    def test_real_store_data_also_hits_the_same_key_bug(
        self, web_app, store_data_willy
    ):
        """Same lowercase/capitalized key mismatch as populate_qc_lines --
        the Stations dropdown is empty against real STORE_DATA too."""
        opts = self._fn(web_app)(None, store_data_willy)
        assert opts == []


# ── 7. run_qc ─────────────────────────────────────────────────────────────────


class TestRunQc:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.IMG_QC_PLOT}.src")

    def _args(self, **overrides):
        args = dict(
            n_run=None,
            fn_name="plot_coverage",
            nav_section="qc",
            store_data=None,
            sel_lines=None,
            sel_stns=None,
            comp=None,
            metric=None,
            method=None,
            thresh=None,
            n_max=None,
            mains_hz=None,
            figsize_key="compact",
            auto_toggle=["auto"],
            theme="dark",
            session_id=None,
            active_lines_store=None,
        )
        args.update(overrides)
        return args

    def _call(self, web_app, trigger, **overrides):
        fn = self._fn(web_app)
        _set_triggered(trigger)
        try:
            return fn(**self._args(**overrides))
        finally:
            _clear_triggered()

    def test_not_on_qc_page_and_not_manual_no_update(self, web_app):
        out = self._call(
            web_app,
            f"{IDs.STORE_DATA}.data",
            nav_section="home",
        )
        assert out == (no_update,) * 5

    def test_on_qc_page_auto_off_no_update(self, web_app):
        out = self._call(
            web_app,
            f"{IDs.QC_PLOT_SELECT}.value",
            nav_section="qc",
            auto_toggle=[],
        )
        assert out == (no_update,) * 5

    def test_manual_run_missing_fn_name(self, web_app):
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name=None,
            session_id="sid",
        )
        assert out[2] == "Select a plot first."
        assert out[3] is True

    def test_auto_trigger_missing_session_no_update(self, web_app):
        out = self._call(
            web_app,
            f"{IDs.STORE_DATA}.data",
            nav_section="qc",
            session_id=None,
        )
        assert out == (no_update,) * 5

    def test_manual_run_no_cached_sites(self, web_app):
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            session_id="no-such-session-xyz",
        )
        assert out[2] == "No data loaded."
        assert out[3] is True
        assert out[0].startswith("data:image")

    def test_auto_trigger_no_cached_sites_no_update(self, web_app):
        out = self._call(
            web_app,
            f"{IDs.STORE_DATA}.data",
            nav_section="qc",
            session_id="no-such-session-xyz",
        )
        assert out == (no_update,) * 5

    def test_manual_run_real_sites_has_ax_plot_dark(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name="plot_coverage",
            session_id=cached_session,
            store_data=store_data_willy,
            theme="dark",
        )
        src, tab, feedback, is_err, toast = out
        assert src.startswith("data:image")
        assert tab == "plot"
        assert feedback == ""
        assert is_err is False

    def test_manual_run_real_sites_light_theme_multipanel(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name="plot_qc_quicklook",
            session_id=cached_session,
            store_data=store_data_willy,
            theme="light",
            figsize_key="wide",
        )
        src = out[0]
        assert src.startswith("data:image")

    def test_manual_run_polar_dispatcher(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name="plot_polar_coverage",
            session_id=cached_session,
            store_data=store_data_willy,
        )
        assert out[0].startswith("data:image")

    def test_active_lines_filter_from_state_store(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._call(
            web_app,
            f"{IDs.QC_PLOT_SELECT}.value",
            nav_section="qc",
            fn_name="plot_coverage",
            session_id=cached_session,
            store_data=store_data_willy,
            active_lines_store={"active": ["L1"], "all": ["L1", "L2"]},
        )
        assert out[0].startswith("data:image")

    def test_active_lines_filter_no_match_shows_warning(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name="plot_coverage",
            session_id=cached_session,
            store_data=store_data_willy,
            sel_lines=["L999-does-not-exist"],
        )
        src, tab, feedback, is_err, toast = out
        assert src.startswith("data:image")
        assert feedback == "No matching lines."
        assert is_err is False

    def test_station_filter_narrows_render(
        self, web_app, cached_session, store_data_willy
    ):
        first_id = store_data_willy["station_records"][0]["ID"]
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name="plot_coverage",
            session_id=cached_session,
            store_data=store_data_willy,
            sel_stns=[first_id],
        )
        assert out[0].startswith("data:image")

    def test_render_exception_is_caught(
        self, monkeypatch, web_app, cached_session, store_data_willy
    ):
        def _boom(*a, **k):
            raise RuntimeError("boom")

        monkeypatch.setattr(qc_mod, "fig_to_src", _boom)
        out = self._call(
            web_app,
            f"{IDs.BTN_QC_RUN}.n_clicks",
            n_run=1,
            fn_name="plot_coverage",
            session_id=cached_session,
            store_data=store_data_willy,
        )
        src, tab, feedback, is_err, toast = out
        assert is_err is True
        assert "plot_coverage" in feedback
        assert "QC error: boom" in toast


# ── 8. run_qc_overview ───────────────────────────────────────────────────────


class TestRunQcOverview:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.IMG_QC_OVERVIEW}.src")

    def test_not_qc_page_no_update(self, web_app):
        out = self._fn(web_app)("home", None, "dark", "sid", None)
        assert out is no_update

    def test_no_cached_sites_no_update(self, web_app):
        out = self._fn(web_app)("qc", None, "dark", "missing-session", None)
        assert out is no_update

    def test_real_sites_dark(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)("qc", store_data_willy, "dark", cached_session, None)
        assert out.startswith("data:image")

    def test_real_sites_light_with_active_lines(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            "qc",
            store_data_willy,
            "light",
            cached_session,
            {"active": ["L1", "L2"], "all": ["L1", "L2"]},
        )
        assert out.startswith("data:image")

    def test_active_lines_no_match_no_update(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            "qc",
            store_data_willy,
            "dark",
            cached_session,
            {"active": ["nope"], "all": ["L1", "L2"]},
        )
        assert out is no_update

    def test_render_exception_no_update(
        self, monkeypatch, web_app, cached_session, store_data_willy
    ):
        def _boom(*a, **k):
            raise RuntimeError("boom")

        monkeypatch.setattr(qc_mod, "fig_to_src", _boom)
        out = self._fn(web_app)("qc", store_data_willy, "dark", cached_session, None)
        assert out is no_update
