# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.interp (Geological Interpretation page).

Strategy
--------
``InterpController.generate()`` never raises -- every ``plot_*`` method
either renders real content or falls back to a placeholder figure (e.g.
"No model loaded") -- so real cached WILLY sites are enough to exercise
the whole callback layer's success path without needing a full inversion
model. Only ``plot_depth_coverage``/``plot_borehole_map``-style "no model
needed" methods are used for the happy path; methods requiring
``state.model`` are exercised via their real "no model loaded" fallback
text, which is itself real production behaviour worth covering.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.interp as interp_mod
from pycsamt.app.web.cache import cache_set
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


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-interp-session"
    cache_set(session_id, willy_sites)
    return session_id


@pytest.fixture(autouse=True)
def _reset_ctrl_state():
    """InterpController is a module-level singleton -- reset between tests."""
    yield
    interp_mod._CTRL.state.model = None
    interp_mod._CTRL.state.sites = None
    interp_mod._CTRL.state.strat_logs = []
    interp_mod._CTRL.state.boreholes = []


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
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id}])
    )


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


# ── T1. switch_interp_cat ────────────────────────────────────────────────────


class TestSwitchInterpCat:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INTERP_ACTIVE_CAT}.data")

    def test_no_clicks_prevents_update(self, web_app):
        n = len(interp_mod._INTERP_CATS)
        with pytest.raises(PreventUpdate):
            _set_triggered("")
            try:
                self._fn(web_app)(*([None] * n))
            finally:
                _clear_triggered()

    def test_matches_triggered_slug(self, web_app):
        _set_triggered("interp-cat-btn-geo.n_clicks")
        try:
            n = len(interp_mod._INTERP_CATS)
            clicks = [0] * n
            clicks[1] = 1  # geo is index 1
            out = self._fn(web_app)(*clicks)
            assert out == "geo"
        finally:
            _clear_triggered()

    def test_unmatched_trigger_prevents_update(self, web_app):
        _set_triggered("something-else.n_clicks")
        try:
            n = len(interp_mod._INTERP_CATS)
            clicks = [0] * n
            clicks[0] = 1
            with pytest.raises(PreventUpdate):
                self._fn(web_app)(*clicks)
        finally:
            _clear_triggered()


# ── T3. update_interp_ctx_bar ─────────────────────────────────────────────


class TestUpdateInterpCtxBar:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.INTERP_INFO_STRIP}.children", IDs.INTERP_ACTIVE_CAT
        )

    def test_no_slug_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)

    def test_known_slug_shows_count(self, web_app):
        out = self._fn(web_app)("geo")
        texts = [
            getattr(c, "children", c) for c in out.children
        ]
        assert any("plot" in str(t) for t in texts)


# ── 1. sync_plots ─────────────────────────────────────────────────────────


class TestSyncPlots:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.INTERP_PLOT}.options")

    def test_known_category(self, web_app):
        opts, value, desc = self._fn(web_app)("geo")
        assert opts
        assert value == opts[0]["value"]
        assert desc

    def test_unknown_category_falls_back_to_first(self, web_app):
        opts, value, desc = self._fn(web_app)("does-not-exist")
        assert opts  # falls back to "Setup & Model"
        assert value is not None


# ── 2. update_desc ────────────────────────────────────────────────────────


class TestUpdateDesc:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.INTERP_DESC}.children", IDs.INTERP_PLOT
        )

    def test_no_fn_name_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "geo")

    def test_known_fn_returns_desc(self, web_app):
        desc = self._fn(web_app)("plot_strat_log", "geo")
        assert "Lithology" in desc

    def test_unknown_fn_no_update(self, web_app):
        from dash import no_update

        out = self._fn(web_app)("not-a-real-fn", "geo")
        assert out is no_update


# ── 3. sync_prereq_card ───────────────────────────────────────────────────


class TestSyncPrereqCard:
    def _fn(self, web_app):
        return _cb_multi(web_app, "interp-prereq-card")

    def test_category_without_prereq(self, web_app):
        style, hint, label = self._fn(web_app)("setup")
        assert style == {"display": "none"}
        assert label == "Run"

    def test_category_with_runnable_prereq(self, web_app):
        style, hint, label = self._fn(web_app)("geo")
        assert style == {"display": "block"}
        assert label == "Run Geological Classification"

    def test_category_with_unavailable_prereq(self, web_app):
        style, hint, label = self._fn(web_app)("field")
        assert style == {"display": "block"}
        assert label == "Not available"


# ── run_prereq ────────────────────────────────────────────────────────────


class TestRunPrereq:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INTERP_PREREQ_OUT}.children")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "geo", "s", None, "dark")

    def test_category_without_prereq_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "setup", "s", None, "dark")

    def test_unavailable_prereq_returns_hint_text(self, web_app):
        out = self._fn(web_app)(1, "field", "s", None, "dark")
        assert "borehole" in out.lower()

    def test_no_data_loaded(self, web_app):
        out = self._fn(web_app)(1, "geo", "no-such-session", None, "dark")
        assert "load a survey first" in out.lower()

    def test_real_run_no_model_loaded(self, web_app, cached_session):
        out = self._fn(web_app)(1, "geo", cached_session, None, "dark")
        assert out == "No model loaded."

    def test_exception_reported(self, web_app, cached_session, monkeypatch):
        def _boom():
            raise RuntimeError("prereq boom")

        monkeypatch.setattr(interp_mod._CTRL, "run_geological", _boom)
        out = self._fn(web_app)(1, "geo", cached_session, None, "dark")
        assert "prereq boom" in out


# ── 4. update_src_hint ────────────────────────────────────────────────────


class TestUpdateSrcHint:
    def _fn(self, web_app):
        return _cb(web_app, "interp-src-hint.children")

    def test_raw_source(self, web_app):
        out = self._fn(web_app)("raw", None)
        assert "raw loaded data" in out

    def test_corrected_no_steps(self, web_app):
        out = self._fn(web_app)("corrected", None)
        assert "No corrections applied yet" in out

    def test_corrected_with_steps(self, web_app):
        out = self._fn(web_app)("corrected", {"steps": [{"fn_name": "x"}]})
        assert "1 step" in out


# ── _build_figure / _get_corrected_sites (direct) ────────────────────────


class TestBuildFigureHelpers:
    def test_get_corrected_sites_no_cache(self):
        out = interp_mod._get_corrected_sites("no-such-session", None, True)
        assert out is None

    def test_get_corrected_sites_no_steps_returns_raw(
        self, cached_session, willy_sites, monkeypatch
    ):
        monkeypatch.setattr(interp_mod, "cache_get", lambda _key: willy_sites)
        out = interp_mod._get_corrected_sites(cached_session, None, True)
        assert out is willy_sites

    def test_get_corrected_sites_applies_steps(self, cached_session):
        steps = [
            {
                "fn_name": "correct_ss_ama",
                "kwargs": {
                    "sort_by": "lon",
                    "half_window": 3,
                    "weights": "tri",
                    "max_skew": 6.0,
                },
                "label": "AMA",
            }
        ]
        out = interp_mod._get_corrected_sites(
            cached_session, {"steps": steps}, True
        )
        assert out is not None

    def test_get_corrected_sites_bad_step_swallowed(
        self, cached_session, willy_sites, monkeypatch
    ):
        monkeypatch.setattr(interp_mod, "cache_get", lambda _key: willy_sites)
        steps = [{"fn_name": "not-a-real-fn", "kwargs": {}, "label": "bad"}]
        out = interp_mod._get_corrected_sites(
            cached_session, {"steps": steps}, True
        )
        assert out is willy_sites

    def test_build_figure_no_data_raises_value_error(self):
        with pytest.raises(ValueError, match="no-data"):
            interp_mod._build_figure(
                "plot_depth_coverage", "setup", "raw", "auto", None,
                None, None, None, "no-such-session", None, None, None,
                "dark",
            )

    def test_build_figure_unknown_method_raises(self, cached_session):
        with pytest.raises(ValueError, match="unknown-method"):
            interp_mod._build_figure(
                "not_a_real_method", "setup", "raw", "auto", None,
                None, None, None, cached_session, None, None, None,
                "dark",
            )

    def test_build_figure_success_real_render(self, cached_session):
        fig, label, cat = interp_mod._build_figure(
            "plot_depth_coverage", "setup", "raw", "auto", None,
            None, None, None, cached_session, None, None, None,
            "dark",
        )
        assert fig is not None
        assert label == "Depth coverage"
        assert cat == "Setup & Model"

    def test_build_figure_applies_line_filter(self, cached_session):
        # Active-lines filter with an empty station_records list can't
        # resolve any ids -> filter_sites_by_lines leaves sites unchanged.
        store_data = {"station_records": []}
        fig, label, cat = interp_mod._build_figure(
            "plot_depth_coverage", "setup", "raw", "auto", None,
            None, None, None, cached_session, {"active": ["L1"]}, store_data,
            None, "dark",
        )
        assert fig is not None

    def test_build_figure_custom_figsize(self, cached_session):
        fig, _label, _cat = interp_mod._build_figure(
            "plot_depth_coverage", "setup", "raw", "wide", None,
            None, None, None, cached_session, None, None, None,
            "dark",
        )
        assert tuple(fig.get_size_inches()) == interp_mod._FIGSIZES["wide"]

    def test_build_figure_no_figure_raises(self, cached_session, monkeypatch):
        monkeypatch.setattr(
            interp_mod._CTRL, "generate", lambda *a, **k: None
        )
        with pytest.raises(ValueError, match="no-figure"):
            interp_mod._build_figure(
                "plot_depth_coverage", "setup", "raw", "auto", None,
                None, None, None, cached_session, None, None, None,
                "dark",
            )

    def test_build_figure_corrected_source(self, cached_session):
        steps = [
            {
                "fn_name": "correct_ss_ama",
                "kwargs": {
                    "sort_by": "lon",
                    "half_window": 3,
                    "weights": "tri",
                    "max_skew": 6.0,
                },
                "label": "AMA",
            }
        ]
        fig, _label, _cat = interp_mod._build_figure(
            "plot_depth_coverage", "setup", "corrected", "auto", None,
            None, None, None, cached_session, None, None, {"steps": steps},
            "dark",
        )
        assert fig is not None

    def test_build_figure_extra_kwargs_forwarded(self, cached_session):
        fig, _label, _cat = interp_mod._build_figure(
            "plot_depth_coverage", "setup", "raw", "auto", "viridis",
            500, 0.1, 100.0, cached_session, None, None, None, "dark",
        )
        assert fig is not None


# ── 5. run_interp ─────────────────────────────────────────────────────────


class TestRunInterp:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "img-interp-setup.src", IDs.BTN_INTERP_RUN
        )

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                None, "plot_depth_coverage", "setup", "raw", "auto", None,
                None, None, None, "s", None, None, None, "dark",
            )

    def test_no_fn_name_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1, None, "setup", "raw", "auto", None, None, None, None,
                "s", None, None, None, "dark",
            )

    def test_no_data_loaded(self, web_app):
        out = self._fn(web_app)(
            1, "plot_depth_coverage", "setup", "raw", "auto", None, None,
            None, None, "no-such-session", None, None, None, "dark",
        )
        assert "No data loaded" in out[10].children

    def test_success_path(self, web_app, cached_session):
        out = self._fn(web_app)(
            1, "plot_depth_coverage", "setup", "raw", "auto", None, None,
            None, None, cached_session, None, None, None, "dark",
        )
        imgs = out[:9]
        assert imgs[0] is not None and imgs[0] != ""
        assert out[9]  # interp-last-src.data
        assert out[12] is False  # toast not open

    def test_unknown_method_reports_error(self, web_app, cached_session):
        out = self._fn(web_app)(
            1, "not_a_real_method", "setup", "raw", "auto", None, None,
            None, None, cached_session, None, None, None, "dark",
        )
        assert out[12] is True
        assert "Unknown method" in out[13]

    def test_no_figure_reports_error(self, web_app, cached_session, monkeypatch):
        monkeypatch.setattr(
            interp_mod._CTRL, "generate", lambda *a, **k: None
        )
        out = self._fn(web_app)(
            1, "plot_depth_coverage", "setup", "raw", "auto", None, None,
            None, None, cached_session, None, None, None, "dark",
        )
        assert out[12] is True
        assert "no figure" in out[13].lower()

    def test_generic_exception_reports_error(
        self, web_app, cached_session, monkeypatch
    ):
        def _boom(*a, **k):
            raise RuntimeError("set_sites boom")

        monkeypatch.setattr(interp_mod._CTRL, "set_sites", _boom)
        out = self._fn(web_app)(
            1, "plot_depth_coverage", "setup", "raw", "auto", None, None,
            None, None, cached_session, None, None, None, "dark",
        )
        assert out[12] is True
        assert "set_sites boom" in out[13]

    def test_corrected_source_with_steps_tags_corrected(
        self, web_app, cached_session
    ):
        steps = [
            {
                "fn_name": "correct_ss_ama",
                "kwargs": {
                    "sort_by": "lon",
                    "half_window": 3,
                    "weights": "tri",
                    "max_skew": 6.0,
                },
                "label": "AMA",
            }
        ]
        out = self._fn(web_app)(
            1, "plot_depth_coverage", "setup", "corrected", "auto", None,
            None, None, None, cached_session, None, None,
            {"steps": steps}, "dark",
        )
        strip_texts = [
            getattr(c, "children", "") for c in out[10].children
        ]
        assert any("corrected" in str(t) for t in strip_texts)


# ── 6. export_plot ────────────────────────────────────────────────────────


class TestExportPlot:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INTERP_DOWNLOAD}.data")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                None, "png", None, "plot_depth_coverage", "setup", "raw",
                "auto", None, None, None, None, "s", None, None, None,
                "dark",
            )

    def test_no_fn_name_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1, "png", None, None, "setup", "raw", "auto", None, None,
                None, None, "s", None, None, None, "dark",
            )

    def test_png_fast_path_uses_last_src(self, web_app):
        last_src = "data:image/png;base64,aGVsbG8="
        out = self._fn(web_app)(
            1, "png", last_src, "plot_depth_coverage", "setup", "raw",
            "auto", None, None, None, None, "s", None, None, None, "dark",
        )
        assert out["filename"] == "plot_depth_coverage.png"

    def test_png_fast_path_falls_through_on_bad_src(
        self, web_app, cached_session
    ):
        out = self._fn(web_app)(
            1, "png", "not-a-valid-data-uri", "plot_depth_coverage",
            "setup", "raw", "auto", None, None, None, None,
            cached_session, None, None, None, "dark",
        )
        assert out["filename"] == "plot_depth_coverage.png"

    def test_svg_regenerates_figure(self, web_app, cached_session):
        out = self._fn(web_app)(
            1, "svg", None, "plot_depth_coverage", "setup", "raw", "auto",
            None, None, None, None, cached_session, None, None, None,
            "dark",
        )
        assert out["filename"] == "plot_depth_coverage.svg"

    def test_regenerate_exception_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1, "svg", None, "plot_depth_coverage", "setup", "raw",
                "auto", None, None, None, None, "no-such-session", None,
                None, None, "dark",
            )
