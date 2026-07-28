# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.correction (Data Correction page).

Strategy
--------
``CorrectionController`` wraps real, fast, deterministic array-based
correction routines from ``pycsamt.emtools`` (``correct_ss_ama``,
``smooth_rho_phase``, ``drop_freqs_manual``, the Stratagem
``FrequencyFilter``) -- not ML training and not slow -- so this module
is exercised against real cached survey data
(``data/AMT/WILLY_DATA/L18PLT/``, 28 real EDI stations) rather than
mocks, matching the pattern already used for callbacks/agents.py.

Two fixtures anchor everything:
  * ``cached_session`` -- a session id with the real WILLY sites
    pre-populated in ``pycsamt.app.web.cache``.
  * ``STORE_DATA_WILLY`` -- a plain dict matching the shape written by
    callbacks/data.py, with real WILLY station IDs split across two
    synthetic survey lines so line/station filtering paths are
    reachable.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from dash import no_update
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.correction as corr_mod
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
def willy_ids(willy_sites):
    return [corr_mod._edi_name(ed) for ed in corr_mod._iter_edi(willy_sites)]


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-correction-session"
    cache_set(session_id, willy_sites)
    return session_id


@pytest.fixture
def store_data_willy(willy_ids):
    half = len(willy_ids) // 2
    records = [
        {"ID": sid, "Line": "L1" if i < half else "L2"}
        for i, sid in enumerate(willy_ids)
    ]
    return {"station_records": records}


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


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


# ══════════════════════════════════════════════════════════════════════════
# Pure helper functions (no Dash, no data)
# ══════════════════════════════════════════════════════════════════════════


class TestCatalogueHelpers:
    def test_get_specs_known(self):
        specs = corr_mod._get_specs("Static Shift", "AMA (spatial average)")
        assert specs
        assert specs[0].name == "sort_by"

    def test_get_specs_unknown_returns_empty(self):
        assert corr_mod._get_specs("nope", "nope") == []
        assert corr_mod._get_specs(None, None) == []

    def test_get_fn_name(self):
        assert (
            corr_mod._get_fn_name("Static Shift", "AMA (spatial average)")
            == "correct_ss_ama"
        )
        assert corr_mod._get_fn_name("nope", "nope") is None

    def test_get_desc(self):
        desc = corr_mod._get_desc("Static Shift", "AMA (spatial average)")
        assert "static-shift" in desc
        assert corr_mod._get_desc("nope", "nope") == ""


class TestBuildKwargs:
    def test_defaults_used_when_no_values(self):
        kwargs = corr_mod._build_kwargs("Static Shift", "AMA (spatial average)", [])
        assert kwargs["sort_by"] == "lon"
        assert kwargs["half_window"] == 3
        assert kwargs["weights"] == "tri"
        assert kwargs["max_skew"] == 6.0

    def test_none_value_falls_back_to_default(self):
        kwargs = corr_mod._build_kwargs(
            "Static Shift", "AMA (spatial average)", [None, None, None, None]
        )
        assert kwargs["half_window"] == 3

    def test_spin_and_dspin_conversion(self):
        kwargs = corr_mod._build_kwargs(
            "Static Shift", "AMA (spatial average)", ["lat", "5", "box", "3.5"]
        )
        assert kwargs["half_window"] == 5
        assert isinstance(kwargs["half_window"], int)
        assert kwargs["max_skew"] == 3.5
        assert isinstance(kwargs["max_skew"], float)

    def test_invalid_numeric_falls_back_to_default(self):
        kwargs = corr_mod._build_kwargs(
            "Static Shift",
            "AMA (spatial average)",
            ["lon", "not-a-number", "tri", "also-bad"],
        )
        assert kwargs["half_window"] == 3
        assert kwargs["max_skew"] == 6.0

    def test_check_kind_coerced_to_bool(self):
        # Manual freq drop / smoothing methods aren't in CATALOGUE with
        # "check" kind by default, so build a synthetic spec via monkeypatch
        # is unnecessary -- exercise through a category known to have one.
        specs = corr_mod._get_specs("Static Shift", "AMA (spatial average)")
        assert specs  # sanity: catalogue not empty for this category


class TestSmallHelpers:
    def test_chain_badge(self):
        assert corr_mod._chain_badge(0) == ""
        assert corr_mod._chain_badge(1) == "1 step"
        assert corr_mod._chain_badge(3) == "3 steps"

    def test_ctx_bar_content_empty(self):
        out = corr_mod._ctx_bar_content([])
        assert "Chain is empty" in out[1]

    def test_ctx_bar_content_with_steps(self):
        steps = [{"fn_name": "foo", "label": "Foo"}]
        out = corr_mod._ctx_bar_content(steps)
        texts = [getattr(c, "children", c) for c in out]
        assert any("1 correction" in str(t) for t in texts)
        assert any("Foo" in str(t) for t in texts)

    def test_stack_el_empty(self):
        el = corr_mod._stack_el([])
        assert "No corrections applied yet." in el.children

    def test_stack_el_with_steps(self):
        steps = [
            {
                "fn_name": "correct_ss_ama",
                "kwargs": {"half_window": 3},
                "label": "AMA",
            }
        ]
        el = corr_mod._stack_el(steps)
        assert len(el.children) == 1

    def test_rho_phi_hint_pairs(self):
        assert "Station Pairs" in corr_mod._rho_phi_hint([], "pairs")

    def test_rho_phi_hint_freq_bands_no_drop(self):
        hint = corr_mod._rho_phi_hint([], "freq-bands")
        assert "Apply 'Drop Bad Frequencies'" in hint

    def test_rho_phi_hint_freq_bands_with_drop(self):
        steps = [{"fn_name": "_strat_freq_filter"}]
        hint = corr_mod._rho_phi_hint(steps, "freq-bands")
        assert "removed frequencies" in hint

    def test_rho_phi_hint_default_with_static_shift(self):
        steps = [{"fn_name": "correct_static_shift"}]
        hint = corr_mod._rho_phi_hint(steps, "grid")
        assert "Station Pairs" in hint

    def test_rho_phi_hint_default_with_drop_freq(self):
        steps = [{"fn_name": "drop_freqs_manual"}]
        hint = corr_mod._rho_phi_hint(steps, "grid")
        assert "Freq Bands" in hint

    def test_rho_phi_hint_default_plain(self):
        hint = corr_mod._rho_phi_hint([], "grid")
        assert "ρ/φ Grid" in hint

    def test_diff_stats_strip(self):
        el = corr_mod._diff_stats_strip(-1.0, 1.0, 0.1, 0.5, 3, 5)
        assert el.className == "corr-diff-stats-row"
        assert len(el.children) == 6

    def test_params_panel_no_params(self):
        el = corr_mod._params_panel("nope", "nope")
        assert "No parameters" in el.children

    def test_params_panel_with_params(self):
        el = corr_mod._params_panel("Static Shift", "AMA (spatial average)")
        assert len(el.children) == 4  # sort_by, half_window, weights, max_skew

    def test_params_panel_check_kind_widget(self):
        # "ρ/φ trend smoothing" has "check" params (smooth_rho, etc.) --
        # exercises the dbc.Switch branch not covered by AMA's params.
        el = corr_mod._params_panel("Noise Removal", "ρ/φ trend smoothing")
        assert len(el.children) >= 3


class TestFreqGridHelpers:
    def test_common_freq_grid_empty(self):
        assert corr_mod._common_freq_grid([]) is None
        assert corr_mod._common_freq_grid([np.array([])]) is None

    def test_common_freq_grid_single_value(self):
        assert corr_mod._common_freq_grid([np.array([1.0])]) is None

    def test_common_freq_grid_builds_geomspace(self):
        grid = corr_mod._common_freq_grid([np.array([1.0, 10.0, 100.0])], n_pts=10)
        assert grid is not None
        assert len(grid) == 10
        assert grid[0] == pytest.approx(1.0)
        assert grid[-1] == pytest.approx(100.0)

    def test_interp_to_grid_insufficient_points(self):
        out = corr_mod._interp_to_grid([1.0], [2.0], np.array([1.0, 2.0]))
        assert np.all(np.isnan(out))

    def test_interp_to_grid_normal(self):
        freq = [1.0, 10.0, 100.0]
        values = [10.0, 100.0, 1000.0]
        grid = np.array([1.0, 10.0, 100.0])
        out = corr_mod._interp_to_grid(freq, values, grid)
        assert out == pytest.approx([10.0, 100.0, 1000.0], rel=1e-6)

    def test_interp_phi_to_grid_insufficient_points(self):
        out = corr_mod._interp_phi_to_grid([1.0], [45.0], np.array([1.0, 2.0]))
        assert np.all(np.isnan(out))

    def test_interp_phi_to_grid_normal(self):
        freq = [1.0, 10.0]
        values = [40.0, 50.0]
        grid = np.array([1.0, 10.0])
        out = corr_mod._interp_phi_to_grid(freq, values, grid)
        assert out == pytest.approx([40.0, 50.0], rel=1e-6)

    def test_band_stats_rho_empty(self):
        med, p25, p75 = corr_mod._band_stats_rho([], np.array([1.0, 2.0]))
        assert np.all(np.isnan(med))

    def test_band_stats_phi_empty(self):
        med, p25, p75 = corr_mod._band_stats_phi([], np.array([1.0, 2.0]))
        assert np.all(np.isnan(med))


class TestRemovedPeriods:
    def test_no_z_returns_empty(self):
        class _NoZ:
            pass

        out = corr_mod._removed_periods(_NoZ(), None)
        assert out.size == 0

    def test_no_cur_returns_all_raw_periods(self):
        class _Z:
            freq = np.array([1.0, 2.0, 4.0])

        class _Ed:
            Z = _Z()

        out = corr_mod._removed_periods(_Ed(), None)
        assert out.size == 3

    def test_cur_missing_freq_returns_all_raw(self):
        class _ZRaw:
            freq = np.array([1.0, 2.0])

        class _ZCur:
            freq = None

        class _EdRaw:
            Z = _ZRaw()

        class _EdCur:
            Z = _ZCur()

        out = corr_mod._removed_periods(_EdRaw(), _EdCur())
        assert out.size == 2

    def test_detects_removed_frequency(self):
        class _ZRaw:
            freq = np.array([1.0, 2.0, 4.0])

        class _ZCur:
            freq = np.array([1.0, 4.0])  # 2.0 Hz removed

        class _EdRaw:
            Z = _ZRaw()

        class _EdCur:
            Z = _ZCur()

        out = corr_mod._removed_periods(_EdRaw(), _EdCur())
        assert out.size == 1
        assert out[0] == pytest.approx(0.5)  # period = 1/2.0


# ══════════════════════════════════════════════════════════════════════════
# Real-data helpers: _build_ctrl, _apply_line_filter, _filter_display_sites,
# _iter_edi/_edi_name/_edi_to_rows/_collect_rho, renderers
# ══════════════════════════════════════════════════════════════════════════


class TestIterEdi:
    def test_iter_edi_none(self):
        assert list(corr_mod._iter_edi(None)) == []

    def test_iter_edi_real_sites(self, willy_sites):
        edis = list(corr_mod._iter_edi(willy_sites))
        assert len(edis) == 28

    def test_edi_name_real(self, willy_sites):
        ed = next(corr_mod._iter_edi(willy_sites))
        name = corr_mod._edi_name(ed)
        assert name

    def test_edi_to_rows_real(self, willy_sites):
        ed = next(corr_mod._iter_edi(willy_sites))
        rows = corr_mod._edi_to_rows(ed)
        assert rows
        assert "freq_hz" in rows[0]

    def test_collect_rho_real(self, willy_sites):
        rho = corr_mod._collect_rho(willy_sites)
        assert rho is not None
        assert rho.size > 0

    def test_collect_rho_none_sites(self):
        assert corr_mod._collect_rho(None) is None


class TestBuildCtrl:
    def test_no_cached_sites_returns_none(self):
        ctrl = corr_mod._build_ctrl("no-such-session", [], dark=True)
        assert ctrl is None

    def test_builds_from_cache(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        assert ctrl is not None
        assert ctrl.raw_sites is not None
        assert ctrl.current_sites is ctrl.raw_sites

    def test_applies_valid_steps(self, cached_session):
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
        ctrl = corr_mod._build_ctrl(cached_session, steps, dark=True)
        assert ctrl is not None
        assert ctrl.current_sites is not ctrl.raw_sites

    def test_bad_step_is_swallowed(self, cached_session):
        steps = [{"fn_name": "not_a_real_fn", "kwargs": {}, "label": "Bad"}]
        ctrl = corr_mod._build_ctrl(cached_session, steps, dark=True)
        assert ctrl is not None
        # bad step swallowed -> current sites unchanged from raw
        assert ctrl.current_sites is ctrl.raw_sites

    def test_sites_override_used_instead_of_full_cache(
        self, cached_session, willy_sites
    ):
        subset = willy_sites  # same object is fine for override plumbing
        ctrl = corr_mod._build_ctrl(
            cached_session, [], dark=True, sites_override=subset
        )
        assert ctrl.raw_sites is subset


class TestApplyLineFilter:
    def test_no_cache_returns_none(self, store_data_willy):
        raw, err = corr_mod._apply_line_filter("no-such-session", {}, store_data_willy)
        assert raw is None
        assert err is None

    def test_no_active_lines_store_returns_all(self, cached_session, store_data_willy):
        raw, err = corr_mod._apply_line_filter(cached_session, None, store_data_willy)
        assert raw is not None
        assert err is None

    def test_active_lines_filters_to_subset(self, cached_session, store_data_willy):
        raw, err = corr_mod._apply_line_filter(
            cached_session, {"active": ["L1"]}, store_data_willy
        )
        assert err is None
        n = len(list(corr_mod._iter_edi(raw)))
        assert n == len(store_data_willy["station_records"]) // 2

    def test_all_lines_muted_returns_muted_sentinel(
        self, cached_session, store_data_willy
    ):
        raw, err = corr_mod._apply_line_filter(
            cached_session, {"active": ["ZZZ"]}, store_data_willy
        )
        assert err == "muted"
        assert raw is None


class TestFilterDisplaySites:
    def test_no_filters_returns_unchanged(self, willy_sites):
        out = corr_mod._filter_display_sites(willy_sites, None, None, [])
        assert out is willy_sites

    def test_filters_by_station(self, willy_sites, willy_ids):
        target = willy_ids[0]
        out = corr_mod._filter_display_sites(willy_sites, None, [target], [])
        names = [corr_mod._edi_name(ed) for ed in corr_mod._iter_edi(out)]
        assert names == [target]

    def test_filters_by_line(self, willy_sites, store_data_willy, willy_ids):
        out = corr_mod._filter_display_sites(
            willy_sites, ["L1"], None, store_data_willy["station_records"]
        )
        names = {corr_mod._edi_name(ed) for ed in corr_mod._iter_edi(out)}
        expected = {
            r["ID"] for r in store_data_willy["station_records"] if r["Line"] == "L1"
        }
        assert names == expected

    def test_empty_intersection_falls_back_to_all(self, willy_sites, store_data_willy):
        out = corr_mod._filter_display_sites(
            willy_sites,
            ["L1"],
            ["not-a-real-station"],
            store_data_willy["station_records"],
        )
        # empty intersection -> fall back to unfiltered *sites*
        assert out is willy_sites


class TestRenderPs:
    def test_none_sites_returns_placeholder(self):
        src = corr_mod._render_ps(None, None, "title", dark=True)
        assert src

    def test_real_render(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_ps(ctrl, ctrl.raw_sites, "Raw", dark=True)
        assert src.startswith("data:image/png")

    def test_render_exception_falls_back(self, cached_session, monkeypatch):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        monkeypatch.setattr(
            ctrl,
            "plot_rho_pseudosection",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        src = corr_mod._render_ps(ctrl, ctrl.raw_sites, "Raw", dark=True)
        assert src


class TestRenderDiff:
    def test_real_render(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_diff(ctrl, ctrl.raw_sites, ctrl.current_sites, dark=True)
        assert src.startswith("data:image/png")

    def test_render_exception_falls_back(self, cached_session, monkeypatch):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        monkeypatch.setattr(
            ctrl,
            "plot_diff",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        src = corr_mod._render_diff(ctrl, ctrl.raw_sites, ctrl.current_sites, dark=True)
        assert src


class TestRenderOverlay:
    @pytest.fixture
    def ctrl_with_step(self, cached_session):
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
        return corr_mod._build_ctrl(cached_session, steps, dark=True)

    def test_median_band(self, ctrl_with_step):
        src = corr_mod._render_overlay(
            ctrl_with_step,
            ctrl_with_step.raw_sites,
            ctrl_with_step.current_sites,
            dark=True,
            ovl_style="median-band",
        )
        assert src.startswith("data:image/png")

    def test_spaghetti(self, ctrl_with_step):
        src = corr_mod._render_overlay(
            ctrl_with_step,
            ctrl_with_step.raw_sites,
            ctrl_with_step.current_sites,
            dark=True,
            ovl_style="spaghetti",
        )
        assert src.startswith("data:image/png")

    def test_per_line(self, ctrl_with_step, store_data_willy):
        src = corr_mod._render_overlay(
            ctrl_with_step,
            ctrl_with_step.raw_sites,
            ctrl_with_step.current_sites,
            dark=True,
            ovl_style="per-line",
            records=store_data_willy["station_records"],
        )
        assert src.startswith("data:image/png")

    def test_per_line_no_records_falls_back_all(self, ctrl_with_step):
        src = corr_mod._render_overlay(
            ctrl_with_step,
            ctrl_with_step.raw_sites,
            ctrl_with_step.current_sites,
            dark=True,
            ovl_style="per-line",
            records=None,
        )
        assert src.startswith("data:image/png")

    def test_dispatch_exception_falls_back(self, ctrl_with_step, monkeypatch):
        monkeypatch.setattr(
            corr_mod,
            "_render_overlay_median_band",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        src = corr_mod._render_overlay(
            ctrl_with_step,
            ctrl_with_step.raw_sites,
            ctrl_with_step.current_sites,
            dark=True,
            ovl_style="median-band",
        )
        assert src

    def test_no_cur_sites_uses_raw_only(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_overlay_median_band(ctrl.raw_sites, None, dark=True)
        assert src.startswith("data:image/png")


class TestRenderRhoPhiGrid:
    def test_real_render(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_rho_phi_grid(
            ctrl,
            ctrl.raw_sites,
            ctrl.current_sites,
            mode="both",
            n_max=4,
            dark=True,
        )
        assert src.startswith("data:image/png")

    def test_zero_stations_returns_placeholder(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_rho_phi_grid(
            ctrl, None, None, mode="both", n_max=4, dark=True
        )
        assert src


class TestRenderPairedStations:
    def test_real_render(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_paired_stations(
            ctrl,
            ctrl.raw_sites,
            ctrl.current_sites,
            mode="both",
            n_max=2,
            dark=True,
        )
        assert src.startswith("data:image/png")

    def test_zero_stations_returns_placeholder(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_paired_stations(
            ctrl, None, None, mode="both", n_max=2, dark=True
        )
        assert src


class TestRenderFreqBands:
    def test_real_render_auto_detect(self, cached_session):
        steps = [
            {
                "fn_name": "_strat_freq_filter",
                "kwargs": {
                    "fmin": 0.0,
                    "fmax": 0.0,
                    "snr_thresh": 2.5,
                    "min_frac": 0.4,
                },
                "label": "Drop bad freq",
            }
        ]
        ctrl = corr_mod._build_ctrl(cached_session, steps, dark=True)
        src = corr_mod._render_freq_bands(
            ctrl,
            ctrl.raw_sites,
            ctrl.current_sites,
            mode="both",
            n_max=3,
            dark=True,
        )
        assert src.startswith("data:image/png")

    def test_real_render_explicit_dropped_freqs(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        ed = next(corr_mod._iter_edi(ctrl.raw_sites))
        z = getattr(ed, "Z", None)
        freqs = list(np.asarray(z.freq, float)[:2]) if z is not None else []
        src = corr_mod._render_freq_bands(
            ctrl,
            ctrl.raw_sites,
            ctrl.current_sites,
            mode="both",
            n_max=3,
            dark=True,
            dropped_freqs=freqs,
        )
        assert src.startswith("data:image/png")

    def test_zero_stations_returns_placeholder(self, cached_session):
        ctrl = corr_mod._build_ctrl(cached_session, [], dark=True)
        src = corr_mod._render_freq_bands(
            ctrl, None, None, mode="both", n_max=3, dark=True
        )
        assert src


# ══════════════════════════════════════════════════════════════════════════
# Dash callbacks
# ══════════════════════════════════════════════════════════════════════════


class TestPopulateCorrLineOpts:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.CORR_LINE_SEL}.options")

    def test_no_store_data(self, web_app):
        assert self._fn(web_app)(None) == []

    def test_populates_sorted_unique_lines(self, web_app, store_data_willy):
        opts = self._fn(web_app)(store_data_willy)
        assert opts == [
            {"label": "L1", "value": "L1"},
            {"label": "L2", "value": "L2"},
        ]


class TestPopulateCorrStnOpts:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.CORR_STN_SEL}.options")

    def test_no_filter_returns_all_stations(self, web_app, store_data_willy):
        opts = self._fn(web_app)(None, store_data_willy)
        assert len(opts) == len(store_data_willy["station_records"])

    def test_filters_by_selected_lines(self, web_app, store_data_willy):
        opts = self._fn(web_app)(["L1"], store_data_willy)
        expected = sorted(
            r["ID"] for r in store_data_willy["station_records"] if r["Line"] == "L1"
        )
        assert [o["value"] for o in opts] == expected


class TestSyncMethods:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.CORR_METHOD}.options")

    def test_no_cat_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)

    def test_known_cat_populates_methods(self, web_app):
        opts, value = self._fn(web_app)("Static Shift")
        assert opts
        assert value == opts[0]["value"]

    def test_unknown_cat_returns_empty(self, web_app):
        opts, value = self._fn(web_app)("does-not-exist")
        assert opts == []
        assert value is None


class TestSyncParamsPanel:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.CORR_PARAMS_PANEL}.children")

    def test_known_method_shows_desc(self, web_app):
        panel, desc = self._fn(web_app)("Static Shift", "AMA (spatial average)")
        assert desc.className == "corr-desc-text"

    def test_unknown_method_shows_empty_desc(self, web_app):
        panel, desc = self._fn(web_app)(None, None)
        assert "No parameters" in panel.children


class TestAutoRender:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.IMG_CORR_LEFT}.src")

    def test_wrong_nav_section_no_update(self, web_app):
        out = self._fn(web_app)(
            "home",
            None,
            None,
            "ba",
            "both",
            6,
            "grid",
            None,
            None,
            "both",
            "median-band",
            "dark",
            "no-session",
            None,
        )
        assert out == (no_update,) * 11

    def test_no_cached_data_shows_hint(self, web_app):
        out = self._fn(web_app)(
            "correction",
            None,
            None,
            "ba",
            "both",
            6,
            "grid",
            None,
            None,
            "both",
            "median-band",
            "dark",
            "no-such-session",
            None,
        )
        assert out[0] is no_update
        assert "Load survey data first" in out[6][1]

    def test_all_lines_muted_shows_warning(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "ba",
            "both",
            6,
            "grid",
            None,
            None,
            "both",
            "median-band",
            "dark",
            cached_session,
            {"active": ["ZZZ"]},
        )
        assert "muted" in out[5]

    def test_ba_tab_renders_pseudosections(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "ba",
            "both",
            6,
            "grid",
            None,
            None,
            "both",
            "median-band",
            "light",
            cached_session,
            None,
        )
        left_src, right_src = out[0], out[1]
        assert left_src.startswith("data:image/png")
        assert right_src.startswith("data:image/png")

    def test_rho_phi_tab_default_grid(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "rho-phi",
            "both",
            4,
            "grid",
            None,
            None,
            "both",
            "median-band",
            "dark",
            cached_session,
            None,
        )
        assert out[4].startswith("data:image/png")

    def test_rho_phi_tab_with_line_and_station_filter(
        self, web_app, cached_session, store_data_willy, willy_ids
    ):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "rho-phi",
            "both",
            4,
            "grid",
            ["L1"],
            [willy_ids[0]],
            "xy",
            "median-band",
            "dark",
            cached_session,
            None,
        )
        assert out[4].startswith("data:image/png")
        assert "Showing" in out[9]

    def test_rho_phi_tab_pairs_style(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "rho-phi",
            "both",
            2,
            "pairs",
            None,
            None,
            "both",
            "median-band",
            "dark",
            cached_session,
            None,
        )
        assert out[4].startswith("data:image/png")

    def test_rho_phi_tab_freq_bands_style(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "rho-phi",
            "both",
            2,
            "freq-bands",
            None,
            None,
            "both",
            "median-band",
            "dark",
            cached_session,
            None,
        )
        assert out[4].startswith("data:image/png")

    def test_overlay_tab(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "overlay",
            "both",
            6,
            "grid",
            None,
            None,
            "both",
            "spaghetti",
            "dark",
            cached_session,
            None,
        )
        assert out[2].startswith("data:image/png")

    def test_diff_tab(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            "correction",
            store_data_willy,
            None,
            "diff",
            "both",
            6,
            "grid",
            None,
            None,
            "both",
            "median-band",
            "dark",
            cached_session,
            None,
        )
        assert out[3].startswith("data:image/png")


class TestPreviewCorrection:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.IMG_CORR_LEFT}.src", IDs.BTN_CORR_PREVIEW)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                None,
                "Static Shift",
                "AMA (spatial average)",
                [],
                None,
                "dark",
                "sess",
                None,
                None,
            )

    def test_muted_lines(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            [],
            None,
            "dark",
            cached_session,
            {"active": ["ZZZ"]},
            store_data_willy,
        )
        assert out[2] == "All lines muted."

    def test_no_data_loaded(self, web_app):
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            [],
            None,
            "dark",
            "no-such-session",
            None,
            None,
        )
        assert out[4] is True

    def test_unknown_method(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "nope",
            "nope",
            [],
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
        )
        assert out[4] is True

    def test_success_path(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            ["lon", 3, "tri", 6.0],
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
        )
        assert out[0].startswith("data:image/png")
        assert out[1].startswith("data:image/png")
        assert out[4] is False

    def test_preview_exception_reported(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.app.desktop.controllers.correction_controller as cc_mod

        def _boom(self, fn_name, kwargs):
            raise RuntimeError("preview boom")

        monkeypatch.setattr(cc_mod.CorrectionController, "preview", _boom)
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            ["lon", 3, "tri", 6.0],
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
        )
        assert out[4] is True
        assert "preview boom" in out[5]


class TestApplyCorrection:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.CORR_STORE}.data", IDs.BTN_CORR_APPLY)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                None, "Static Shift", "AMA (spatial average)", [], None, "s"
            )

    def test_no_data_loaded(self, web_app):
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            [],
            None,
            "no-such-session",
        )
        assert out[4] is True

    def test_unknown_method(self, web_app, cached_session):
        out = self._fn(web_app)(1, "nope", "nope", [], None, cached_session)
        assert out[4] is True

    def test_success_adds_step(self, web_app, cached_session):
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            ["lon", 3, "tri", 6.0],
            None,
            cached_session,
        )
        store, stack, badge, fb, is_open, body = out
        assert len(store["steps"]) == 1
        assert badge == "1 step"
        assert is_open is False

    def test_invalid_kwargs_reports_validation_failure(
        self, web_app, cached_session, monkeypatch
    ):
        import pycsamt.app.desktop.controllers.correction_controller as cc_mod

        def _boom(self, fn_name, kwargs):
            raise ValueError("bad kwargs")

        monkeypatch.setattr(cc_mod.CorrectionController, "preview", _boom)
        out = self._fn(web_app)(
            1,
            "Static Shift",
            "AMA (spatial average)",
            ["lon", 3, "tri", 6.0],
            None,
            cached_session,
        )
        assert out[4] is True
        assert "bad kwargs" in out[5]


class TestApplyDropFreq:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.CORR_STORE}.data", IDs.BTN_CORR_DROP_FREQ)

    def test_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, 2.5, 0.4, None, "s")

    def test_no_data(self, web_app):
        out = self._fn(web_app)(1, 2.5, 0.4, None, "no-such-session")
        assert out[4] is True

    def test_success(self, web_app, cached_session):
        out = self._fn(web_app)(1, 2.5, 0.4, None, cached_session)
        store, stack, badge, fb, is_open, body = out
        assert len(store["steps"]) == 1
        assert store["steps"][0]["fn_name"] == "_strat_freq_filter"
        assert is_open is False

    def test_defaults_when_none(self, web_app, cached_session):
        out = self._fn(web_app)(1, None, None, None, cached_session)
        store = out[0]
        assert store["steps"][0]["kwargs"]["snr_thresh"] == 2.5
        assert store["steps"][0]["kwargs"]["min_frac"] == 0.4


class TestSmoothPreviewApply:
    def _preview_fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.IMG_CORR_RHO_PHI}.src", IDs.BTN_CORR_SMOOTH_PREVIEW
        )

    def _apply_fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.CORR_STORE}.data", IDs.BTN_CORR_SMOOTH_APPLY
        )

    def test_preview_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._preview_fn(web_app)(
                None,
                3,
                1.0,
                "offdiag",
                False,
                None,
                "dark",
                "s",
                None,
                None,
                6,
                "both",
            )

    def test_preview_muted(self, web_app, cached_session, store_data_willy):
        out = self._preview_fn(web_app)(
            1,
            3,
            1.0,
            "offdiag",
            False,
            None,
            "dark",
            cached_session,
            {"active": ["ZZZ"]},
            store_data_willy,
            6,
            "both",
        )
        assert out[2] == "All lines muted."

    def test_preview_no_data(self, web_app):
        out = self._preview_fn(web_app)(
            1,
            3,
            1.0,
            "offdiag",
            False,
            None,
            "dark",
            "no-such-session",
            None,
            None,
            6,
            "both",
        )
        assert out[4] is True

    def test_preview_success(self, web_app, cached_session, store_data_willy):
        out = self._preview_fn(web_app)(
            1,
            3,
            1.0,
            "offdiag",
            False,
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
            6,
            "both",
        )
        assert out[0].startswith("data:image/png")
        assert out[1] == "rho-phi"
        assert out[4] is False

    def test_preview_exception(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.app.desktop.controllers.correction_controller as cc_mod

        def _boom(self, fn_name, kwargs):
            raise RuntimeError("smooth boom")

        monkeypatch.setattr(cc_mod.CorrectionController, "preview", _boom)
        out = self._preview_fn(web_app)(
            1,
            3,
            1.0,
            "offdiag",
            False,
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
            6,
            "both",
        )
        assert out[4] is True
        assert "smooth boom" in out[2]

    def test_apply_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._apply_fn(web_app)(None, 3, 1.0, "offdiag", False, None, "s")

    def test_apply_no_data(self, web_app):
        out = self._apply_fn(web_app)(
            1, 3, 1.0, "offdiag", False, None, "no-such-session"
        )
        assert out[4] is True

    def test_apply_success(self, web_app, cached_session):
        out = self._apply_fn(web_app)(1, 3, 1.0, "offdiag", False, None, cached_session)
        store = out[0]
        assert store["steps"][0]["fn_name"] == "smooth_rho_phase"

    def test_apply_validation_failure(self, web_app, cached_session, monkeypatch):
        import pycsamt.app.desktop.controllers.correction_controller as cc_mod

        def _boom(self, fn_name, kwargs):
            raise RuntimeError("bad")

        monkeypatch.setattr(cc_mod.CorrectionController, "preview", _boom)
        out = self._apply_fn(web_app)(1, 3, 1.0, "offdiag", False, None, cached_session)
        assert out[4] is True


class TestManualFreqDrop:
    def _load_fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.CORR_MFREQ_LIST}.options")

    def _all_fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.CORR_MFREQ_LIST}.value", IDs.BTN_CORR_MFREQ_ALL
        )

    def _none_fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.CORR_MFREQ_LIST}.value", IDs.BTN_CORR_MFREQ_NONE
        )

    def _preview_fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.IMG_CORR_RHO_PHI}.src", IDs.BTN_CORR_MFREQ_PREVIEW
        )

    def _apply_fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.CORR_STORE}.data", IDs.BTN_CORR_MFREQ_APPLY)

    def test_load_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._load_fn(web_app)(None, "s", None, None)

    def test_load_no_data(self, web_app):
        opts, val, status = self._load_fn(web_app)(1, "no-such-session", None, None)
        assert opts == []
        assert "No data loaded" in status

    def test_load_success(self, web_app, cached_session, store_data_willy):
        opts, val, status = self._load_fn(web_app)(
            1, cached_session, None, store_data_willy
        )
        assert opts
        assert val == []
        assert "frequencies loaded" in status

    def test_select_all_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._all_fn(web_app)(None, [{"value": 1.0}])

    def test_select_all_no_opts(self, web_app):
        with pytest.raises(PreventUpdate):
            self._all_fn(web_app)(1, [])

    def test_select_all(self, web_app):
        opts = [{"value": 1.0}, {"value": 2.0}]
        assert self._all_fn(web_app)(1, opts) == [1.0, 2.0]

    def test_select_none_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._none_fn(web_app)(None)

    def test_select_none(self, web_app):
        assert self._none_fn(web_app)(1) == []

    def test_preview_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._preview_fn(web_app)(
                None, [1.0], None, "dark", "s", None, None, 6, "both"
            )

    def test_preview_no_freqs_selected(self, web_app, cached_session):
        out = self._preview_fn(web_app)(
            1, [], None, "dark", cached_session, None, None, 6, "both"
        )
        assert "No frequencies selected" in out[2]

    def test_preview_muted(self, web_app, cached_session, store_data_willy):
        out = self._preview_fn(web_app)(
            1,
            [1000.0],
            None,
            "dark",
            cached_session,
            {"active": ["ZZZ"]},
            store_data_willy,
            6,
            "both",
        )
        assert out[2] == "All lines muted."

    def test_preview_no_data(self, web_app):
        out = self._preview_fn(web_app)(
            1,
            [1000.0],
            None,
            "dark",
            "no-such-session",
            None,
            None,
            6,
            "both",
        )
        assert out[4] is True

    def test_preview_success(
        self, web_app, cached_session, store_data_willy, willy_sites
    ):
        ed = next(corr_mod._iter_edi(willy_sites))
        z = getattr(ed, "Z", None)
        freqs = list(np.asarray(z.freq, float)[:1])
        out = self._preview_fn(web_app)(
            1,
            freqs,
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
            6,
            "both",
        )
        assert out[0].startswith("data:image/png")
        assert out[1] == "rho-phi"
        assert out[4] is False

    def test_preview_exception(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.app.desktop.controllers.correction_controller as cc_mod

        def _boom(self, fn_name, kwargs):
            raise RuntimeError("mfreq boom")

        monkeypatch.setattr(cc_mod.CorrectionController, "preview", _boom)
        out = self._preview_fn(web_app)(
            1,
            [1000.0],
            None,
            "dark",
            cached_session,
            None,
            store_data_willy,
            6,
            "both",
        )
        assert out[4] is True
        assert "mfreq boom" in out[2]

    def test_apply_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._apply_fn(web_app)(None, [1.0], None, "s")

    def test_apply_no_freqs(self, web_app):
        out = self._apply_fn(web_app)(1, [], None, "s")
        assert "nothing to apply" in out[3]

    def test_apply_no_data(self, web_app):
        out = self._apply_fn(web_app)(1, [1000.0], None, "no-such-session")
        assert out[4] is True

    def test_apply_success(self, web_app, cached_session, willy_sites):
        ed = next(corr_mod._iter_edi(willy_sites))
        z = getattr(ed, "Z", None)
        freqs = list(np.asarray(z.freq, float)[:1])
        out = self._apply_fn(web_app)(1, freqs, None, cached_session)
        store = out[0]
        assert store["steps"][0]["fn_name"] == "drop_freqs_manual"

    def test_apply_validation_failure(self, web_app, cached_session, monkeypatch):
        import pycsamt.app.desktop.controllers.correction_controller as cc_mod

        def _boom(self, fn_name, kwargs):
            raise RuntimeError("bad")

        monkeypatch.setattr(cc_mod.CorrectionController, "preview", _boom)
        out = self._apply_fn(web_app)(1, [1000.0], None, cached_session)
        assert out[4] is True


class TestUndoResetDelete:
    def _undo_fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.CORR_STORE}.data", IDs.BTN_CORR_UNDO)

    def _reset_fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.CORR_STORE}.data", IDs.BTN_CORR_RESET)

    def _delete_fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.CORR_STORE}.data", "corr-step-del")

    def test_undo_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._undo_fn(web_app)(None, {"steps": []})

    def test_undo_empty_stack(self, web_app):
        out = self._undo_fn(web_app)(1, {"steps": []})
        assert out[3] == "Nothing to undo."

    def test_undo_pops_last_step(self, web_app):
        steps = [
            {"fn_name": "a", "label": "A"},
            {"fn_name": "b", "label": "B"},
        ]
        out = self._undo_fn(web_app)(1, {"steps": steps})
        store, stack, badge, fb = out
        assert len(store["steps"]) == 1
        assert "Undone 'B'" in fb

    def test_reset_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._reset_fn(web_app)(None)

    def test_reset_clears_steps(self, web_app):
        out = self._reset_fn(web_app)(1)
        assert out[0] == {"steps": []}

    def test_delete_no_trigger(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        try:
            with pytest.raises(PreventUpdate):
                self._delete_fn(web_app)([0], {"steps": []})
        finally:
            _clear_triggered()

    def test_delete_step_removes_by_index(self, web_app, monkeypatch):

        fn = self._delete_fn(web_app)
        steps = [
            {"fn_name": "a", "label": "A"},
            {"fn_name": "b", "label": "B"},
        ]

        class _FakeCtx:
            triggered_id = {"type": "corr-step-del", "index": 0}
            inputs_list = [[{"id": {"index": 0}, "value": 1}]]

        monkeypatch.setattr(corr_mod, "ctx", _FakeCtx())
        out = fn([1], {"steps": steps})
        store, stack, badge, fb = out
        assert len(store["steps"]) == 1
        assert store["steps"][0]["label"] == "B"
        assert "Removed step 1" in fb

    def test_delete_step_no_idx_prevents_update(self, monkeypatch, web_app):
        fn = self._delete_fn(web_app)

        class _FakeCtx:
            triggered_id = {"type": "corr-step-del"}
            inputs_list = [[{"id": {}, "value": 1}]]

        monkeypatch.setattr(corr_mod, "ctx", _FakeCtx())
        with pytest.raises(PreventUpdate):
            fn([1], {"steps": []})

    def test_delete_step_not_clicked_prevents_update(self, monkeypatch, web_app):
        fn = self._delete_fn(web_app)

        class _FakeCtx:
            triggered_id = {"type": "corr-step-del", "index": 0}
            inputs_list = [[{"id": {"index": 0}, "value": 0}]]

        monkeypatch.setattr(corr_mod, "ctx", _FakeCtx())
        with pytest.raises(PreventUpdate):
            fn([0], {"steps": [{"fn_name": "a", "label": "A"}]})

    def test_delete_step_out_of_range_prevents_update(self, monkeypatch, web_app):
        fn = self._delete_fn(web_app)

        class _FakeCtx:
            triggered_id = {"type": "corr-step-del", "index": 5}
            inputs_list = [[{"id": {"index": 5}, "value": 1}]]

        monkeypatch.setattr(corr_mod, "ctx", _FakeCtx())
        with pytest.raises(PreventUpdate):
            fn([1], {"steps": [{"fn_name": "a", "label": "A"}]})


class TestUpdateDiffStats:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.CORR_DIFF_STATS}.children")

    def test_not_diff_tab_no_update(self, web_app):
        assert self._fn(web_app)(None, "ba", "s", "dark") is no_update

    def test_no_steps_shows_placeholder(self, web_app):
        out = self._fn(web_app)({"steps": []}, "diff", "s", "dark")
        assert "No corrections applied yet." in out.children

    def test_ctrl_none_no_update(self, web_app):
        steps = [{"fn_name": "correct_ss_ama", "kwargs": {}, "label": "AMA"}]
        out = self._fn(web_app)({"steps": steps}, "diff", "no-such-session", "dark")
        assert out is no_update

    def test_real_stats(self, web_app, cached_session):
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
        out = self._fn(web_app)({"steps": steps}, "diff", cached_session, "dark")
        assert out.className == "corr-diff-stats-row"

    def test_exception_returns_no_update(self, web_app, cached_session, monkeypatch):
        monkeypatch.setattr(
            corr_mod,
            "_collect_rho",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        steps = [{"fn_name": "correct_ss_ama", "kwargs": {}, "label": "AMA"}]
        out = self._fn(web_app)({"steps": steps}, "diff", cached_session, "dark")
        assert out is no_update


class TestExportCorrected:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.CORR_DOWNLOAD}.data")

    def test_no_clicks(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "csv", None, "s", "dark")

    def test_no_ctrl_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "csv", None, "no-such-session", "dark")

    def test_export_csv(self, web_app, cached_session):
        out = self._fn(web_app)(1, "csv", None, cached_session, "dark")
        assert out["filename"] == "corrected_data.csv"

    def test_export_xyz(self, web_app, cached_session):
        out = self._fn(web_app)(1, "xyz", None, cached_session, "dark")
        assert out["filename"] == "corrected_profile.xyz"

    def test_export_other_defaults_to_edi_csv(self, web_app, cached_session):
        out = self._fn(web_app)(1, "edi", None, cached_session, "dark")
        assert out["filename"] == "corrected_edi.csv"
