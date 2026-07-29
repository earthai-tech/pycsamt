# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.advanced (Advanced Plots page).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — cached via ``cache_set`` and fed to the
module-level ``_CTRL`` singleton so ``run_advanced``/``train_atom`` and
the ``_filter_sites``/``_pt_grid_profiles`` helpers exercise the real
emtools dispatch (fast rose/fingerprint plots only -- heavy ML training
in ``train_atom`` is mocked at the controller method).

Since ``_CTRL`` is a process-wide singleton (module attribute, not
per-app-instance), every test that depends on its ``_sites`` state sets
it explicitly rather than relying on ordering.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from dash import no_update

import pycsamt.app.web.callbacks.advanced as adv_mod
from pycsamt.app.web.callbacks.advanced import _filter_sites, _pt_grid_profiles
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
    """station_records-shaped list matching the real willy_sites, split
    into two synthetic "lines" so line-filter tests have something to
    group on."""
    names = [s.name for s in willy_sites]
    return [
        {"ID": name, "Line": "L1" if i % 2 == 0 else "L2"}
        for i, name in enumerate(names)
    ]


@pytest.fixture
def cached_session(willy_sites):
    from pycsamt.app.web.cache import cache_set

    session_id = "test-advanced-session"
    cache_set(session_id, willy_sites)
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
    """Disambiguate callbacks sharing an (allow_duplicate) output by the
    real Input component id from the callback_map's own metadata, since
    duplicate-output keys are hash-suffixed and not otherwise distinguishable
    by substring alone."""
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(i.get("id") == input_id for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


# ── _filter_sites ────────────────────────────────────────────────────────────


class TestFilterSites:
    def test_none_sites_returns_none(self):
        assert _filter_sites(None, None, None, None) is None

    def test_no_filter_no_active_lines_returns_all(self, willy_sites):
        assert _filter_sites(willy_sites, None, None, None) is willy_sites

    def test_explicit_station_filter(self, willy_sites, willy_records):
        wanted = [willy_records[0]["ID"]]
        result = _filter_sites(
            willy_sites, wanted, {"station_records": willy_records}, None
        )
        names = [s.name for s in result]
        assert names == wanted

    def test_active_lines_filter(self, willy_sites, willy_records):
        result = _filter_sites(
            willy_sites,
            None,
            {"station_records": willy_records},
            {"active": ["L1"]},
        )
        names = {s.name for s in result}
        expected = {r["ID"] for r in willy_records if r["Line"] == "L1"}
        assert names == expected

    def test_filter_yielding_nothing_falls_back_to_original(
        self, willy_sites, willy_records
    ):
        result = _filter_sites(
            willy_sites,
            ["station-that-does-not-exist"],
            {"station_records": willy_records},
            None,
        )
        assert result is willy_sites

    def test_sites_import_failure_falls_back_to_original(
        self, willy_sites, willy_records, monkeypatch
    ):
        import pycsamt.site.base as site_base_mod

        monkeypatch.setattr(
            site_base_mod,
            "Sites",
            property(lambda self: (_ for _ in ()).throw(RuntimeError())),
        )
        result = _filter_sites(
            willy_sites,
            [willy_records[0]["ID"]],
            {"station_records": willy_records},
            None,
        )
        assert result is willy_sites


# ── _pt_grid_profiles ────────────────────────────────────────────────────────


class TestPtGridProfiles:
    def test_groups_by_line(self, willy_records):
        result = _pt_grid_profiles({"station_records": willy_records}, None, per_line=2)
        assert set(result.keys()) == {"L1", "L2"}
        assert all(len(v) <= 2 for v in result.values())

    def test_restricted_to_active_lines(self, willy_records):
        result = _pt_grid_profiles(
            {"station_records": willy_records}, {"active": ["L1"]}, per_line=2
        )
        assert set(result.keys()) == {"L1"}

    def test_empty_records_returns_empty_dict(self):
        assert _pt_grid_profiles({}, None, per_line=2) == {}

    def test_records_without_line_are_skipped(self):
        records = [{"ID": "STA1", "Line": ""}]
        assert _pt_grid_profiles({"station_records": records}, None, per_line=2) == {}


# ── switch_adv_tab ───────────────────────────────────────────────────────────


class TestSwitchAdvTab:
    def test_no_trigger_returns_no_update(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = _cb(web_app, f"{IDs.ADV_ACTIVE_TAB}.data")
        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        try:
            result = fn(*([0] * 6))
        finally:
            cc.context_value.set({})
        assert result is no_update

    def test_trigger_from_tab_button_switches(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = _cb(web_app, f"{IDs.ADV_ACTIVE_TAB}.data")
        cc.context_value.set(
            AttributeDict(triggered_inputs=[{"prop_id": "adv-tab-pt.n_clicks"}])
        )
        try:
            result = fn(*([1] * 6))
        finally:
            cc.context_value.set({})
        assert result == "pt"


# ── sync_plots ───────────────────────────────────────────────────────────────


class TestSyncPlots:
    def test_no_active_tab_raises_prevent_update(self, web_app):
        from dash.exceptions import PreventUpdate

        fn = _cb_multi(web_app, f"{IDs.ADV_PLOT}.options", f"{IDs.ADV_PLOT}.value")
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_valid_tab_returns_options_and_first_value(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.ADV_PLOT}.options", f"{IDs.ADV_PLOT}.value")
        opts, value = fn("strike")
        assert len(opts) > 0
        assert value == opts[0]["value"]

    def test_unknown_tab_falls_back_to_strike_plots(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.ADV_PLOT}.options", f"{IDs.ADV_PLOT}.value")
        opts, _value = fn("bogus-tab")
        from pycsamt.app.desktop.controllers.advanced_controller import (
            STRIKE_PLOTS,
        )

        assert len(opts) == len(STRIKE_PLOTS)


# ── update_param_panel ───────────────────────────────────────────────────────


class TestUpdateParamPanel:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.ADV_DESCRIPTION}.children", "adv-atom-section")

    def test_no_fn_name(self, web_app):
        desc, period, atom, pt_sta, pt_grid = self._fn(web_app)(None)
        assert desc == ""
        assert period["display"] == "none"

    def test_period_fn_shows_period_section(self, web_app):
        desc, period, atom, pt_sta, pt_grid = self._fn(web_app)("plot_phase_tensor_map")
        assert period["display"] == "block"
        assert atom["display"] == "none"

    def test_atom_fn_shows_atom_section(self, web_app):
        _desc, _period, atom, _pt_sta, _pt_grid = self._fn(web_app)(
            "plot_atom_psection"
        )
        assert atom["display"] == "block"

    def test_pt_strip_fn_shows_pt_station_section(self, web_app):
        *_rest, pt_sta, pt_grid = self._fn(web_app)("plot_phase_tensor_strip")
        assert pt_sta["display"] == "block"
        assert pt_grid["display"] == "none"

    def test_pt_strip_grid_fn_shows_pt_grid_section(self, web_app):
        *_rest, pt_sta, pt_grid = self._fn(web_app)("plot_phase_tensor_strip_grid")
        assert pt_grid["display"] == "block"

    def test_other_fn_hides_all_sections(self, web_app):
        _desc, period, atom, pt_sta, pt_grid = self._fn(web_app)("plot_strike_rose")
        assert period["display"] == "none"
        assert atom["display"] == "none"
        assert pt_sta["display"] == "none"
        assert pt_grid["display"] == "none"


# ── update_info_bar ──────────────────────────────────────────────────────────


class TestUpdateInfoBar:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.ADV_INFO_BAR}.children")

    def test_defaults_when_nothing_selected(self, web_app):
        parts = self._fn(web_app)(None, None, None)
        assert len(parts) >= 1

    def test_includes_label_and_station_count(self, web_app):
        parts = self._fn(web_app)("strike", "plot_strike_rose", {"n_stations": 5})
        text_dump = str(parts)
        assert "5 station" in text_dump


# ── update_filters ───────────────────────────────────────────────────────────


class TestUpdateFilters:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            f"{IDs.ADV_LINE_FILTER}.options",
            f"{IDs.ADV_DATA_BAR}.children",
        )

    def test_no_store_data_shows_empty_bar(self, web_app):
        line_opts, stn_opts, bar = self._fn(web_app)(None, None)
        assert line_opts == []
        assert stn_opts == []

    def test_with_records_and_no_active_lines_derives_all(self, web_app, willy_records):
        line_opts, stn_opts, bar = self._fn(web_app)(
            {
                "n_stations": len(willy_records),
                "station_records": willy_records,
            },
            None,
        )
        assert {o["value"] for o in line_opts} == {"L1", "L2"}
        assert len(stn_opts) == len(willy_records)

    def test_with_explicit_active_lines(self, web_app, willy_records):
        line_opts, stn_opts, bar = self._fn(web_app)(
            {
                "n_stations": len(willy_records),
                "station_records": willy_records,
            },
            {"active": ["L1"]},
        )
        assert len(stn_opts) == sum(1 for r in willy_records if r["Line"] == "L1")


# ── line_to_stations ─────────────────────────────────────────────────────────


class TestLineToStations:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            f"{IDs.ADV_STN_FILTER}.value",
            f"{IDs.ADV_STN_FILTER}.options",
        )

    def test_no_store_data_returns_no_update(self, web_app):
        value, opts = self._fn(web_app)(["L1"], None, None)
        assert value is no_update
        assert opts is no_update

    def test_no_selection_returns_none_with_all_options(self, web_app, willy_records):
        value, opts = self._fn(web_app)(None, {"station_records": willy_records}, None)
        assert value is None
        assert len(opts) == len(willy_records)

    def test_selection_filters_stations(self, web_app, willy_records):
        value, opts = self._fn(web_app)(
            ["L1"], {"station_records": willy_records}, None
        )
        expected = {r["ID"] for r in willy_records if r["Line"] == "L1"}
        assert set(value) == expected


# ── update_pt_station_options ────────────────────────────────────────────────


class TestUpdatePtStationOptions:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.ADV_PT_STATION}.options")

    def test_no_store_data_returns_empty(self, web_app):
        assert self._fn(web_app)(None, None) == []

    def test_with_records_all_lines(self, web_app, willy_records):
        opts = self._fn(web_app)({"station_records": willy_records}, None)
        assert len(opts) == len(willy_records)

    def test_with_active_lines_filter(self, web_app, willy_records):
        opts = self._fn(web_app)({"station_records": willy_records}, {"active": ["L2"]})
        expected = {r["ID"] for r in willy_records if r["Line"] == "L2"}
        assert {o["value"] for o in opts} == expected


# ── Station All / None ───────────────────────────────────────────────────────


class TestStationAllNone:
    def _select_all_fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.ADV_STN_FILTER}.value", IDs.ADV_STN_ALL)

    def _deselect_all_fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.ADV_STN_FILTER}.value", IDs.ADV_STN_NONE)

    def test_select_all_no_click_returns_no_update(self, web_app):
        fn = self._select_all_fn(web_app)
        assert fn(0, [{"value": "a"}]) is no_update

    def test_select_all_returns_every_value(self, web_app):
        fn = self._select_all_fn(web_app)
        result = fn(1, [{"value": "a"}, {"value": "b"}])
        assert result == ["a", "b"]

    def test_select_all_no_options_returns_no_update(self, web_app):
        fn = self._select_all_fn(web_app)
        assert fn(1, None) is no_update

    def test_deselect_all_no_click_returns_no_update(self, web_app):
        fn = self._deselect_all_fn(web_app)
        assert fn(0) is no_update

    def test_deselect_all_returns_empty_list(self, web_app):
        fn = self._deselect_all_fn(web_app)
        assert fn(1) == []


# ── train_atom ───────────────────────────────────────────────────────────────


class TestTrainAtom:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            f"{IDs.ADV_TRAIN_SPINNER}.children",
            f"{IDs.ADV_TRAIN_STATUS}.children",
        )

    def test_no_clicks_raises_prevent_update(self, web_app):
        from dash.exceptions import PreventUpdate

        with pytest.raises(PreventUpdate):
            self._fn(web_app)(0, "s1", 6, 40, None, None, None)

    def test_no_data_shows_warning(self, web_app):
        _spinner, status = self._fn(web_app)(
            1, "no-such-session", 6, 40, None, None, None
        )
        assert "No data loaded" in str(status)

    def test_success_trains_and_reports(self, web_app, cached_session, monkeypatch):
        monkeypatch.setattr(
            adv_mod._CTRL,
            "train_dim_model",
            lambda n_atoms, n_iter: {"D": None},
        )
        _spinner, status = self._fn(web_app)(
            1, cached_session, 6, 40, None, {"station_records": []}, None
        )
        assert "Trained" in str(status)

    def test_exception_reports_warning(self, web_app, cached_session, monkeypatch):
        def _boom(n_atoms, n_iter):
            raise RuntimeError("training failed")

        monkeypatch.setattr(adv_mod._CTRL, "train_dim_model", _boom)
        _spinner, status = self._fn(web_app)(
            1, cached_session, 6, 40, None, {"station_records": []}, None
        )
        assert "training failed" in str(status)


# ── run_advanced ─────────────────────────────────────────────────────────────


class TestRunAdvanced:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            f"{IDs.IMG_ADV_STRIKE}.src",
            f"{IDs.ADV_SPINNER}.children",
        )

    def test_no_clicks_raises_prevent_update(self, web_app):
        from dash.exceptions import PreventUpdate

        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                0,
                "strike",
                "plot_strike_rose",
                "s",
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )

    def test_no_fn_name_raises_prevent_update(self, web_app):
        from dash.exceptions import PreventUpdate

        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1,
                "strike",
                None,
                "s",
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )

    def test_generic_plot_success_writes_active_tab_only(self, web_app, cached_session):
        result = self._fn(web_app)(
            1,
            "strike",
            "plot_strike_rose",
            cached_session,
            None,
            None,
            None,
            None,
            "10x6",
            None,
            None,
        )
        *imgs, spinner, is_open, body = result
        assert is_open is False
        assert imgs[0] != no_update  # strike tab (index 0) written
        assert all(v is no_update for v in imgs[1:])

    def test_generic_plot_exception_opens_toast(
        self, web_app, cached_session, monkeypatch
    ):
        monkeypatch.setattr(
            adv_mod._CTRL,
            "draw",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("draw boom")),
        )
        result = self._fn(web_app)(
            1,
            "strike",
            "plot_strike_rose",
            cached_session,
            None,
            None,
            None,
            None,
            "10x6",
            None,
            None,
        )
        *imgs, spinner, is_open, body = result
        assert is_open is True
        assert "draw boom" in body

    def test_pt_strip_without_station_returns_warning(self, web_app, cached_session):
        result = self._fn(web_app)(
            1,
            "pt",
            "plot_phase_tensor_strip",
            cached_session,
            None,
            None,
            None,
            None,
            "10x6",
            None,
            None,
        )
        *imgs, spinner, is_open, body = result
        assert is_open is True
        assert "Pick a station" in body

    def test_pt_strip_with_station_success(
        self, web_app, cached_session, willy_records
    ):
        station = willy_records[0]["ID"]
        result = self._fn(web_app)(
            1,
            "pt",
            "plot_phase_tensor_strip",
            cached_session,
            {"station_records": willy_records},
            None,
            None,
            None,
            "10x6",
            station,
            None,
        )
        *imgs, spinner, is_open, body = result
        assert is_open is False
        assert imgs[_adv_tab_index("pt")] != no_update

    def test_pt_strip_grid_without_profiles_returns_warning(
        self, web_app, cached_session
    ):
        result = self._fn(web_app)(
            1,
            "pt",
            "plot_phase_tensor_strip_grid",
            cached_session,
            {"station_records": []},
            None,
            None,
            None,
            "10x6",
            None,
            4,
        )
        *imgs, spinner, is_open, body = result
        assert is_open is True
        assert "No survey lines" in body

    def test_pt_strip_grid_with_profiles_success(
        self, web_app, cached_session, willy_records
    ):
        result = self._fn(web_app)(
            1,
            "pt",
            "plot_phase_tensor_strip_grid",
            cached_session,
            {"station_records": willy_records},
            None,
            None,
            None,
            "10x6",
            None,
            4,
        )
        *imgs, spinner, is_open, body = result
        assert is_open is False


def _adv_tab_index(tab: str) -> int:
    return ["strike", "pt", "induction", "impedance", "depth", "survey"].index(tab)
