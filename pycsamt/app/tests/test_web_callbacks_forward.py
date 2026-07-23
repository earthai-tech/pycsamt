# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.forward (Forward Modelling page).

Strategy
--------
* The 1D/2D/3D solvers (``MT1DForward``/``MT2DForward``/``MT3DForward``)
  are real, fast, analytic/numeric solvers (not ML training), so
  ``run_forward`` is exercised for real with small grids/point counts
  rather than mocked.
* ``ForwardController`` (the module-level ``_CTRL`` singleton in
  forward.py) persists its model library to a real file at
  ``~/.pycsamt/forward_models.json``; ``save_model`` tests redirect the
  module's ``_LIBRARY_PATH`` to a tmp_path first (``_CTRL`` itself is
  already constructed by the time a test runs, but ``save_library()``
  re-reads the module-global path fresh on every call).
"""

from __future__ import annotations

import pytest
from dash import no_update
from dash.exceptions import PreventUpdate

import pycsamt.app.desktop.controllers.forward_controller as fc_mod
import pycsamt.app.web.callbacks.forward as fwd_mod
from pycsamt.app.web.layout import IDs


@pytest.fixture(autouse=True)
def _isolated_library(tmp_path, monkeypatch):
    monkeypatch.setattr(fc_mod, "_LIBRARY_PATH", tmp_path / "models.json")


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


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id}])
    )


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


# ── Tab switching ────────────────────────────────────────────────────────────


class TestTabSwitching:
    def test_switch_dim_no_trigger_returns_no_update(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.FWD_ACTIVE_TAB}.data", "fwd-dim-btn-1d"
        )
        _set_triggered("")
        try:
            assert fn(1, 0, 0) is no_update
        finally:
            _clear_triggered()

    def test_switch_dim_returns_tab_key(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.FWD_ACTIVE_TAB}.data", "fwd-dim-btn-2d"
        )
        _set_triggered("fwd-dim-btn-2d.n_clicks")
        try:
            assert fn(0, 1, 0) == "fwd-tab-2d"
        finally:
            _clear_triggered()

    def test_switch_view_tab_returns_tab_key(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.FWD_ACTIVE_TAB}.data", "fwd-view-tab-3d"
        )
        _set_triggered("fwd-view-tab-3d.n_clicks")
        try:
            assert fn(0, 0, 1) == "fwd-tab-3d"
        finally:
            _clear_triggered()

    def test_switch_view_tab_no_trigger_returns_no_update(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.FWD_ACTIVE_TAB}.data", "fwd-view-tab-1d"
        )
        _set_triggered("")
        try:
            assert fn(0, 0, 0) is no_update
        finally:
            _clear_triggered()


# ── Context info bar ─────────────────────────────────────────────────────────


class TestUpdateCtxBar:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.FWD_CTX_BAR}.children")

    def test_defaults_to_1d(self, web_app):
        parts = self._fn(web_app)(None, None)
        assert "1-D" in str(parts)

    def test_1d_tab_shows_method(self, web_app):
        parts = self._fn(web_app)("fwd-tab-1d", "CSAMT1D")
        assert "CSAMT1D" in str(parts)

    def test_non_1d_tab_hides_method(self, web_app):
        parts = self._fn(web_app)("fwd-tab-2d", "CSAMT1D")
        assert "CSAMT1D" not in str(parts)
        assert "2-D MT" in str(parts)


# ── sync_method_ui / sync_2d_model / sync_3d_model ───────────────────────────


class TestSyncMethodUi:
    def _fn(self, web_app):
        return _cb_multi(
            web_app, f"{IDs.FWD_CSAMT_CARD}.style", f"{IDs.FWD_DIM}.data"
        )

    def test_default_mt1d(self, web_app):
        csamt, tem, label, fmin, fmax, dim = self._fn(web_app)(None)
        assert csamt["display"] == "none"
        assert tem["display"] == "none"
        assert "Frequency" in label
        assert dim == "MT1D"

    def test_csamt1d_shows_csamt_card(self, web_app):
        csamt, tem, *_rest = self._fn(web_app)("CSAMT1D")
        assert csamt["display"] == "block"
        assert tem["display"] == "none"

    def test_tem1d_shows_tem_card_and_time_label(self, web_app):
        csamt, tem, label, fmin, fmax, dim = self._fn(web_app)("TEM1D")
        assert tem["display"] == "block"
        assert "Time range" in label
        assert (fmin, fmax) == (-6, -2)


class TestSync2D3DModel:
    def test_2d_anomaly(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.FWD2_ANOMALY_CARD}.style")
        anomaly, random_ = fn("anomaly")
        assert anomaly["display"] == "block"
        assert random_["display"] == "none"

    def test_2d_random(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.FWD2_ANOMALY_CARD}.style")
        anomaly, random_ = fn("random")
        assert anomaly["display"] == "none"
        assert random_["display"] == "block"

    def test_3d_block(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.FWD3_BLOCK_CARD}.style")
        block, random_ = fn("block")
        assert block["display"] == "block"
        assert random_["display"] == "none"


# ── load_preset ──────────────────────────────────────────────────────────────


class TestLoadPreset:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            f"{IDs.FWD_LAYER_TABLE}.data",
            f"{IDs.FWD_HALFSPACE_RHO}.value",
        )

    def test_no_clicks_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(0, "sedimentary")

    def test_no_preset_name_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None)

    def test_success_loads_real_preset(self, web_app):
        rows, halfspace, msg = self._fn(web_app)(1, "sedimentary")
        assert rows and all("resistivity" in r for r in rows)
        assert isinstance(halfspace, float)
        assert "loaded" in msg

    def test_unknown_preset_reports_error(self, web_app):
        rows, halfspace, msg = self._fn(web_app)(1, "not-a-real-geology")
        assert rows is no_update
        assert halfspace is no_update
        assert msg.startswith("Error:")


# ── add_layer / renumber_layers ──────────────────────────────────────────────


class TestAddLayer:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.FWD_LAYER_TABLE}.data", IDs.BTN_FWD_ADD_LAYER
        )

    def test_no_clicks_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(0, [])

    def test_appends_row_with_default_when_empty(self, web_app):
        rows = self._fn(web_app)(1, None)
        assert rows == [{"layer": 1, "resistivity": 100.0, "thickness": 500}]

    def test_appends_row_reusing_last_resistivity(self, web_app):
        existing = [{"layer": 1, "resistivity": 250.0, "thickness": 300}]
        rows = self._fn(web_app)(1, existing)
        assert rows[-1] == {
            "layer": 2,
            "resistivity": 250.0,
            "thickness": 500,
        }


class TestRenumberLayers:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.FWD_LAYER_TABLE}.data", IDs.FWD_LAYER_TABLE
        )

    def test_no_data_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)

    def test_renumbers_sequentially(self, web_app):
        data = [
            {"layer": 5, "resistivity": 1, "thickness": 1},
            {"layer": 9, "resistivity": 2, "thickness": 2},
        ]
        result = self._fn(web_app)(data)
        assert [r["layer"] for r in result] == [1, 2]


# ── run_forward ──────────────────────────────────────────────────────────────


def _run_forward_args(**overrides):
    """Full positional-arg tuple for run_forward with light/fast defaults."""
    args = dict(
        n_clicks=1,
        dimension="fwd-tab-1d",
        table_data=[{"layer": 1, "resistivity": 100.0, "thickness": 300.0}],
        halfspace_rho=1000.0,
        method="MT1D",
        freq_min=-1,
        freq_max=1,
        n_freq=5,
        offset=5000,
        loop_r=50,
        m2d_type="halfspace",
        bg2=100.0,
        nx2=4,
        nz2=4,
        x2max=2000,
        z2max=2000,
        n2st=3,
        a2rho=5,
        ax2lo=500,
        ax2hi=1000,
        az2lo=100,
        az2hi=300,
        n2lay=3,
        seed2=1,
        plot2="model",
        m3d_type="halfspace",
        bg3=100.0,
        nx3=3,
        ny3=3,
        nz3=3,
        x3max=2000,
        y3max=2000,
        z3max=2000,
        nx3st=2,
        ny3st=2,
        a3rho=5,
        ax3lo=500,
        ax3hi=1000,
        ay3lo=500,
        ay3hi=1000,
        az3lo=100,
        az3hi=300,
        n3lay=3,
        seed3=1,
        plot3="model",
        comp3="xy",
        freq3idx=0,
        theme="dark",
        session_id="test-fwd-session",
    )
    args.update(overrides)
    return tuple(args.values())


class TestRunForward:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            f"{IDs.IMG_FWD_1D}.src",
            f"{IDs.FWD_SPINNER}.children",
        )

    def test_no_clicks_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(*_run_forward_args(n_clicks=0))

    def test_run_1d_mt1d_success(self, web_app):
        result = self._fn(web_app)(*_run_forward_args())
        img1d, img2d, img3d, spinner, feedback, is_open, body = result
        assert is_open is False
        assert img1d != no_update
        assert img2d is no_update
        assert img3d is no_update
        assert "MT1D done" in feedback

    def test_run_1d_csamt_success(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(method="CSAMT1D")
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False
        assert "CSAMT1D done" in feedback

    def test_run_1d_tem_success(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(
                method="TEM1D", freq_min=-6, freq_max=-4, n_freq=5
            )
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False
        assert "TEM1D done" in feedback

    def test_run_1d_missing_layers_reports_error(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(table_data=[])
        )
        img1d, img2d, img3d, spinner, feedback, is_open, body = result
        assert is_open is True
        assert "at least one finite layer" in body

    def test_run_2d_halfspace_success(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-2d")
        )
        img1d, img2d, img3d, spinner, feedback, is_open, body = result
        assert is_open is False
        assert img1d is no_update
        assert img2d != no_update
        assert img3d is no_update
        assert "MT2D done" in feedback

    def test_run_2d_anomaly_model(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-2d", m2d_type="anomaly")
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_2d_random_model(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-2d", m2d_type="random")
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_2d_from_layers_model(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-2d", m2d_type="from_layers")
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    @pytest.mark.parametrize(
        "plot2",
        [
            "psd_te_rho",
            "psd_tm_rho",
            "psd_te_phi",
            "psd_tm_phi",
            "profiles_te",
            "profiles_tm",
        ],
    )
    def test_run_2d_pseudosection_and_profile_plots(self, web_app, plot2):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-2d", plot2=plot2)
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_3d_halfspace_success(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-3d")
        )
        img1d, img2d, img3d, spinner, feedback, is_open, body = result
        assert is_open is False
        assert img3d != no_update
        assert "MT3D" in feedback

    def test_run_3d_block_model_map_plot(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(
                dimension="fwd-tab-3d", m3d_type="block", plot3="map"
            )
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_3d_random_model_tensor_plot(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(
                dimension="fwd-tab-3d", m3d_type="random", plot3="tensor"
            )
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_3d_section_plot(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-3d", plot3="section")
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_3d_unknown_plot_type_falls_back_to_model(self, web_app):
        result = self._fn(web_app)(
            *_run_forward_args(dimension="fwd-tab-3d", plot3="bogus")
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False

    def test_run_1d_tem_with_no_valid_response(self, web_app, monkeypatch):
        """Force TEM1DForward to return an all-NaN response so
        _plot_tem_response's "No valid TEM response" text branch runs."""
        import numpy as np

        class _FakeTemResp:
            def __init__(self, times):
                self.times = times
                self.dBz_dt = np.full_like(times, np.nan)

        class _FakeTEM1DForward:
            def __init__(self, times, loop_radius):
                self._times = times

            def run(self, model):
                return _FakeTemResp(self._times)

        import pycsamt.forward.em1d as em1d_mod

        monkeypatch.setattr(em1d_mod, "TEM1DForward", _FakeTEM1DForward)
        result = self._fn(web_app)(
            *_run_forward_args(
                dimension="fwd-tab-1d",
                method="TEM1D",
                freq_min=-6,
                freq_max=-4,
                n_freq=5,
            )
        )
        *_imgs, feedback, is_open, _body = result
        assert is_open is False
        assert "TEM1D done" in feedback

    def test_run_light_theme(self, web_app):
        result = self._fn(web_app)(*_run_forward_args(theme="light"))
        *_imgs, feedback, is_open, _body = result
        assert is_open is False


# ── save_model ───────────────────────────────────────────────────────────────


class TestSaveModel:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.FWD_FEEDBACK}.children", IDs.BTN_FWD_SAVE
        )

    def test_no_clicks_raises_prevent_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(0, "M1", [], 1000.0, "MT1D")

    def test_no_name_prompts_user(self, web_app):
        msg = self._fn(web_app)(1, "", [], 1000.0, "MT1D")
        assert "Enter a model name" in msg

    def test_success_saves_and_reports(self, web_app):
        table_data = [
            {"layer": 1, "resistivity": 100.0, "thickness": 300.0}
        ]
        msg = self._fn(web_app)(1, "MyModel", table_data, 1000.0, "MT1D")
        assert "Saved 'MyModel'" in msg
        assert "MyModel" in fwd_mod._CTRL.model_names

    def test_failure_reports_error(self, web_app):
        msg = self._fn(web_app)(1, "Bad", [], 1000.0, "MT1D")
        assert msg.startswith("Save error:")
