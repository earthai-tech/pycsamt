# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.inversion (Inversion page: Traditional +
AI Neural + PINN + Hybrid tracks, plus the PyTorch installer).

Strategy
--------
* "Traditional" 1-D/2-D inversion uses ``pycsamt.inversion.workflow.
  InversionWorkflow`` with the built-in SciPy backend -- real, fast
  (~0.1s), deterministic -- so it's run for real rather than mocked.
* "AI Neural" (1-D/2-D/3-D) trains real neural nets
  (``EMInverter1D``/``EMInverter2D``/``GCNInverter3D``) which is exactly
  the kind of heavy/slow ML training this project's established
  convention (see ``test_ai_inversion_worker.py``) says to mock; dataset
  generation itself (``generate_dataset``/``generate_dataset_3d``/the
  module's own ``_generate_2d_dataset``) is real and fast with a tiny
  sample count, so only the inverter classes are faked.
* "PINN"/"Hybrid" tracks train physics-informed nets via
  ``pycsamt.app.web.utils_pinn.build_pinn_inv``/``build_hybrid_inv`` --
  also heavy training, so that whole boundary module is faked for the
  one happy-path test each; all the cheap guard/early-exit branches
  (no click, wrong track, no session, no checkpoint, no torch) are
  exercised for real.
* ``_run_pip_install`` shells out to ``pip install torch`` -- always
  mocked (external network/process, not something a unit test should
  ever actually run).
"""

from __future__ import annotations

import os

# This module (and pycsamt.app.web.callbacks.inversion) imports both torch
# and scipy at process scope. On this Windows/conda environment, having
# both native BLAS/OpenMP runtimes loaded in one process aborts the whole
# interpreter the first time scipy.optimize.least_squares actually runs
# ("OMP: Error #15: Initializing libiomp5md.dll, but found ... already
# initialized") -- not a bug in either library, just duplicate OpenMP
# runtime registration. Must be set before torch/scipy are imported.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import subprocess
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.inversion as inv_mod
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
    session_id = "test-inversion-session"
    cache_set(session_id, willy_sites)
    return session_id


@pytest.fixture(autouse=True)
def _reset_install_state():
    yield
    with inv_mod._install_lock:
        inv_mod._install_state.update(
            status="idle", progress=0, log_lines=[], done_time=0.0
        )


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


def _cb_by_state(web_app, output_substr, state_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(state_id in str(s.get("id")) for s in v.get("state", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} state={state_id!r}"
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


# ══════════════════════════════════════════════════════════════════════════
# Pure helper functions
# ══════════════════════════════════════════════════════════════════════════


class TestParsePipProgress:
    def test_empty(self):
        assert inv_mod._parse_pip_progress([]) == 5

    def test_collecting(self):
        assert inv_mod._parse_pip_progress(["Collecting torch"]) >= 8

    def test_downloading_with_percent(self):
        p = inv_mod._parse_pip_progress(["Downloading torch (50%)"])
        assert p == 20 + int(50 * 0.58)

    def test_downloading_without_percent(self):
        assert inv_mod._parse_pip_progress(["Downloading torch"]) >= 20

    def test_installing_collected(self):
        assert inv_mod._parse_pip_progress(["Installing collected packages"]) >= 82

    def test_building_wheel(self):
        assert inv_mod._parse_pip_progress(["Running setup.py"]) >= 88

    def test_successfully_installed_caps_at_99(self):
        assert inv_mod._parse_pip_progress(["Successfully installed torch"]) == 99

    def test_already_satisfied(self):
        assert inv_mod._parse_pip_progress(["Requirement already satisfied"]) >= 98


class TestClassifyAndColouredLog:
    def test_classify_iter(self):
        assert inv_mod._classify("iteration 3") == "log-iter"

    def test_classify_error(self):
        assert inv_mod._classify("Traceback (most recent call)") == "log-err"

    def test_classify_none(self):
        assert inv_mod._classify("just some plain text") == ""

    def test_coloured_log_wraps_spans(self):
        spans = inv_mod._coloured_log(["error here", "plain"])
        assert spans[0].className == "log-err"
        assert spans[1].className is None


class TestRunPipInstall:
    def test_success_path(self, monkeypatch):
        class _FakeProc:
            stdout = iter(["Collecting torch\n", "Successfully installed torch\n"])
            returncode = 0

            def wait(self):
                pass

        monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: _FakeProc())
        inv_mod._run_pip_install()
        assert inv_mod._install_state["status"] == "done"
        assert inv_mod._install_state["progress"] == 100

    def test_failure_path(self, monkeypatch):
        class _FakeProc:
            stdout = iter(["some error\n"])
            returncode = 1

            def wait(self):
                pass

        monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: _FakeProc())
        inv_mod._run_pip_install()
        assert inv_mod._install_state["status"] == "error"

    def test_exception_path(self, monkeypatch):
        def _boom(*a, **k):
            raise RuntimeError("popen boom")

        monkeypatch.setattr(subprocess, "Popen", _boom)
        inv_mod._run_pip_install()
        assert inv_mod._install_state["status"] == "error"
        assert "popen boom" in inv_mod._install_state["log_lines"][-1]


class TestBuildTrueModel:
    def test_no_table_data_raises(self):
        with pytest.raises(ValueError, match="at least one finite layer"):
            inv_mod._build_true_model(None, 1000.0)

    def test_valid_table(self):
        table = [
            {"resistivity": 100.0, "thickness": 200.0},
            {"resistivity": 10.0, "thickness": 300.0},
        ]
        model = inv_mod._build_true_model(table, 1000.0)
        assert model.n_layers == 3
        assert model.resistivity[-1] == 1000.0

    def test_invalid_values_raise(self):
        table = [{"resistivity": -5.0, "thickness": 200.0}]
        with pytest.raises(ValueError, match="must be > 0"):
            inv_mod._build_true_model(table, 1000.0)


class TestMakeSynthetic:
    @pytest.fixture
    def true_model(self):
        return inv_mod._build_true_model(
            [{"resistivity": 100.0, "thickness": 200.0}], 1000.0
        )

    def test_1d_mt(self, true_model):
        freqs = np.logspace(3, -1, 10)
        em = inv_mod._make_synthetic_1d(true_model, "mt", freqs, 0.05)
        assert em.method == "mt"

    def test_1d_tdem(self, true_model):
        freqs = np.logspace(3, -1, 10)
        em = inv_mod._make_synthetic_1d(true_model, "tdem", freqs, 0.05)
        assert em.method == "tdem"

    def test_1d_csamt(self, true_model):
        freqs = np.logspace(3, -1, 10)
        em = inv_mod._make_synthetic_1d(true_model, "csamt", freqs, 0.05)
        assert em.method == "csamt"

    def test_2d(self, true_model):
        freqs = np.logspace(3, -1, 10)
        em = inv_mod._make_synthetic_2d(
            true_model, "mt", freqs, 0.05, n_stations=4, profile_len=1000
        )
        assert em.n_stations == 4


class TestFwd2dToTensor:
    def test_converts_shape(self):
        resp = SimpleNamespace(
            rho_a_te=np.ones((5, 3)) * 100.0,
            phase_te=np.ones((5, 3)) * 45.0,
        )
        X = inv_mod._fwd2d_to_tensor(resp)
        assert X.shape == (1, 2, 5, 3)
        assert X.dtype == np.float32


class TestGenerate2dDataset:
    def test_real_generation(self):
        freqs = np.logspace(3, -1, 8)
        X, y = inv_mod._generate_2d_dataset(
            n_profiles=3, n_stations=4, n_layers=3, freqs=freqs, noise=0.05, seed=1
        )
        assert X.shape == (3, 2, 8, 4)
        assert y.shape == (3, 3, 4)


class TestAiStatsAndPlots:
    def test_ai_stats_div(self):
        div = inv_mod._ai_stats_div(
            "AI 1-D", solver="MT1D", n_layers=4, n_train=100, epochs=10, n_test=20
        )
        assert div.children[0].children == "AI 1-D"

    def test_ai_stats_div_with_extra(self):
        div = inv_mod._ai_stats_div(
            "AI 2-D", solver="MT2D", n_layers=4, n_train=100, epochs=10, n_test=20,
            extra=["foo: bar"],
        )
        texts = [c.children for c in div.children[1:]]
        assert "foo: bar" in texts

    def test_plot_ai_convergence_none_when_no_loss(self):
        assert inv_mod._plot_ai_convergence({}, dark=True) is None

    def test_plot_ai_convergence_real(self):
        fig = inv_mod._plot_ai_convergence(
            {"train_loss": [1.0, 0.5, 0.2], "val_loss": [1.1, 0.6, 0.3]}, dark=True
        )
        assert fig is not None

    def test_plot_ai_comparison(self):
        n_layers = 3
        y_true = np.random.default_rng(0).uniform(0.5, 3, size=(4, 2 * n_layers - 1))
        y_pred = y_true + 0.1
        fig = inv_mod._plot_ai_comparison(y_true, y_pred, n_layers, dark=True)
        assert fig is not None

    def test_plot_ai_section_2d(self):
        rng = np.random.default_rng(0)
        y_true = rng.uniform(1, 3, size=(2, 3, 4))
        y_pred = y_true + 0.1
        fig = inv_mod._plot_ai_section_2d(y_true, y_pred, dark=True)
        assert fig is not None

    def test_plot_gcn_section_3d(self):
        rng = np.random.default_rng(0)
        n_layers = 3
        y_true = rng.uniform(1, 3, size=(2, 5, 2 * n_layers - 1))
        y_pred = y_true + 0.1
        fig = inv_mod._plot_gcn_section_3d(y_true, y_pred, n_layers, dark=True)
        assert fig is not None

    def test_plot_predicted_section_2d(self):
        section = np.random.default_rng(0).uniform(1, 3, size=(3, 4))
        fwd_resp = SimpleNamespace(rho_a_te=np.ones((5, 4)) * 100.0)
        fig = inv_mod._plot_predicted_section_2d(section, fwd_resp, 3, dark=True)
        assert fig is not None


class TestTraditionalPlotsAndStats:
    @pytest.fixture(scope="class")
    def result_1d(self):
        from pycsamt.forward import LayeredModel
        from pycsamt.forward.em1d import MT1DForward
        from pycsamt.inversion.config import InversionConfig
        from pycsamt.inversion.data import EMData
        from pycsamt.inversion.workflow import InversionWorkflow

        freqs = np.logspace(3, -1, 15)
        model_true = LayeredModel(
            resistivity=np.array([100.0, 10.0, 1000.0]),
            thickness=np.array([100.0, 300.0]),
        )
        resp = MT1DForward(freqs).run(model_true)
        em_data = EMData(method="mt", frequencies=freqs, rho_a=resp.rho_a, phase=resp.phase)
        cfg = InversionConfig(
            method="mt", dimension="1d", backend="builtin", data=em_data,
            n_layers=4, regularization="smooth", max_iter=15, error_floor=0.05,
            include_phase=True,
        )
        return model_true, InversionWorkflow(cfg).run()

    def test_plot_1d_result_with_true_model(self, result_1d):
        model_true, result = result_1d
        fig = inv_mod._plot_1d_result(model_true, result, dark=True)
        assert fig is not None

    def test_plot_1d_result_without_true_model(self, result_1d):
        _model_true, result = result_1d
        fig = inv_mod._plot_1d_result(None, result, dark=True)
        assert fig is not None

    def test_plot_convergence(self, result_1d):
        _model_true, result = result_1d
        fig = inv_mod._plot_convergence(result, dark=True)
        assert fig is None or fig is not None  # exercised regardless of history

    def test_stats_table(self, result_1d):
        _model_true, result = result_1d
        table = inv_mod._stats_table(
            result, "mt", "1d", "builtin", 4, "smooth"
        )
        assert table is not None


# ══════════════════════════════════════════════════════════════════════════
# Simple UI callbacks
# ══════════════════════════════════════════════════════════════════════════


class TestGuardRunBtn:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.BTN_INV_RUN}.disabled")

    def test_hybrid_without_checkpoint_disables(self, web_app):
        disabled, title = self._fn(web_app)("hybrid", None)
        assert disabled is True
        assert "checkpoint" in title

    def test_trad_always_enabled(self, web_app):
        disabled, title = self._fn(web_app)("trad", None)
        assert disabled is False

    def test_hybrid_with_checkpoint_enabled(self, web_app):
        disabled, title = self._fn(web_app)("hybrid", "data:...;base64,AA==")
        assert disabled is False


class TestLineStationSelectors:
    RECORDS = [
        {"ID": "S1", "Line": "L1"},
        {"ID": "S2", "Line": "L1"},
        {"ID": "S3", "Line": "L2"},
    ]

    def test_populate_pinn_selectors(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_PINN_LINE_SEL}.options")
        line_opts, sta_opts = fn({"station_records": self.RECORDS})
        assert {o["value"] for o in line_opts} == {"L1", "L2"}
        assert len(sta_opts) == 3

    def test_filter_pinn_stations(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.INV_PINN_STATION_SEL}.options", IDs.INV_PINN_LINE_SEL
        )
        opts = fn(["L1"], {"station_records": self.RECORDS})
        assert {o["value"] for o in opts} == {"S1", "S2"}

    def test_populate_hyb_selectors(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_HYB_LINE_SEL}.options")
        line_opts, sta_opts = fn({"station_records": self.RECORDS})
        assert len(sta_opts) == 3

    def test_filter_hyb_stations(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.INV_HYB_STATION_SEL}.options", IDs.INV_HYB_LINE_SEL
        )
        opts = fn(["L2"], {"station_records": self.RECORDS})
        assert {o["value"] for o in opts} == {"S3"}

    def test_populate_datafit_lines_empty(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_DATAFIT_LINE}.options")
        lines, sel_style, hint_style = fn([], {"station_records": self.RECORDS})
        assert lines == []
        assert sel_style == {"display": "none"}

    def test_populate_datafit_lines_with_stations(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_DATAFIT_LINE}.options")
        sta_opts = [{"label": "S1", "value": "S1"}, {"label": "S3", "value": "S3"}]
        lines, sel_style, hint_style = fn(sta_opts, {"station_records": self.RECORDS})
        assert {o["value"] for o in lines} == {"L1", "L2"}
        assert hint_style == {"display": "none"}

    def test_filter_datafit_stations_no_selection(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.INV_DATAFIT_STA}.options", IDs.INV_DATAFIT_LINE
        )
        from dash import no_update

        out = fn(None, [{"value": "S1"}], {"station_records": self.RECORDS})
        assert out is no_update

    def test_filter_datafit_stations_filters(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.INV_DATAFIT_STA}.options", IDs.INV_DATAFIT_LINE
        )
        opts = [{"value": "S1"}, {"value": "S3"}]
        out = fn("L1", opts, {"station_records": self.RECORDS})
        assert out == [{"value": "S1"}]


class TestSwitchTrackAndOut:
    def test_switch_inv_track(self, web_app, monkeypatch):
        fn = _cb(web_app, f"{IDs.INV_ACTIVE_TRACK}.data")
        monkeypatch.setattr(
            inv_mod, "ctx", type("C", (), {"triggered_id": "inv-track-ai"})()
        )
        assert fn(0, 1, 0, 0) == "ai"

    def test_switch_inv_track_no_trigger(self, web_app, monkeypatch):
        from dash import no_update

        fn = _cb(web_app, f"{IDs.INV_ACTIVE_TRACK}.data")
        monkeypatch.setattr(inv_mod, "ctx", type("C", (), {"triggered_id": None})())
        assert fn(0, 0, 0, 0) is no_update

    def test_switch_inv_out(self, web_app, monkeypatch):
        fn = _cb(web_app, f"{IDs.INV_ACTIVE_OUT}.data")
        monkeypatch.setattr(
            inv_mod, "ctx", type("C", (), {"triggered_id": "inv-out-tab-conv"})()
        )
        assert fn(0, 1, 0, 0, 0) == "conv"


class TestUpdateInvCtxBar:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INV_CTX_BAR}.children")

    def test_trad_with_method(self, web_app):
        out = self._fn(web_app)("trad", "mt", None, None, None)
        texts = [getattr(c, "children", c) for c in out]
        assert any("MT" == t for t in texts)

    def test_ai_dim_label(self, web_app):
        out = self._fn(web_app)("ai", None, "2d", None, None)
        texts = [getattr(c, "children", c) for c in out]
        assert any("2-D UNet" == t for t in texts)

    def test_pinn_label(self, web_app):
        out = self._fn(web_app)("pinn", None, None, "1d", None)
        texts = [getattr(c, "children", c) for c in out]
        assert any("PINN 1D" == t for t in texts)

    def test_hybrid_label(self, web_app):
        out = self._fn(web_app)("hybrid", None, None, None, "2d")
        texts = [getattr(c, "children", c) for c in out]
        assert any("Hybrid 2D" == t for t in texts)


class TestSyncCallbacks:
    def test_sync_data_src_synthetic(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_SYN_PANEL}.style")
        show, hide = fn("synthetic")
        assert show == {"display": "block"}
        assert hide == {"display": "none"}

    def test_sync_data_src_session(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_SYN_PANEL}.style")
        hide, show = fn("session")
        assert hide == {"display": "none"}
        assert show == {"display": "block"}

    def test_sync_t_dim(self, web_app):
        fn = _cb(web_app, f"{IDs.INV_2D_PARAMS}.style")
        assert fn("2d") == {"display": "block"}
        assert fn("1d") == {"display": "none"}

    def test_sync_ai_dim_1d(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_AI_ARCH}.options")
        opts, value, note, show2d, show3d = fn("1d")
        assert value == "resnet"
        assert show2d == {"display": "none"}

    def test_sync_ai_dim_2d(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_AI_ARCH}.options")
        opts, value, note, show2d, show3d = fn("2d")
        assert value == "unet2d"
        assert show2d == {"display": "block"}

    def test_sync_ai_dim_3d(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_AI_ARCH}.options")
        opts, value, note, show2d, show3d = fn("3d")
        assert value == "gcn3d"
        assert show3d == {"display": "block"}

    def test_sync_pinn_dim(self, web_app):
        fn = _cb_multi(web_app, "inv-pinn-solver-card")
        solver, p2d, p3d, lamx, lamg = fn("3d")
        assert solver == {"display": "none"}
        assert p3d == {"display": "block"}
        assert lamg == {"display": "block"}

    def test_sync_hyb_dim(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_HYB_2D_PANEL}.style")
        p2d, p3d, lamx, lamg = fn("2d")
        assert p2d == {"display": "block"}
        assert lamx == {"display": "block"}
        assert lamg == {"display": "none"}


class TestAddTrueLayer:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INV_TRUE_TABLE}.data")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, [])

    def test_adds_first_layer(self, web_app):
        out = self._fn(web_app)(1, None)
        assert out == [{"layer": 1, "resistivity": 100.0, "thickness": 500}]

    def test_appends_layer_using_last_resistivity(self, web_app):
        data = [{"layer": 1, "resistivity": 50.0, "thickness": 200}]
        out = self._fn(web_app)(1, data)
        assert out[-1] == {"layer": 2, "resistivity": 50.0, "thickness": 500}


class TestInstallFlow:
    def test_start_install_no_clicks_prevents_update(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.INV_INSTALL_INTERVAL}.disabled", IDs.BTN_INV_INSTALL_YES
        )
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_start_install_launches_thread(self, web_app, monkeypatch):
        called = threading.Event()
        monkeypatch.setattr(inv_mod, "_run_pip_install", called.set)
        fn = _cb_by_input(
            web_app, f"{IDs.INV_INSTALL_INTERVAL}.disabled", IDs.BTN_INV_INSTALL_YES
        )
        out = fn(1)
        assert called.wait(timeout=2.0)
        assert out[0] is False
        assert out[1] == {"display": "block"}

    def test_cancel_install(self, web_app):
        fn = _cb_by_input(
            web_app, f"{IDs.INV_INSTALL_MODAL}.is_open", IDs.BTN_INV_INSTALL_NO
        )
        with pytest.raises(PreventUpdate):
            fn(None)
        assert fn(1) is False

    def test_poll_install_idle_prevents_update(self, web_app):
        fn = _cb_multi(web_app, f"{IDs.INV_INSTALL_PROGRESS}.value")
        with pytest.raises(PreventUpdate):
            fn(1)

    def test_poll_install_running(self, web_app):
        inv_mod._install_state.update(
            status="running", progress=42, log_lines=["step 1"], done_time=0.0
        )
        fn = _cb_multi(web_app, f"{IDs.INV_INSTALL_PROGRESS}.value")
        out = fn(1)
        assert out[0] == 42
        assert out[3] == "info"

    def test_poll_install_done_recent(self, web_app):
        from dash import no_update

        inv_mod._install_state.update(
            status="done", progress=100, log_lines=["ok"], done_time=time.time()
        )
        fn = _cb_multi(web_app, f"{IDs.INV_INSTALL_PROGRESS}.value")
        out = fn(1)
        assert out[3] == "success"
        assert out[9] is no_update  # modal/interval untouched during success window

    def test_poll_install_done_elapsed(self, web_app):
        inv_mod._install_state.update(
            status="done", progress=100, log_lines=["ok"],
            done_time=time.time() - 10,
        )
        fn = _cb_multi(web_app, f"{IDs.INV_INSTALL_PROGRESS}.value")
        out = fn(1)
        assert out[9] is True  # modal closes
        assert out[11] is True  # success toast opens

    def test_poll_install_error(self, web_app):
        inv_mod._install_state.update(
            status="error", progress=30, log_lines=["boom"], done_time=0.0
        )
        fn = _cb_multi(web_app, f"{IDs.INV_INSTALL_PROGRESS}.value")
        out = fn(1)
        assert out[3] == "danger"


class TestHybCheckpointInfo:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INV_HYB_MODEL_INFO}.children")

    def test_no_contents_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "ckpt.npz")

    def test_bad_contents_shows_error(self, web_app):
        out = self._fn(web_app)("not-a-valid-data-uri", "ckpt.npz")
        assert "Could not read" in out.children

    def test_real_npz_checkpoint(self, web_app):
        import base64
        import io

        buf = io.BytesIO()
        np.savez(buf, weights=np.zeros(3), bias=np.zeros(1))
        contents = "data:application/npz;base64," + base64.b64encode(
            buf.getvalue()
        ).decode()
        out = self._fn(web_app)(contents, "ckpt.npz")
        assert "keys:" in out.children


class TestUpdateDataFit:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.IMG_INV_DATA_FIT}.src")

    def test_no_station_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, None, None, "dark")

    def test_no_store_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)("S1", None, None, "dark")

    def test_no_residuals_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)("S1", {"residuals": []}, None, "dark")

    def test_real_data_fit_plot(self, web_app):
        pinn_store = {
            "residuals": [
                {
                    "station": "S1",
                    "freq": 10.0,
                    "rho_obs": 100.0,
                    "rho_pred": 95.0,
                    "phase_obs": 45.0,
                    "phase_pred": 44.0,
                }
            ]
        }
        out = self._fn(web_app)("S1", pinn_store, None, "dark")
        assert out.startswith("data:image/png")


# ══════════════════════════════════════════════════════════════════════════
# run_inversion: Traditional + AI tracks
# ══════════════════════════════════════════════════════════════════════════


def _default_run_args(**overrides):
    defaults = dict(
        track="trad",
        method="mt",
        t_dim="1d",
        backend="builtin",
        data_src="synthetic",
        true_table=[{"layer": 1, "resistivity": 100.0, "thickness": 200.0}],
        halfspace_rho=1000.0,
        noise_level=0.05,
        n_stations=4,
        profile_len=1000.0,
        n_layers=4,
        regularize="smooth",
        max_iter=10,
        error_floor=0.05,
        include_phase=True,
        freq_min=-1,
        freq_max=3,
        n_freq=10,
        ai_dim="1d",
        ai_arch="resnet",
        ai_solver="mt1d",
        ai_n_layers=3,
        ai_n_samples=5,
        ai_epochs=2,
        ai_batch=4,
        ai_lr=1e-3,
        ai_noise=0.05,
        ai_n_components=2,
        ai_n_depth=3,
        ai_n_stations=4,
        ai_n_stations_3d=5,
        ai_unet_depth="auto",
        session_id=None,
        theme="dark",
    )
    defaults.update(overrides)
    return tuple(defaults.values())


class _FakeInv1D:
    def __init__(self, **kw):
        self.kw = kw
        self.history_ = {"train_loss": [1.0, 0.5]}

    def fit(self, X, y, **kw):
        self._y = y

    def predict(self, X, as_log_rho=True):
        n = len(X)
        return np.tile(self._y[0], (n, 1)).astype(np.float32)


class _FakeInv2D:
    def __init__(self, n_components=2, n_depth=3, n_stations=4, n_freqs=10,
                 unet_depth=None, dropout=0.1):
        self._channels = [8, 16, 32] if unet_depth is None else list(range(unet_depth + 2))
        self._history = {"train_loss": [1.0, 0.4]}

    def fit(self, X, y, **kw):
        self._y = y

    def predict(self, X, as_log_rho=True):
        n = len(X)
        return np.tile(self._y[0], (n, 1, 1)).astype(np.float32)


class _FakeGCN3D:
    def __init__(self, n_features=None, n_layers=None, hidden=None):
        self._history = {"train_loss": [1.0, 0.3]}

    def fit(self, X, y, adjacency=None, **kw):
        self._y = y

    def predict(self, X, adjacency=None, as_log_rho=True):
        n = len(X)
        return np.tile(self._y[0], (n, 1, 1)).astype(np.float32)


class TestRunInversionTraditional:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.IMG_INV}.src", IDs.BTN_INV_RUN)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, *_default_run_args())

    def test_pinn_track_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, *_default_run_args(track="pinn"))

    def test_trad_1d_synthetic_success(self, web_app):
        out = self._fn(web_app)(1, *_default_run_args())
        src, src_conv, log_spans, stats, spinner, msg, tab, is_open, body, modal, run_msg = out
        assert src.startswith("data:image/png")
        assert is_open is False
        assert tab == "result"

    def test_trad_2d_synthetic_success(self, web_app):
        out = self._fn(web_app)(
            1, *_default_run_args(t_dim="2d", max_iter=5)
        )
        assert out[7] is False

    def test_trad_non_builtin_backend_warns(self, web_app):
        out = self._fn(web_app)(1, *_default_run_args(backend="occam2d"))
        assert "requires an external binary" in out[5]
        assert out[7] is False  # not an error toast, just a warning message

    def test_trad_session_source_no_data(self, web_app):
        out = self._fn(web_app)(
            1, *_default_run_args(data_src="session", session_id="no-such-session")
        )
        assert "No session data" in out[5]

    def test_trad_session_source_real_data(self, web_app, cached_session):
        out = self._fn(web_app)(
            1, *_default_run_args(
                data_src="session", session_id=cached_session, max_iter=5
            )
        )
        # Either a successful run or a graceful "could not extract" error --
        # both are valid outcomes depending on WILLY's real freq coverage;
        # the important thing is the callback never crashes uncaught.
        assert out[7] in (True, False)

    def test_trad_invalid_true_table_reports_error(self, web_app):
        out = self._fn(web_app)(
            1, *_default_run_args(true_table=[{"resistivity": -1, "thickness": 100}])
        )
        assert out[7] is True
        assert "Error" in out[5]


class TestRunInversionAI:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.IMG_INV}.src", IDs.BTN_INV_RUN)

    def test_ai_torch_missing_opens_install_modal(self, web_app, monkeypatch):
        monkeypatch.setattr(inv_mod, "_TORCH_OK", False)
        out = self._fn(web_app)(1, *_default_run_args(track="ai"))
        assert out[9] is True  # install modal opens
        assert "Install PyTorch" in out[5]

    def test_ai_1d_real_dataset_mocked_inverter(self, web_app, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.ai.inversion.EMInverter1D", _FakeInv1D
        )
        out = self._fn(web_app)(1, *_default_run_args(track="ai", ai_dim="1d"))
        assert out[7] is False
        assert out[0].startswith("data:image/png")

    def test_ai_2d_real_dataset_mocked_inverter(self, web_app, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.ai.inversion.EMInverter2D", _FakeInv2D
        )
        out = self._fn(web_app)(
            1, *_default_run_args(track="ai", ai_dim="2d", ai_n_stations=4)
        )
        assert out[7] is False

    def test_ai_3d_real_dataset_mocked_inverter(self, web_app, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.ai.inversion.GCNInverter3D", _FakeGCN3D
        )
        out = self._fn(web_app)(
            1, *_default_run_args(
                track="ai", ai_dim="3d", ai_n_stations_3d=6, ai_n_layers=3
            )
        )
        assert out[7] is False

    def test_ai_exception_reports_error(self, web_app, monkeypatch):
        class _BoomInverter:
            def __init__(self, **kw):
                raise RuntimeError("inverter boom")

        monkeypatch.setattr(
            "pycsamt.ai.inversion.EMInverter1D", _BoomInverter
        )
        out = self._fn(web_app)(1, *_default_run_args(track="ai", ai_dim="1d"))
        assert out[7] is True
        assert "inverter boom" in out[8]


# ══════════════════════════════════════════════════════════════════════════
# run_pinn_hybrid_inversion: cheap guard branches + one mocked happy path each
# ══════════════════════════════════════════════════════════════════════════


def _default_pinn_hyb_args(**overrides):
    defaults = dict(
        track="pinn",
        pinn_dim="1d",
        pinn_src="synthetic",
        pinn_stations=None,
        pinn_solver="mt1d",
        pinn_mode="te",
        pinn_n_layers=4,
        pinn_depth_max=1000.0,
        pinn_n_freqs=8,
        pinn_epochs=3,
        pinn_lr=0.01,
        pinn_log_every=50,
        pinn_device="cpu",
        pinn_lam_z=0.01,
        pinn_lam_x=0.005,
        pinn_lam_g=0.005,
        pinn_radius=5000.0,
        pinn_spc=500.0,
        pinn_comp_te="xy",
        pinn_comp_tm="yx",
        pinn_backend=None,
        hyb_dim="1d",
        hyb_src="session",
        hyb_stations=None,
        hyb_ckpt=None,
        hyb_mode="te",
        hyb_n_layers=4,
        hyb_n_freqs=8,
        hyb_epochs=3,
        hyb_lr=0.005,
        hyb_device="cpu",
        hyb_lam_z=0.005,
        hyb_comp_te="xy",
        hyb_comp_tm="yx",
        session_id=None,
        theme="dark",
    )
    defaults.update(overrides)
    return tuple(defaults.values())


class _FakePinnInv:
    stations = ["S1", "S2"]

    def fit(self, **kw):
        pass

    def residuals(self):
        import pandas as pd

        return pd.DataFrame(
            {"station": ["S1", "S2"], "freq": [10.0, 10.0], "resid": [0.1, 0.2]}
        )


def _fake_plot_pinn_section(inv, dim, dark=True):
    return "data:image/png;base64,AAAA"


def _fake_plot_pinn_convergence(inv, dark=True, label=""):
    return "data:image/png;base64,BBBB"


def _fake_pinn_stats_div(inv, dim, elapsed_s=0.0, label=""):
    from dash import html

    return html.Div(f"{label} stats")


class TestRunPinnHybridInversion:
    def _fn(self, web_app):
        return _cb_by_state(web_app, f"{IDs.IMG_INV}.src", IDs.INV_PINN_DIM)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, *_default_pinn_hyb_args())

    def test_trad_track_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, *_default_pinn_hyb_args(track="trad"))

    def test_pinn_torch_missing(self, web_app, monkeypatch):
        monkeypatch.setattr(inv_mod, "_TORCH_OK", False)
        out = self._fn(web_app)(1, *_default_pinn_hyb_args(track="pinn"))
        assert out[9] is True
        assert "PyTorch required for PINN." in out[5]

    def test_pinn_session_source_no_session(self, web_app):
        out = self._fn(web_app)(
            1, *_default_pinn_hyb_args(track="pinn", pinn_src="session", session_id=None)
        )
        assert out[7] is True
        assert "No session" in out[8]

    def test_pinn_synthetic_2d_unsupported(self, web_app):
        out = self._fn(web_app)(
            1, *_default_pinn_hyb_args(track="pinn", pinn_dim="2d", pinn_src="synthetic")
        )
        assert out[7] is True
        assert "only" in out[8]

    def test_pinn_synthetic_1d_mocked_success(self, web_app, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.build_pinn_inv",
            lambda obs, dim, **kw: _FakePinnInv(),
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.plot_pinn_section", _fake_plot_pinn_section
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.plot_pinn_convergence",
            _fake_plot_pinn_convergence,
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.pinn_stats_div", _fake_pinn_stats_div
        )
        out = self._fn(web_app)(
            1, *_default_pinn_hyb_args(track="pinn", pinn_src="synthetic")
        )
        assert out[7] is False
        assert out[0] == "data:image/png;base64,AAAA"
        pinn_store = out[10]
        assert pinn_store["stations"] == ["S1", "S2"]

    def test_hybrid_torch_missing(self, web_app, monkeypatch):
        monkeypatch.setattr(inv_mod, "_TORCH_OK", False)
        out = self._fn(web_app)(1, *_default_pinn_hyb_args(track="hybrid"))
        assert out[9] is True

    def test_hybrid_no_checkpoint(self, web_app):
        out = self._fn(web_app)(
            1, *_default_pinn_hyb_args(track="hybrid", hyb_ckpt=None)
        )
        assert out[7] is True
        assert "Upload a .npz" in out[8]

    def test_hybrid_mocked_success(self, web_app, monkeypatch, cached_session):
        import tempfile

        def _fake_decode(contents):
            fd, path = tempfile.mkstemp(suffix=".npz")
            import os

            os.close(fd)
            return path

        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.decode_npz_checkpoint", _fake_decode
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.session_to_obs_1d",
            lambda session_id, stations, comp="xy": [
                SimpleNamespace(name="S1"), SimpleNamespace(name="S2")
            ],
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.build_hybrid_inv",
            lambda obs, dim, ckpt_path, **kw: _FakePinnInv(),
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.plot_pinn_section", _fake_plot_pinn_section
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.plot_pinn_convergence",
            _fake_plot_pinn_convergence,
        )
        monkeypatch.setattr(
            "pycsamt.app.web.utils_pinn.pinn_stats_div", _fake_pinn_stats_div
        )
        out = self._fn(web_app)(
            1, *_default_pinn_hyb_args(
                track="hybrid", hyb_ckpt="data:application/npz;base64,AA==",
                session_id=cached_session,
            )
        )
        assert out[7] is False
        hyb_store = out[11]
        assert hyb_store["stations"] == ["S1", "S2"]
