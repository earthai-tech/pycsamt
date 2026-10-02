# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for AIInversionWorker (pycsamt.app.desktop.workers.ai_inversion_worker).

Strategy
--------
``_run_1d`` still builds its own dataset/training pipeline directly, so
real dataset generation (``generate_dataset``) and real neural-net
training (``EMInverter1D``) are monkeypatched at their source modules
(imported locally inside the method, so patching the source module is
enough) with lightweight fakes that preserve the real shape contracts:

  * fake ``generate_dataset(n_samples=..., freqs=..., n_layers=..., ...)``
    returns an object with ``.X`` shape ``(n_samples, 2*len(freqs))`` and
    ``.y`` shape ``(n_samples, 2*n_layers_upper - 1)`` -- matching what
    the real dataset object provides.
  * fake inverter classes accept the same constructor/``fit``/``predict``
    keyword contract and return a plausibly-shaped ``y_pred``.

``_run_2d``/``_run_3d`` instead delegate entirely to
``pycsamt.agents.inv2d_agent.Inv2DAgent``/``inv3d_agent.Inv3DAgent``
(already-validated classes with their own dataset generation and
mt1d/mt2d/mt3d physics dispatch), so those tests fake the *agent classes*
themselves (``__init__(**kw)`` + ``execute(input_data) -> AgentResult``)
instead of the dataset/inverter primitives.

This exercises 100% of the worker's own control flow (progress/log
signals, delegation, error/status handling, result dict shape) without
depending on the real ML stack's runtime behavior.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.workers.ai_inversion_worker import AIInversionWorker


def _fake_generate_dataset(*, n_samples, freqs, n_layers, **kw):
    n_freq = len(freqs)
    nl = n_layers[1] if isinstance(n_layers, (tuple, list)) else n_layers
    X = np.random.rand(n_samples, 2 * n_freq)
    y = np.random.rand(n_samples, 2 * nl - 1)
    return SimpleNamespace(X=X, y=y)


class _FakeInverter:
    instances = []

    def __init__(self, **kw):
        self.init_kw = kw
        self._y = None
        type(self).instances.append(self)

    def fit(self, X, y, **kw):
        self.fit_kw = kw
        self._X = X
        self._y = y

    def predict(self, X_obs, **kw):
        self.predict_kw = kw
        n = len(X_obs)
        last_dim = self._y.shape[-1] if self._y.ndim > 1 else 1
        return np.zeros((n, last_dim))


@pytest.fixture
def fake_backend(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.forward.batch.generate_dataset", _fake_generate_dataset
    )
    monkeypatch.setattr("pycsamt.ai.inversion.inv1d.EMInverter1D", _FakeInverter)
    monkeypatch.setattr("pycsamt.ai.inversion.inv2d.EMInverter2D", _FakeInverter)
    monkeypatch.setattr("pycsamt.ai.inversion.inv3d.GCNInverter3D", _FakeInverter)
    _FakeInverter.instances = []
    return _FakeInverter


# ── 1D ────────────────────────────────────────────────────────────────────


class Test1D:
    def test_run_1d_success(self, qapp, fake_backend):
        params = {
            "dim": "1D",
            "arch": "resnet",
            "n_layers": 3,
            "solver": "mt1d",
            "epochs": 2,
            "batch_size": 4,
            "lr": 1e-3,
            "n_samples": 5,
            "n_freq": 4,
            "f_min": 1e-2,
            "f_max": 1e2,
            "noise_level": 0.05,
            "geology": None,
        }
        w = AIInversionWorker(params)
        results = []
        logs = []
        progresses = []
        w.finished.connect(results.append)
        w.log_line.connect(logs.append)
        w.progress.connect(progresses.append)
        w.run()

        assert len(results) == 1
        res = results[0]
        assert res["dim"] == "1D"
        assert res["y_pred"].shape[0] == 5  # default: first 5 training rows
        assert res["n_layers"] == 3
        assert progresses[-1] == 100
        assert any("complete" in line.lower() for line in logs)

    def test_run_1d_with_explicit_x_obs(self, qapp, fake_backend):
        x_obs = np.random.rand(2, 8)  # 2 stations, 2*n_freq=8
        params = {
            "dim": "1D",
            "n_layers": 3,
            "n_samples": 5,
            "n_freq": 4,
            "X_obs": x_obs,
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert results[0]["y_pred"].shape[0] == 2
        assert results[0]["X_obs"].shape == (2, 8)

    def test_run_1d_exception_reports_error(self, qapp, monkeypatch):
        def _boom(**kw):
            raise RuntimeError("dataset generation failed")

        monkeypatch.setattr("pycsamt.forward.batch.generate_dataset", _boom)
        w = AIInversionWorker({"dim": "1D"})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["dataset generation failed"]


# ── 2D / 3D ──────────────────────────────────────────────────────────────
#
# _run_2d/_run_3d delegate to pycsamt.agents.inv2d_agent.Inv2DAgent /
# inv3d_agent.Inv3DAgent (already-validated classes that own their own
# dataset generation, mt1d/mt2d(/mt3d) physics dispatch, and training) --
# not a hand-rolled dataset/training pipeline of the worker's own. Fakes
# here stand in for the *agent classes themselves* (constructor kwargs +
# ``execute(input_data) -> AgentResult``-shaped object), exercising the
# worker's own delegation, error/status handling, and result-dict shape
# without depending on real ML training.


class _FakeAgentResult:
    def __init__(self, *, status="success", data=None, warnings=None,
                 error=None, summary=""):
        self.status = status
        self.data = data or {}
        self.warnings = warnings or []
        self.error = error
        self.summary = summary

    def get(self, key, default=None):
        return self.data.get(key, default)


class _FakeAgent:
    instances = []

    def __init__(self, **kw):
        self.init_kw = kw
        type(self).instances.append(self)

    def execute(self, input_data):
        self.input_data = input_data
        return _FakeAgentResult(
            data={
                "figures": {"section": object()},
                "inverter": object(),
                "rms_global": 0.42,
                "physics": self.init_kw.get("physics", "mt1d"),
            },
            summary="fake agent run complete",
        )


class _FakeFailingAgent:
    def __init__(self, **kw):
        self.init_kw = kw

    def execute(self, input_data):
        return _FakeAgentResult(
            status="failed", error="synthetic failure", warnings=["heads up"]
        )


@pytest.fixture
def fake_agent(monkeypatch):
    _FakeAgent.instances = []
    monkeypatch.setattr(
        "pycsamt.agents.inv2d_agent.Inv2DAgent", _FakeAgent
    )
    monkeypatch.setattr(
        "pycsamt.agents.inv3d_agent.Inv3DAgent", _FakeAgent
    )
    return _FakeAgent


class Test2D:
    def test_run_2d_success(self, qapp, fake_agent):
        params = {
            "dim": "2D",
            "physics": "mt2d",
            "n_components": 2,
            "n_depth": 10,
            "n_stations": 3,
            "n_freq": 4,
            "epochs": 2,
            "n_samples": 2,
            "f_min": 1e-2,
            "f_max": 1e2,
            "sites": object(),
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert len(results) == 1
        res = results[0]
        assert res["dim"] == "2D"
        assert res["agent_result"].status == "success"
        assert fake_agent.instances[0].init_kw["physics"] == "mt2d"
        assert "sites" in fake_agent.instances[0].input_data

    def test_run_2d_without_sites_reports_error(self, qapp, monkeypatch):
        w = AIInversionWorker({"dim": "2D"})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors and "No sites loaded" in errors[0]

    def test_run_2d_agent_failure_reports_error(self, qapp, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.agents.inv2d_agent.Inv2DAgent", _FakeFailingAgent
        )
        w = AIInversionWorker({"dim": "2D", "sites": object()})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["synthetic failure"]


class Test3D:
    def test_run_3d_success(self, qapp, fake_agent):
        params = {
            "dim": "3D",
            "physics": "mt3d",
            "n_layers": 3,
            "hidden": [16, 8],
            "dropout": 0.1,
            "epochs": 2,
            "n_samples": 2,
            "radius": 1000.0,
            "f_min": 1e-2,
            "f_max": 1e1,
            "n_freq": 4,
            "sites": object(),
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert len(results) == 1
        res = results[0]
        assert res["dim"] == "3D"
        assert res["agent_result"].status == "success"
        assert fake_agent.instances[0].init_kw["physics"] == "mt3d"
        assert fake_agent.instances[0].init_kw["n_mc"] == 0

    def test_run_3d_without_sites_reports_error(self, qapp, monkeypatch):
        w = AIInversionWorker({"dim": "3D"})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors and "No sites loaded" in errors[0]

    def test_run_3d_agent_failure_reports_error(self, qapp, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.agents.inv3d_agent.Inv3DAgent", _FakeFailingAgent
        )
        w = AIInversionWorker({"dim": "3D", "sites": object()})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["synthetic failure"]


class TestDimDispatch:
    def test_default_dim_is_1d(self, qapp, fake_backend):
        params = {"n_samples": 5, "n_layers": 3, "n_freq": 4}  # no "dim" key
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert results[0]["dim"] == "1D"
