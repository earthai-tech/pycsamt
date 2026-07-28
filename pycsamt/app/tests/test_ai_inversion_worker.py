# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for AIInversionWorker (pycsamt.app.desktop.workers.ai_inversion_worker).

Strategy
--------
Real dataset generation (``generate_dataset``) and real neural-net
training (``EMInverter1D`` / ``EMInverter2D`` / ``GCNInverter3D``) are both
too slow and non-deterministic for a unit test, so both are monkeypatched
at their source modules (they're imported locally inside each
``_run_Nd`` method, so patching the source module is enough) with
lightweight fakes that preserve the real shape contracts the worker's own
reshape/transpose arithmetic depends on:

  * fake ``generate_dataset(n_samples=..., freqs=..., n_layers=..., ...)``
    returns an object with ``.X`` shape ``(n_samples, 2*len(freqs))`` and
    ``.y`` shape ``(n_samples, 2*n_layers_upper - 1)`` -- matching what
    the real dataset object provides.
  * fake inverter classes accept the same constructor/``fit``/``predict``
    keyword contract and return a plausibly-shaped ``y_pred``.

This exercises 100% of the worker's own control flow (progress/log
signals, X_obs default-vs-provided branching, reshape math, result dict
shape) without depending on the real ML stack's runtime behavior.
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


# ── 2D ────────────────────────────────────────────────────────────────────


class Test2D:
    def test_run_2d_success(self, qapp, fake_backend):
        params = {
            "dim": "2D",
            "n_components": 2,
            "n_depth": 10,
            "n_stations": 3,
            "n_freq": 4,
            "epochs": 2,
            "batch_size": 2,
            "lr": 1e-3,
            "n_samples": 2,
            "f_min": 1e-2,
            "f_max": 1e2,
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert len(results) == 1
        res = results[0]
        assert res["dim"] == "2D"
        assert res["n_stations"] == 3
        assert res["y_pred"].shape[0] == 2  # default X_2d[:2]

    def test_run_2d_with_explicit_x_obs(self, qapp, fake_backend):
        x_obs = np.random.rand(1, 2, 4, 3)  # (n, 2, n_freqs, n_stations)
        params = {
            "dim": "2D",
            "n_stations": 3,
            "n_freq": 4,
            "n_samples": 2,
            "X_obs": x_obs,
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert results[0]["y_pred"].shape[0] == 1

    def test_run_2d_exception_reports_error(self, qapp, monkeypatch):
        def _boom(**kw):
            raise RuntimeError("2D dataset failed")

        monkeypatch.setattr("pycsamt.forward.batch.generate_dataset", _boom)
        w = AIInversionWorker({"dim": "2D"})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["2D dataset failed"]


# ── 3D ────────────────────────────────────────────────────────────────────


class Test3D:
    def test_run_3d_success(self, qapp, fake_backend):
        params = {
            "dim": "3D",
            "n_features": 10,
            "n_layers": 3,
            "hidden": [16, 8],
            "dropout": 0.1,
            "epochs": 2,
            "batch_size": 2,
            "n_samples": 2,
            "n_sta": 4,
            "radius": 1000.0,
            "f_min": 1e-2,
            "f_max": 1e1,
            "n_freq": 4,
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert len(results) == 1
        res = results[0]
        assert res["dim"] == "3D"
        assert res["coords"] is not None
        assert res["y_pred"].shape[0] == 4  # default X_3d[:n_sta]

    def test_run_3d_with_explicit_x_obs(self, qapp, fake_backend):
        x_obs = np.random.rand(2, 8)
        coords_obs = np.random.rand(2, 2)
        params = {
            "dim": "3D",
            "n_sta": 4,
            "n_samples": 2,
            "n_freq": 4,
            "X_obs": x_obs,
            "coords_obs": coords_obs,
        }
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert results[0]["y_pred"].shape[0] == 2

    def test_run_3d_exception_reports_error(self, qapp, monkeypatch):
        def _boom(**kw):
            raise RuntimeError("3D dataset failed")

        monkeypatch.setattr("pycsamt.forward.batch.generate_dataset", _boom)
        w = AIInversionWorker({"dim": "3D"})
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["3D dataset failed"]


class TestDimDispatch:
    def test_default_dim_is_1d(self, qapp, fake_backend):
        params = {"n_samples": 5, "n_layers": 3, "n_freq": 4}  # no "dim" key
        w = AIInversionWorker(params)
        results = []
        w.finished.connect(results.append)
        w.run()
        assert results[0]["dim"] == "1D"
