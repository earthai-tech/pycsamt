# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.inv3d_agent`.

Cheap guard-clause tests use lightweight fake sites (bypassing
``ensure_sites``). The heavier ``mt1d``-physics happy path stubs
``generate_dataset``, ``build_adjacency``, and ``GCNInverter3D`` (all
imported locally inside ``execute()``) with fast, deterministic fakes so
the real training/prediction/figure/topography/LLM code in ``execute()``
runs for real without needing a DL backend or real Maxwell solves. The
``physics="mt3d"`` branch (real 3-D Maxwell training data) is left to the
existing, much slower inversion-battery tests.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.agents.inv3d_agent import Inv3DAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


class _FakeZ:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _FakeSite:
    def __init__(self, z=None, freq=None, station="st"):
        self.station = station
        if z is not None:
            self.Z = _FakeZ(z, freq)


def _zblock(n=8, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    freq = np.linspace(1.0, 100.0, n)
    return z, freq


def _passthrough_ensure_sites(monkeypatch, sites):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def _make_sites(n=6, seed=0):
    sites = []
    for i in range(n):
        z, fr = _zblock(n=8, seed=seed + i)
        sites.append(_FakeSite(z, fr, station=f"s{i}"))
    return sites


class _FakeDataset:
    def __init__(self, n, n_freqs, n_layers, seed=0):
        rng = np.random.default_rng(seed)
        self.X = rng.normal(size=(n, n_freqs, 4)).astype(np.float32)
        self.y = rng.normal(size=(n, n_layers)).astype(np.float32)


class _FakeInverter:
    def __init__(self, n_features, n_layers, hidden=None, dropout=0.1):
        self.n_features = n_features
        self.n_layers = n_layers
        self._history = {"val_loss": [1.0, 0.5, 0.3]}
        self._meta = {"best_val_loss": 0.3}
        self._n_out = 2 * n_layers - 1

    def fit(self, X, y, *, adjacency, epochs, batch_size, patience, verbose):
        return None

    def predict(self, X_obs, *, adjacency):
        n_sta = X_obs.shape[0]
        rng = np.random.default_rng(1)
        return rng.normal(
            loc=2.0, scale=0.3, size=(n_sta, self._n_out)
        ).astype(np.float32)

    def predict_with_uncertainty(self, X_obs, *, adjacency, n_mc):
        n_sta = X_obs.shape[0]
        rng = np.random.default_rng(2)
        mu = rng.normal(size=(n_sta, self._n_out)).astype(np.float32)
        sigma = np.abs(rng.normal(size=(n_sta, self._n_out))).astype(
            np.float32
        )
        return mu, sigma


def _mock_heavy_deps(monkeypatch, *, backend_available=True, n_layers=3):
    import pycsamt.ai.inversion.inv3d as inv3d_mod
    import pycsamt.ai.nets.gcn as gcn_mod
    import pycsamt.backends as backends
    import pycsamt.forward.batch as batch_mod

    monkeypatch.setattr(
        batch_mod,
        "generate_dataset",
        lambda **k: _FakeDataset(
            n=k["n_samples"], n_freqs=len(k["freqs"]), n_layers=n_layers
        ),
    )
    monkeypatch.setattr(
        gcn_mod,
        "build_adjacency",
        lambda coords, radius: np.eye(len(coords), dtype=np.float32),
    )
    monkeypatch.setattr(inv3d_mod, "GCNInverter3D", _FakeInverter)
    monkeypatch.setattr(
        backends,
        "get_backend_instance",
        lambda: object() if backend_available else None,
    )


# ── cheap guard-clause tests ─────────────────────────────────────────────────


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = Inv3DAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_invalid_freqs_fails(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv3DAgent()
    result = agent.execute({"sites": sites, "freqs": [1.0]})
    assert result.status == "failed"
    assert "freqs" in result.error


def test_invalid_depth_max_fails(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv3DAgent()
    result = agent.execute({"sites": sites, "depth_max": -5.0})
    assert result.status == "failed"
    assert "depth_max" in result.error


def test_no_z_and_bad_data_stations_are_skipped(monkeypatch):
    no_z = _FakeSite(z=None, station="no_z")
    z_bad = np.ones((1, 2, 2), dtype=complex)
    fr_bad = np.array([1.0])
    bad_data = _FakeSite(z_bad, fr_bad, station="bad_data")
    sites = [no_z, bad_data]
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv3DAgent()
    result = agent.execute(
        {"sites": sites, "freqs": np.logspace(-2, 2, 8)}
    )
    assert result.status == "failed"
    assert any("no Z data" in w for w in result.warnings)
    assert any("feature extraction failed" in w for w in result.warnings)


def test_fewer_than_two_usable_stations_fails(monkeypatch):
    sites = _make_sites(1)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv3DAgent()
    result = agent.execute(
        {"sites": sites, "freqs": np.logspace(-2, 2, 8)}
    )
    assert result.status == "failed"
    assert "usable station" in result.error


def test_adjacency_shape_mismatch_warns_and_rebuilds(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch)
    agent = Inv3DAgent(n_layers=3, n_train_profiles=2, epochs=1, n_mc=0)
    bad_adjacency = np.ones((2, 2))  # wrong shape for 4 stations
    result = agent.execute(
        {
            "sites": sites,
            "freqs": np.logspace(-2, 2, 8),
            "adjacency": bad_adjacency,
            "topography": False,
        }
    )
    assert any("rebuilding from coordinates" in w for w in result.warnings)


def test_no_dl_backend_fails_after_dataset_generation(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch, backend_available=False)
    agent = Inv3DAgent(n_layers=3, n_train_profiles=2, epochs=1)
    result = agent.execute(
        {"sites": sites, "freqs": np.logspace(-2, 2, 8)}
    )
    assert result.status == "failed"
    assert "requires PyTorch or TensorFlow" in result.error


def test_training_exception_fails(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch)
    import pycsamt.ai.inversion.inv3d as inv3d_mod

    class _BoomInverter(_FakeInverter):
        def fit(self, *a, **k):
            raise RuntimeError("gcn fit boom")

    monkeypatch.setattr(inv3d_mod, "GCNInverter3D", _BoomInverter)
    agent = Inv3DAgent(n_layers=3, n_train_profiles=2, epochs=1)
    result = agent.execute(
        {"sites": sites, "freqs": np.logspace(-2, 2, 8)}
    )
    assert result.status == "failed"
    assert "GCNInverter3D training failed" in result.error


def test_prediction_exception_fails(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch)
    import pycsamt.ai.inversion.inv3d as inv3d_mod

    class _BoomInverter(_FakeInverter):
        def predict(self, *a, **k):
            raise RuntimeError("predict boom")

    monkeypatch.setattr(inv3d_mod, "GCNInverter3D", _BoomInverter)
    agent = Inv3DAgent(n_layers=3, n_train_profiles=2, epochs=1)
    result = agent.execute(
        {"sites": sites, "freqs": np.logspace(-2, 2, 8)}
    )
    assert result.status == "failed"
    assert "3-D prediction failed" in result.error


# ── full mocked happy path ───────────────────────────────────────────────────


def test_happy_path_with_uncertainty_and_topography(monkeypatch, tmp_output):
    sites = _make_sites(6)
    coords = np.array(
        [[0.0, 0.0], [1000.0, 0.0], [2000.0, 500.0],
         [500.0, 1500.0], [1500.0, 1800.0], [2500.0, 900.0]]
    )
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch, n_layers=3)
    agent = Inv3DAgent(
        api_key="fake-key",
        n_layers=3,
        n_train_profiles=2,
        epochs=1,
        n_mc=5,
    )
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute(
        {
            "sites": sites,
            "freqs": np.logspace(-2, 2, 8),
            "coords": coords,
            "output_dir": str(tmp_output),
            "topography": {
                "elevation_m": [10.0, 12.0, 20.0, 15.0, 30.0, 25.0],
                "chainage_km": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            },
        }
    )
    assert result.status == "success"
    assert result["pred_uncertainty"] is not None
    assert result["topography"]["applied"] is True
    assert result.llm_interpretation == "mocked interpretation"
    assert "depth_slices" in result["figures"]
    assert "resistivity_section" in result["figures"]
    assert "uncertainty_map" in result["figures"]
    assert result["figure_paths"]


def test_mc_dropout_uncertainty_exception_is_recorded(monkeypatch):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch, n_layers=3)
    import pycsamt.ai.inversion.inv3d as inv3d_mod

    class _BoomUncertainty(_FakeInverter):
        def predict_with_uncertainty(self, *a, **k):
            raise RuntimeError("mc boom")

    monkeypatch.setattr(inv3d_mod, "GCNInverter3D", _BoomUncertainty)
    agent = Inv3DAgent(n_layers=3, n_train_profiles=2, epochs=1, n_mc=5)
    result = agent.execute(
        {
            "sites": sites,
            "freqs": np.logspace(-2, 2, 8),
            "topography": False,
        }
    )
    assert any("MC-dropout uncertainty failed" in w for w in result.warnings)
    assert result["pred_uncertainty"] is None
