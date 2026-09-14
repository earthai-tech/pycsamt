# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.joint_agent`.

Cheap guard-clause tests use lightweight fake sites (bypassing
``ensure_sites``). The heavy DRCNN training path stubs ``generate_dataset``
and ``JointInverter`` (both imported locally inside ``execute()``) with
fast, deterministic fakes, so training/predict/figures/LLM in ``execute()``
run for real without a DL backend. The private helpers
(``_collect_secondary_features``, ``_extract_sec_features``,
``_forward_rms_joint``, ``_plot_joint_section``) are also tested directly.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.agents.joint_agent import (
    JointInversionAgent,
    _collect_secondary_features,
    _extract_sec_features,
    _forward_rms_joint,
    _plot_joint_section,
)

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


def _passthrough_ensure_sites(monkeypatch, sites, secondary=None):
    import pycsamt.emtools._core as core

    def _fake(x, **k):
        if secondary is not None and x is secondary:
            return secondary
        return sites

    monkeypatch.setattr(core, "ensure_sites", _fake)


def _make_sites(n=4, seed=0):
    return [
        _FakeSite(*_zblock(n=8, seed=seed + i), station=f"s{i}")
        for i in range(n)
    ]


class _FakeDataset:
    def __init__(self, n, n_freqs, n_layers, seed=0):
        rng = np.random.default_rng(seed)
        self.X = rng.normal(size=(n, n_freqs, 4)).astype(np.float32)
        self.y = rng.normal(size=(n, n_layers)).astype(np.float32)


class _FakeJointInverter:
    def __init__(self, n_features_list, n_layers, growth_rate=32):
        self.n_layers = n_layers

    def fit(self, Xs, y, *, epochs, batch_size, patience, verbose):
        return None

    def predict(self, Xs):
        n_sta = Xs[0].shape[0]
        rng = np.random.default_rng(3)
        return rng.normal(
            loc=2.0, scale=0.3, size=(n_sta, self.n_layers)
        ).astype(np.float32)


def _mock_heavy_deps(monkeypatch, *, backend_available=True, n_layers=3):
    import pycsamt.ai.inversion.joint as joint_mod
    import pycsamt.backends as backends
    import pycsamt.forward.batch as batch_mod

    monkeypatch.setattr(
        batch_mod,
        "generate_dataset",
        lambda **k: _FakeDataset(
            n=k["n_samples"], n_freqs=len(k["freqs"]), n_layers=n_layers
        ),
    )
    monkeypatch.setattr(joint_mod, "JointInverter", _FakeJointInverter)
    monkeypatch.setattr(
        backends,
        "get_backend_instance",
        lambda: object() if backend_available else None,
    )


# ── cheap guard-clause tests ─────────────────────────────────────────────────


def test_no_dl_backend_fails(monkeypatch):
    import pycsamt.backends as backends

    monkeypatch.setattr(backends, "get_backend_instance", lambda: None)
    agent = JointInversionAgent()
    result = agent.execute({"sites": [object()]})
    assert result.status == "failed"
    assert "requires PyTorch or TensorFlow" in result.error


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.backends as backends
    import pycsamt.emtools._core as core

    monkeypatch.setattr(backends, "get_backend_instance", lambda: object())
    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = JointInversionAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_secondary_path_exception_warns_and_falls_back(monkeypatch):
    sites = _make_sites(3)
    import pycsamt.backends as backends
    import pycsamt.emtools._core as core

    monkeypatch.setattr(backends, "get_backend_instance", lambda: object())

    def _fake_ensure(x, **k):
        if x == "/bad/secondary":
            raise RuntimeError("cannot read secondary")
        return sites

    monkeypatch.setattr(core, "ensure_sites", _fake_ensure)
    _mock_heavy_deps(monkeypatch, n_layers=3)
    agent = JointInversionAgent(n_layers=3, n_train_samples=8, epochs=1)
    result = agent.execute(
        {
            "sites": sites,
            "secondary_path": "/bad/secondary",
            "n_layers": 3,
        }
    )
    assert any(
        "Could not load secondary sites" in w for w in result.warnings
    )


def test_no_z_and_bad_mt_data_are_skipped(monkeypatch):
    no_z = _FakeSite(z=None, station="no_z")
    z_bad = np.ones((1, 2, 2), dtype=complex)
    fr_bad = np.array([1.0])
    bad_data = _FakeSite(z_bad, fr_bad, station="bad_data")
    sites = [no_z, bad_data]
    import pycsamt.backends as backends

    monkeypatch.setattr(backends, "get_backend_instance", lambda: object())
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = JointInversionAgent()
    result = agent.execute({"sites": sites})
    assert result.status == "failed"
    assert any("skipped (bad MT data)" in w for w in result.warnings)


def test_fewer_than_two_usable_stations_fails(monkeypatch):
    sites = _make_sites(1)
    import pycsamt.backends as backends

    monkeypatch.setattr(backends, "get_backend_instance", lambda: object())
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = JointInversionAgent()
    result = agent.execute({"sites": sites})
    assert result.status == "failed"
    assert "usable stations" in result.error


# ── full mocked happy path ───────────────────────────────────────────────────


def test_happy_path_success(monkeypatch, tmp_output):
    sites = _make_sites(4)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch, n_layers=3)
    agent = JointInversionAgent(
        api_key="fake-key",
        n_layers=3,
        n_freqs_primary=8,
        n_freqs_secondary=4,
        n_train_samples=8,
        epochs=1,
    )
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute(
        {"sites": sites, "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert result.llm_interpretation == "mocked interpretation"
    assert "joint_section" in result["figures"]
    assert result["figure_paths"]
    assert len(result["predictions"]) == 4


def test_dataset_generation_exception_fails(monkeypatch):
    sites = _make_sites(3)
    _passthrough_ensure_sites(monkeypatch, sites)
    import pycsamt.backends as backends
    import pycsamt.forward.batch as batch_mod

    monkeypatch.setattr(backends, "get_backend_instance", lambda: object())
    monkeypatch.setattr(
        batch_mod,
        "generate_dataset",
        lambda **k: (_ for _ in ()).throw(RuntimeError("gen boom")),
    )
    agent = JointInversionAgent(n_train_samples=8, epochs=1)
    result = agent.execute({"sites": sites})
    assert result.status == "failed"
    assert "Synthetic dataset generation failed" in result.error


def test_training_exception_fails(monkeypatch):
    sites = _make_sites(3)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch, n_layers=3)
    import pycsamt.ai.inversion.joint as joint_mod

    class _BoomInverter(_FakeJointInverter):
        def fit(self, *a, **k):
            raise RuntimeError("fit boom")

    monkeypatch.setattr(joint_mod, "JointInverter", _BoomInverter)
    agent = JointInversionAgent(n_layers=3, n_train_samples=8, epochs=1)
    result = agent.execute({"sites": sites})
    assert result.status == "failed"
    assert "JointInverter training failed" in result.error


def test_prediction_exception_fails(monkeypatch):
    sites = _make_sites(3)
    _passthrough_ensure_sites(monkeypatch, sites)
    _mock_heavy_deps(monkeypatch, n_layers=3)
    import pycsamt.ai.inversion.joint as joint_mod

    class _BoomInverter(_FakeJointInverter):
        def predict(self, *a, **k):
            raise RuntimeError("predict boom")

    monkeypatch.setattr(joint_mod, "JointInverter", _BoomInverter)
    agent = JointInversionAgent(n_layers=3, n_train_samples=8, epochs=1)
    result = agent.execute({"sites": sites})
    assert result.status == "failed"
    assert "Joint prediction failed" in result.error


# ── private helper direct tests ──────────────────────────────────────────────


def test_collect_secondary_features_with_secondary_sites():
    z, fr = _zblock(n=6, seed=5)
    sec_sites = [_FakeSite(z, fr, "sec1"), _FakeSite(z=None, station="sec2")]
    freqs_sec = np.logspace(-2, 1, 4)
    X_sec = _collect_secondary_features(
        ["s0", "s1"], primary_sites=[], secondary_sites=sec_sites,
        freqs_mt=np.logspace(-2, 2, 8), freqs_sec=freqs_sec,
        n_feat_sec=8, warnings=[],
    )
    assert X_sec.shape == (2, 8)


def test_collect_secondary_features_fallback_to_primary():
    z, fr = _zblock(n=6, seed=6)
    primary = [_FakeSite(z, fr, "p0"), _FakeSite(z=None, station="p1")]
    warnings: list[str] = []
    X_sec = _collect_secondary_features(
        ["p0", "p1"], primary_sites=primary, secondary_sites=None,
        freqs_mt=np.logspace(-2, 2, 8), freqs_sec=np.logspace(-2, 1, 4),
        n_feat_sec=8, warnings=warnings,
    )
    assert X_sec.shape == (2, 8)
    assert any("No secondary dataset provided" in w for w in warnings)


def test_extract_sec_features_insufficient_samples_returns_none():
    z = np.ones((1, 2, 2), dtype=complex)
    fr = np.array([1.0])
    feat = _extract_sec_features(z, fr, np.logspace(-2, 1, 4))
    assert feat is None


def test_forward_rms_joint_name_mismatch_and_no_z():
    z, fr = _zblock(n=6, seed=7)
    sites = [_FakeSite(z, fr, "match"), _FakeSite(z=None, station="other")]
    # station_idx=99 forces name-based matching only; "no_such" never matches
    rms = _forward_rms_joint(
        sites, "no_such_station", 99, np.zeros(3), np.logspace(-2, 2, 6), 4
    )
    assert rms is None


def test_forward_rms_joint_success():
    from pycsamt.forward import LayeredModel, MT1DForward

    freqs = np.logspace(-2, 2, 12)
    lm = LayeredModel(
        resistivity=np.array([100.0, 500.0, 50.0]),
        thickness=np.array([200.0, 800.0]),
    )
    resp = MT1DForward(freqs=freqs).run(lm)
    rho_xy = np.asarray(resp.rho_a)
    rho_xy = rho_xy[:, 0, 1] if rho_xy.ndim == 3 else rho_xy

    z = np.zeros((len(freqs), 2, 2), dtype=complex)
    # z[:,0,1] chosen so 0.2/f * |z|^2 reproduces rho_xy
    z[:, 0, 1] = np.sqrt(rho_xy / (0.2 / freqs))
    site = _FakeSite(z, freqs, station="s0")

    rms = _forward_rms_joint(
        [site], "s0", 0, np.log10(np.array([100.0, 500.0, 50.0])), freqs, 3
    )
    assert rms is not None
    assert rms < 0.5


def test_plot_joint_section_skips_missing_predictions():
    fig = _plot_joint_section(
        predictions={"s0": np.array([1.0, 2.0, 3.0])},
        n_layers=3,
        freqs=np.logspace(-2, 2, 8),
        station_names=["s0", "s1"],
        modalities=["mt", "tem"],
    )
    assert fig is not None
    assert fig.axes


def test_plot_joint_section_no_stations_returns_none():
    fig = _plot_joint_section(
        predictions={}, n_layers=3, freqs=np.logspace(-2, 2, 8),
        station_names=[], modalities=["mt", "tem"],
    )
    assert fig is None
