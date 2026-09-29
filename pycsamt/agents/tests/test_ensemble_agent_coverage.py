# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.ensemble_agent`.

Mocks the whole DL stack (``get_backend_instance``, ``generate_dataset``,
``EMInverter1D``, ``EnsembleInverter``) plus the ``emtools._core``
helpers and the module-level ``_z_to_features``/``_forward_rms``/
``_plot_uncertainty_section`` names (all bound directly on
``pycsamt.agents.ensemble_agent`` at import time, per the
"monkeypatch module-level import gotcha" — patching the defining
module would not be seen by ``execute()``) so the whole pipeline runs
in milliseconds and every branch — backend guard, dataset/training/
calibration/coverage exceptions, the per-station
skip/warn/exception loop, the ensemble-save branch, both figure
blocks, and the LLM branch — can be exercised deterministically. The
private ``_forward_rms``/``_plot_uncertainty_section`` helpers are
also unit-tested directly for their own edge branches.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pycsamt.agents.ensemble_agent as ens_mod
from pycsamt.agents.ensemble_agent import (
    EnsembleAgent,
    _forward_rms,
    _plot_uncertainty_section,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

N_FEATS = 10


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


class _FakeED:
    def __init__(self, name, z="ok", fr=None):
        self.name = name
        self.z = z
        self.fr = fr if fr is not None else np.logspace(-2, 2, 5)


class _FakeDataset:
    def __init__(self, n=200, n_feats=N_FEATS):
        rng = np.random.default_rng(0)
        self.X = rng.normal(size=(n, n_feats))
        self.y = rng.normal(size=(n, n_feats))


class _FakeEnsemble:
    def __init__(
        self,
        *,
        fit_raises=False,
        calibrate_raises=False,
        coverage_raises=False,
        raise_markers=frozenset(),
        save_raises=False,
        profile_return="fig",
    ):
        self.fit_raises = fit_raises
        self.calibrate_raises = calibrate_raises
        self.coverage_raises = coverage_raises
        self.raise_markers = raise_markers
        self.save_raises = save_raises
        self.profile_return = profile_return
        self.saved_to = None

    def fit(self, X, y, **kw):
        if self.fit_raises:
            raise RuntimeError("fit boom")

    def calibrate(self, X, y, alpha=0.10):
        if self.calibrate_raises:
            raise RuntimeError("calibrate boom")

    def coverage(self, X, y):
        if self.coverage_raises:
            raise RuntimeError("coverage boom")
        return 0.9

    def predict_with_uncertainty(self, X_in):
        marker = float(X_in[0, 0])
        if marker in self.raise_markers:
            raise RuntimeError("predict boom")
        return np.full((1, N_FEATS), 2.0), np.full((1, N_FEATS), 0.1)

    def predict_intervals(self, X_in):
        return (
            np.full((1, N_FEATS), 2.0),
            np.full((1, N_FEATS), 1.5),
            np.full((1, N_FEATS), 2.5),
        )

    def save(self, path):
        if self.save_raises:
            raise RuntimeError("save boom")
        self.saved_to = path

    def plot_uncertainty_profile(self, X_in, sample_idx=0):
        if self.profile_return == "none":
            return None
        if self.profile_return == "bad":
            return object()
        return plt.figure()


def _mock_stack(
    monkeypatch,
    *,
    backend_ok=True,
    dataset_raises=False,
    ensemble_kwargs=None,
    z_features=None,
    forward_rms_value=1.0,
    plot_unc_section="fig",
):
    """Patch the full DL/emtools stack for EnsembleAgent.execute()."""
    import pycsamt.ai.inversion.ensemble as ensemble_mod
    import pycsamt.ai.inversion.inv1d as inv1d_mod
    import pycsamt.backends as backends_mod
    import pycsamt.emtools._core as core
    import pycsamt.forward.batch as fwd_batch_mod

    monkeypatch.setattr(
        backends_mod,
        "get_backend_instance",
        lambda: (object() if backend_ok else None),
    )

    def _gen_dataset(**kw):
        if dataset_raises:
            raise RuntimeError("dataset boom")
        return _FakeDataset()

    monkeypatch.setattr(fwd_batch_mod, "generate_dataset", _gen_dataset)
    monkeypatch.setattr(inv1d_mod, "EMInverter1D", lambda **kw: object())

    ensemble = _FakeEnsemble(**(ensemble_kwargs or {}))
    monkeypatch.setattr(
        ensemble_mod, "EnsembleInverter", lambda **kw: ensemble
    )

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)
    monkeypatch.setattr(core, "_iter_items", lambda sites: list(sites))
    monkeypatch.setattr(core, "_name", lambda ed, i: ed.name)
    monkeypatch.setattr(core, "_get_z_block", lambda ed, **k: (None, ed.z, ed.fr))

    z_features = z_features or {}

    def _fake_z2f(ed, z, fr, freqs, include_phase=True):
        if ed.name in z_features:
            return z_features[ed.name]
        marker = float(abs(hash(ed.name)) % 97)
        arr = np.zeros(N_FEATS)
        arr[0] = marker
        return arr

    monkeypatch.setattr(ens_mod, "_z_to_features", _fake_z2f)

    if forward_rms_value is None:
        monkeypatch.setattr(ens_mod, "_forward_rms", lambda *a, **k: None)
    else:
        monkeypatch.setattr(
            ens_mod, "_forward_rms", lambda *a, **k: forward_rms_value
        )

    if plot_unc_section == "fig":
        monkeypatch.setattr(
            ens_mod, "_plot_uncertainty_section", lambda *a, **k: plt.figure()
        )
    elif plot_unc_section == "none":
        monkeypatch.setattr(
            ens_mod, "_plot_uncertainty_section", lambda *a, **k: None
        )
    elif plot_unc_section == "raise":
        def _raise(*a, **k):
            raise RuntimeError("unc section boom")

        monkeypatch.setattr(ens_mod, "_plot_uncertainty_section", _raise)

    return ensemble


# ── import / backend guard ──────────────────────────────────────────────────


def test_no_backend_fails(monkeypatch):
    _mock_stack(monkeypatch, backend_ok=False)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert result.status == "failed"
    assert "requires PyTorch or TensorFlow" in result.error


def test_no_sites_fails(monkeypatch):
    _mock_stack(monkeypatch)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute({})
    assert result.status == "failed"


def test_ensure_sites_exception_fails(monkeypatch):
    _mock_stack(monkeypatch)
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_dataset_generation_exception_fails(monkeypatch):
    _mock_stack(monkeypatch, dataset_raises=True)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert result.status == "failed"
    assert "Dataset generation failed" in result.error


def test_training_exception_fails(monkeypatch):
    _mock_stack(monkeypatch, ensemble_kwargs={"fit_raises": True})
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert result.status == "failed"
    assert "Ensemble training failed" in result.error


# ── calibration / coverage branches ──────────────────────────────────────


def test_calibrate_disabled_skips_block(monkeypatch):
    _mock_stack(monkeypatch)
    result = EnsembleAgent(
        n_train_samples=10, epochs=1, calibrate=False
    ).execute({"sites": [_FakeED("S1")]})
    assert result.status in {"success", "needs_review"}
    assert not any("Calibration failed" in w for w in result.warnings)


def test_calibrate_exception_warns(monkeypatch):
    _mock_stack(monkeypatch, ensemble_kwargs={"calibrate_raises": True})
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert any("Calibration failed" in w for w in result.warnings)


def test_coverage_exception_warns(monkeypatch):
    _mock_stack(monkeypatch, ensemble_kwargs={"coverage_raises": True})
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert any("Coverage computation" in w for w in result.warnings)
    import math

    assert math.isnan(result["coverage"])


# ── per-station prediction loop branches ─────────────────────────────────


def test_per_station_branches(monkeypatch):
    ed_ok = _FakeED("S_OK")
    ed_zerofail = _FakeED("S_ZNONE", z=None)
    ed_nofeat = _FakeED("S_NOFEAT")
    ed_predfail = _FakeED("S_PREDFAIL")

    z_features = {
        "S_OK": np.concatenate([[1.0], np.zeros(N_FEATS - 1)]),
        "S_NOFEAT": None,
        "S_PREDFAIL": np.concatenate([[2.0], np.zeros(N_FEATS - 1)]),
    }
    ensemble = _mock_stack(
        monkeypatch,
        z_features=z_features,
        ensemble_kwargs={"raise_markers": {2.0}},
    )
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [ed_ok, ed_zerofail, ed_nofeat, ed_predfail]}
    )
    assert result.status == "success"
    assert "S_OK" in result["pred_mean"]
    assert "S_ZNONE" not in result["pred_mean"]
    assert "S_NOFEAT" not in result["pred_mean"]
    assert "S_PREDFAIL" not in result["pred_mean"]
    assert any("could not build feature vector" in w for w in result.warnings)
    assert any("Prediction for S_PREDFAIL" in w for w in result.warnings)


def test_forward_rms_none_skips_rms_list(monkeypatch):
    _mock_stack(monkeypatch, forward_rms_value=None)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert result.status == "success"
    import math

    assert math.isnan(result["rms_global"])


def test_no_predictions_skips_figures_and_llm(monkeypatch):
    ed = _FakeED("S1", z=None)
    _mock_stack(monkeypatch)
    agent = EnsembleAgent(n_train_samples=10, epochs=1, api_key="fake-key")
    agent.query_llm = lambda *a, **k: "should not be called"
    result = agent.execute({"sites": [ed]})
    assert result.status == "needs_review"
    assert result["figures"] == {}
    assert result.llm_interpretation is None


# ── ensemble-save branch ──────────────────────────────────────────────────


def test_save_ensemble_success(monkeypatch, tmp_output):
    ensemble = _mock_stack(monkeypatch)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")], "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert ensemble.saved_to is not None


def test_save_ensemble_exception_warns(monkeypatch, tmp_output):
    _mock_stack(monkeypatch, ensemble_kwargs={"save_raises": True})
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")], "output_dir": str(tmp_output)}
    )
    assert any("Could not save ensemble" in w for w in result.warnings)


# ── uncertainty-section figure branches ──────────────────────────────────


def test_uncertainty_section_none_skips_block(monkeypatch):
    _mock_stack(monkeypatch, plot_unc_section="none")
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert "uncertainty_section" not in result["figures"]


def test_uncertainty_section_exception_warns(monkeypatch):
    _mock_stack(monkeypatch, plot_unc_section="raise")
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert any("Uncertainty section figure" in w for w in result.warnings)


def test_uncertainty_section_saved_with_output_dir(monkeypatch, tmp_output):
    _mock_stack(monkeypatch)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")], "output_dir": str(tmp_output)}
    )
    assert result["figure_paths"].get("uncertainty_section")


def test_uncertainty_section_no_output_dir_skips_path(monkeypatch):
    _mock_stack(monkeypatch)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert "uncertainty_section" in result["figures"]
    assert "uncertainty_section" not in result["figure_paths"]


# ── uncertainty-profile figure branches ──────────────────────────────────


def test_uncertainty_profile_none_and_bad_skip(monkeypatch):
    _mock_stack(monkeypatch, ensemble_kwargs={"profile_return": "none"})
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert "uncertainty_profile" not in result["figures"]

    _mock_stack(monkeypatch, ensemble_kwargs={"profile_return": "bad"})
    result2 = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert "uncertainty_profile" not in result2["figures"]


def test_uncertainty_profile_saved_with_output_dir(monkeypatch, tmp_output):
    _mock_stack(monkeypatch)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")], "output_dir": str(tmp_output)}
    )
    assert result["figure_paths"].get("uncertainty_profile")


def test_uncertainty_profile_first_station_not_found(monkeypatch):
    """When the second `_iter_items` pass can't re-locate the first
    predicted station (e.g. a non-deterministic naming scheme), X_first
    stays None and the profile figure is skipped without error."""
    _mock_stack(monkeypatch)
    import pycsamt.emtools._core as core

    call = {"n": 0}

    def _flaky_name(ed, i):
        call["n"] += 1
        # first pass (prediction loop) resolves real names; second pass
        # (profile lookup) never matches any of them.
        return ed.name if call["n"] <= 1 else "NEVER_MATCHES"

    monkeypatch.setattr(core, "_name", _flaky_name)
    result = EnsembleAgent(n_train_samples=10, epochs=1).execute(
        {"sites": [_FakeED("S1")]}
    )
    assert result.status == "success"
    assert "uncertainty_profile" not in result["figures"]


# ── LLM interpretation branch ──────────────────────────────────────────────


def test_llm_interpretation_branch(monkeypatch):
    _mock_stack(monkeypatch)
    agent = EnsembleAgent(
        n_train_samples=10, epochs=1, api_key="fake-key"
    )
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute({"sites": [_FakeED("S1")]})
    assert result.llm_interpretation == "mocked interpretation"


# ── private helper: _forward_rms ────────────────────────────────────────


def test_forward_rms_returns_none_when_z_missing():
    ed = _FakeED("S1", z=None)
    assert _forward_rms(ed, np.zeros(3), np.logspace(-2, 2, 5), 3) is None


def test_forward_rms_returns_none_when_too_few_finite_points():
    freqs = np.logspace(-2, 2, 5)
    z = np.zeros((5, 2, 2), dtype=complex)
    z[:, 0, 1] = 0.0  # rho_obs all zero -> mask always False
    ed = _FakeED("S1", z=z, fr=freqs)
    assert _forward_rms(ed, np.array([2.0, 2.0]), freqs, 3) is None


def test_forward_rms_swallows_exception(monkeypatch):
    import pycsamt.forward as fwd_mod

    class _Boom:
        def __init__(self, *a, **k):
            raise RuntimeError("layered model boom")

    monkeypatch.setattr(fwd_mod, "LayeredModel", _Boom)
    freqs = np.logspace(-2, 2, 5)
    z = np.ones((5, 2, 2), dtype=complex)
    ed = _FakeED("S1", z=z, fr=freqs)
    assert _forward_rms(ed, np.array([2.0, 2.0, 2.0]), freqs, 3) is None


# ── private helper: _plot_uncertainty_section ───────────────────────────


def test_plot_uncertainty_section_empty_returns_none():
    assert _plot_uncertainty_section({}, {}, {}, {}, 3, np.logspace(-2, 2, 5)) is None


def test_plot_uncertainty_section_missing_std_for_station():
    freqs = np.logspace(-2, 2, 5)
    pred_mean = {"S1": np.array([1.0, 2.0, 3.0]), "S2": np.array([1.0, 2.0, 3.0])}
    pred_std = {"S1": np.array([0.1, 0.1, 0.1])}  # S2 missing -> stays NaN
    fig = _plot_uncertainty_section(pred_mean, pred_std, {}, {}, 3, freqs)
    assert fig is not None
