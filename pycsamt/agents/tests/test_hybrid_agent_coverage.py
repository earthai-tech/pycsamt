# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.hybrid_agent`.

Bypasses the (heavy, DL-backend-dependent) real ``HybridInverterXD``
fitting by stubbing ``HybridInversionAgent._run`` directly, so the
``execute()``-level guard clauses, figure/LLM/predict branches can be
exercised deterministically and cheaply: the invalid-``dim`` guard, the
no-DL-backend guard, the missing-sites guard, the ``ensure_sites``
exception, the ``_run`` exception, the three figure exception guards, the
LLM interpretation path, and the ``predict()``/``stage1_models()``
exception guards.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.agents.hybrid_agent import HybridInversionAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


class _FakeInv:
    def __init__(
        self,
        n_sites=2,
        n_layers=4,
        predict_ok=True,
        stage1_ok=True,
    ):
        self.stations = [f"s{i}" for i in range(n_sites)]
        self.n_sites = n_sites
        self.n_layers = n_layers
        self.depth_max = 1000.0
        self._predict_ok = predict_ok
        self._stage1_ok = stage1_ok

    def predict(self):
        if not self._predict_ok:
            raise RuntimeError("predict boom")
        return ["model"] * self.n_sites

    def stage1_models(self):
        if not self._stage1_ok:
            raise RuntimeError("stage1_models boom")
        return ["s1_model"] * self.n_sites


def _mock_backend(monkeypatch, available=True):
    import pycsamt.backends as backends

    monkeypatch.setattr(
        backends,
        "get_backend_instance",
        lambda: object() if available else None,
    )


def _stub_run(monkeypatch, agent, *, mat=None, s1_mat=None, conv_df=None, inv=None):
    inv = inv or _FakeInv()
    monkeypatch.setattr(
        agent, "_run", lambda *a, **k: (inv, mat, s1_mat, conv_df, None, None)
    )
    return inv


# ── __init__ / guard-clause tests ────────────────────────────────────────────


def test_invalid_dim_raises():
    with pytest.raises(ValueError, match="dim must be 1, 2, or 3"):
        HybridInversionAgent(dim=4)


def test_no_dl_backend_fails(monkeypatch):
    _mock_backend(monkeypatch, available=False)
    agent = HybridInversionAgent()
    result = agent.execute(
        {"ai_inverter": object(), "sites": object()}
    )
    assert result.status == "failed"
    assert "requires PyTorch or TensorFlow" in result.error


def test_no_sites_fails(monkeypatch):
    _mock_backend(monkeypatch, available=True)
    agent = HybridInversionAgent()
    result = agent.execute({"ai_inverter": object()})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = HybridInversionAgent()
    result = agent.execute(
        {"ai_inverter": object(), "sites": object()}
    )
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_run_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)
    agent = HybridInversionAgent()
    monkeypatch.setattr(
        agent,
        "_run",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("fit boom")),
    )
    result = agent.execute(
        {"ai_inverter": object(), "sites": object()}
    )
    assert result.status == "failed"
    assert "Hybrid fitting failed" in result.error


# ── figure / LLM / predict branch tests (via a stubbed _run) ────────────────


def _base_execute(monkeypatch, agent, **stub_kwargs):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)
    inv = _stub_run(monkeypatch, agent, **stub_kwargs)
    result = agent.execute(
        {"ai_inverter": object(), "sites": object()}
    )
    return result, inv


def test_stage2_section_plot_exception(monkeypatch):
    import pycsamt.agents.hybrid_agent as ha

    monkeypatch.setattr(
        ha,
        "_plot_pinn_section",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("s2 plot boom")),
    )
    agent = HybridInversionAgent()
    mat = np.random.default_rng(0).normal(size=(4, 2))
    result, _ = _base_execute(monkeypatch, agent, mat=mat)
    assert any("Stage-2 section plot" in w for w in result.warnings)


def test_stage1_section_plot_exception(monkeypatch):
    import pycsamt.agents.hybrid_agent as ha

    calls = {"n": 0}

    def _maybe_boom(*a, **k):
        calls["n"] += 1
        if calls["n"] == 1:
            return None  # stage-2 plot: skip cleanly
        raise RuntimeError("s1 plot boom")

    monkeypatch.setattr(ha, "_plot_pinn_section", _maybe_boom)
    agent = HybridInversionAgent()
    rng = np.random.default_rng(1)
    mat = rng.normal(size=(4, 2))
    s1_mat = rng.normal(size=(4, 2))
    result, _ = _base_execute(monkeypatch, agent, mat=mat, s1_mat=s1_mat)
    assert any("Stage-1 section plot" in w for w in result.warnings)


def test_convergence_plot_exception(monkeypatch):
    import pandas as pd

    import pycsamt.agents.hybrid_agent as ha

    monkeypatch.setattr(
        ha,
        "_plot_loss_curves",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("conv plot boom")),
    )
    agent = HybridInversionAgent()
    conv_df = pd.DataFrame({"iter": [1, 2], "loss": [1.0, 0.5]})
    result, _ = _base_execute(monkeypatch, agent, conv_df=conv_df)
    assert any("Convergence plot" in w for w in result.warnings)


def test_llm_interpretation_called_when_api_key_set(monkeypatch):
    agent = HybridInversionAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    mat = np.random.default_rng(2).normal(size=(4, 2))
    result, _ = _base_execute(monkeypatch, agent, mat=mat)
    assert result.llm_interpretation == "mocked interpretation"


def test_predict_exception_is_recorded(monkeypatch):
    agent = HybridInversionAgent(dim=1)
    inv = _FakeInv(predict_ok=False)
    result, _ = _base_execute(monkeypatch, agent, inv=inv)
    assert any("predict():" in w for w in result.warnings)


def test_stage1_models_exception_is_recorded(monkeypatch):
    agent = HybridInversionAgent(dim=1)
    inv = _FakeInv(stage1_ok=False)
    result, _ = _base_execute(monkeypatch, agent, inv=inv)
    assert any("stage1_models():" in w for w in result.warnings)
