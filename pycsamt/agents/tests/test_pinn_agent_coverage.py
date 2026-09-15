# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.pinn_agent`.

Mirrors the ``hybrid_agent`` coverage approach: stubs
``PINNInversionAgent._run`` (and, for the dim dispatcher, ``_run_2d`` /
``_run_3d``) to avoid the heavy, DL-backend-dependent real
``PINNInverterXD`` fitting, so the ``execute()``-level guard clauses,
figure/LLM/predict branches, and the dim dispatcher can be exercised
deterministically. Also directly tests the ``_plot_pinn_section`` /
``_plot_loss_curves`` no-data early returns.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.agents.pinn_agent import (
    PINNInversionAgent,
    _plot_loss_curves,
    _plot_pinn_section,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


class _FakeInv:
    def __init__(self, n_sites=2, predict_ok=True):
        self.stations = [f"s{i}" for i in range(n_sites)]
        self.n_sites = n_sites
        self._predict_ok = predict_ok

    def predict(self):
        if not self._predict_ok:
            raise RuntimeError("predict boom")
        return ["model"] * self.n_sites


def _mock_backend(monkeypatch, available=True):
    import pycsamt.backends as backends

    monkeypatch.setattr(
        backends,
        "get_backend_instance",
        lambda: object() if available else None,
    )


def _stub_run(monkeypatch, agent, *, mat=None, loss_df=None, inv=None):
    inv = inv or _FakeInv()
    monkeypatch.setattr(
        agent, "_run", lambda *a, **k: (inv, mat, loss_df, None)
    )
    return inv


def _base_execute(monkeypatch, agent, **stub_kwargs):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)
    inv = _stub_run(monkeypatch, agent, **stub_kwargs)
    result = agent.execute({"sites": object()})
    return result, inv


# ── __init__ / guard-clause tests ────────────────────────────────────────────


def test_invalid_dim_raises():
    with pytest.raises(ValueError, match="dim must be 1, 2, or 3"):
        PINNInversionAgent(dim=5)


def test_no_dl_backend_fails(monkeypatch):
    _mock_backend(monkeypatch, available=False)
    agent = PINNInversionAgent()
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "requires PyTorch or TensorFlow" in result.error


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = PINNInversionAgent()
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_run_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)
    agent = PINNInversionAgent()
    monkeypatch.setattr(
        agent,
        "_run",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("fit boom")),
    )
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "PINN fitting failed" in result.error


def test_run_dispatches_to_2d_and_3d(monkeypatch):
    agent = PINNInversionAgent()
    monkeypatch.setattr(
        agent, "_run_2d", lambda *a, **k: ("ran_2d",)
    )
    monkeypatch.setattr(
        agent, "_run_3d", lambda *a, **k: ("ran_3d",)
    )
    assert agent._run(2, None, 10, 2000.0, 100, []) == ("ran_2d",)
    assert agent._run(3, None, 10, 2000.0, 100, []) == ("ran_3d",)


# ── figure / LLM / predict branch tests (via a stubbed _run) ────────────────


def test_section_plot_exception(monkeypatch):
    import pycsamt.agents.pinn_agent as pa

    monkeypatch.setattr(
        pa,
        "_plot_pinn_section",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("plot boom")),
    )
    agent = PINNInversionAgent(n_layers=4)
    mat = np.random.default_rng(0).normal(size=(4, 2))
    result, _ = _base_execute(monkeypatch, agent, mat=mat)
    assert any("Section plot" in w for w in result.warnings)


def test_convergence_plot_exception(monkeypatch):
    import pandas as pd

    import pycsamt.agents.pinn_agent as pa

    monkeypatch.setattr(
        pa,
        "_plot_loss_curves",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("conv boom")),
    )
    agent = PINNInversionAgent()
    loss_df = pd.DataFrame({"epoch": [1, 2], "loss": [1.0, 0.5]})
    result, _ = _base_execute(monkeypatch, agent, loss_df=loss_df)
    assert any("Convergence plot" in w for w in result.warnings)


def test_llm_interpretation_called_when_api_key_set(monkeypatch):
    agent = PINNInversionAgent(api_key="fake-key", n_layers=4)
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    mat = np.random.default_rng(1).normal(size=(4, 2))
    result, _ = _base_execute(monkeypatch, agent, mat=mat)
    assert result.llm_interpretation == "mocked interpretation"


def test_predict_exception_is_recorded(monkeypatch):
    agent = PINNInversionAgent(dim=1)
    inv = _FakeInv(predict_ok=False)
    result, _ = _base_execute(monkeypatch, agent, inv=inv)
    assert any("predict():" in w for w in result.warnings)


# ── _plot_pinn_section / _plot_loss_curves direct tests ─────────────────────


def test_plot_pinn_section_no_stations_returns_none():
    fig = _plot_pinn_section(
        mat=np.zeros((3, 0)),
        station_names=[],
        n_layers=3,
        depths_km=np.array([0.0, 1.0, 2.0, 3.0]),
    )
    assert fig is None


def test_plot_loss_curves_no_data_returns_none():
    assert _plot_loss_curves(None) is None

    import pandas as pd

    assert _plot_loss_curves(pd.DataFrame()) is None
