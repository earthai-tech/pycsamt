# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.inversion_backend`.

Mocks ``InversionConfig``/``run_inversion``/``available_backends``/
``get_backend`` (all imported locally inside ``execute()``) so the whole
pipeline runs deterministically without a real physics backend: the
``pycsamt.inversion`` import guard, ``ensure_sites`` exception, the
backend-fallback branches (unregistered name / registered-but-unavailable),
the ``InversionConfig`` construction exception, the ``log_rho_section`` /
``station_names`` extraction exception guards, the two figure exception
guards, and the convergence-figure block.
"""

from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from pycsamt.agents.inversion_backend import InversionBackendAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


class _FakeInvResult:
    def __init__(
        self,
        rms=0.05,
        n_iter=12,
        model=None,
        station_names=None,
        history=None,
    ):
        self.rms = rms
        self.n_iter = n_iter
        self.model = model
        self.station_names = station_names
        self.history = history


def _mock_inversion(
    monkeypatch,
    *,
    result=None,
    config_raises=False,
    run_raises=False,
):
    import pycsamt.inversion as inv_pkg

    monkeypatch.setattr(inv_pkg, "available_backends", lambda: ["builtin"])

    if config_raises:
        monkeypatch.setattr(
            inv_pkg,
            "InversionConfig",
            lambda **k: (_ for _ in ()).throw(ValueError("bad config")),
        )
    else:
        monkeypatch.setattr(inv_pkg, "InversionConfig", lambda **k: k)

    if run_raises:
        monkeypatch.setattr(
            inv_pkg,
            "run_inversion",
            lambda cfg: (_ for _ in ()).throw(RuntimeError("solve boom")),
        )
    else:
        monkeypatch.setattr(
            inv_pkg, "run_inversion", lambda cfg: result or _FakeInvResult()
        )


def _passthrough_ensure_sites(monkeypatch, sites=object()):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)
    return sites


def test_inversion_import_error_fails(monkeypatch):
    fake_mod = types.ModuleType("pycsamt.inversion")
    monkeypatch.setitem(sys.modules, "pycsamt.inversion", fake_mod)
    agent = InversionBackendAgent()
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "pycsamt.inversion not available" in result.error


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = InversionBackendAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_unregistered_backend_falls_back_to_builtin(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(monkeypatch)
    agent = InversionBackendAgent(backend="not_a_real_backend")
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("not in registry" in w for w in result.warnings)
    assert result["backend"] == "builtin"


def test_registered_but_unavailable_backend_falls_back(
    monkeypatch, tmp_output
):
    import pycsamt.inversion as inv_pkg
    import pycsamt.inversion.backends as backends_mod

    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(monkeypatch)
    monkeypatch.setattr(
        inv_pkg, "available_backends", lambda: ["builtin", "simpeg"]
    )
    monkeypatch.setattr(
        backends_mod,
        "get_backend",
        lambda name: (_ for _ in ()).throw(ImportError("no simpeg")),
    )
    agent = InversionBackendAgent(backend="simpeg")
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("registered but unavailable" in w for w in result.warnings)
    assert result["backend"] == "builtin"


def test_inversion_config_exception_fails(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(monkeypatch, config_raises=True)
    agent = InversionBackendAgent()
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "InversionConfig construction failed" in result.error


def test_run_inversion_exception_fails(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(monkeypatch, run_raises=True)
    agent = InversionBackendAgent()
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "Inversion failed" in result.error


def test_log_rho_section_extraction_exception_is_recorded(
    monkeypatch, tmp_output
):
    class _BoomModel:
        def get(self, *a, **k):
            raise RuntimeError("model boom")

    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(
        monkeypatch,
        result=_FakeInvResult(model=_BoomModel(), station_names=["s1"]),
    )
    agent = InversionBackendAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any(
        "Could not extract log_rho_section" in w for w in result.warnings
    )


def test_station_names_extraction_exception_is_swallowed(
    monkeypatch, tmp_output
):
    class _BoomResult(_FakeInvResult):
        @property
        def station_names(self):
            raise RuntimeError("names boom")

        @station_names.setter
        def station_names(self, value):
            pass

    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(monkeypatch, result=_BoomResult())
    agent = InversionBackendAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result["station_names"] is None


def test_inversion_section_figure_exception_is_recorded(
    monkeypatch, tmp_output
):
    import pycsamt.agents.inversion_backend as ib

    monkeypatch.setattr(
        ib,
        "_plot_inversion_section",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("plot boom")),
    )
    log_rho = np.random.default_rng(0).normal(size=(4, 3))
    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(
        monkeypatch,
        result=_FakeInvResult(model={"log_rho_section": log_rho}),
    )
    agent = InversionBackendAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("Inversion section figure" in w for w in result.warnings)


def test_convergence_figure_success_and_llm(monkeypatch, tmp_output):
    import matplotlib.pyplot as plt

    import pycsamt.ai.plot.convergence as conv_mod

    fig_marker = plt.figure()
    monkeypatch.setattr(
        conv_mod, "plot_convergence", lambda history: fig_marker
    )
    log_rho = np.random.default_rng(1).normal(size=(4, 3))
    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(
        monkeypatch,
        result=_FakeInvResult(
            model={"log_rho_section": log_rho},
            history={"rms": [1.0, 0.5, 0.1]},
        ),
    )
    agent = InversionBackendAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result["figures"]["convergence"] is fig_marker
    assert result.llm_interpretation == "mocked interpretation"


def test_convergence_figure_exception_is_recorded(monkeypatch, tmp_output):
    import pycsamt.ai.plot.convergence as conv_mod

    monkeypatch.setattr(
        conv_mod,
        "plot_convergence",
        lambda history: (_ for _ in ()).throw(RuntimeError("conv boom")),
    )
    log_rho = np.random.default_rng(2).normal(size=(4, 3))
    _passthrough_ensure_sites(monkeypatch)
    _mock_inversion(
        monkeypatch,
        result=_FakeInvResult(
            model={"log_rho_section": log_rho}, history={"rms": [1.0]}
        ),
    )
    agent = InversionBackendAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("Convergence figure" in w for w in result.warnings)
