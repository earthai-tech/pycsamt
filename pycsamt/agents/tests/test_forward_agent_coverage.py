# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.forward`.

``test_forward.py`` only exercises the 1-D happy path with default
figures. This file adds: the ``pycsamt.forward`` import guard, the
unsupported-``dim`` guard, real (small, fast) 2-D and 3-D forward runs
with the LLM interpretation path, 1-D exception branches (solve
failure, ρa/phase extraction, RMS computation success/failure), the
``_build_layered_model`` / ``_build_grid_2d`` / ``_build_grid_3d``
construction-failure and model-type branches, the 1-D/2-D/3-D figure
exception guards, and ``_compute_rms_1d``'s direct branches.
"""

from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from pycsamt.agents.forward import ForwardModelAgent, _compute_rms_1d

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _small_2d_model(**overrides):
    cfg = {
        "type": "halfspace",
        "bg_rho": 100.0,
        "nx": 8,
        "nz": 6,
        "x_max": 4000.0,
        "z_max": 2000.0,
        "n_stations": 4,
    }
    cfg.update(overrides)
    return cfg


def _small_3d_model(**overrides):
    cfg = {
        "type": "halfspace",
        "bg_rho": 100.0,
        "nx": 6,
        "ny": 6,
        "nz": 5,
        "x_max": 3000.0,
        "y_max": 3000.0,
        "z_max": 1500.0,
        "nx_stations": 3,
        "ny_stations": 3,
    }
    cfg.update(overrides)
    return cfg


_FREQS = list(np.logspace(-2, 2, 4))


# ── top-level guards ──────────────────────────────────────────────────────────


def test_forward_import_error_fails(monkeypatch):
    fake_mod = types.ModuleType("pycsamt.forward")
    monkeypatch.setitem(sys.modules, "pycsamt.forward", fake_mod)
    agent = ForwardModelAgent()
    result = agent.execute({})
    assert result.status == "failed"
    assert "pycsamt.forward is not available" in result.error


def test_unsupported_dim_fails():
    agent = ForwardModelAgent(dim=4)
    result = agent.execute({})
    assert result.status == "failed"
    assert "not supported" in result.error


# ── real 2-D and 3-D runs ─────────────────────────────────────────────────────


def test_2d_forward_real_run_with_llm(tmp_output):
    agent = ForwardModelAgent(dim=2, api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked 2-D interpretation"
    result = agent.execute(
        {
            "model": _small_2d_model(),
            "freqs": _FREQS,
            "output_dir": str(tmp_output),
        }
    )
    assert result.status == "success"
    assert result["dim"] == 2
    assert result.llm_interpretation == "mocked 2-D interpretation"
    assert result["figures"]
    assert result["figure_paths"]


def test_3d_forward_real_run_with_llm(tmp_output):
    agent = ForwardModelAgent(dim=3, api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked 3-D interpretation"
    result = agent.execute(
        {
            "model": _small_3d_model(),
            "freqs": _FREQS,
            "output_dir": str(tmp_output),
        }
    )
    assert result.status == "success"
    assert result["dim"] == 3
    assert result.llm_interpretation == "mocked 3-D interpretation"
    assert result["figures"]
    assert result["figure_paths"]


def test_2d_grid_anomaly_type():
    agent = ForwardModelAgent(dim=2)
    result = agent.execute(
        {
            "model": _small_2d_model(type="anomaly", anomaly_rho=1.0),
            "freqs": _FREQS,
        }
    )
    assert result.status == "success"


def test_2d_grid_from_1d_layers():
    agent = ForwardModelAgent(dim=2)
    result = agent.execute(
        {
            "model": {
                "resistivities": [200.0, 20.0, 5000.0],
                "thicknesses": [300.0, 800.0],
                "nx": 8,
                "x_max": 3000.0,
                "n_stations": 4,
            },
            "freqs": _FREQS,
        }
    )
    assert result.status == "success"


def test_2d_grid_unknown_type_warns_and_falls_back():
    agent = ForwardModelAgent(dim=2)
    result = agent.execute(
        {
            "model": _small_2d_model(type="not_a_real_type"),
            "freqs": _FREQS,
        }
    )
    assert result.status == "success"
    assert any("Unknown model type" in w for w in result.warnings)


def test_2d_prebuilt_grid_passthrough():
    from pycsamt.forward import Grid2D

    grid = Grid2D.halfspace(**{
        k: v
        for k, v in _small_2d_model().items()
        if k in ("nx", "nz", "x_max", "z_max", "n_stations")
    } | {"rho": 100.0})
    agent = ForwardModelAgent(dim=2)
    result = agent.execute({"grid": grid, "freqs": _FREQS})
    assert result.status == "success"
    assert result["layered_model"] is None


def test_3d_grid_block_anomaly_type():
    agent = ForwardModelAgent(dim=3)
    result = agent.execute(
        {
            "model": _small_3d_model(type="block_anomaly", anomaly_rho=1.0),
            "freqs": _FREQS,
        }
    )
    assert result.status == "success"


def test_3d_grid_unknown_type_warns_and_falls_back():
    agent = ForwardModelAgent(dim=3)
    result = agent.execute(
        {
            "model": _small_3d_model(type="not_a_real_type"),
            "freqs": _FREQS,
        }
    )
    assert result.status == "success"
    assert any("Unknown 3-D model type" in w for w in result.warnings)


def test_2d_grid_construction_exception_fails():
    agent = ForwardModelAgent(dim=2)
    result = agent.execute(
        {"model": _small_2d_model(nx=0), "freqs": _FREQS}
    )
    assert result.status == "failed"
    assert "Grid2D construction failed" in result.error


def test_3d_grid_construction_exception_fails():
    agent = ForwardModelAgent(dim=3)
    result = agent.execute(
        {"model": _small_3d_model(nx=0), "freqs": _FREQS}
    )
    assert result.status == "failed"
    assert "Grid3D construction failed" in result.error


def test_2d_solve_exception_fails(monkeypatch):
    import pycsamt.forward as fwd_pkg

    monkeypatch.setattr(
        fwd_pkg,
        "MT2DForward",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("2d solve boom")),
    )
    agent = ForwardModelAgent(dim=2)
    result = agent.execute(
        {"model": _small_2d_model(), "freqs": _FREQS}
    )
    assert result.status == "failed"
    assert "2-D forward failed" in result.error


def test_3d_solve_exception_fails(monkeypatch):
    import pycsamt.forward as fwd_pkg

    monkeypatch.setattr(
        fwd_pkg,
        "MT3DForward",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("3d solve boom")),
    )
    agent = ForwardModelAgent(dim=3)
    result = agent.execute(
        {"model": _small_3d_model(), "freqs": _FREQS}
    )
    assert result.status == "failed"
    assert "3-D forward failed" in result.error


# ── 1-D exception branches ───────────────────────────────────────────────────


def test_1d_model_construction_exception_fails():
    agent = ForwardModelAgent(dim=1)
    result = agent.execute(
        {"model": {"resistivity": [100.0, -20.0], "thickness": [500.0]}}
    )
    assert result.status == "failed"
    assert "Could not build LayeredModel" in result.error


def test_1d_solve_exception_fails(monkeypatch):
    import pycsamt.forward as fwd_pkg

    monkeypatch.setattr(
        fwd_pkg,
        "MT1DForward",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("1d solve boom")),
    )
    agent = ForwardModelAgent(dim=1)
    result = agent.execute(
        {"model": {"resistivity": [100.0, 20.0], "thickness": [500.0]}}
    )
    assert result.status == "failed"
    assert "1-D forward failed" in result.error


def test_1d_rho_phase_extraction_exception_is_recorded(monkeypatch):
    import pycsamt.forward as fwd_pkg

    class _BadResponse:
        @property
        def rho_a(self):
            raise RuntimeError("rho_a boom")

        phase = np.zeros(4)

    class _FakeSolver:
        def __init__(self, *a, **k):
            pass

        def run(self, layered):
            return _BadResponse()

    monkeypatch.setattr(fwd_pkg, "MT1DForward", _FakeSolver)
    agent = ForwardModelAgent(dim=1)
    result = agent.execute(
        {"model": {"resistivity": [100.0, 20.0], "thickness": [500.0]}}
    )
    assert result.status == "success"
    assert any("Could not extract" in w for w in result.warnings)
    assert result["rho_a"] is None


def test_1d_rms_computed_with_real_sites(edi_dir):
    agent = ForwardModelAgent(dim=1)
    result = agent.execute(
        {
            "model": {"resistivity": [100.0, 20.0], "thickness": [500.0]},
            "sites": str(edi_dir),
        }
    )
    assert result.status == "success"
    assert "RMS" in result.summary


def test_1d_rms_exception_is_recorded(monkeypatch, edi_dir):
    import pycsamt.agents.forward as fwd_agent

    monkeypatch.setattr(
        fwd_agent,
        "_compute_rms_1d",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("rms boom")),
    )
    agent = ForwardModelAgent(dim=1)
    result = agent.execute(
        {
            "model": {"resistivity": [100.0, 20.0], "thickness": [500.0]},
            "sites": str(edi_dir),
        }
    )
    assert result.status == "success"
    assert any("RMS computation failed" in w for w in result.warnings)
    assert result["rms"] is None


def test_1d_figure_fallback_after_combined_plot_fails(monkeypatch):
    import pycsamt.forward as fwd_pkg

    monkeypatch.setattr(
        fwd_pkg,
        "plot_response_and_model_1d",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("combined boom")),
    )
    agent = ForwardModelAgent(dim=1)
    result = agent.execute(
        {"model": {"resistivity": [100.0, 20.0], "thickness": [500.0]}}
    )
    assert result.status == "success"
    assert any("plot_response_and_model_1d" in w for w in result.warnings)
    # the individual response/model plots still get produced as a fallback
    assert "response" in result["figures"] or "model" in result["figures"]


def test_2d_figure_exception_is_recorded(monkeypatch):
    import pycsamt.forward as fwd_pkg

    def _boom_plot_model_2d(*a, **k):
        raise RuntimeError("2d plot boom")

    _boom_plot_model_2d.__name__ = "plot_model_2d"
    monkeypatch.setattr(fwd_pkg, "plot_model_2d", _boom_plot_model_2d)
    agent = ForwardModelAgent(dim=2)
    result = agent.execute(
        {"model": _small_2d_model(), "freqs": _FREQS}
    )
    assert any(
        "plot_model_2d(model_2d)" in w and "2d plot boom" in w
        for w in result.warnings
    )


def test_3d_figure_exception_is_recorded(monkeypatch):
    import pycsamt.forward as fwd_pkg

    def _boom_plot_model_3d(*a, **k):
        raise RuntimeError("3d plot boom")

    _boom_plot_model_3d.__name__ = "plot_model_3d"
    monkeypatch.setattr(fwd_pkg, "plot_model_3d", _boom_plot_model_3d)
    agent = ForwardModelAgent(dim=3)
    result = agent.execute(
        {"model": _small_3d_model(), "freqs": _FREQS}
    )
    assert any(
        "plot_model_3d(model_3d)" in w and "3d plot boom" in w
        for w in result.warnings
    )


# ── _compute_rms_1d direct tests ─────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _FakeSite:
    def __init__(self, z=None, freq=None, station="st"):
        self.station = station
        if z is not None:
            self.Z = _FakeZ(z, freq)


def _fake_response(n=6, seed=0):
    rng = np.random.default_rng(seed)
    return types.SimpleNamespace(
        rho_a=np.abs(rng.normal(loc=100.0, scale=10.0, size=n))
    )


def _passthrough_ensure_sites(monkeypatch, sites):
    # `_compute_rms_1d` calls `ensure_sites` on its `sites_raw` argument;
    # the real implementation converts an unrecognised list into an empty
    # Sites() container, silently dropping the fake site doubles below.
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)


def test_compute_rms_1d_skips_station_with_no_z(monkeypatch):
    freqs = np.logspace(-2, 2, 6)
    site = _FakeSite(z=None, station="no_z")
    _passthrough_ensure_sites(monkeypatch, [site])
    with pytest.raises(ValueError, match="No valid observed data"):
        _compute_rms_1d([site], _fake_response(6), freqs, 0, 1)


def test_compute_rms_1d_skips_station_with_no_finite_mask(monkeypatch):
    freqs = np.logspace(-2, 2, 6)
    z = np.zeros((6, 2, 2), dtype=complex)  # |z|=0 -> rho_obs == 0 -> mask False
    site = _FakeSite(z, freqs, station="zeroz")
    _passthrough_ensure_sites(monkeypatch, [site])
    with pytest.raises(ValueError, match="No valid observed data"):
        _compute_rms_1d([site], _fake_response(6), freqs, 0, 1)


def test_compute_rms_1d_success(monkeypatch):
    freqs = np.logspace(-2, 2, 6)
    rng = np.random.default_rng(3)
    z = rng.normal(size=(6, 2, 2)) + 1j * rng.normal(size=(6, 2, 2))
    site = _FakeSite(z, freqs, station="s1")
    _passthrough_ensure_sites(monkeypatch, [site])
    rms = _compute_rms_1d([site], _fake_response(6, seed=1), freqs, 0, 1)
    assert rms >= 0.0
