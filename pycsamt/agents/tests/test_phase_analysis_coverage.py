# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.phase_analysis`.

``PhaseAnalysisAgent.execute`` chains eight independent
try/except-guarded steps (PT table, dimensionality classification,
strike consensus, and five figures), each with its own
None/empty/exception branch. This module mocks all eight functions
(imported locally inside ``execute()``) with small, fast fakes so
every branch can be exercised deterministically, without needing real
EDI geometry or matplotlib-heavy real plotting logic.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pycsamt.agents.phase_analysis import PhaseAnalysisAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _passthrough_ensure_sites(monkeypatch, sites=object()):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def _mock_all(
    monkeypatch,
    *,
    pt_raises=False,
    pt_empty=False,
    pt_missing_cols=False,
    dim_raises=False,
    strike_result="ok",
    psection_raises=False,
    psection_none=False,
    rose_raises=False,
    rose_none=False,
    rose_bad=False,
    strike_fig_raises=False,
    strike_fig_none=False,
    strike_fig_bad=False,
    dimgrid_raises=False,
    dimgrid_none=False,
    fingerprint_raises=False,
    mohr_raises=False,
):
    import pycsamt.emtools.advanced as advanced_mod
    import pycsamt.emtools.dimensionality as dim_mod
    import pycsamt.emtools.strike as strike_mod
    import pycsamt.emtools.tensor as tensor_mod

    def _pt_table(sites, verbose=0):
        if pt_raises:
            raise RuntimeError("pt boom")
        if pt_empty:
            return pd.DataFrame(columns=["station", "period", "beta", "ellipt"])
        if pt_missing_cols:
            return pd.DataFrame({"station": ["S1"], "period": [1.0]})
        return pd.DataFrame(
            {
                "station": ["S1", "S1"],
                "period": [1.0, 10.0],
                "beta": [1.0, 6.0],
                "ellipt": [0.05, 0.5],
            }
        )

    monkeypatch.setattr(tensor_mod, "build_phase_tensor_table", _pt_table)

    def _dim_classify(sites, skew_th=None, ellipt_th=None, verbose=0):
        if dim_raises:
            raise RuntimeError("dim boom")
        return pd.DataFrame({"station": ["S1"], "cls": ["1D"]})

    monkeypatch.setattr(dim_mod, "classify_dimensionality", _dim_classify)

    def _strike(sites, band=None, verbose=0):
        if strike_result == "raise":
            raise RuntimeError("strike boom")
        if strike_result == "no_ang":
            return {"ang": [1.0]}  # plain dict: hasattr(obj, "ang") is False
        if strike_result == "empty_ang":
            return pd.DataFrame(
                {"ang": pd.Series(dtype=float), "iqr": pd.Series(dtype=float)}
            )
        return pd.DataFrame({"ang": [10.0, 20.0], "iqr": [2.0, 3.0]})

    monkeypatch.setattr(strike_mod, "estimate_strike_consensus", _strike)

    def _plot_psection(sites, period_range=None, figsize=None, verbose=0):
        if psection_raises:
            raise RuntimeError("psection boom")
        if psection_none:
            return None
        return plt.figure()

    monkeypatch.setattr(tensor_mod, "plot_phase_tensor_psection", _plot_psection)

    def _plot_rose(sites, band=None, verbose=0):
        if rose_raises:
            raise RuntimeError("rose boom")
        if rose_none:
            return None
        if rose_bad:
            return object()
        return plt.figure()

    monkeypatch.setattr(tensor_mod, "plot_phase_tensor_rose", _plot_rose)

    def _plot_strike_analysis(sites, band=None, verbose=0):
        if strike_fig_raises:
            raise RuntimeError("strike fig boom")
        if strike_fig_none:
            return None
        if strike_fig_bad:
            return object()
        return plt.figure()

    monkeypatch.setattr(strike_mod, "plot_strike_analysis", _plot_strike_analysis)

    def _plot_dimgrid(sites, skew_th=None, ellipt_th=None, verbose=0):
        if dimgrid_raises:
            raise RuntimeError("dimgrid boom")
        if dimgrid_none:
            return None
        return plt.figure()

    monkeypatch.setattr(dim_mod, "plot_dim_confidence_grid", _plot_dimgrid)

    def _plot_fp(sites, quantities=None, period_range=None):
        if fingerprint_raises:
            raise RuntimeError("fp boom")
        return plt.figure()

    monkeypatch.setattr(advanced_mod, "plot_survey_fingerprint", _plot_fp)

    def _plot_mohr(sites, n_periods=8, verbose=0):
        if mohr_raises:
            raise RuntimeError("mohr boom")
        return plt.figure()

    monkeypatch.setattr(advanced_mod, "plot_impedance_mohr_circles", _plot_mohr)


# ── ensure_sites exception ──────────────────────────────────────────────────


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    result = PhaseAnalysisAgent().execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


# ── PT table / dimensionality classification branches ──────────────────────


def test_pt_and_dim_exceptions_skip_class_counts(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, pt_raises=True, dim_raises=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result.status == "success"
    assert any("build_phase_tensor_table" in w for w in result.warnings)
    assert any("classify_dimensionality" in w for w in result.warnings)
    assert result["pt_table"] is None
    assert (result["n_1d"], result["n_2d"], result["n_3d"]) == (0, 0, 0)


def test_pt_table_missing_columns_swallowed(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, pt_missing_cols=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result.status == "success"
    assert (result["n_1d"], result["n_2d"], result["n_3d"]) == (0, 0, 0)


def test_pt_table_empty_skips_class_counts(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, pt_empty=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result.status == "success"
    assert (result["n_1d"], result["n_2d"], result["n_3d"]) == (0, 0, 0)


def test_class_counts_computed_from_real_table(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result.status == "success"
    assert sum((result["n_1d"], result["n_2d"], result["n_3d"])) == 2


# ── strike consensus branches ───────────────────────────────────────────────


def test_strike_no_ang_attribute_skips_consensus(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, strike_result="no_ang")
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result.status == "success"
    import math

    assert math.isnan(result["strike_consensus"])


def test_strike_empty_ang_skips_consensus(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, strike_result="empty_ang")
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result.status == "success"
    import math

    assert math.isnan(result["strike_consensus"])


def test_strike_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, strike_result="raise")
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert any("estimate_strike_consensus" in w for w in result.warnings)


def test_strike_consensus_computed(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert result["strike_consensus"] == pytest.approx(15.0, abs=1.0)


# ── figure branches: pt_psection, pt_rose, strike_analysis, dim_confidence ──


def test_figures_no_output_dir_skip_paths(monkeypatch):
    """Without output_dir, _save_figure returns None for every figure ->
    every "if p:" branch takes its False arm."""
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    result = PhaseAnalysisAgent().execute(
        {"sites": object(), "run_mohr": True}
    )
    assert result.status == "success"
    assert result["figure_paths"] == {}
    # psection, rose, strike, dim, fingerprint (default on) + mohr (forced on)
    assert len(result["figures"]) == 6


def test_psection_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, psection_raises=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert any("plot_phase_tensor_psection" in w for w in result.warnings)


def test_rose_none_and_bad_object_skip_figure(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, rose_none=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert "pt_rose" not in result["figures"]

    _mock_all(monkeypatch, rose_bad=True)
    result2 = PhaseAnalysisAgent().execute({"sites": object()})
    assert "pt_rose" not in result2["figures"]


def test_rose_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, rose_raises=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert any("plot_phase_tensor_rose" in w for w in result.warnings)


def test_strike_analysis_none_and_bad_object_skip_figure(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, strike_fig_none=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert "strike_analysis" not in result["figures"]

    _mock_all(monkeypatch, strike_fig_bad=True)
    result2 = PhaseAnalysisAgent().execute({"sites": object()})
    assert "strike_analysis" not in result2["figures"]


def test_strike_analysis_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, strike_fig_raises=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert any("plot_strike_analysis" in w for w in result.warnings)


def test_dim_confidence_none_skips_figure(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, dimgrid_none=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert "dim_confidence" not in result["figures"]


def test_dim_confidence_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, dimgrid_raises=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert any("plot_dim_confidence_grid" in w for w in result.warnings)


# ── fingerprint / Mohr optional branches ────────────────────────────────────


def test_fingerprint_disabled_is_skipped(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    result = PhaseAnalysisAgent().execute(
        {"sites": object(), "run_fingerprint": False}
    )
    assert "survey_fingerprint" not in result["figures"]


def test_fingerprint_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, fingerprint_raises=True)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert any("plot_survey_fingerprint" in w for w in result.warnings)


def test_mohr_disabled_by_default(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    result = PhaseAnalysisAgent().execute({"sites": object()})
    assert "mohr_circles" not in result["figures"]


def test_mohr_enabled_and_saved(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    result = PhaseAnalysisAgent().execute(
        {
            "sites": object(),
            "run_mohr": True,
            "output_dir": str(tmp_output),
        }
    )
    assert result["figure_paths"].get("mohr_circles")


def test_mohr_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch, mohr_raises=True)
    result = PhaseAnalysisAgent().execute(
        {"sites": object(), "run_mohr": True}
    )
    assert any("plot_impedance_mohr_circles" in w for w in result.warnings)


# ── LLM interpretation branch ────────────────────────────────────────────────


def test_llm_interpretation_branch(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_all(monkeypatch)
    agent = PhaseAnalysisAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute({"sites": object()})
    assert result.llm_interpretation == "mocked interpretation"
