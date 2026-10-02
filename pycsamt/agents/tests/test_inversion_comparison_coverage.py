# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.inversion_comparison`.

Exercises branches not reached by the synthetic happy-path test in
``test_agents_offline_battery.py``: the ``predictions``-dict extraction
path, the "no recognisable section key" fallback, explicit
``depths_km``/``station_names`` overrides, the tiny-section /
``corrcoef``-failure correlation branches, the figure-saving success
and failure paths, the LLM interpretation branch, and the private
``_extract_section`` / ``_get_depths`` / ``_plot_comparison`` helpers
directly.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.agents.inversion_comparison import (
    InversionComparisonAgent,
    _extract_section,
    _get_depths,
    _plot_comparison,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _synthetic_pair(n_layers=3, n_sta=8, seed=0):
    rng = np.random.default_rng(seed)
    depths = np.linspace(0.05, 2.0, n_layers)
    a = {
        "pred_rho": 10 ** rng.normal(2, 0.3, (n_layers, n_sta)),
        "depths_km": depths,
    }
    b = {
        "pred_rho": 10 ** rng.normal(2, 0.3, (n_layers, n_sta)),
        "depths_km": depths,
    }
    return a, b


# ── _extract_section direct-call branches ──────────────────────────────────


def test_extract_section_none_result_returns_none_tuple():
    mat, stations, depths = _extract_section(None, [], "X")
    assert mat is None and stations is None and depths is None


def test_extract_section_skips_non_matching_and_1d_keys():
    warnings: list[str] = []
    # "pred_rho" present but 1-D -> ndim check fails -> loop continues;
    # no other 2-D key or "predictions" dict -> falls through to warning.
    result = {"pred_rho": np.array([1.0, 2.0, 3.0])}
    mat, stations, depths = _extract_section(result, warnings, "Model A")
    assert mat is None
    assert any("no recognisable section key found" in w for w in warnings)


def test_extract_section_predictions_dict_path():
    preds = {"S1": [1.0, 2.0, 3.0], "S2": [1.0, 2.0]}
    result = {"predictions": preds, "depths_km": [0.1, 0.5, 1.0]}
    mat, stations, depths = _extract_section(result, [], "Model A")
    assert mat.shape == (3, 2)
    assert stations == ["S1", "S2"]
    assert np.isnan(mat[2, 1])  # S2 only has 2 layers
    assert depths is not None


# ── _get_depths ──────────────────────────────────────────────────────────


def test_get_depths_returns_none_when_no_key_matches():
    assert _get_depths(lambda k, d=None: None) is None


def test_get_depths_finds_alternate_key():
    assert _get_depths(lambda k, d=None: [1, 2] if k == "depths" else None) is not None


# ── agent-level branch coverage ─────────────────────────────────────────────


def test_missing_extractable_section_fails():
    agent = InversionComparisonAgent()
    result = agent.execute({"result_a": {}, "result_b": {"pred_rho": [[1.0]]}})
    assert result.status == "failed"
    assert "Model A" in result.error

    result_b_bad = agent.execute({"result_a": {"pred_rho": [[1.0]]}, "result_b": {}})
    assert result_b_bad.status == "failed"
    assert "Model B" in result_b_bad.error


def test_explicit_depths_and_station_names_override(tmp_path):
    a, b = _synthetic_pair()
    agent = InversionComparisonAgent()
    result = agent.execute(
        {
            "result_a": a,
            "result_b": b,
            "depths_km": [0.0, 1.0, 2.0],
            "station_names": [f"ST{i}" for i in range(8)],
            "output_dir": str(tmp_path),
        }
    )
    assert result.status == "success"
    assert result["station_names"][0] == "ST0"


def test_default_depths_when_neither_side_has_them():
    a = {"pred_rho": np.full((3, 4), 10.0)}
    b = {"pred_rho": np.full((3, 4), 12.0)}
    agent = InversionComparisonAgent()
    result = agent.execute({"result_a": a, "result_b": b})
    assert result.status == "success"
    assert len(result["depths_km"]) == 4  # n_layers + 1 default linspace


def test_tiny_section_skips_correlation_block():
    a = {"pred_rho": np.array([[1.0]]), "depths_km": [0.0, 1.0]}
    b = {"pred_rho": np.array([[1.0]]), "depths_km": [0.0, 1.0]}
    agent = InversionComparisonAgent()
    result = agent.execute({"result_a": a, "result_b": b})
    assert result.status == "success"
    assert np.isnan(result["correlation"])
    assert np.isnan(result["rmse"])


def test_corrcoef_exception_is_swallowed(monkeypatch):
    import pycsamt.agents.inversion_comparison as mod

    a, b = _synthetic_pair()

    def _raise(*a_, **k_):
        raise RuntimeError("corrcoef boom")

    monkeypatch.setattr(mod.np, "corrcoef", _raise)
    agent = InversionComparisonAgent()
    result = agent.execute({"result_a": a, "result_b": b})
    assert result.status == "success"
    assert np.isnan(result["correlation"])


def test_figure_saved_when_output_dir_given(tmp_path):
    a, b = _synthetic_pair()
    agent = InversionComparisonAgent()
    result = agent.execute(
        {"result_a": a, "result_b": b, "output_dir": str(tmp_path)}
    )
    assert result.status == "success"
    assert result["figure_paths"].get("comparison")


def test_plot_comparison_none_return_skips_figure_block(monkeypatch):
    import pycsamt.agents.inversion_comparison as mod

    a, b = _synthetic_pair()
    monkeypatch.setattr(mod, "_plot_comparison", lambda *a_, **k_: None)
    agent = InversionComparisonAgent()
    result = agent.execute({"result_a": a, "result_b": b})
    assert result.status == "success"
    assert result["figures"] == {}


def test_plot_comparison_exception_is_recorded_as_warning(monkeypatch):
    import pycsamt.agents.inversion_comparison as mod

    a, b = _synthetic_pair()

    def _raise(*a_, **k_):
        raise RuntimeError("plot boom")

    monkeypatch.setattr(mod, "_plot_comparison", _raise)
    agent = InversionComparisonAgent()
    result = agent.execute({"result_a": a, "result_b": b})
    assert result.status == "success"
    assert any("Comparison figure" in w for w in result.warnings)


def test_llm_interpretation_branch():
    a, b = _synthetic_pair()
    agent = InversionComparisonAgent(api_key="fake-key")
    agent.query_llm = lambda *a_, **k_: "mocked interpretation"
    result = agent.execute({"result_a": a, "result_b": b})
    assert result.llm_interpretation == "mocked interpretation"


# ── _plot_comparison direct call: short depths axis ─────────────────────────


def test_plot_comparison_pads_short_depths_axis():
    n_layers, n_sta = 4, 3
    mat_a = np.random.default_rng(0).normal(2, 0.2, (n_layers, n_sta))
    mat_b = np.random.default_rng(1).normal(2, 0.2, (n_layers, n_sta))
    diff = mat_a - mat_b
    fig = _plot_comparison(
        mat_a,
        mat_b,
        diff,
        depths_km=np.array([0.0, 1.0]),  # shorter than n_layers + 1
        station_names=[f"S{i}" for i in range(n_sta)],
        label_a="A",
        label_b="B",
        corr=0.5,
        rmse=0.1,
    )
    assert fig is not None
