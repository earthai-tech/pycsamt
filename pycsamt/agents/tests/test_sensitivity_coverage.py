# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.sensitivity`.

Mocks ``vertical_resolution``/``plot_sensitivity_depth_section``
(imported locally inside ``execute()``) with fast, deterministic fakes
so the DOI/groupby branches, both figure blocks, and the LLM
interpretation path can all be exercised without real EDI geometry:
the ``emtools.csumt``/``emtools.advanced`` import guard, the
``ensure_sites`` exception, the ``vertical_resolution`` exception, the
groupby-missing / groupby-exception DOI branches, the
``depth_max``/``period_range`` kwarg branches, the figure
None/exception/unsaved branches for both plots, and the LLM branch.
"""

from __future__ import annotations

import sys
import types

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pycsamt.agents.sensitivity import SensitivityAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _passthrough_ensure_sites(monkeypatch, sites=object()):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def _mock_emtools(
    monkeypatch,
    *,
    res_table="ok",
    vres_raises=False,
    plot_sens_raises=False,
    plot_sens_none=False,
):
    import pycsamt.emtools.advanced as advanced_mod
    import pycsamt.emtools.csumt as csumt_mod

    def _vertical_resolution(sites, rho_override=None, verbose=0):
        if vres_raises:
            raise RuntimeError("vres boom")
        if res_table == "no_groupby":
            return object()
        if res_table == "bad_columns":
            return pd.DataFrame({"foo": [1, 2, 3]})
        return pd.DataFrame(
            {
                "station": ["S1", "S1", "S2"],
                "depth_lo_m": [100.0, 300.0, 500.0],
            }
        )

    monkeypatch.setattr(csumt_mod, "vertical_resolution", _vertical_resolution)

    def _plot_sens(sites, component="xy", depth_unit="km", **kw):
        if plot_sens_raises:
            raise RuntimeError("plot_sens boom")
        if plot_sens_none:
            return None
        return plt.figure()

    monkeypatch.setattr(
        advanced_mod, "plot_sensitivity_depth_section", _plot_sens
    )


# ── import guard ─────────────────────────────────────────────────────────


def test_import_error_fails(monkeypatch):
    fake_mod = types.ModuleType("pycsamt.emtools.advanced")
    monkeypatch.setitem(sys.modules, "pycsamt.emtools.advanced", fake_mod)
    result = SensitivityAgent().execute({"path": "/nope"})
    assert result.status == "failed"
    assert "not available" in result.error


# ── ensure_sites / vertical_resolution exceptions ───────────────────────────


def test_ensure_sites_exception_fails(monkeypatch):
    _mock_emtools(monkeypatch)
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    result = SensitivityAgent().execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_vertical_resolution_exception_fails(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch, vres_raises=True)
    result = SensitivityAgent().execute({"sites": object()})
    assert result.status == "failed"
    assert "vertical_resolution failed" in result.error


# ── DOI / groupby branches ──────────────────────────────────────────────────


def test_no_groupby_skips_doi_and_mean(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch, res_table="no_groupby")
    result = SensitivityAgent().execute({"sites": object()})
    assert result.status == "success"
    assert result["doi_per_station"] == {}
    import math

    assert math.isnan(result["mean_doi_km"])


def test_groupby_exception_swallowed(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch, res_table="bad_columns")
    result = SensitivityAgent().execute({"sites": object()})
    assert result.status == "success"
    assert result["doi_per_station"] == {}


def test_doi_computed_from_real_table(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch)
    result = SensitivityAgent().execute({"sites": object()})
    assert result["doi_per_station"] == {"S1": 300.0, "S2": 500.0}
    assert result["mean_doi_km"] == pytest.approx(0.4)


# ── depth_max / period_range kwarg branches ─────────────────────────────────


def test_depth_max_and_period_range_kwargs_forwarded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    seen = {}

    import pycsamt.emtools.advanced as advanced_mod
    import pycsamt.emtools.csumt as csumt_mod

    monkeypatch.setattr(
        csumt_mod,
        "vertical_resolution",
        lambda sites, rho_override=None, verbose=0: pd.DataFrame(
            {"station": ["S1"], "depth_lo_m": [200.0]}
        ),
    )

    def _plot_sens(sites, component="xy", depth_unit="km", **kw):
        seen.update(kw)
        return plt.figure()

    monkeypatch.setattr(
        advanced_mod, "plot_sensitivity_depth_section", _plot_sens
    )
    result = SensitivityAgent().execute(
        {
            "sites": object(),
            "depth_max": 5.0,
            "period_range": [0.01, 100.0],
        }
    )
    assert result.status == "success"
    assert seen["depth_max"] == 5.0
    assert seen["period_range"] == (0.01, 100.0)


# ── sensitivity-section figure branches ─────────────────────────────────────


def test_sensitivity_figure_none_skips_block(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch, plot_sens_none=True)
    result = SensitivityAgent().execute({"sites": object()})
    assert "sensitivity_section" not in result["figures"]


def test_sensitivity_figure_no_output_dir_skips_path(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch)
    result = SensitivityAgent().execute({"sites": object()})
    assert "sensitivity_section" in result["figures"]
    assert "sensitivity_section" not in result["figure_paths"]


def test_sensitivity_figure_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch, plot_sens_raises=True)
    result = SensitivityAgent().execute({"sites": object()})
    assert any(
        "plot_sensitivity_depth_section" in w for w in result.warnings
    )


# ── DOI-bar figure branches ──────────────────────────────────────────────────


def test_doi_bar_skipped_when_no_doi(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch, res_table="no_groupby")
    result = SensitivityAgent().execute({"sites": object()})
    assert "doi_bar" not in result["figures"]


def test_doi_bar_no_output_dir_skips_path(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch)
    result = SensitivityAgent().execute({"sites": object()})
    assert "doi_bar" in result["figures"]
    assert "doi_bar" not in result["figure_paths"]


def test_doi_bar_saved_with_output_dir(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch)
    result = SensitivityAgent().execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result["figure_paths"].get("doi_bar")
    assert result["figure_paths"].get("sensitivity_section")


def test_doi_bar_exception_recorded(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch)
    import pycsamt.agents.sensitivity as sens_mod

    def _raise(doi):
        raise RuntimeError("doi bar boom")

    monkeypatch.setattr(sens_mod, "_plot_doi_bar", _raise)
    result = SensitivityAgent().execute({"sites": object()})
    assert any("DOI bar figure" in w for w in result.warnings)


# ── LLM interpretation branch ────────────────────────────────────────────────


def test_llm_interpretation_branch(monkeypatch):
    _passthrough_ensure_sites(monkeypatch)
    _mock_emtools(monkeypatch)
    agent = SensitivityAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute({"sites": object()})
    assert result.llm_interpretation == "mocked interpretation"
