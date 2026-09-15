# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.static_shift`.

Targets the branches left uncovered by the broader agent battery tests:
missing 'sites'/'path', ``ensure_sites`` failure, each correction method
(loess / bilateral / refmedian / unknown-fallback), the correction
exception fallback, the delta-stats exception guard, the three figure
exception guards, the "no significant shift" narrative branch, and the
``_collect_rho`` helper's continue / matching-length / empty-result
branches.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.agents.static_shift import StaticShiftAgent, _collect_rho

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ── fake site helpers for direct _collect_rho tests ─────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _FakeSite:
    def __init__(self, z=None, freq=None, station="st"):
        self.station = station
        if z is not None:
            self.Z = _FakeZ(z, freq)


def _zblock(n=4, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    freq = np.linspace(1.0, 100.0, n)
    return z, freq


def _passthrough_ensure_sites(monkeypatch, sites):
    """Make ``ensure_sites`` return *sites* unchanged (identity-preserving).

    ``StaticShiftAgent.execute`` calls ``ensure_sites`` before every
    correction path, which may return a converted/copied object. Tests
    that need to assert identity on ``corrected_sites`` bypass that
    conversion so the object flowing through the agent is exactly the
    fixture instance.
    """
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


# ── _collect_rho direct tests ────────────────────────────────────────────────


def test_collect_rho_skips_sites_with_no_z_and_returns_empty():
    site = _FakeSite(z=None, station="no_z")
    mat, freqs, labels = _collect_rho([site])
    assert mat is None
    assert freqs is None
    assert labels == []


def test_collect_rho_matching_length_second_station_appends():
    z1, fr = _zblock(4, seed=1)
    z2, _ = _zblock(4, seed=2)
    s1 = _FakeSite(z1, fr, station="s1")
    s2 = _FakeSite(z2, fr, station="s2")
    mat, freqs, labels = _collect_rho([s1, s2])
    assert mat is not None
    assert mat.shape == (4, 2)
    assert labels == ["s1", "s2"]
    assert np.array_equal(freqs, fr)


# ── agent-level tests ─────────────────────────────────────────────────────────


def test_no_sites_or_path(no_llm_kw):
    agent = StaticShiftAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_ensure_sites_raises(no_llm_kw, monkeypatch):
    import pycsamt.emtools._core as core

    def _boom(*a, **k):
        raise ValueError("bad path")

    monkeypatch.setattr(core, "ensure_sites", _boom)
    agent = StaticShiftAgent(**no_llm_kw)
    result = agent.execute({"path": "/does/not/matter"})
    assert result.status == "failed"
    assert "bad path" in result.error


def test_method_none_skips_correction(no_llm_kw, loaded_sites):
    agent = StaticShiftAgent(**no_llm_kw, method="none")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("skipped" in w for w in result.warnings)
    assert result["delta_stats"]["n_shifted"] == 0


def test_ama_no_shift_reports_narrative_note(
    no_llm_kw, loaded_sites, monkeypatch
):
    # method != "none"/"skip" but the correction is a no-op -> rho_after
    # == rho_before, n_shifted == 0 -> the "no significant shift" note
    # branch (as opposed to the "skipped" one above).
    import pycsamt.emtools.ss as ss

    monkeypatch.setattr(ss, "correct_ss_ama", lambda sites, **k: sites)
    agent = StaticShiftAgent(**no_llm_kw, method="ama")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert "no significant shift" in result.summary
    assert any("strongly" in w or "3-D" in w for w in result.warnings)
    assert result["delta_stats"]["n_shifted"] == 0


def test_method_loess(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.ss as ss

    _passthrough_ensure_sites(monkeypatch, loaded_sites)
    monkeypatch.setattr(
        ss, "estimate_ss_loess", lambda sites, **k: "loess_result"
    )
    monkeypatch.setattr(
        ss, "apply_ss_factors", lambda sites, res, **k: sites
    )
    agent = StaticShiftAgent(**no_llm_kw, method="loess")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["corrected_sites"] is loaded_sites


def test_method_bilateral(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.ss as ss

    _passthrough_ensure_sites(monkeypatch, loaded_sites)
    monkeypatch.setattr(
        ss, "estimate_ss_bilateral", lambda sites, **k: "bilateral_result"
    )
    monkeypatch.setattr(
        ss, "apply_ss_factors", lambda sites, res, **k: sites
    )
    agent = StaticShiftAgent(**no_llm_kw, method="bilateral")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["corrected_sites"] is loaded_sites


def test_method_refmedian_dataframe(no_llm_kw, loaded_sites, monkeypatch):
    import pandas as pd
    import pycsamt.emtools.ss as ss

    df = pd.DataFrame(
        {"station": ["s1", "s2"], "fac_z": [1.1, 0.9]}
    )
    monkeypatch.setattr(ss, "estimate_ss_refmedian", lambda sites, **k: df)
    monkeypatch.setattr(
        ss, "apply_ss_factors", lambda sites, res, **k: sites
    )
    agent = StaticShiftAgent(**no_llm_kw, method="refmedian")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["shift_factors"] == {"s1": 1.1, "s2": 0.9}


def test_method_refmedian_unexpected_type_falls_back_to_ama(
    no_llm_kw, loaded_sites, monkeypatch
):
    import pycsamt.emtools.ss as ss

    monkeypatch.setattr(
        ss, "estimate_ss_refmedian", lambda sites, **k: {"not": "a df"}
    )
    monkeypatch.setattr(
        ss, "correct_ss_ama", lambda sites, **k: sites
    )
    agent = StaticShiftAgent(**no_llm_kw, method="refmedian")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("unexpected type" in w for w in result.warnings)


def test_unknown_method_falls_back_to_ama(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.ss as ss

    monkeypatch.setattr(ss, "correct_ss_ama", lambda sites, **k: sites)
    agent = StaticShiftAgent(**no_llm_kw, method="not_a_real_method")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("falling back to AMA" in w for w in result.warnings)


def test_correction_exception_falls_back_to_raw_data(
    no_llm_kw, loaded_sites, monkeypatch
):
    import pycsamt.emtools.ss as ss

    def _boom(sites, **k):
        raise RuntimeError("ama exploded")

    _passthrough_ensure_sites(monkeypatch, loaded_sites)
    monkeypatch.setattr(ss, "correct_ss_ama", _boom)
    agent = StaticShiftAgent(**no_llm_kw, method="ama")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("ama exploded" in w for w in result.warnings)
    assert result["corrected_sites"] is loaded_sites


def test_delta_stats_exception_is_swallowed(
    no_llm_kw, loaded_sites, monkeypatch
):
    import pycsamt.agents.static_shift as ss_mod

    calls = {"n": 0}

    def _fake_collect_rho(sites):
        calls["n"] += 1
        if calls["n"] == 1:
            return np.ones((5, 3)), np.linspace(1, 5, 5), ["a", "b", "c"]
        # second call (post-correction) returns a shape that cannot
        # broadcast against the first -> `rho_after - rho_before` raises.
        return np.ones((5, 4)), np.linspace(1, 5, 5), ["a", "b", "c", "d"]

    monkeypatch.setattr(ss_mod, "_collect_rho", _fake_collect_rho)
    agent = StaticShiftAgent(**no_llm_kw, method="none")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["delta_stats"] == {}


def test_figure_exceptions_are_captured_as_warnings(
    no_llm_kw, loaded_sites, monkeypatch
):
    import pycsamt.emtools.ss as ss

    def _boom(*a, **k):
        raise RuntimeError("plot failed")

    monkeypatch.setattr(ss, "plot_ss_summary", _boom)
    monkeypatch.setattr(ss, "plot_ss_1d_curves", _boom)
    monkeypatch.setattr(ss, "ss_comparison_psection", _boom)

    agent = StaticShiftAgent(**no_llm_kw, method="none")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    joined = " ".join(result.warnings)
    assert "plot_ss_summary" in joined
    assert "plot_ss_1d_curves" in joined
    assert "ss_comparison_psection" in joined
    assert result["figures"] == {}


def test_ama_happy_path_saves_figures(no_llm_kw, loaded_sites, tmp_output):
    agent = StaticShiftAgent(**no_llm_kw, method="ama")
    result = agent.execute(
        {"sites": loaded_sites, "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert result["figures"]
    assert result["figure_paths"]
    for p in result["figure_paths"].values():
        assert Path(p).exists()
