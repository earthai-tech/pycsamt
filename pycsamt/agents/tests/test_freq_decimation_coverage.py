# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.freq_decimation`.

Targets: ``ensure_sites`` failure, the (dead-code-protected) ``qc_result``
parsing try/except, the per-station no-Z skip, the period-range filter,
the SNR-proxy exception guard, the no-good-frequency warning, the
log-spaced decimation branch, the figure exception guard, the LLM
interpretation path, and ``_plot_selection_summary``'s missing-frequency
break.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.agents.freq_decimation import (
    FrequencyDecimationAgent,
    _plot_selection_summary,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _passthrough_ensure_sites(monkeypatch, sites):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def test_no_sites_or_path(no_llm_kw):
    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_ensure_sites_raises(no_llm_kw, monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad path")),
    )
    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad path" in result.error


def test_qc_result_dict_snr_section_is_parsed(no_llm_kw, loaded_sites):
    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute(
        {
            "sites": loaded_sites,
            "qc_result": {"snr_section": {"a": 1}},
        }
    )
    assert result.status in ("success", "needs_review")


def test_qc_result_get_raising_is_swallowed(no_llm_kw, loaded_sites):
    class _BoomQC:
        def get(self, *a, **k):
            raise RuntimeError("boom")

    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute(
        {"sites": loaded_sites, "qc_result": _BoomQC()}
    )
    assert result.status in ("success", "needs_review")


def test_station_with_no_z_or_freq_is_skipped(no_llm_kw, monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core, "_get_z_block", lambda ed, **k: (None, None, None)
    )
    _passthrough_ensure_sites(monkeypatch, [object()])
    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute({"sites": [object()]})
    assert result.status == "needs_review"
    assert result["n_original"] == 0


def test_period_range_filter_applied(no_llm_kw, loaded_sites):
    agent = FrequencyDecimationAgent(**no_llm_kw, snr_threshold=0.0)
    result_full = agent.execute({"sites": loaded_sites})
    result_narrow = agent.execute(
        {"sites": loaded_sites, "period_range": [1e-3, 1e-2]}
    )
    assert result_narrow["n_original"] <= result_full["n_original"]


def test_snr_proxy_exception_is_swallowed(no_llm_kw, monkeypatch):
    import pycsamt.emtools._core as core

    bad_z = np.zeros((4, 2), dtype=complex)  # wrong ndim -> IndexError
    fr = np.array([1.0, 2.0, 3.0, 4.0])
    monkeypatch.setattr(
        core, "_get_z_block", lambda ed, **k: (None, bad_z, fr)
    )
    _passthrough_ensure_sites(monkeypatch, [object()])
    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute({"sites": [object()]})
    assert result.status in ("success", "needs_review")


def test_no_good_frequencies_warns_and_skips(no_llm_kw, loaded_sites):
    agent = FrequencyDecimationAgent(**no_llm_kw, snr_threshold=1e12)
    result = agent.execute({"sites": loaded_sites})
    assert any("no frequencies pass SNR threshold" in w for w in result.warnings)
    assert all(p.size == 0 for p in result["selected_periods"].values())


def test_decimation_reduces_period_count(no_llm_kw, loaded_sites):
    agent = FrequencyDecimationAgent(
        **no_llm_kw, n_per_decade=1, snr_threshold=0.0
    )
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["n_selected"] < result["n_original"]


def test_figure_exception_is_captured(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.agents.freq_decimation as fd

    monkeypatch.setattr(
        fd,
        "_plot_selection_summary",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("plot boom")),
    )
    agent = FrequencyDecimationAgent(**no_llm_kw)
    result = agent.execute({"sites": loaded_sites})
    assert any("Selection summary figure" in w for w in result.warnings)


def test_llm_interpretation_called_when_api_key_set(loaded_sites):
    agent = FrequencyDecimationAgent(api_key="fake-key", snr_threshold=0.0)
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute({"sites": loaded_sites})
    assert result.llm_interpretation == "mocked interpretation"


def test_plot_selection_summary_breaks_when_freq_missing():
    class _NoZSite:
        station = "s1"

    fig = _plot_selection_summary(
        selected_periods={"s1": np.array([1.0])},
        dead_band_mask={"s1": np.array([False])},
        sites=[_NoZSite()],
        n_per_decade=6,
        snr_threshold=3.0,
    )
    assert fig is not None
    assert fig.axes
