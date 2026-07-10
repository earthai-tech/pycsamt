# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for Mare2DEMAgent (no LLM key, no binary required)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.agents import Mare2DEMAgent
from pycsamt.agents._base import AgentResult

# Data directory shared with mare2dem model tests
_DATA_DIR = Path(__file__).parents[3] / "data" / "mare2dem"
_HILL_DIR  = _DATA_DIR / "hill"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def agent():
    return Mare2DEMAgent(n_procs=4, initial_rho=1.0, target_rms=1.0)


@pytest.fixture(scope="module")
def hill_emdata():
    if not (_HILL_DIR / "hill.emdata").exists():
        pytest.skip("Hill example data not found.")
    return str(_HILL_DIR / "hill.emdata")


# ---------------------------------------------------------------------------
# Instantiation
# ---------------------------------------------------------------------------

class TestMare2DEMAgentInit:

    def test_default_params(self):
        a = Mare2DEMAgent()
        assert a.n_procs == 4
        assert a.use_mpi is True
        assert a.initial_rho == 1.0
        assert a.target_rms  == 1.0
        assert a.max_iterations == 150

    def test_custom_params(self):
        a = Mare2DEMAgent(n_procs=16, initial_rho=10.0, target_rms=0.95)
        assert a.n_procs == 16
        assert a.initial_rho == 10.0

    def test_invalid_provider_raises(self):
        with pytest.raises(ValueError):
            Mare2DEMAgent(llm_provider="invalid_llm")

    def test_lazy_import_from_agents(self):
        from pycsamt.agents import Mare2DEMAgent as M2D
        assert M2D is Mare2DEMAgent


# ---------------------------------------------------------------------------
# Mode: prepare — from existing .emdata
# ---------------------------------------------------------------------------

class TestPrepareFromExisting:

    def test_prepare_succeeds(self, agent, hill_emdata, tmp_path):
        result = agent.execute({
            "emdata":     hill_emdata,
            "output_dir": str(tmp_path),
            "mode":       "prepare",
        })
        assert result.status in ("success", "needs_review")
        assert result["data_path"] is not None

    def test_prepare_writes_three_files(self, agent, hill_emdata, tmp_path):
        result = agent.execute({
            "emdata":     hill_emdata,
            "output_dir": str(tmp_path),
        })
        # data, model, settings
        assert result["data_path"].exists()
        assert result["resistivity_path"].exists()
        assert result["settings_path"].exists()

    def test_output_dir_created(self, agent, hill_emdata, tmp_path):
        out = tmp_path / "subdir" / "run"
        agent.execute({
            "emdata":     hill_emdata,
            "output_dir": str(out),
        })
        assert out.exists()

    def test_data_statistics_populated(self, agent, hill_emdata, tmp_path):
        result = agent.execute({
            "emdata":     hill_emdata,
            "output_dir": str(tmp_path),
        })
        assert result["n_mt_receivers"] == 228
        assert result["n_data"] > 0

    def test_output_dir_key_present(self, agent, hill_emdata, tmp_path):
        result = agent.execute({
            "emdata":     hill_emdata,
            "output_dir": str(tmp_path),
        })
        assert result["output_dir"] == str(tmp_path)

    def test_binary_found_key_present(self, agent, hill_emdata, tmp_path):
        result = agent.execute({
            "emdata":     hill_emdata,
            "output_dir": str(tmp_path),
        })
        assert isinstance(result["binary_found"], bool)


# ---------------------------------------------------------------------------
# Mode: prepare — from inline MT survey config
# ---------------------------------------------------------------------------

class TestPrepareFromSurveyConfig:

    def test_mt_survey_basic(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [0.01, 0.1, 1.0, 10.0, 100.0],
                "rx_y":        list(np.linspace(-3000, 3000, 10)),
                "rx_type":     "marine",
                "lTE":         True,
                "lTM":         True,
            },
            "topo":       -800.0,
            "output_dir": str(tmp_path),
        })
        assert result.status in ("success", "needs_review")
        assert result["n_mt_receivers"] == 10
        assert result["n_data"] > 0

    def test_csem_survey_basic(self, agent, tmp_path):
        result = agent.execute({
            "csem": {
                "frequencies": [0.25, 1.0],
                "tx_y":        list(np.linspace(-2000, 2000, 5)),
                "rx_y":        list(np.linspace(-4000, 4000, 8)),
                "rx_type":     "marine",
                "lEx":         True,
            },
            "topo":       -1200.0,
            "output_dir": str(tmp_path),
        })
        assert result["n_csem_transmitters"] == 5
        assert result["n_data"] > 0

    def test_combined_mt_csem(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0, 10.0],
                "rx_y":        list(np.linspace(-2000, 2000, 5)),
                "lTE": True,
            },
            "csem": {
                "frequencies": [0.5],
                "tx_y":        [0.0, 1000.0],
                "rx_y":        list(np.linspace(-1000, 1000, 4)),
                "lEx":         True,
            },
            "topo":       -500.0,
            "output_dir": str(tmp_path),
        })
        assert result["n_mt_receivers"] >= 5
        assert result["n_csem_transmitters"] == 2
        assert result["n_data"] > 0

    def test_no_source_fails(self, agent, tmp_path):
        result = agent.execute({
            "output_dir": str(tmp_path),
        })
        assert result.status == "failed"
        assert "data source" in result.error.lower()

    def test_result_contains_all_keys(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
        })
        for key in (
            "data_path", "resistivity_path", "settings_path",
            "binary_found", "source_downloaded",
            "n_mt_receivers", "n_csem_transmitters", "n_data",
            "final_rms", "n_iterations", "converged", "result",
            "output_dir",
        ):
            assert key in result, f"Missing key: {key!r}"


# ---------------------------------------------------------------------------
# Mode: run — no binary, should warn and not crash
# ---------------------------------------------------------------------------

class TestRunModeNoBinary:

    def test_run_no_binary_warns(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
            "mode":       "run",
        })
        # should not fail outright; binary warning expected
        assert result.status in ("success", "needs_review")
        if not result["binary_found"]:
            assert any("binary" in w.lower() for w in result.warnings)

    def test_run_result_keys_null(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
            "mode":       "run",
        })
        if not result["binary_found"]:
            assert result["final_rms"]  is None
            assert result["n_iterations"] == 0
            assert result["converged"]  is False


# ---------------------------------------------------------------------------
# Mode: report — scan existing run directory
# ---------------------------------------------------------------------------

class TestReportMode:

    @pytest.fixture(scope="class")
    def inv_dir(self):
        d = _DATA_DIR / "demo_mt_inversion"
        if not d.exists():
            pytest.skip("demo_mt_inversion data not found.")
        return d

    def test_report_basic(self, agent, inv_dir):
        result = agent.execute({
            "output_dir": str(inv_dir),
            "mode":       "report",
        })
        assert result.status in ("success", "needs_review")

    def test_report_log_parsed(self, agent, inv_dir):
        result = agent.execute({
            "output_dir": str(inv_dir),
            "mode":       "report",
        })
        assert result["n_iterations"] == 6

    def test_report_final_rms(self, agent, inv_dir):
        result = agent.execute({
            "output_dir": str(inv_dir),
            "mode":       "report",
        })
        assert result["final_rms"] is not None
        assert result["final_rms"] < 1.5

    def test_report_no_crash_empty_dir(self, agent, tmp_path):
        result = agent.execute({
            "output_dir": str(tmp_path),
            "mode":       "report",
        })
        # empty dir → result loaded but no log; should not crash
        assert result.status in ("success", "needs_review", "failed")


# ---------------------------------------------------------------------------
# Invalid mode
# ---------------------------------------------------------------------------

def test_invalid_mode(tmp_path):
    agent = Mare2DEMAgent()
    result = agent.execute({
        "emdata":     "some.emdata",
        "output_dir": str(tmp_path),
        "mode":       "invalid_mode",
    })
    assert result.status == "failed"


# ---------------------------------------------------------------------------
# AgentResult contract
# ---------------------------------------------------------------------------

class TestAgentResultContract:

    def test_result_bool_success(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
        })
        assert bool(result)  # success → truthy

    def test_result_failed_bool(self):
        r = AgentResult.failed("test error")
        assert not bool(r)

    def test_result_getitem(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
        })
        assert result["output_dir"] == str(tmp_path)

    def test_result_get_default(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
        })
        assert result.get("nonexistent_key", 42) == 42

    def test_result_repr(self, agent, tmp_path):
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
        })
        r = repr(result)
        assert "AgentResult" in r
        assert "status=" in r

    def test_no_llm_interp_without_key(self, agent, tmp_path):
        # agent has no api_key → llm_interpretation should be None
        result = agent.execute({
            "mt": {
                "frequencies": [1.0],
                "rx_y": [0.0],
                "lTE": True,
            },
            "output_dir": str(tmp_path),
        })
        assert result.llm_interpretation is None
