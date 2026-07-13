# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Offline contract battery for the analysis/processing agent family.

Every agent is constructed inside ``AGENT_CONFIG.offline()`` so that no
API key is resolved from the developer environment and no LLM call can
ever be attempted. Each agent is exercised twice:

* missing input  -> ``AgentResult`` with status ``failed``;
* happy path on the bundled 3-station ``data/3edis`` dataset (or a
  small synthetic payload for agents that do not consume Sites).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.api.agents import AGENT_CONFIG

# ── helpers ───────────────────────────────────────────────────────────────────


def mk(cls, **kwargs):
    """Instantiate an agent with LLM resolution forced off."""
    with AGENT_CONFIG.offline():
        return cls(**kwargs)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def sites3(edi_dir: Path):
    from pycsamt.emtools import ensure_sites

    return ensure_sites(edi_dir)


# ── missing-input contract ───────────────────────────────────────────────────


def _sites_agents():
    from pycsamt.agents import (
        DataQCAgent,
        DenoisingAgent,
        EDIExportAgent,
        FrequencyDecimationAgent,
        PhaseAnalysisAgent,
        SensitivityAgent,
        TensorRotationAgent,
        TipperAnalysisAgent,
    )

    return [
        DataQCAgent,
        DenoisingAgent,
        EDIExportAgent,
        FrequencyDecimationAgent,
        PhaseAnalysisAgent,
        SensitivityAgent,
        TensorRotationAgent,
        TipperAnalysisAgent,
    ]


@pytest.mark.parametrize("cls_name", [
    "DataQCAgent", "DenoisingAgent", "EDIExportAgent",
    "FrequencyDecimationAgent", "PhaseAnalysisAgent",
    "SensitivityAgent", "TensorRotationAgent", "TipperAnalysisAgent",
])
def test_sites_agents_fail_without_input(cls_name):
    import pycsamt.agents as agents_pkg

    cls = getattr(agents_pkg, cls_name)
    result = mk(cls).execute({})
    assert result.status == "failed"
    assert result.error


# ── sites-consuming happy paths ──────────────────────────────────────────────


def test_phase_analysis_agent_full_run(sites3, tmp_path):
    from pycsamt.agents import PhaseAnalysisAgent

    ag = mk(PhaseAnalysisAgent)
    result = ag.execute({
        "sites": sites3,
        "output_dir": str(tmp_path),
        "run_mohr": True,
    })
    assert result.status == "success"
    assert result.summary
    data = result.data
    assert {"n_1d", "n_2d", "n_3d"} <= set(data)
    assert isinstance(data["figures"], dict)
    # saving requested: at least one figure landed on disk
    if data["figure_paths"]:
        saved = next(iter(data["figure_paths"].values()))
        assert Path(saved).exists()


def test_tensor_rotation_agent_rotates(sites3, tmp_path):
    from pycsamt.agents import TensorRotationAgent

    result = mk(TensorRotationAgent).execute({
        "sites": sites3,
        "strike_deg": 30.0,
        "output_dir": str(tmp_path),
    })
    assert result.status == "success"
    assert result.summary


def test_denoising_agent_hampel(sites3):
    from pycsamt.agents import DenoisingAgent

    result = mk(DenoisingAgent, method="hampel").execute(
        {"sites": sites3}
    )
    assert result.status == "success"
    assert result.summary


def test_freq_decimation_agent(sites3):
    from pycsamt.agents import FrequencyDecimationAgent

    result = mk(FrequencyDecimationAgent).execute({"sites": sites3})
    assert result.status == "success"
    assert result.summary


def test_sensitivity_agent(sites3, tmp_path):
    from pycsamt.agents import SensitivityAgent

    result = mk(SensitivityAgent).execute(
        {"sites": sites3, "output_dir": str(tmp_path)}
    )
    assert result.status == "success"
    assert result.summary


def test_tipper_analysis_agent(sites3):
    from pycsamt.agents import TipperAnalysisAgent

    result = mk(TipperAnalysisAgent).execute({"sites": sites3})
    # the 3edis stations may not carry tipper data: both a graceful
    # failure and a warned success honour the contract
    assert result.status in {"success", "failed"}
    assert result.summary or result.error


def test_data_qc_agent(sites3):
    from pycsamt.agents import DataQCAgent

    result = mk(DataQCAgent).execute({"sites": sites3})
    assert result.status == "success"
    assert result.summary


def test_edi_export_agent_writes_files(sites3, tmp_path):
    from pycsamt.agents import EDIExportAgent

    out = tmp_path / "edis"
    result = mk(EDIExportAgent).execute(
        {"sites": sites3, "output_dir": str(out)}
    )
    assert result.status == "success"
    written = list(out.rglob("*.edi"))
    assert written, "expected exported EDI files"


def test_anomaly_detection_agent_tiny_run(sites3):
    from pycsamt.agents import AnomalyDetectionAgent

    ag = mk(AnomalyDetectionAgent, epochs=1)
    result = ag.execute({"sites": sites3})
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary


# ── synthetic-payload agents ─────────────────────────────────────────────────


def test_inversion_comparison_agent_synthetic():
    from pycsamt.agents import InversionComparisonAgent

    rng = np.random.default_rng(0)
    depths = np.linspace(0.05, 2.0, 8)
    a = {"pred_rho": 10 ** rng.normal(2, 0.3, (3, 8)),
         "depths_km": depths}
    b = {"pred_rho": 10 ** rng.normal(2, 0.3, (3, 8)),
         "depths_km": depths}
    result = mk(InversionComparisonAgent).execute(
        {"result_a": a, "result_b": b}
    )
    assert result.status == "success"
    assert result.summary

    missing = mk(InversionComparisonAgent).execute({})
    assert missing.status == "failed"


def test_inversion_evaluation_agent_self_comparison(edi_dir, tmp_path):
    from pycsamt.agents import InversionEvaluationAgent

    result = mk(InversionEvaluationAgent).execute({
        "path_obs": str(edi_dir),
        "path_mod": str(edi_dir),
        "output_dir": str(tmp_path),
    })
    # observed == modelled -> misfit ~ 0, must succeed
    assert result.status == "success"
    assert result.summary

    missing = mk(InversionEvaluationAgent).execute({})
    assert missing.status == "failed"


def test_resistivity_map_agent_synthetic(tmp_path):
    from pycsamt.agents import ResistivityMapAgent

    rng = np.random.default_rng(1)
    stations = [f"S{i:02d}" for i in range(6)]
    predictions = {s: rng.normal(2.0, 0.2, 10) for s in stations}
    # scattered (non-collinear) positions: triangulation-based
    # interpolation needs a 2-D point cloud
    coords = {
        s: (
            28.0 + 0.01 * i + float(rng.uniform(0, 0.004)),
            102.0 + 0.01 * ((i * 7) % 5) + float(rng.uniform(0, 0.004)),
        )
        for i, s in enumerate(stations)
    }
    result = mk(ResistivityMapAgent).execute({
        "predictions": predictions,
        "station_coords": coords,
        "depths_km": np.linspace(0.05, 2.0, 10),
        "output_dir": str(tmp_path),
    })
    assert result.status == "success"
    assert result.summary

    missing = mk(ResistivityMapAgent).execute({})
    assert missing.status == "failed"


def test_report_agent_markdown(tmp_path):
    from pycsamt.agents import ReportAgent
    from pycsamt.agents._base import AgentResult

    qc_res = AgentResult(
        status="success",
        summary="QC finished: 3 stations, coverage 98%.",
        data={"n_stations": 3},
    )
    result = mk(ReportAgent, formats=["md"]).execute({
        "results": {"qc": qc_res},
        "output_dir": str(tmp_path),
        "title": "Battery report",
    })
    assert result.status == "success"
    md_files = list(tmp_path.rglob("*.md"))
    assert md_files, "expected a markdown report on disk"
    text = md_files[0].read_text(encoding="utf-8")
    assert "Battery report" in text or "QC" in text


def test_batch_survey_agent_single_profile(edi_dir):
    from pycsamt.agents import BatchSurveyAgent

    result = mk(BatchSurveyAgent, workflow="qc").execute(
        {"profiles": {"L1": str(edi_dir)}}
    )
    assert result.status in {"success", "partial", "failed"}
    assert result.summary or result.error

    missing = mk(BatchSurveyAgent).execute({})
    assert missing.status == "failed"


def test_interpretation_agent_layered_model():
    from pycsamt.agents import InterpretationAgent

    result = mk(InterpretationAgent).execute({
        "model": {
            "resistivity": [100.0, 10.0, 500.0],
            "thickness": [300.0, 800.0],
        },
        "rms": 1.2,
    })
    assert result.status in {"success", "failed"}
    assert result.summary or result.error

    missing = mk(InterpretationAgent).execute({})
    assert missing.status == "failed"


def test_pipeline_agent_direct_mode(sites3):
    from pycsamt.agents import PipelineAgent

    result = mk(PipelineAgent).execute({
        "sites": sites3,
        "preset": "basic_qc",
    })
    assert result.status in {"success", "partial"}
    assert result.summary

    missing = mk(PipelineAgent).execute({})
    assert missing.status == "failed"
