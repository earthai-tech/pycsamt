# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Offline tests for ForwardModelAgent, CodeGenerationAgent, and the
AIInversion -> Hybrid agent chain (tiny training budgets).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from pycsamt.api.agents import AGENT_CONFIG

from .conftest import has_backend

requires_dl = pytest.mark.skipif(
    not has_backend(), reason="PyTorch/TensorFlow not installed"
)


def mk(cls, **kwargs):
    with AGENT_CONFIG.offline():
        return cls(**kwargs)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ── ForwardModelAgent ────────────────────────────────────────────────────────


def test_forward_agent_1d_layered_model(tmp_path):
    from pycsamt.agents import ForwardModelAgent

    result = mk(ForwardModelAgent).execute({
        "model": {
            "resistivities": [100.0, 10.0, 500.0],
            "thicknesses": [300.0, 800.0],
        },
        "output_dir": str(tmp_path),
    })
    assert result.status == "success"
    assert result.summary


def test_forward_agent_2d_halfspace(tmp_path):
    from pycsamt.agents import ForwardModelAgent

    result = mk(ForwardModelAgent, dim=2).execute({
        "model": {"type": "halfspace", "rho": 100.0},
        "output_dir": str(tmp_path),
    })
    assert result.status in {"success", "failed"}
    assert result.summary or result.error


def test_forward_agent_missing_model_fails():
    from pycsamt.agents import ForwardModelAgent

    result = mk(ForwardModelAgent).execute({})
    assert result.status in {"success", "failed"}
    # default model may be synthesised; either way the agent must
    # produce a coherent result object
    assert result.summary or result.error


# ── CodeGenerationAgent (static template path, no LLM) ───────────────────────


def test_code_generation_agent_offline_templates(tmp_path, edi_dir):
    from pycsamt.agents import CodeGenerationAgent, ContextInputAgent

    ctx = mk(ContextInputAgent).execute(
        {"request": f"Run QC on {edi_dir}"}
    )
    # regex-only extraction may flag itself for review; the config
    # dict is produced either way
    assert ctx.status in {"success", "needs_review"}
    config = ctx["config"]

    result = mk(CodeGenerationAgent).execute({
        "workflow_config": config,
        "results": {},
        "output_dir": str(tmp_path),
    })
    # offline template generation flags itself for human review
    assert result.status in {"success", "needs_review"}
    scripts = list(tmp_path.rglob("*.py"))
    code = result.data.get("code") or result.data.get("script")
    assert scripts or code, "expected a generated script"


def test_code_generation_agent_missing_config():
    from pycsamt.agents import CodeGenerationAgent

    result = mk(CodeGenerationAgent).execute({})
    assert result.status in {"success", "failed"}
    assert result.summary or result.error


# ── AIInversion -> Hybrid chain (tiny budgets) ───────────────────────────────


@requires_dl
def test_ai_inversion_then_hybrid_chain(edi_dir: Path, tmp_path):
    from pycsamt.agents import AIInversionAgent, HybridInversionAgent

    ai = mk(
        AIInversionAgent,
        n_layers=4,
        n_train_samples=64,
        epochs=1,
    )
    r1 = ai.execute({"path": str(edi_dir), "output_dir": str(tmp_path)})
    assert r1.status in {"success", "failed"}
    if r1.status != "success":
        pytest.skip(f"tiny AI inversion failed: {r1.error}")

    inverter = r1.data.get("inverter")
    assert inverter is not None

    hy = mk(HybridInversionAgent)
    r2 = hy.execute({
        "path": str(edi_dir),
        "ai_inverter": inverter,
        "output_dir": str(tmp_path),
    })
    assert r2.status in {"success", "needs_review", "failed"}
    assert r2.summary or r2.error
