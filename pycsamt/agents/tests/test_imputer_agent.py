# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.agents.imputer_agent`.

Real data: ``data/MT/broken-hill/edis`` -- unlike the bundled 3-EDI
dataset, several of its stations have genuinely missing frequency rows
(see ``docs/scripts/generate_user_guide_ai_processing_imputer_figures.py``),
so this exercises the real gap-filling path, not a synthetic stand-in.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from pycsamt.agents.imputer_agent import ImputerAgent
from pycsamt.agents.tests.conftest import has_backend

_ROOT = Path(__file__).resolve().parents[3]
_BROKEN_HILL = _ROOT / "data" / "MT" / "broken-hill" / "edis"
_HAS_BH = _BROKEN_HILL.exists() and any(_BROKEN_HILL.glob("*.edi"))
_L18PLT = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18PLT.exists() and any(_L18PLT.glob("*.edi"))

pytestmark = pytest.mark.skipif(
    not has_backend(), reason="requires PyTorch or TensorFlow"
)


def test_no_sites_or_path(no_llm_kw):
    agent = ImputerAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_requires_a_dl_backend(no_llm_kw, monkeypatch):
    # get_backend_instance is imported function-locally inside execute();
    # patching the source module (not this one) is what actually takes
    # effect at call time.
    import pycsamt.backends as backends_mod

    monkeypatch.setattr(backends_mod, "get_backend_instance", lambda: None)
    agent = ImputerAgent(**no_llm_kw)
    result = agent.execute({"path": "does/not/matter"})
    assert result.status == "failed"
    assert "PyTorch or TensorFlow" in result.error


@pytest.mark.skipif(not _HAS_BH, reason="Broken Hill dataset not bundled")
def test_fills_real_gaps_on_broken_hill(no_llm_kw, tmp_output):
    agent = ImputerAgent(**no_llm_kw, epochs=12, channels=(16, 32, 16))
    result = agent.execute(
        {"path": str(_BROKEN_HILL), "output_dir": str(tmp_output)}
    )

    assert result.status == "success"
    assert result["n_missing_rows"] > 0
    assert result["n_gapped_stations"] > 0
    assert len(result["gapped_stations"]) == result["n_gapped_stations"]
    assert result["gap_table"] is not None
    assert "imputer_summary" in result["figures"]
    assert result.llm_interpretation is None  # no api_key
    assert result.cost_estimate_usd == 0.0


@pytest.mark.skipif(not _HAS_L18, reason="WILLY L18PLT dataset not bundled")
def test_no_missing_cells_short_circuits(no_llm_kw):
    """WILLY L18PLT has zero genuinely missing cells -- real short-circuit
    path, not a mocked one."""
    agent = ImputerAgent(**no_llm_kw)
    result = agent.execute({"path": str(_L18PLT)})

    assert result.status == "success"
    assert result["n_missing_rows"] == 0
    assert "nothing to fill" in result.summary.lower()
    assert result["figures"] == {}
