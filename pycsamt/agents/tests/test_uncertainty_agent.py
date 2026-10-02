# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.agents.uncertainty_agent`.

Real data: ``data/AMT/WILLY_DATA/L18PLT`` (clean, no NaN feature rows)
and ``data/MT/broken-hill/edis`` (53 non-finite z_err_frac rows out of
2006 -- the regression case for the NaN-RMSE bug this test pins: an
early version dropped no rows before computing the held-out RMSE/R^2,
so any NaN landing in the validation split propagated straight through).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from pycsamt.agents.uncertainty_agent import UncertaintyCalibrationAgent
from pycsamt.agents.tests.conftest import has_backend

_ROOT = Path(__file__).resolve().parents[3]
_L18PLT = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18PLT.exists() and any(_L18PLT.glob("*.edi"))
_BROKEN_HILL = _ROOT / "data" / "MT" / "broken-hill" / "edis"
_HAS_BH = _BROKEN_HILL.exists() and any(_BROKEN_HILL.glob("*.edi"))

pytestmark = pytest.mark.skipif(
    not has_backend(), reason="requires PyTorch or TensorFlow"
)


def test_no_sites_or_path(no_llm_kw):
    agent = UncertaintyCalibrationAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_requires_a_dl_backend(no_llm_kw, monkeypatch):
    import pycsamt.backends as backends_mod

    monkeypatch.setattr(backends_mod, "get_backend_instance", lambda: None)
    agent = UncertaintyCalibrationAgent(**no_llm_kw)
    result = agent.execute({"path": "does/not/matter"})
    assert result.status == "failed"
    assert "PyTorch or TensorFlow" in result.error


@pytest.mark.skipif(not _HAS_L18, reason="WILLY L18PLT dataset not bundled")
def test_calibrates_clean_l18plt_data(no_llm_kw, tmp_output):
    agent = UncertaintyCalibrationAgent(
        **no_llm_kw, epochs=20, hidden=(16, 8)
    )
    result = agent.execute(
        {"path": str(_L18PLT), "output_dir": str(tmp_output)}
    )

    assert result.status == "success"
    assert result.warnings == []  # no non-finite rows on this dataset
    assert not math.isnan(result["rmse"])
    assert not math.isnan(result["r2"])
    assert result["rmse"] >= 0.0
    assert len(result["table"]) > 0
    assert "z_err_frac_calibrated" in result["table"].columns
    assert "uncertainty_summary" in result["figures"]


@pytest.mark.skipif(not _HAS_BH, reason="Broken Hill dataset not bundled")
def test_drops_non_finite_rows_before_computing_metrics(no_llm_kw):
    """Regression: Broken Hill has 53 non-finite z_err_frac rows; an
    earlier version let them leak into the held-out RMSE/R^2 computation,
    producing NaN metrics instead of a real number."""
    agent = UncertaintyCalibrationAgent(**no_llm_kw, epochs=15)
    result = agent.execute({"path": str(_BROKEN_HILL)})

    assert result.status == "success"
    assert any("non-finite" in w for w in result.warnings)
    assert not math.isnan(result["rmse"])
    assert not math.isnan(result["r2"])


def test_too_few_rows_fails_cleanly(no_llm_kw, monkeypatch):
    import pandas as pd

    import pycsamt.agents.uncertainty_agent as mod

    monkeypatch.setattr(
        mod, "ensure_sites", lambda *a, **k: object(), raising=False
    )
    monkeypatch.setattr(
        "pycsamt.emtools._core.ensure_sites", lambda *a, **k: object()
    )
    monkeypatch.setattr(
        "pycsamt.ai.processing.build_uncertainty_features_table",
        lambda sites: pd.DataFrame({"z_err_frac": [0.01, 0.02]}),
    )
    agent = UncertaintyCalibrationAgent(**no_llm_kw)
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert ">= 5" in result.error
