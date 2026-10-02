# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.agents.distortion_agent`.

Real data: ``data/AMT/WILLY_DATA/L18PLT`` (28 stations, the same line
used throughout ``pycsamt.ai.processing.distortion``'s own user guide).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from pycsamt.agents.distortion_agent import (
    DistortionClassificationAgent,
    _SELF_TRAINING_WARNING,
)
from pycsamt.agents.tests.conftest import has_backend

_ROOT = Path(__file__).resolve().parents[3]
_L18PLT = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18PLT.exists() and any(_L18PLT.glob("*.edi"))

pytestmark = pytest.mark.skipif(
    not has_backend(), reason="requires PyTorch or TensorFlow"
)


def test_no_sites_or_path(no_llm_kw):
    agent = DistortionClassificationAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error
    # A hard input-validation failure never reaches the self-training
    # code path, so no disclaimer is attached here.
    assert result.warnings == []


def test_requires_a_dl_backend(no_llm_kw, monkeypatch):
    import pycsamt.backends as backends_mod

    monkeypatch.setattr(backends_mod, "get_backend_instance", lambda: None)
    agent = DistortionClassificationAgent(**no_llm_kw)
    result = agent.execute({"path": "does/not/matter"})
    assert result.status == "failed"
    assert "PyTorch or TensorFlow" in result.error


@pytest.mark.skipif(not _HAS_L18, reason="WILLY L18PLT dataset not bundled")
def test_triages_real_l18plt_stations_and_always_warns(no_llm_kw, tmp_output):
    agent = DistortionClassificationAgent(**no_llm_kw, epochs=40)
    result = agent.execute(
        {"path": str(_L18PLT), "output_dir": str(tmp_output)}
    )

    assert result.status == "success"
    assert set(result["regime_counts"]) <= {
        "clean", "static_shift_only", "distorted",
    }
    assert sum(result["regime_counts"].values()) == len(result["table"])
    assert "regime_label" in result["table"].columns
    assert "distortion_summary" in result["figures"]
    # The self-training-instability disclaimer is attached on every
    # successful run, not only when something looks unusual.
    assert _SELF_TRAINING_WARNING in result.warnings


@pytest.mark.skipif(not _HAS_L18, reason="WILLY L18PLT dataset not bundled")
def test_routes_to_real_correction_tools(no_llm_kw):
    agent = DistortionClassificationAgent(
        **no_llm_kw, epochs=40, route_corrections=True
    )
    result = agent.execute({"path": str(_L18PLT)})

    assert result.status == "success"
    counts = result["regime_counts"]
    if counts.get("static_shift_only", 0) > 0:
        assert result["ss_corrected_sites"] is not None
    else:
        assert result["ss_corrected_sites"] is None
    if counts.get("distorted", 0) > 0:
        assert result["distorted_corrected_sites"] is not None
    else:
        assert result["distorted_corrected_sites"] is None


@pytest.mark.skipif(not _HAS_L18, reason="WILLY L18PLT dataset not bundled")
def test_route_corrections_false_skips_routing(no_llm_kw):
    agent = DistortionClassificationAgent(
        **no_llm_kw, epochs=40, route_corrections=False
    )
    result = agent.execute({"path": str(_L18PLT)})

    assert result.status == "success"
    assert result["ss_corrected_sites"] is None
    assert result["distorted_corrected_sites"] is None


def test_too_few_stations_fails_cleanly(no_llm_kw, monkeypatch):
    import pandas as pd

    monkeypatch.setattr(
        "pycsamt.emtools._core.ensure_sites", lambda *a, **k: object()
    )
    monkeypatch.setattr(
        "pycsamt.ai.processing.build_distortion_features_table",
        lambda sites: pd.DataFrame(
            {"delta_log10_rho": [0.1], "twist_deg": [1.0], "shear": [0.0]}
        ),
    )
    agent = DistortionClassificationAgent(**no_llm_kw)
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert ">= 5" in result.error
