"""Focused coverage for orchestration and joint-inversion agent helpers."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import numpy as np

matplotlib.use("Agg")


def test_inversion_prep_occam_modem_and_unknown(tmp_path, monkeypatch):
    from pycsamt.agents.inversion_prep import InversionPrepAgent
    import pycsamt.emtools._core as core
    import pycsamt.models.occam2d as occam

    ed = object()
    monkeypatch.setattr(core, "ensure_sites", lambda value, verbose=0: [ed])
    monkeypatch.setattr(core, "_iter_items", iter)
    monkeypatch.setattr(core, "_get_z_block", lambda item: (None, None, np.array([1., 10., 0.])))
    monkeypatch.setattr(
        occam, "write_occam2d_data", lambda *a, **k: None, raising=False
    )
    agent = InversionPrepAgent(api_key="key")
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: "recommendations")
    result = agent.execute({"sites": [ed], "output_dir": str(tmp_path), "period_range": [.01, 2]})
    assert result.status == "success" and result.data["n_periods"] == 2
    assert result.llm_interpretation == "recommendations"
    assert InversionPrepAgent(code="modem").execute({"sites": [ed], "output_dir": str(tmp_path)}).status == "needs_review"
    unknown = InversionPrepAgent(code="other").execute({"sites": [ed], "output_dir": str(tmp_path)})
    assert unknown.status == "needs_review" and unknown.warnings


def test_joint_secondary_features_and_plot(monkeypatch):
    from pycsamt.agents.joint_agent import (
        _collect_secondary_features, _extract_sec_features,
        _pad_or_trim, _plot_joint_section,
    )
    import pycsamt.emtools._core as core

    z = np.ones((4, 2, 2), complex)
    z[:, 0, 1] = np.array([1+1j, 2+1j, 3+1j, 4+1j])
    freqs = np.array([1., 2., 4., 8.])
    assert _extract_sec_features(z, freqs, np.array([1., 4.])).shape == (4,)
    assert _extract_sec_features(z[:1], freqs[:1], np.array([1.])) is None
    monkeypatch.setattr(core, "_iter_items", iter)
    monkeypatch.setattr(core, "_get_z_block", lambda item: (None, z, freqs))
    warnings = []
    primary = [object(), object()]
    arr = _collect_secondary_features(
        ["A", "B"], primary, None, freqs, np.array([1., 4.]), 4, warnings,
    )
    assert arr.shape == (2, 4) and warnings
    arr2 = _collect_secondary_features(
        ["A", "B"], primary, [object()], freqs, np.array([1., 4.]), 4, [],
    )
    assert arr2.shape == (2, 4)
    assert _pad_or_trim(np.ones((2, 2)), 2).shape == (2, 2)
    assert _pad_or_trim(np.ones((2, 3)), 2).shape == (2, 2)
    assert _pad_or_trim(np.ones((2, 2)), 4).shape == (2, 4)
    predictions = {"A": np.array([1., 2., 3.]), "B": np.array([2., 3., 4.])}
    assert _plot_joint_section(predictions, 3, freqs, ["A", "B"], ["mt", "tem"]) is not None
    assert _plot_joint_section({}, 3, freqs, [], ["mt"]) is None


def test_joint_forward_rms(monkeypatch):
    from pycsamt.agents.joint_agent import _forward_rms_joint
    import pycsamt.emtools._core as core
    import pycsamt.forward as forward

    z = np.ones((3, 2, 2), complex)
    ed = SimpleNamespace(rho=np.full((3, 2, 2), 10.))
    monkeypatch.setattr(core, "_iter_items", iter)
    monkeypatch.setattr(core, "_name", lambda item, i: "A")
    monkeypatch.setattr(core, "_get_z_block", lambda item: (None, z, np.array([1., 10., 100.])))
    monkeypatch.setattr(forward, "LayeredModel", lambda **kwargs: kwargs)
    class Fwd:
        def __init__(self, **kwargs): pass
        def run(self, model): return SimpleNamespace(rho_a=np.full((3, 2, 2), 10.))
    monkeypatch.setattr(forward, "MT1DForward", Fwd)
    rms = _forward_rms_joint([ed], "A", 0, np.array([1., 2., 3.]), np.array([1., 10., 100.]), 3)
    assert rms < 1e-12
