"""Synthetic tipper tests independent of optional EDI fixtures."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.agents.tipper_analysis import (
    TipperAnalysisAgent,
    _plot_induction_arrows,
    _plot_tipper_amplitude,
    _plot_tipper_pseudosection,
)


def _station(tipper, freqs=(10.0, 1.0)):
    return SimpleNamespace(
        station="S01",
        freqs=np.asarray(freqs),
        Tip=SimpleNamespace(tipper=tipper),
    )


def test_synthetic_tipper_execute(monkeypatch, tmp_path):
    import pycsamt.emtools._core as core

    tipper = np.array([[[0.1 + 0.2j, 0.3 + 0.1j]], [[0.2 + 0.1j, 0.1 + 0.4j]]])
    station = _station(tipper)
    monkeypatch.setattr(core, "ensure_sites", lambda value, **_: value)
    monkeypatch.setattr(core, "_iter_items", lambda sites: iter(sites))
    monkeypatch.setattr(core, "_name", lambda ed, _i: ed.station)
    monkeypatch.setattr(
        core,
        "_get_z_block",
        lambda ed: (None, np.ones((2, 2, 2), dtype=complex), ed.freqs),
    )

    agent = TipperAnalysisAgent(period_ref=0.1)
    result = agent.execute({"sites": [station], "output_dir": str(tmp_path)})
    assert result.status == "success"
    assert result["n_stations_with_tipper"] == 1
    assert len(result["figures"]) == 3
    assert result["figure_paths"]


def test_tipper_skips_bad_stations(monkeypatch):
    import pycsamt.emtools._core as core

    stations = [
        _station(None),
        _station(np.ones((3, 2), dtype=complex)),
        _station(np.ones((2, 2), dtype=complex), freqs=()),
    ]
    monkeypatch.setattr(core, "ensure_sites", lambda value, **_: value)
    monkeypatch.setattr(core, "_iter_items", lambda sites: iter(sites))
    monkeypatch.setattr(core, "_name", lambda ed, i: f"S{i}")
    monkeypatch.setattr(
        core,
        "_get_z_block",
        lambda ed: (None, np.empty((len(ed.freqs), 2, 2)), ed.freqs),
    )
    result = TipperAnalysisAgent().execute({"sites": stations})
    assert result.status == "failed"
    assert "No tipper" in result.error


def test_plot_helpers_empty_and_populated():
    assert _plot_tipper_amplitude([], []) is None
    assert _plot_induction_arrows([], None, "parkinson") is None
    assert _plot_tipper_pseudosection([], []) is None
    rows = [
        {
            "station": "A",
            "period_s": 1.0,
            "amplitude": 0.2,
            "arrow_x": 0.1,
            "arrow_y": -0.1,
        },
        {
            "station": "A",
            "period_s": 10.0,
            "amplitude": 0.4,
            "arrow_x": 0.2,
            "arrow_y": 0.1,
        },
    ]
    figs = [
        _plot_tipper_amplitude(rows, ["A", "missing"]),
        _plot_induction_arrows(rows[:1], 1.0, "weise"),
        _plot_tipper_pseudosection(rows, ["A"]),
    ]
    assert all(fig is not None for fig in figs)
    for fig in figs:
        plt.close(fig)
