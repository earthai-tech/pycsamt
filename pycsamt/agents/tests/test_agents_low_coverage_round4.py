"""Additional edge coverage for analysis and mapping agents."""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")


def test_resistivity_map_interpolation_and_plot():
    from pycsamt.agents.resistivity_map import _interpolate_map, _plot_depth_maps

    cx = np.array([0., 1000., 0., 1000.])
    cy = np.array([0., 0., 1000., 1000.])
    vals = np.array([1., 2., 3., 4.])
    axis = np.linspace(0, 1000, 4)
    Xg, Yg = np.meshgrid(axis, axis)
    assert _interpolate_map(cx, cy, vals, Xg, Yg, "idw").shape == Xg.shape
    linear = _interpolate_map(cx, cy, vals, Xg, Yg, "linear")
    depth_map = {
        "depth_km": 1.0, "grid_x": axis, "grid_y": axis,
        "grid_rho": linear, "station_x": cx, "station_y": cy,
        "station_rho": vals,
    }
    assert _plot_depth_maps([], []) is None
    assert _plot_depth_maps([depth_map], ["A", "B", "C", "D"]) is not None


def test_sensitivity_bar_empty_and_populated():
    from pycsamt.agents.sensitivity import _plot_doi_bar

    assert _plot_doi_bar({}) is None
    assert _plot_doi_bar({"A": 1000., "B": 2500.}) is not None


def test_qc_gracefully_collects_all_tool_failures(monkeypatch):
    from pycsamt.agents.qc import DataQCAgent
    import pycsamt.emtools._core as core
    import pycsamt.emtools.qc as qc

    monkeypatch.setattr(core, "ensure_sites", lambda value, verbose=0: value)

    def fail(*args, **kwargs):
        raise RuntimeError("unavailable")

    for name in (
        "build_qc_table", "qc_flags", "station_confidence_table",
        "frequency_confidence_table", "plot_frequency_confidence_psection",
        "plot_confidence_profile",
    ):
        monkeypatch.setattr(qc, name, fail)
    result = DataQCAgent().execute({"sites": object()})
    assert result.status == "success"
    assert len(result.warnings) == 6 and not result.data["figures"]


def test_qc_bad_sites(monkeypatch):
    from pycsamt.agents.qc import DataQCAgent
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core, "ensure_sites",
        lambda value, verbose=0: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    result = DataQCAgent().execute({"sites": object()})
    assert result.status == "failed" and "bad sites" in result.error
