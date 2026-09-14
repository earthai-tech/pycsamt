from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.api.station import (
    PYCSAMT_STATION_RENDERING,
    PyCSAMTStationRendering,
    StationAxisStyle,
    StationMarkerStyle,
    configure_station_rendering,
    reset_station_rendering,
)


# ─────────────────────────────────────────────────────────────────────────
# StationAxisStyle.compute_every / label_indices
# ─────────────────────────────────────────────────────────────────────────


def test_compute_every_explicit_int_overrides_auto():
    style = StationAxisStyle(every=3)
    assert style.compute_every(100) == 3


def test_compute_every_long_labels_increase_needed_step():
    style = StationAxisStyle(max_labels=5, every="auto")
    short = style.compute_every(200, figwidth_in=10.0, max_label_len=2)
    long = style.compute_every(200, figwidth_in=10.0, max_label_len=20)
    assert long >= short


def test_compute_every_falls_back_to_ceil_beyond_nice_steps():
    style = StationAxisStyle(max_labels=1, every="auto")
    step = style.compute_every(100000, figwidth_in=0.01, max_label_len=2)
    assert step > 100


def test_label_indices_empty_when_no_labels():
    style = StationAxisStyle()
    idx = style.label_indices([])
    assert idx.size == 0


# ─────────────────────────────────────────────────────────────────────────
# StationAxisStyle.apply: validation and empty input
# ─────────────────────────────────────────────────────────────────────────


def test_apply_rejects_invalid_side():
    fig, ax = plt.subplots()
    style = StationAxisStyle(side="left")
    with pytest.raises(ValueError):
        style.apply(ax, [1, 2, 3])
    plt.close(fig)


def test_apply_returns_empty_for_no_positions():
    fig, ax = plt.subplots()
    style = StationAxisStyle()
    idx = style.apply(ax, [])
    assert idx.size == 0
    plt.close(fig)


def test_apply_sets_xlim_when_given():
    fig, ax = plt.subplots()
    style = StationAxisStyle()
    style.apply(ax, np.arange(5), xlim=(-1.0, 10.0))
    assert ax.get_xlim() == (-1.0, 10.0)
    plt.close(fig)


def test_apply_hides_labels_when_show_labels_false():
    fig, ax = plt.subplots()
    style = StationAxisStyle(show_labels=False)
    idx = style.apply(ax, np.arange(5), [f"S{i}" for i in range(5)])
    ticklabels = [t.get_text() for t in ax.get_xticklabels()]
    assert all(text == "" for text in ticklabels)
    assert len(idx) == 5
    plt.close(fig)


def test_apply_bottom_side_configures_bottom_axis():
    fig, ax = plt.subplots()
    style = StationAxisStyle(side="bottom")
    style.apply(ax, np.arange(5), [f"S{i}" for i in range(5)])
    assert ax.xaxis.get_label_position() == "bottom"
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
# StationAxisStyle.apply: topo-surface mode
# ─────────────────────────────────────────────────────────────────────────


def test_apply_topo_elev_mode_draws_markers_and_labels():
    fig, ax = plt.subplots()
    ax.set_ylim(0, 100)
    style = StationAxisStyle(side="top")
    positions = np.arange(5, dtype=float)
    elev = np.array([10.0, 20.0, 15.0, 25.0, 30.0])
    idx = style.apply(
        ax, positions, [f"S{i}" for i in range(5)], topo_elev=elev,
    )
    assert len(idx) == 5
    assert ax.collections  # markers drawn via scatter
    assert ax.texts  # labels drawn inline
    plt.close(fig)


def test_apply_topo_elev_mode_inverted_axis_and_no_markers_labels():
    fig, ax = plt.subplots()
    ax.set_ylim(100, 0)  # inverted y-axis (toward_top < 0 branch)
    style = StationAxisStyle(
        side="top", show_markers=False, show_labels=False,
    )
    positions = np.arange(3, dtype=float)
    elev = np.array([10.0, 20.0, 15.0])
    idx = style.apply(ax, positions, topo_elev=elev)
    assert len(idx) == 3
    assert not ax.collections
    assert not ax.texts
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
# PyCSAMTStationRendering
# ─────────────────────────────────────────────────────────────────────────


def test_style_for_rejects_unknown_preset():
    renderer = PyCSAMTStationRendering()
    with pytest.raises(ValueError):
        renderer.style_for("bogus")


def test_renderer_apply_delegates_to_style_for():
    fig, ax = plt.subplots()
    renderer = PyCSAMTStationRendering()
    idx = renderer.apply(ax, np.arange(4), preset="survey")
    assert len(idx) == 4
    plt.close(fig)


def test_use_preset_copies_into_pseudosection_slot():
    renderer = PyCSAMTStationRendering()
    renderer.use_preset("survey")
    assert renderer.pseudosection.side == renderer.survey.side
    assert renderer.pseudosection is not renderer.survey
    renderer.reset()


def test_context_with_preset_uses_and_restores():
    renderer = PyCSAMTStationRendering()
    original_side = renderer.pseudosection.side
    with renderer.context(preset="survey"):
        assert renderer.pseudosection.side == renderer.survey.side
    assert renderer.pseudosection.side == original_side


def test_reset_restores_defaults_after_configure():
    renderer = PyCSAMTStationRendering()
    renderer.configure(pseudosection__max_labels=1)
    assert renderer.pseudosection.max_labels == 1
    renderer.reset()
    assert renderer.pseudosection.max_labels == 14


def test_summary_and_repr_contain_preset_names():
    renderer = PyCSAMTStationRendering()
    text = renderer.summary()
    assert "pseudosection" in text
    assert "inversion" in text
    assert "survey" in text
    assert repr(renderer) == text


def test_module_level_configure_and_reset_helpers():
    old = PYCSAMT_STATION_RENDERING.pseudosection.max_labels
    try:
        configure_station_rendering(pseudosection__max_labels=2)
        assert PYCSAMT_STATION_RENDERING.pseudosection.max_labels == 2
    finally:
        reset_station_rendering()
    assert PYCSAMT_STATION_RENDERING.pseudosection.max_labels == old
