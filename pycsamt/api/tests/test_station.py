from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pycsamt.api.station import (
    PYCSAMT_STATION_RENDERING,
    StationAxisStyle,
    StationMarkerStyle,
)


def test_station_axis_auto_decimates_labels():
    style = StationAxisStyle(max_labels=5)
    idx = style.label_indices([f"S{i:02d}" for i in range(20)], 8.0)

    assert 0 in idx
    assert 19 in idx
    assert len(idx) <= 7


def test_station_axis_applies_top_ticks_and_markers():
    fig, ax = plt.subplots()
    style = StationAxisStyle(
        side="top",
        label_pad=11.0,
        marker=StationMarkerStyle(marker="v", facecolor="white"),
    )

    idx = style.apply(ax, np.arange(5), [f"S{i}" for i in range(5)])

    assert len(idx) == 5
    assert ax.xaxis.get_label_position() == "top"
    assert ax.collections
    assert style.label_pad == 11.0
    plt.close(fig)


def test_station_rendering_context_restores_value():
    old = PYCSAMT_STATION_RENDERING.pseudosection.max_labels

    with PYCSAMT_STATION_RENDERING.context(
        pseudosection__max_labels=3,
    ):
        assert PYCSAMT_STATION_RENDERING.pseudosection.max_labels == 3

    assert PYCSAMT_STATION_RENDERING.pseudosection.max_labels == old


def test_inversion_preset_uses_filled_black_triangle():
    marker = PYCSAMT_STATION_RENDERING.inversion.marker

    assert marker.marker == "v"
    assert marker.facecolor == "black"
