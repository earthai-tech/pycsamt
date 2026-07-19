# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for profile maps."""

from __future__ import annotations

from pycsamt.map import ProfileMapOptions
from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.profile import (
    build_profile_map,
    build_pseudosection,
)


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        import numpy as np

        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


class _Sites:
    def as_list(self):
        return [_Edi("S00"), _Edi("S01")]


def test_profile_map_options_defaults() -> None:
    opts = ProfileMapOptions()
    assert opts.quantity == "rho"
    assert opts.components == ("xy", "yx")


def test_profile_renderers_build_figures() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = ProfileMapOptions(components=("xy",))
    assert build_profile_map(data, opts).data
    assert build_pseudosection(data, opts).data


def test_profile_options_filter_and_distance_axis() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = ProfileMapOptions(
        components=("xy",),
        period_range=(0.01, 2.0),
        value_range=(1.0, 500.0),
        x_axis="distance",
    )
    fig = build_pseudosection(data, opts)
    assert fig.data


def test_profile_matplotlib_backend_builds_figure() -> None:
    data = MapData(sites=_Sites())
    opts = ProfileMapOptions(backend="matplotlib")
    fig = build_pseudosection(data, opts)
    assert hasattr(fig, "savefig")


class _MultiLineSites:
    def as_list(self):
        return [
            _Edi("A00"),
            _Edi("A01"),
            _Edi("B00"),
            _Edi("B01"),
        ]


def _multi_line_data() -> MapData:
    return MapData(
        sites=_MultiLineSites(),
        stations=(
            StationRecord("A00", 1.0, 2.0, 10.0, "LA", 0),
            StationRecord("A01", 1.1, 2.1, 20.0, "LA", 1),
            StationRecord("B00", 1.2, 2.2, 30.0, "LB", 2),
            StationRecord("B01", 1.3, 2.3, 40.0, "LB", 3),
        ),
    )


def test_by_line_pseudosection_builds_one_panel_per_line() -> None:
    data = _multi_line_data()
    opts = ProfileMapOptions(component="xy", by_line=True)
    fig = build_pseudosection(data, opts)
    assert fig.data
    # One heatmap trace per line, each keeping its own station axis.
    assert len(fig.data) == 2
    assert {tuple(trace.x) for trace in fig.data} == {
        ("A00", "A01"),
        ("B00", "B01"),
    }


def test_by_line_pseudosection_matplotlib_backend() -> None:
    data = _multi_line_data()
    opts = ProfileMapOptions(
        component="xy",
        by_line=True,
        backend="matplotlib",
    )
    fig = build_pseudosection(data, opts)
    assert hasattr(fig, "savefig")


def test_by_line_pseudosection_single_line_falls_back() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = ProfileMapOptions(component="xy", by_line=True)
    fig = build_pseudosection(data, opts)
    assert len(fig.data) == 1
