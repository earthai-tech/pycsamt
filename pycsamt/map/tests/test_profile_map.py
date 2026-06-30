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
