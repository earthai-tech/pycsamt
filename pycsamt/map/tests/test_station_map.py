# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for station maps."""

from __future__ import annotations

from pycsamt.map import StationMapOptions, build_station_map
from pycsamt.map._core import MapData, StationRecord


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
        return [_Edi("S00"), _Edi("S01"), _Edi("S02")]


def test_station_map_options_defaults() -> None:
    opts = StationMapOptions()
    assert opts.overlay == "index"
    assert opts.backend == "plotly"
    assert opts.show_profiles is True
    assert opts.elevation_mode == "markers"


def test_station_map_builds_plotly_figure() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.2, 30.0, "L1", 2),
        ),
    )
    fig = build_station_map(data, StationMapOptions())
    assert fig.data


def test_station_map_depth_and_density_overlay() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.2, 30.0, "L1", 2),
        ),
    )
    opts = StationMapOptions(
        overlay="skin_depth",
        frequency=10.0,
        show_contours=True,
        log_color=True,
    )
    fig = build_station_map(data, opts)
    assert len(fig.data) >= 2


def test_station_map_matplotlib_backend() -> None:
    data = MapData(sites=None)
    opts = StationMapOptions(backend="matplotlib")
    fig = build_station_map(data, opts)
    assert hasattr(fig, "savefig")


def test_station_elevation_contour_mode() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.2, 30.0, "L1", 2),
        ),
    )
    opts = StationMapOptions(
        overlay="elevation",
        elevation_mode="contours",
        backend="matplotlib",
        contour_interp="linear",
    )
    fig = build_station_map(data, opts)
    assert len(fig.axes[0].collections) > 1

    plotly = build_station_map(
        data,
        StationMapOptions(
            overlay="elevation",
            elevation_mode="contours",
            contour_interp="linear",
        ),
    )
    assert len(plotly.layout.map.layers) == 1


def test_station_elevation_mode_rejects_unknown_value() -> None:
    import pytest

    with pytest.raises(ValueError, match="elevation_mode"):
        build_station_map(
            MapData(sites=None),
            StationMapOptions(elevation_mode="raster"),
        )
