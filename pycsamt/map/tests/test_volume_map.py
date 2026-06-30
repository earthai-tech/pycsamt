# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for volume maps."""

from __future__ import annotations

import numpy as np

from pycsamt.map import VolumeMapOptions
from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.volume import build_3d_map


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


class _VarZ:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.array(
            [
                [[100.0, 100.0], [100.0, 100.0]],
                [[100.0, 10000.0], [100.0, 100.0]],
            ]
        )
        self.phase = np.ones((2, 2, 2)) * 45.0


class _VarEdi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _VarZ()


class _VarSites:
    def as_list(self):
        return [_VarEdi("S00"), _VarEdi("S01")]


def test_volume_map_options_defaults() -> None:
    opts = VolumeMapOptions()
    assert opts.mode == "fence"
    assert opts.quantity == "resistivity"


def test_volume_renderers_build_figures() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    assert build_3d_map(data, VolumeMapOptions()).data
    block = build_3d_map(data, VolumeMapOptions(mode="block"))
    depth = build_3d_map(data, VolumeMapOptions(mode="depth"))
    assert block.data
    assert depth.data


def test_volume_slice_and_surface_controls() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = VolumeMapOptions(
        mode="depth",
        n_slices=3,
        surface_count=4,
        show_contours=True,
        azimuth=30.0,
    )
    fig = build_3d_map(data, opts)
    assert len(fig.data) == 3


def test_volume_surface_and_period_filter() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = VolumeMapOptions(
        mode="surface",
        period_range=(0.01, 0.2),
        surface_count=2,
    )
    fig = build_3d_map(data, opts)
    assert fig.data


def test_volume_can_color_by_phase() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = VolumeMapOptions(
        mode="fence",
        quantity="phase",
        show_labels=False,
    )
    fig = build_3d_map(data, opts)
    assert fig.data


def test_depth_mode_honors_rho_range() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    opts = VolumeMapOptions(
        mode="depth",
        n_slices=2,
        rho_range=(1.0, 500.0),
        show_contours=False,
    )
    fig = build_3d_map(data, opts)
    colors = np.asarray(
        fig.data[-1].surfacecolor,
        dtype=float,
    )
    assert np.isnan(colors).any()


def test_block_mode_honors_azimuth() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L2", 1),
        ),
    )
    opts = VolumeMapOptions(mode="block", azimuth=90.0)
    fig = build_3d_map(data, opts)
    assert max(fig.data[0].x) > 1.0
    assert max(abs(y) for y in fig.data[0].y) < 1e-6


def _var_data() -> MapData:
    return MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 100.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 200.0, "L1", 1),
        ),
    )


def test_log_vs_linear_colorbar_range() -> None:
    log = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", log_color=True,
                         value_range=(10.0, 1000.0)),
    )
    lin = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", log_color=False,
                         value_range=(10.0, 1000.0)),
    )
    assert np.isclose(log.data[0].cmin, 1.0)
    assert np.isclose(log.data[0].cmax, 3.0)
    assert np.isclose(lin.data[0].cmin, 10.0)
    assert np.isclose(lin.data[0].cmax, 1000.0)


def test_show_stations_adds_marker_trace() -> None:
    fig = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", show_stations=True,
                         station_symbol="circle", station_size=6),
    )
    markers = [t for t in fig.data if getattr(t, "name", "") == "stations"]
    assert len(markers) == 1
    assert markers[0].marker.symbol == "circle"
    assert len(markers[0].x) == 2  # two stations


def test_topography_shifts_surface_and_adds_terrain() -> None:
    flat = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", topography=False,
                         show_terrain=False),
    )
    topo = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", topography=True,
                         show_terrain=True),
    )
    flat_mean = float(np.nanmean(np.asarray(flat.data[0].z)))
    topo_mean = float(np.nanmean(np.asarray(topo.data[0].z)))
    # mean elevation is 150 m -> the draped surface sits ~150 m higher
    assert topo_mean - flat_mean > 100.0
    assert any(
        "terrain" in str(getattr(t, "name", "")) for t in topo.data
    )
    assert "Elevation" in topo.layout.scene.zaxis.title.text
