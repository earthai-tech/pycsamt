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
        topography=False,
        show_terrain=False,
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
        # Keep both available periods (0.1 s, 1.0 s) — go.Isosurface
        # needs at least two depth samples to build anything.
        period_range=(0.05, 2.0),
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


class _Sites4:
    def as_list(self):
        return [_Edi("S00"), _Edi("S01"), _Edi("S02"), _Edi("S03")]


def test_block_mode_ignores_azimuth() -> None:
    """go.Volume needs an axis-aligned rectilinear (x, y, z) grid to
    reconstruct a smooth block; rotating it by azimuth would shear
    that grid and break the reconstruction, so — unlike fence/depth
    mode — block mode deliberately ignores azimuth, matching
    ``pycsamt.app.web``'s block-mode 3-D map.
    """
    data = MapData(
        sites=_Sites4(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.0, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.1, 10.0, "L2", 2),
            StationRecord("S03", 1.1, 2.1, 20.0, "L2", 3),
        ),
    )
    base = build_3d_map(data, VolumeMapOptions(mode="block"))
    rotated = build_3d_map(data, VolumeMapOptions(mode="block", azimuth=90.0))
    assert base.data and rotated.data
    np.testing.assert_allclose(
        np.asarray(base.data[0].x), np.asarray(rotated.data[0].x)
    )
    np.testing.assert_allclose(
        np.asarray(base.data[0].y), np.asarray(rotated.data[0].y)
    )


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
        VolumeMapOptions(
            mode="fence", log_color=True, value_range=(10.0, 1000.0)
        ),
    )
    lin = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence", log_color=False, value_range=(10.0, 1000.0)
        ),
    )
    assert np.isclose(log.data[0].cmin, 1.0)
    assert np.isclose(log.data[0].cmax, 3.0)
    assert np.isclose(lin.data[0].cmin, 10.0)
    assert np.isclose(lin.data[0].cmax, 1000.0)


def test_show_stations_adds_marker_trace() -> None:
    fig = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence",
            show_stations=True,
            station_symbol="circle",
            station_size=6,
        ),
    )
    markers = [t for t in fig.data if getattr(t, "name", "") == "stations"]
    assert len(markers) == 1
    assert markers[0].marker.symbol == "circle"
    assert len(markers[0].x) == 2  # two stations


def test_axis_units_scale_independently() -> None:
    import numpy as np

    def xspan(fig):
        xs = [
            v
            for t in fig.data
            if type(t).__name__ == "Surface"
            for v in (np.nanmin(t.x), np.nanmax(t.x))
        ]
        return max(xs) - min(xs)

    def zspan(fig):
        zs = [
            v
            for t in fig.data
            if type(t).__name__ == "Surface"
            for v in (np.nanmin(t.z), np.nanmax(t.z))
        ]
        return max(zs) - min(zs)

    m = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", x_unit="m", depth_unit="m"),
    )
    km = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", x_unit="km", depth_unit="km"),
    )
    # km view is 1000x smaller numerically than the metre view
    assert np.isclose(xspan(m), xspan(km) * 1000.0, rtol=1e-3)
    assert np.isclose(zspan(m), zspan(km) * 1000.0, rtol=1e-3)
    assert "(m)" in m.layout.scene.xaxis.title.text
    assert "(km)" in km.layout.scene.zaxis.title.text


def test_topography_shifts_surface_and_adds_terrain() -> None:
    flat = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence", topography=False, show_terrain=False, depth_unit="m"
        ),
    )
    topo = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence", topography=True, show_terrain=True, depth_unit="m"
        ),
    )
    flat_mean = float(np.nanmean(np.asarray(flat.data[0].z)))
    topo_mean = float(np.nanmean(np.asarray(topo.data[0].z)))
    # mean elevation is 150 m -> the draped surface sits ~150 m higher
    assert topo_mean - flat_mean > 100.0
    assert any("terrain" in str(getattr(t, "name", "")) for t in topo.data)
    assert "Elevation" in topo.layout.scene.zaxis.title.text


class _WideZ:
    freq = [1000.0, 100.0, 10.0, 1.0]

    def __init__(self, base: float) -> None:
        rng = np.random.RandomState(int(base))
        self.resistivity = np.abs(rng.rand(4, 2, 2)) * 400.0 + base
        self.phase = np.ones((4, 2, 2)) * 45.0


class _WideEdi:
    def __init__(self, station: str, base: float) -> None:
        self.station = station
        self.Z = _WideZ(base)


def _wide_data() -> MapData:
    """A single line with 4 stations / 4 periods — big enough to
    trigger the fence-panel smoothing/resampling path."""
    stations = tuple(
        StationRecord(f"S0{i}", 1.0, 2.0 + i * 0.01, 10.0 * i, "L1", i)
        for i in range(4)
    )
    edis = [_WideEdi(f"S0{i}", 10.0 + i) for i in range(4)]
    return MapData(sites=edis, stations=stations)


def test_smooth_sections_upsamples_and_stays_in_range() -> None:
    raw = build_3d_map(
        _wide_data(),
        VolumeMapOptions(mode="fence", smooth_sections=False),
    )
    smooth = build_3d_map(
        _wide_data(),
        VolumeMapOptions(mode="fence", smooth_sections=True, section_res=40),
    )
    raw_x = np.asarray(raw.data[0].x)
    smooth_x = np.asarray(smooth.data[0].x)
    assert smooth_x.shape[1] > raw_x.shape[1]
    assert smooth_x.shape[0] > raw_x.shape[0]

    raw_v = np.asarray(raw.data[0].surfacecolor)
    smooth_v = np.asarray(smooth.data[0].surfacecolor)
    # cubic-spline smoothing must not overshoot the source data range
    assert np.nanmin(smooth_v) >= np.nanmin(raw_v) - 1e-6
    assert np.nanmax(smooth_v) <= np.nanmax(raw_v) + 1e-6


def test_smooth_sections_disabled_keeps_raw_grid() -> None:
    fig = build_3d_map(
        _wide_data(),
        VolumeMapOptions(mode="fence", smooth_sections=False),
    )
    # 4 stations, 4 periods -> the raw (unsmoothed) grid shape
    assert np.asarray(fig.data[0].x).shape == (4, 4)
