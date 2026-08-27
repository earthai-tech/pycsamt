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


class _ThreeLineZ:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _ThreeLineEdi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _ThreeLineZ()


class _ThreeLineSites:
    def as_list(self):
        return [_ThreeLineEdi(s) for s in ("A0", "A1", "B0", "B1", "C0", "C1")]


def test_depth_mode_stacks_lines_in_real_offset_order() -> None:
    """``go.Surface`` interpolates between *adjacent* rows of the
    stacked (x, y, z) grid -- lines must be vstacked in real
    cross-strike (y-offset) order, not dict/insertion order, or the
    reconstructed surface connects unrelated lines and twists on
    itself. Regression for a real bug: with real geometry where line
    "B" sits farthest out and "C" sits in between (offset order
    A < C < B), the old code stacked rows in insertion order
    (A, B, C) instead.
    """
    data = MapData(
        sites=_ThreeLineSites(),
        stations=(
            StationRecord("A0", 0.0, 0.0, 10.0, "A", 0),
            StationRecord("A1", 0.001, 0.0, 10.0, "A", 1),
            StationRecord("B0", 0.0, 0.02, 10.0, "B", 2),
            StationRecord("B1", 0.001, 0.02, 10.0, "B", 3),
            StationRecord("C0", 0.0, 0.01, 10.0, "C", 4),
            StationRecord("C1", 0.001, 0.01, 10.0, "C", 5),
        ),
    )
    opts = VolumeMapOptions(
        mode="depth", n_slices=1, show_terrain=False, show_contours=False
    )
    fig = build_3d_map(data, opts)
    y_rows = np.asarray(fig.data[0].y, dtype=float)[:, 0]
    assert np.all(np.diff(y_rows) >= 0) or np.all(np.diff(y_rows) <= 0)


class _UnorderedZ:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _UnorderedEdi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _UnorderedZ()


class _UnorderedSites:
    def as_list(self):
        return [_UnorderedEdi(s) for s in ("A", "B", "C")]


def test_depth_mode_sorts_columns_by_real_chainage_within_a_line() -> None:
    """``go.Surface`` also interpolates between *adjacent columns* of
    one line's own row -- those must be in real along-profile (x)
    order too, not whatever order the source pivot table's columns
    happen to come in (e.g. lexical by station name). Regression for
    a real bug: station names sort alphabetically A, B, C, but their
    real coordinates place B farthest along the line and C in the
    middle (real spatial order A, C, B) -- the old code left the row
    in the pivot table's A, B, C column order instead.
    """
    data = MapData(
        sites=_UnorderedSites(),
        stations=(
            StationRecord("A", 0.000, 0.0, 10.0, "L22", 0),
            StationRecord("B", 0.003, 0.0, 10.0, "L22", 1),
            StationRecord("C", 0.001, 0.0, 10.0, "L22", 2),
        ),
    )
    opts = VolumeMapOptions(
        mode="depth", n_slices=1, show_terrain=False, show_contours=False
    )
    fig = build_3d_map(data, opts)
    x_row = np.asarray(fig.data[0].x, dtype=float)[0]
    assert np.all(np.diff(x_row) >= 0) or np.all(np.diff(x_row) <= 0)


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


def test_fence_and_depth_colorbar_range_is_stable_under_rho_range_filter() -> None:
    """Narrowing ``rho_range`` must only change which cells
    :func:`_filtered_values` masks to NaN (what is shown), never
    cmin/cmax (the colour scale they are shown on). Regression for a
    real bug: a narrow selection (e.g. a "conductive" 1-100 ohm.m
    slice) used to get stretched across the whole colourscale and
    render as one near-solid colour instead of the same hues it had
    before filtering.
    """
    data = _var_data()
    unfiltered_fence = build_3d_map(data, VolumeMapOptions(mode="fence"))
    filtered_fence = build_3d_map(
        data, VolumeMapOptions(mode="fence", rho_range=(1.0, 500.0))
    )
    assert np.isclose(unfiltered_fence.data[0].cmin, filtered_fence.data[0].cmin)
    assert np.isclose(unfiltered_fence.data[0].cmax, filtered_fence.data[0].cmax)
    # the filter must still actually do something: some cells masked out
    assert np.isnan(
        np.asarray(filtered_fence.data[0].surfacecolor, dtype=float)
    ).any()

    unfiltered_depth = build_3d_map(
        data,
        VolumeMapOptions(
            mode="depth", n_slices=2, show_contours=False, show_terrain=False
        ),
    )
    filtered_depth = build_3d_map(
        data,
        VolumeMapOptions(
            mode="depth",
            n_slices=2,
            rho_range=(1.0, 500.0),
            show_contours=False,
            show_terrain=False,
        ),
    )
    assert np.isclose(unfiltered_depth.data[-1].cmin, filtered_depth.data[-1].cmin)
    assert np.isclose(unfiltered_depth.data[-1].cmax, filtered_depth.data[-1].cmax)


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


def test_block_mode_drapes_with_topography() -> None:
    """The dense volume grid must follow real terrain when topography
    is enabled (elevation varies station to station) -- regression for
    a real bug where block mode never applied any topography drape at
    all, unlike fence/depth mode (which already shift their panel by
    real station elevation).
    """
    data = MapData(
        sites=_Sites4(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 100.0, "L1", 0),
            StationRecord("S01", 1.1, 2.0, 300.0, "L1", 1),
            StationRecord("S02", 1.0, 2.1, 100.0, "L2", 2),
            StationRecord("S03", 1.1, 2.1, 300.0, "L2", 3),
        ),
    )
    draped = build_3d_map(data, VolumeMapOptions(mode="block", topography=True))
    flat = build_3d_map(data, VolumeMapOptions(mode="block", topography=False))
    assert draped.data and flat.data

    draped_z = np.asarray(draped.data[0].z, dtype=float)
    flat_z = np.asarray(flat.data[0].z, dtype=float)

    # The flat (no-topography) datum never goes above 0 (the surface).
    assert np.nanmax(flat_z) <= 0.0
    # The draped grid's top must reach real terrain highs (up to
    # 300 m), not stop at the flat datum.
    assert np.nanmax(draped_z) > 100.0

    # Draping introduces genuine air gaps above the local terrain. The
    # grid is kept whole for go.Volume, so those no-data cells are
    # pushed below isomin (rendered transparent) rather than dropped --
    # the draped block must carry more of them than the flat one.
    draped_values = np.asarray(draped.data[0].value, dtype=float)
    flat_values = np.asarray(flat.data[0].value, dtype=float)
    draped_hidden = (draped_values < float(draped.data[0].isomin)).sum()
    flat_hidden = (flat_values < float(flat.data[0].isomin)).sum()
    assert draped_hidden > flat_hidden


class _SpreadEdi:
    """4 stations / 2 lines with a wide resistivity spread (~30-3000
    ohm.m) so resistivity-band filters have something to clamp to."""

    def __init__(self, station: str, base: float) -> None:
        self.station = station

        class _Z:
            freq = [1000.0, 100.0, 10.0, 1.0]

            def __init__(self) -> None:
                self.resistivity = base * np.logspace(0, 2, 4)[:, None, None] * np.ones((4, 2, 2))
                self.phase = np.ones((4, 2, 2)) * 45.0

        self.Z = _Z()


def _spread_data() -> MapData:
    edis = [_SpreadEdi(f"L{li}S{si}", 30.0 + 20 * si)
            for li in range(2) for si in range(2)]

    class _S:
        def as_list(self):
            return edis

    stations = tuple(
        StationRecord(f"L{li}S{si}", 1.0 + si * 0.02, 2.0 + li * 0.02,
                      50.0, f"L{li}", li * 2 + si)
        for li in range(2) for si in range(2)
    )
    return MapData(sites=_S(), stations=stations)


def test_block_mode_visible_with_topography_and_no_nan_in_value() -> None:
    """Regression: block mode used to hand go.Volume a value array
    still carrying the NaN cells from the topography drape, and
    go.Volume then rendered nothing. The grid is now kept whole with
    no-data cells pushed below isomin instead.
    """
    fig = build_3d_map(_spread_data(), VolumeMapOptions(mode="block", topography=True))
    assert fig.data and type(fig.data[0]).__name__ == "Volume"
    val = np.asarray(fig.data[0].value, dtype=float)
    assert val.size > 0
    assert np.isfinite(val).all()  # no NaN reaches go.Volume
    t = fig.data[0]
    # isomin/isomax sit inside the colour span, never a degenerate point
    assert float(t.isomin) < float(t.isomax)
    assert float(t.cmin) < float(t.cmax)


def test_block_mode_rho_range_isolates_the_band_not_the_whole_model() -> None:
    """A resistivity-band selection renders one closed iso-surface body
    -- NOT a go.Volume of the whole model tinted by the band (which is
    what ``isomin``/``isomax`` alone would give, since the isomax surface
    wraps everything above it). One that misses the model entirely draws
    a note, not a blank scene.
    """
    data = _spread_data()
    full = build_3d_map(data, VolumeMapOptions(mode="block"))
    overlap = build_3d_map(
        data, VolumeMapOptions(mode="block", rho_range=(100.0, 1000.0))
    )
    assert overlap.data and type(overlap.data[0]).__name__ == "Isosurface"
    t = overlap.data[0]
    # the colour scale is untouched by the filter -- same as unfiltered
    assert np.isclose(t.cmin, full.data[0].cmin)
    assert np.isclose(t.cmax, full.data[0].cmax)
    # every renderable cell is inside the band; the rest sit at a
    # far-below sentinel so the isomax surface has nothing to wrap
    val = np.asarray(t.value, dtype=float)
    keep = val > float(t.isomin) - 1.0
    assert keep.any()
    assert val[keep].min() >= float(t.isomin) - 1e-6
    assert val[keep].max() <= float(t.isomax) + 1e-6
    assert (val <= float(t.isomin) - 1.0).any()  # sentinel cells exist

    disjoint = build_3d_map(
        data, VolumeMapOptions(mode="block", rho_range=(1e6, 1e7))
    )
    assert not disjoint.data  # no trace
    assert disjoint.layout.annotations  # ...but a note explaining why


def test_block_mode_station_markers_stay_axis_aligned() -> None:
    """Station markers must sit in the same unrotated ``(x, y)`` frame
    as block mode's own axis-aligned volume grid -- regression for a
    real bug where markers kept rotating by azimuth even though the
    block itself never does (``test_block_mode_ignores_azimuth``), so
    markers drifted outside the block's own footprint whenever
    ``azimuth != 0``.
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
    opts = VolumeMapOptions(
        mode="block", azimuth=90.0, show_stations=True, show_labels=False
    )
    fig = build_3d_map(data, opts)
    volume_trace, marker_trace = fig.data[0], fig.data[-1]

    x_lo, x_hi = float(np.min(volume_trace.x)), float(np.max(volume_trace.x))
    y_lo, y_hi = float(np.min(volume_trace.y)), float(np.max(volume_trace.y))
    margin = 1e-6 + 0.05 * max(x_hi - x_lo, 1.0)
    mx = np.asarray(marker_trace.x, dtype=float)
    my = np.asarray(marker_trace.y, dtype=float)
    assert mx.min() >= x_lo - margin
    assert mx.max() <= x_hi + margin
    assert my.min() >= y_lo - margin
    assert my.max() <= y_hi + margin


def test_surface_mode_inherits_topography_drape_and_axis_aligned_markers() -> None:
    """Isosurface (``mode="surface"``) shares block mode's dense
    ``(x, y, z)`` grid builder (:func:`_dense_volume_grid`), so both
    the topography-drape fix and the axis-aligned-marker fix must
    apply there too, not just to block mode.
    """
    data = MapData(
        sites=_Sites4(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 100.0, "L1", 0),
            StationRecord("S01", 1.1, 2.0, 300.0, "L1", 1),
            StationRecord("S02", 1.0, 2.1, 100.0, "L2", 2),
            StationRecord("S03", 1.1, 2.1, 300.0, "L2", 3),
        ),
    )
    draped = build_3d_map(
        data,
        VolumeMapOptions(mode="surface", topography=True, period_range=(0.05, 2.0)),
    )
    flat = build_3d_map(
        data,
        VolumeMapOptions(mode="surface", topography=False, period_range=(0.05, 2.0)),
    )
    assert draped.data and flat.data
    assert np.nanmax(np.asarray(flat.data[0].z, dtype=float)) <= 0.0
    assert np.nanmax(np.asarray(draped.data[0].z, dtype=float)) > 100.0

    opts = VolumeMapOptions(
        mode="surface",
        azimuth=90.0,
        show_stations=True,
        show_labels=False,
        period_range=(0.05, 2.0),
    )
    fig = build_3d_map(data, opts)
    iso_trace, marker_trace = fig.data[0], fig.data[-1]
    x_lo, x_hi = float(np.min(iso_trace.x)), float(np.max(iso_trace.x))
    y_lo, y_hi = float(np.min(iso_trace.y)), float(np.max(iso_trace.y))
    margin = 1e-6 + 0.05 * max(x_hi - x_lo, 1.0)
    mx = np.asarray(marker_trace.x, dtype=float)
    my = np.asarray(marker_trace.y, dtype=float)
    assert mx.min() >= x_lo - margin
    assert mx.max() <= x_hi + margin
    assert my.min() >= y_lo - margin
    assert my.max() <= y_hi + margin


class _ManyZ:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _ManyEdi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _ManyZ()


class _ManySites:
    def as_list(self):
        return [_ManyEdi(f"S{i:02d}") for i in range(12)]


def _many_station_data() -> MapData:
    return MapData(
        sites=_ManySites(),
        stations=tuple(
            StationRecord(f"S{i:02d}", 1.0 + i * 0.01, 2.0, 10.0, "L1", i)
            for i in range(12)
        ),
    )


def test_thin_indices_keeps_first_and_last_and_respects_cap() -> None:
    from pycsamt.map.volume import _thin_indices

    assert list(_thin_indices(12, None)) == list(range(12))
    assert list(_thin_indices(12, 20)) == list(range(12))
    assert list(_thin_indices(0, 5)) == []

    thinned = _thin_indices(12, 5)
    assert len(thinned) <= 5
    assert thinned[0] == 0
    assert thinned[-1] == 11

    assert list(_thin_indices(12, 1)) == [0]


def test_max_stations_thins_markers_per_line() -> None:
    """New option: cap how many station markers/labels are drawn per
    line, evenly spaced, so a crowded fence doesn't overload with
    every single station. Unset (default) must keep showing every
    station, unchanged from before this option existed.
    """
    data = _many_station_data()
    full = build_3d_map(data, VolumeMapOptions(mode="fence", show_stations=True))
    thinned = build_3d_map(
        data,
        VolumeMapOptions(mode="fence", show_stations=True, max_stations=5),
    )
    full_markers = [t for t in full.data if getattr(t, "name", "") == "stations"]
    thinned_markers = [
        t for t in thinned.data if getattr(t, "name", "") == "stations"
    ]
    assert len(full_markers[0].x) == 12
    assert len(thinned_markers[0].x) <= 5
    # First and last real station must survive thinning.
    assert thinned_markers[0].hovertext[0].endswith("S00")
    assert thinned_markers[0].hovertext[-1].endswith("S11")


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
        VolumeMapOptions(mode="fence", log_color=True, value_range=(10.0, 1000.0)),
    )
    lin = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", log_color=False, value_range=(10.0, 1000.0)),
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


def test_triangle_down_filled_uses_a_mesh3d_trace() -> None:
    # Plotly's Scatter3d marker enum has no triangle at all -- pycsamt's
    # own "triangle-down" extension must render as real geometry
    # instead of an (invalid) sprite symbol.
    fig = build_3d_map(
        _var_data(),
        VolumeMapOptions(mode="fence", show_stations=True, station_symbol="triangle-down"),
    )
    stations = [t for t in fig.data if getattr(t, "name", "") == "stations"]
    assert len(stations) == 1
    mesh = stations[0]
    assert type(mesh).__name__ == "Mesh3d"
    # 2 stations x 3 vertices/triangle = 6 vertices, 1 face each.
    assert len(mesh.x) == 6
    assert len(mesh.i) == 2


def test_triangle_down_open_uses_a_scatter3d_line_loop() -> None:
    fig = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence", show_stations=True, station_symbol="triangle-down-open"
        ),
    )
    stations = [t for t in fig.data if getattr(t, "name", "") == "stations"]
    assert len(stations) == 1
    outline = stations[0]
    assert type(outline).__name__ == "Scatter3d"
    assert outline.mode == "lines"
    # 2 stations x (4 loop points + 1 None separator) = 10 entries.
    assert len(outline.x) == 10
    # A None gap must separate the two triangles' loops.
    assert outline.x[4] is None

def test_triangle_down_with_labels_adds_scene_annotations() -> None:
    fig = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence",
            show_stations=True,
            station_symbol="triangle-down",
            station_labels=True,
        ),
    )
    # the marker geometry itself is one Mesh3d trace, no text trace
    stations = [t for t in fig.data if getattr(t, "name", "") == "stations"]
    assert [type(t).__name__ for t in stations] == ["Mesh3d"]
    # labels ride on scene annotations (so they can be rotated)
    notes = fig.layout.scene.annotations
    assert [n.text for n in notes] == ["S00", "S01"]
    assert all(n.textangle == 0 for n in notes)


def test_station_label_angle_rotates_scene_annotations() -> None:
    """New option: station labels can be tilted (like a 2-D section's
    tick labels) since they are drawn as scene annotations."""
    fig = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence",
            show_stations=True,
            station_labels=True,
            station_label_angle=90.0,
        ),
    )
    notes = fig.layout.scene.annotations
    assert len(notes) == 2
    assert all(n.textangle == 90 for n in notes)
    # each label is bottom-anchored and lifted clear of the marker so a
    # rotated label sits on top of the glyph, not through it
    assert all(n.yanchor == "bottom" and n.yshift > 0 for n in notes)
    # markers stay a plain marker trace (no text mode)
    marker = next(t for t in fig.data if getattr(t, "name", "") == "stations")
    assert marker.mode == "markers"


def test_station_label_fraction_thins_labels_not_markers() -> None:
    """Label density (``station_label_fraction``) drops labels only --
    every station marker still shows so no station is hidden."""
    data = _many_station_data()  # 12 stations, one line
    full = build_3d_map(
        data,
        VolumeMapOptions(mode="fence", show_stations=True, station_labels=True),
    )
    quarter = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence",
            show_stations=True,
            station_labels=True,
            station_label_fraction=0.25,
        ),
    )
    m_full = next(t for t in full.data if getattr(t, "name", "") == "stations")
    m_qtr = next(t for t in quarter.data if getattr(t, "name", "") == "stations")
    assert len(m_full.x) == len(m_qtr.x) == 12  # markers untouched
    assert len(full.layout.scene.annotations) == 12
    assert 2 <= len(quarter.layout.scene.annotations) <= 4


def test_station_label_names_labels_only_the_listed_stations() -> None:
    data = _many_station_data()
    fig = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence",
            show_stations=True,
            station_labels=True,
            station_label_names=("S02", "S07"),
        ),
    )
    marker = next(t for t in fig.data if getattr(t, "name", "") == "stations")
    assert len(marker.x) == 12  # all markers
    labelled = {n.text for n in fig.layout.scene.annotations}
    assert labelled == {"S02", "S07"}


def test_triangle_down_size_scales_with_station_size() -> None:
    small = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence", show_stations=True, station_symbol="triangle-down", station_size=2,
        ),
    )
    big = build_3d_map(
        _var_data(),
        VolumeMapOptions(
            mode="fence", show_stations=True, station_symbol="triangle-down", station_size=12,
        ),
    )
    small_mesh = [t for t in small.data if getattr(t, "name", "") == "stations"][0]
    big_mesh = [t for t in big.data if getattr(t, "name", "") == "stations"][0]
    small_width = max(small_mesh.x) - min(small_mesh.x)
    big_width = max(big_mesh.x) - min(big_mesh.x)
    assert big_width > small_width


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


def test_smoothed_fence_rho_range_opens_holes_not_recolour() -> None:
    """Regression: with section smoothing on, a resistivity-range
    selection must open real holes in the fence (NaN in both the
    colour and the mesh), not leave a full surface that merely
    recolours. ``_fill_nan_2d`` used to fill the masked cells so the
    spline could fit and never punched them back through.
    """
    base = VolumeMapOptions(mode="fence", smooth_sections=True, section_res=40)
    unfiltered = build_3d_map(_wide_data(), base)
    filtered = build_3d_map(
        _wide_data(),
        VolumeMapOptions(
            mode="fence",
            smooth_sections=True,
            section_res=40,
            rho_range=(10.0, 100.0),
        ),
    )
    u_color = np.asarray(unfiltered.data[0].surfacecolor, dtype=float)
    f_color = np.asarray(filtered.data[0].surfacecolor, dtype=float)
    f_z = np.asarray(filtered.data[0].z, dtype=float)
    # the unfiltered smoothed panel is gap-free; the filtered one masks
    # cells out of range in BOTH the colour and the geometry
    assert not np.isnan(u_color).any()
    assert np.isnan(f_color).any()
    assert np.isnan(f_z).any()
    # and some cells inside 10-100 ohm.m survive -- it is not all masked
    assert np.isfinite(f_color).any()
    # the colour scale is still the full unfiltered range
    assert np.isclose(unfiltered.data[0].cmin, filtered.data[0].cmin)
    assert np.isclose(unfiltered.data[0].cmax, filtered.data[0].cmax)
