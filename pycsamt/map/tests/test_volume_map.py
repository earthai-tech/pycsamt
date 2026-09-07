# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for volume maps."""

from __future__ import annotations

import numpy as np
import pytest

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


def _layered_section_data() -> MapData:
    """3 lines, a real (z, rho) inversion section each: conductive
    overburden over a resistive layer over a conductive base, plus a
    resistive body down one flank -- so a horizontal slice is never
    laterally uniform."""
    z = np.linspace(20.0, 900.0, 24)
    sections, recs = {}, []
    for li in range(3):
        nsta = 10 + li * 2
        sta = np.array([f"L{li}S{i:02d}" for i in range(nsta)], dtype=object)
        zz, xx = np.meshgrid(z, np.arange(nsta), indexing="ij")
        rho = np.where(zz < 200, 30.0, np.where(zz < 550, 500.0, 12.0))
        rho = np.where(xx > nsta * 0.7, rho * 5.0, rho)
        elev = 100.0 + 30.0 * np.sin(np.arange(nsta) * 0.5) + li * 10
        sections[f"L{li}"] = {
            "stations": sta,
            "rho": rho,
            "z": z,
            "elev": elev,
        }
        for i, s in enumerate(sta):
            recs.append(
                StationRecord(str(s), float("nan"), float("nan"),
                              float(elev[i]), f"L{li}", len(recs))
            )
    return MapData(
        sites=None, stations=tuple(recs), metadata={"sections": sections}
    )


def test_depth_slices_stay_inside_the_data_and_show_lateral_structure() -> None:
    """Regression: slices were generated across the raw depth *window*
    (e.g. 0-800 m), so the endpoints landed outside the model's real
    z-extent and rendered as all-NaN ghost surfaces. Every drawn slice
    must now carry real, laterally-varying resistivity.
    """
    from pycsamt.map.volume import _profile_grids, _slice_depths

    data = _layered_section_data()
    opts = VolumeMapOptions(mode="depth", n_slices=5, depth_range=(0.0, 5000.0))
    grids = _profile_grids(data, opts)
    z_all = np.concatenate([g["z"] for g in grids.values()])
    depths = _slice_depths(grids, opts)
    assert depths.min() >= z_all.min() - 1e-6
    assert depths.max() <= z_all.max() + 1e-6

    fig = build_3d_map(data, opts)
    surfaces = [t for t in fig.data if type(t).__name__ == "Surface"]
    assert surfaces
    for s in surfaces:
        sc = np.asarray(s.surfacecolor, dtype=float)
        assert np.isfinite(sc).any()  # no all-NaN ghost slice
        assert np.nanstd(sc) > 1e-3  # laterally varying, not uniform


def test_depth_slice_colorbar_is_fixed_across_depth_windows() -> None:
    """The colour scale must not move when the depth window / slice
    depth changes -- same principle as ignoring ``rho_range``."""
    data = _layered_section_data()
    full = build_3d_map(data, VolumeMapOptions(mode="depth", n_slices=4))
    shallow = build_3d_map(
        data, VolumeMapOptions(mode="depth", n_slices=1, depth_range=(0.0, 150.0))
    )
    f0, s0 = full.data[0], shallow.data[0]
    assert np.isclose(f0.cmin, s0.cmin) and np.isclose(f0.cmax, s0.cmax)


def test_depth_single_slice_cuts_through_the_window_middle() -> None:
    from pycsamt.map.volume import _profile_grids, _slice_depths

    data = _layered_section_data()
    opts = VolumeMapOptions(mode="depth", n_slices=1, depth_range=(0.0, 200.0))
    depths = _slice_depths(_profile_grids(data, opts), opts)
    assert len(depths) == 1
    assert 80.0 <= depths[0] <= 120.0  # ~100 m, not 0 m


def test_depth_slice_rho_range_filters_after_interpolating_the_depth() -> None:
    """Regression: the resistivity filter used to NaN whole cells
    *before* the depth interpolation, so ``np.interp`` bridged across an
    out-of-band layer and painted the slice a solid colour where it
    should be masked. It must now filter against the resistivity
    actually interpolated at the slice depth.
    """
    from pycsamt.map.volume import _profile_grids, _values_at_depth

    data = _layered_section_data()  # 30 / 500 / 12 ohm.m layers, split at 200/550 m
    cond = VolumeMapOptions(mode="depth", rho_range=(1.0, 100.0))
    grids = _profile_grids(data, cond)
    grid = next(iter(grids.values()))

    # a slice right in the 500 ohm.m resistive layer -> every left-flank
    # column masked (the right flank is 5x higher, also masked)
    mid = _values_at_depth(grid, 350.0, cond)
    assert np.isnan(mid).all()

    # a slice in the conductive overburden -> left flank (30 ohm.m) shows,
    # right flank (150 ohm.m) is masked
    top = _values_at_depth(grid, 90.0, cond)
    assert np.isfinite(top).any() and np.isnan(top).any()


def test_depth_slice_rho_range_opens_holes_not_recolour() -> None:
    """Regression: a resistivity-range selection in depth-slice mode
    must open real holes in the mesh (NaN in both the colour and the
    ``z`` geometry), not leave a full, solid slice tinted by a NaN
    colour. ``go.Surface`` does not treat a NaN *surfacecolor* cell as
    transparent on its own -- only NaN in ``z`` opens a real gap (see
    ``_fence_figure``'s own ``hole``/``zz`` step, applied here to
    ``_depth_figure`` too). Before the fix, every cell kept rendering
    at its own true colour regardless of the selected range -- the
    filter had no visible effect at all in this one mode.
    """
    data = _layered_section_data()  # 30 / 500 / 12 ohm.m layers, split at 200/550 m
    # A shallow slice (~90 m, inside the 30/150 ohm.m conductive
    # overburden layer, per test_depth_slice_rho_range_filters_after_
    # interpolating_the_depth above): the left flank (30 ohm.m) is
    # inside 1-100 ohm.m, the right flank (150 ohm.m, x > 0.7*nsta) is
    # not -- a genuine partial mask, not all-or-nothing.
    window = dict(mode="depth", n_slices=1, depth_range=(20.0, 150.0))
    unfiltered = build_3d_map(data, VolumeMapOptions(**window))
    filtered = build_3d_map(
        data,
        VolumeMapOptions(**window, rho_range=(1.0, 100.0)),
    )
    u_color = np.asarray(unfiltered.data[0].surfacecolor, dtype=float)
    f_color = np.asarray(filtered.data[0].surfacecolor, dtype=float)
    f_z = np.asarray(filtered.data[0].z, dtype=float)
    # Filtering must open *more* holes than the unfiltered panel's own
    # (line-padding) gaps -- the resistivity range, not just the
    # ragged station counts, is what is being tested here.
    assert np.isnan(f_color).sum() > np.isnan(u_color).sum()
    assert np.isnan(f_z).any()
    assert np.array_equal(np.isnan(f_color), np.isnan(f_z))
    assert np.isfinite(f_color).any()
    # the colour scale is still the full unfiltered range
    assert np.isclose(unfiltered.data[0].cmin, filtered.data[0].cmin)
    assert np.isclose(unfiltered.data[0].cmax, filtered.data[0].cmax)


def test_volume_smoothing_refines_and_softens_without_moving_the_footprint() -> None:
    """Opt-in ``volume_smoothing`` reconstructs the block / iso-surface
    volume on a *finer* lattice (so Plotly rounds the surfaces instead
    of faceting -- the Geosoft-voxel look) and then gently blurs it.
    ``0.0`` (default) leaves the raw lattice untouched. Neither the data
    footprint (bounding box) nor the colour scale may move.
    """
    from pycsamt.map.volume import _dense_volume_grid, _profile_grids

    data = _layered_section_data()
    raw_opts = VolumeMapOptions(mode="block", volume_smoothing=0.0)
    sm_opts = VolumeMapOptions(mode="block", volume_smoothing=1.5)
    xr, _, zr, raw = _dense_volume_grid(_profile_grids(data, raw_opts), raw_opts)
    xs, _, zs, sm = _dense_volume_grid(_profile_grids(data, sm_opts), sm_opts)

    # finer lattice, same spatial extent (identical x/z span, just
    # sampled denser -- the blur must never bleed the volume outward)
    assert sm.size > raw.size
    assert np.isclose(xs.min(), xr.min()) and np.isclose(xs.max(), xr.max())
    assert np.nanmin(zs) >= np.nanmin(zr) - 1e-6
    assert np.nanmax(zs) <= np.nanmax(zr) + 1e-6

    # same along-profile footprint: finite data starts/ends at the same x
    def _x_span(vol, xa):
        fin = np.isfinite(vol).any(axis=(1, 2))
        xi = np.where(fin)[0]
        return xa[xi[0]], xa[xi[-1]]

    np.testing.assert_allclose(
        _x_span(sm, xs), _x_span(raw, xr), rtol=0.02, atol=1e-6
    )

    # genuinely smoother: smaller normalised cell-to-cell gradients
    g_raw = np.nanstd(np.diff(raw, axis=2)) / (np.nanstd(raw) or 1.0)
    g_sm = np.nanstd(np.diff(sm, axis=2)) / (np.nanstd(sm) or 1.0)
    assert g_sm < g_raw

    raw_fig = build_3d_map(data, raw_opts)
    sm_fig = build_3d_map(data, sm_opts)
    assert np.isclose(raw_fig.data[0].cmin, sm_fig.data[0].cmin)
    assert np.isclose(raw_fig.data[0].cmax, sm_fig.data[0].cmax)


def test_volume_smoothing_levels_scale_lattice_density() -> None:
    """Each stronger preset (Light -> Medium -> Strong -> Very strong)
    reconstructs the volume on a denser lattice, so raising the control
    keeps removing visible facets."""
    from pycsamt.map.volume import (
        _dense_volume_grid,
        _profile_grids,
        _volume_smoothing_params,
    )

    data = _layered_section_data()
    sizes = []
    for strength in (0.8, 1.5, 2.5, 4.0):
        opts = VolumeMapOptions(mode="block", volume_smoothing=strength)
        _, _, _, vol = _dense_volume_grid(_profile_grids(data, opts), opts)
        sizes.append(vol.size)
    assert sizes == sorted(sizes)
    assert sizes[-1] > 2 * sizes[0]

    # 4.0 resolves to the "very strong" preset (the densest factor)
    assert _volume_smoothing_params(4.0)[0] >= _volume_smoothing_params(2.5)[0]
    assert _volume_smoothing_params(0.0) == (1.0, 0, 0, 0.0)


def test_volume_smoothing_feathers_the_banded_isosurface_cliff() -> None:
    """A resistivity-band body used to render as a hard sentinel/in-band
    step (voxel staircase). With smoothing on, the composed ``value``
    field carries a smooth ramp between the two, so ``go.Isosurface``
    cuts a rounded envelope."""
    from pycsamt.map.volume import _volume_point_cloud, _profile_grids

    data = _layered_section_data()
    band = dict(mode="block", rho_range=(20.0, 120.0))
    raw = _volume_point_cloud(
        _profile_grids(data, VolumeMapOptions(**band, volume_smoothing=0.0)),
        VolumeMapOptions(**band, volume_smoothing=0.0),
    )
    sm = _volume_point_cloud(
        _profile_grids(data, VolumeMapOptions(**band, volume_smoothing=2.5)),
        VolumeMapOptions(**band, volume_smoothing=2.5),
    )
    assert raw is not None and sm is not None
    raw_vals = np.unique(np.round(raw[3], 6))
    sm_vals = np.unique(np.round(sm[3], 6))
    # raw: essentially two populations (sentinel + the in-band hues);
    # feathered: many intermediate values bridging them
    assert sm_vals.size > raw_vals.size * 3
    # colour scale still the full-model range, not the band
    assert np.isclose(raw[6], sm[6]) and np.isclose(raw[7], sm[7])


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


def _air_fill_section_data(*, air_value: float = float("nan")) -> MapData:
    """One line, a real (z, rho) section whose shallowest rows are a
    ModEM-style air / overburden fill over a normal 12-500 ohm.m earth,
    plus a handful of genuine 5e4 ohm.m resistive outliers -- the
    Baohuashan fig08-vs-MapView discrepancy, reduced to a single curtain.

    ``air_value`` is ``nan`` by default (the state after the adapter /
    loader defensive mask), or a finite ~1e12 to exercise the in-view
    ``rho_display_max`` cutoff.
    """
    z = np.linspace(20.0, 900.0, 24)
    nsta = 12
    sta = np.array([f"S{i:02d}" for i in range(nsta)], dtype=object)
    zz, xx = np.meshgrid(z, np.arange(nsta), indexing="ij")
    rho = np.where(zz < 200, 30.0, np.where(zz < 550, 500.0, 12.0))
    rho[(zz > 800) & (xx >= nsta - 1)] = 5e4  # ~2 real resistive outliers
    rho[:3, :] = air_value  # three air-fill rows on top
    elev = np.full(nsta, 100.0)
    recs = [
        StationRecord(str(s), float("nan"), float("nan"), 100.0, "L0", i)
        for i, s in enumerate(sta)
    ]
    return MapData(
        sites=None,
        stations=tuple(recs),
        metadata={"sections": {"L0": {"stations": sta, "rho": rho, "z": z, "elev": elev}}},
    )


def test_crange_percentile_robust_default_clips_outliers() -> None:
    """With air already masked to NaN (the real pipeline), the 2nd-98th
    percentile default keeps the auto colour scale on the bulk earth
    instead of a handful of 5e4 ohm.m outliers."""
    data = _air_fill_section_data()
    fig = build_3d_map(data, VolumeMapOptions(mode="fence", log_color=True))
    # bulk earth tops out at 500 -> log10 < 3; 5e4 outliers -> log10 ~ 4.7
    assert float(fig.data[0].cmax) < 4.0


def test_crange_percentile_none_is_raw_minmax() -> None:
    data = _air_fill_section_data()
    fig = build_3d_map(
        data,
        VolumeMapOptions(mode="fence", log_color=True, crange_percentile=None),
    )
    # raw max is 5e4 -> log10 ~ 4.7, not clipped
    assert float(fig.data[0].cmax) > 4.5


def test_rho_display_max_masks_cells_and_rescales_colour() -> None:
    data = _air_fill_section_data(air_value=1e12)
    base = VolumeMapOptions(mode="fence", log_color=True, crange_percentile=None)
    no_cut = build_3d_map(data, base)
    cut = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence",
            log_color=True,
            crange_percentile=None,
            rho_display_max=1e6,
        ),
    )
    # cells above 1e6 ohm.m are punched out of the panel geometry ...
    assert np.isnan(np.asarray(cut.data[0].z, dtype=float)).any()
    # ... and no longer inflate the colour scale.
    assert float(cut.data[0].cmax) < float(no_cut.data[0].cmax)
    # no_cut colourbar reaches the 1e12 air (log10 ~ 12); the cut one
    # only the 5e4 real earth outliers (log10 ~ 4.7).
    assert float(no_cut.data[0].cmax) > 10.0
    assert float(cut.data[0].cmax) < 5.0


def test_rho_display_max_clamps_block_iso_band() -> None:
    data = _air_fill_section_data(air_value=1e12)
    fig = build_3d_map(
        data,
        VolumeMapOptions(
            mode="block", log_color=True, rho_display_max=1e6
        ),
    )
    iso = [t for t in fig.data if type(t).__name__ in ("Isosurface", "Volume")]
    assert iso
    # the iso band's upper bound is clamped to log10(1e6) = 6, well
    # below the 1e12 air fill.
    assert float(iso[0].isomax) <= 6.0 + 1e-9


def _dup_x_section_data() -> MapData:
    """Two lines whose station projection collapses a pair onto one
    along-strike position each -- a real ModEM multi-line geometry that
    used to make RegularGridInterpolator raise and leave the whole
    block/iso volume empty."""
    z = np.linspace(20.0, 900.0, 20)
    sections, recs = {}, []
    for li in range(2):
        nsta = 10
        xs = np.arange(nsta, dtype=float) * 100.0
        xs[3] = xs[2]  # two stations at the same along-strike x
        sta = np.array([f"L{li}S{i:02d}" for i in range(nsta)], dtype=object)
        zz, xx = np.meshgrid(z, xs, indexing="ij")
        rho = np.where(zz < 400, 40.0, 800.0) * (1.0 + 0.1 * (xx / 900.0))
        sections[f"L{li}"] = {
            "stations": sta, "rho": rho, "z": z,
            "elev": np.full(nsta, 100.0 + li * 20.0),
        }
        for i, s in enumerate(sta):
            recs.append(
                StationRecord(str(s), float("nan"), float("nan"),
                              100.0 + li * 20.0, f"L{li}", len(recs))
            )
    return MapData(
        sites=None, stations=tuple(recs), metadata={"sections": sections}
    )


def test_block_and_iso_render_with_duplicate_station_x() -> None:
    data = _dup_x_section_data()
    for mode, kind in (("block", "Volume"), ("surface", "Isosurface")):
        fig = build_3d_map(data, VolumeMapOptions(mode=mode))
        bodies = [t for t in fig.data if type(t).__name__ == kind]
        assert bodies, f"{mode}: no {kind} trace was drawn"
        val = np.asarray(bodies[0].value, dtype=float)
        assert np.isfinite(val).any()
        # no "empty block" annotation
        assert not getattr(fig.layout, "annotations", None)


# ---------------------------------------------------------------------------
# Interpretation overlay: ``VolumeMapOptions.geology`` discrete colour bands
# ---------------------------------------------------------------------------

_GEOLOGY_BANDS = (
    (10.0, 200.0, "#E9C46A"),
    (1000.0, 20000.0, "#8D99AE"),
)


def test_geology_unset_keeps_the_continuous_named_colorscale() -> None:
    """Default behaviour (no Interpretation overlay applied) must be
    bit-for-bit unchanged: a continuous, named Plotly colorscale."""
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    fig = build_3d_map(data, VolumeMapOptions(mode="fence"))
    scale = fig.data[0].colorscale
    # A continuous scale has many more than 6 stops and no exact
    # duplicate-value pair (the hard-step signature of a banded scale).
    values = [s[0] for s in scale]
    assert len(scale) > 6
    assert len(set(values)) == len(values)


def test_geology_bands_render_as_hard_stepped_colorscale_on_every_mode() -> None:
    """``options.geology`` must take over the colour scale (and cmin/cmax)
    on all four 3-D modes -- block/fence/depth/iso -- while every other
    render behaviour (trace kind, isomin/isomax, masking) stays as-is."""
    from pycsamt.map.styles import geology_crange

    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    lo, hi = geology_crange(_GEOLOGY_BANDS)
    for mode in ("fence", "block", "depth", "surface"):
        fig = build_3d_map(
            data, VolumeMapOptions(mode=mode, geology=_GEOLOGY_BANDS)
        )
        colored = [
            t for t in fig.data if getattr(t, "colorscale", None) is not None
        ]
        assert colored, f"{mode}: no coloured trace found"
        trace = colored[0]
        scale = trace.colorscale
        values = [s[0] for s in scale]
        # Hard-step signature: a repeated stop value at every band edge.
        assert len(values) != len(set(values))
        assert scale[0][0] == 0.0 and scale[-1][0] == 1.0
        assert trace.cmin == pytest.approx(lo)
        assert trace.cmax == pytest.approx(hi)


def test_geology_bands_do_not_change_which_cells_are_masked() -> None:
    """The Interpretation overlay only recolours -- it must not interact
    with the independent ``rho_range`` visibility-band masking."""
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    plain = build_3d_map(
        data, VolumeMapOptions(mode="depth", rho_range=(1.0, 500.0))
    )
    with_geology = build_3d_map(
        data,
        VolumeMapOptions(
            mode="depth", rho_range=(1.0, 500.0), geology=_GEOLOGY_BANDS
        ),
    )
    plain_nan = np.isnan(np.asarray(plain.data[-1].surfacecolor, dtype=float))
    geo_nan = np.isnan(
        np.asarray(with_geology.data[-1].surfacecolor, dtype=float)
    )
    assert np.array_equal(plain_nan, geo_nan)


_NAMED_GEOLOGY_BANDS = (
    (10.0, 200.0, "#E9C46A", "Sand"),
    (1000.0, 20000.0, "#8D99AE", "Granodiorite"),
)


def test_geology_named_bands_replace_the_colorbar_ticks_with_rock_names() -> None:
    """Applying a legend whose bands carry a name masks the numeric
    Ω·m colourbar with the geology legend itself, on every 3-D mode."""
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    for mode in ("fence", "block", "depth", "surface"):
        fig = build_3d_map(
            data, VolumeMapOptions(mode=mode, geology=_NAMED_GEOLOGY_BANDS)
        )
        colored = [
            t for t in fig.data if getattr(t, "colorscale", None) is not None
        ]
        assert colored, f"{mode}: no coloured trace found"
        # The native colourbar is suppressed in favour of the
        # non-overlapping on-canvas legend chips.
        assert colored[0].showscale is False
        chip_texts = {a.text.strip() for a in fig.layout.annotations}
        assert {"Sand", "Granodiorite"} <= chip_texts
        assert fig.layout.margin.r > 0


def test_geology_legend_toggle_frees_the_plot_width() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    shown = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence", geology=_NAMED_GEOLOGY_BANDS, geology_legend=True,
        ),
    )
    hidden = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence", geology=_NAMED_GEOLOGY_BANDS, geology_legend=False,
        ),
    )
    assert shown.layout.margin.r > 0
    assert hidden.layout.margin.r == 0
    assert not any(
        "Sand" in (a.text or "") for a in hidden.layout.annotations
    )
    # Legend hidden still means the classification/colouring itself
    # stays applied -- only the on-canvas key disappears.
    assert hidden.data[0].showscale is False


def test_geology_unnamed_bands_still_hide_the_native_colorbar() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    fig = build_3d_map(
        data, VolumeMapOptions(mode="fence", geology=_GEOLOGY_BANDS)
    )
    assert fig.data[0].showscale is False
    # Unnamed bands fall back to a "rho_min-rho_max" chip label instead
    # of a rock name, but the legend itself still renders.
    assert any("–" in (a.text or "") for a in fig.layout.annotations)


# ---------------------------------------------------------------------------
# geology_legend_style: "swatch" (default) vs. the previous "colorbar"
# ---------------------------------------------------------------------------


def test_geology_legend_style_colorbar_shows_the_native_colorbar_with_named_ticks() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    for mode in ("fence", "block", "depth", "surface"):
        fig = build_3d_map(
            data,
            VolumeMapOptions(
                mode=mode, geology=_NAMED_GEOLOGY_BANDS,
                geology_legend_style="colorbar",
            ),
        )
        colored = [
            t for t in fig.data if getattr(t, "colorscale", None) is not None
        ]
        assert colored, f"{mode}: no coloured trace found"
        assert colored[0].showscale is True
        assert list(colored[0].colorbar.ticktext) == ["Sand", "Granodiorite"]
        # the "colorbar" style relies on Plotly's own colourbar margin,
        # not the custom on-canvas legend's reserved margin.
        assert fig.layout.margin.r == 0
        assert not any(
            "Sand" in (a.text or "") for a in fig.layout.annotations
        )


def test_geology_legend_style_swatch_is_the_default() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    default_style = build_3d_map(
        data, VolumeMapOptions(mode="fence", geology=_NAMED_GEOLOGY_BANDS)
    )
    explicit_swatch = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence", geology=_NAMED_GEOLOGY_BANDS,
            geology_legend_style="swatch",
        ),
    )
    assert default_style.data[0].showscale is False
    assert explicit_swatch.data[0].showscale is False
    assert default_style.layout.margin.r == explicit_swatch.layout.margin.r > 0


def test_geology_legend_toggle_also_hides_the_colorbar_style() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    fig = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence", geology=_NAMED_GEOLOGY_BANDS,
            geology_legend_style="colorbar", geology_legend=False,
        ),
    )
    assert fig.data[0].showscale is False


# ---------------------------------------------------------------------------
# geology_fill="pattern": true pattern-texture rendering (fence + depth)
# ---------------------------------------------------------------------------


def _checkerboard_stencil(size: int = 8) -> np.ndarray:
    stencil = np.zeros((size, size))
    stencil[:, size // 2 :] = 1.0
    return stencil


def test_pattern_fill_fence_varies_within_the_textured_band_only() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    patterns = {"Granodiorite": _checkerboard_stencil()}
    fig_pattern = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence", geology=_NAMED_GEOLOGY_BANDS,
            geology_fill="pattern", geology_patterns=patterns,
            pattern_tile_size_m=5.0,
        ),
    )
    fig_solid = build_3d_map(
        data, VolumeMapOptions(mode="fence", geology=_NAMED_GEOLOGY_BANDS)
    )
    sc_pattern = np.asarray(fig_pattern.data[0].surfacecolor, dtype=float)
    sc_solid = np.asarray(fig_solid.data[0].surfacecolor, dtype=float)
    # The untextured "Sand" band's cells (value == log10(100)) are
    # identical either way -- only the textured band's cells changed.
    sand_mask = np.isclose(sc_solid, np.log10(100.0))
    assert np.allclose(sc_pattern[sand_mask], sc_solid[sand_mask])
    assert not np.allclose(
        np.nan_to_num(sc_pattern), np.nan_to_num(sc_solid)
    )
    # more colorscale stops than the flat/solid case -- an actual
    # gradient, not one colour.
    assert len(fig_pattern.data[0].colorscale) > len(fig_solid.data[0].colorscale)


def test_pattern_fill_depth_slice_varies_within_the_textured_band_only() -> None:
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    patterns = {"Granodiorite": _checkerboard_stencil()}
    fig = build_3d_map(
        data,
        VolumeMapOptions(
            mode="depth", n_slices=2, geology=_NAMED_GEOLOGY_BANDS,
            geology_fill="pattern", geology_patterns=patterns,
            pattern_tile_size_m=5.0,
        ),
    )
    surfaces = [t for t in fig.data if getattr(t, "surfacecolor", None) is not None]
    assert surfaces
    values = np.concatenate(
        [np.asarray(t.surfacecolor, dtype=float).ravel() for t in surfaces]
    )
    values = values[np.isfinite(values)]
    # both the checkerboard's two density levels show up within the
    # textured band's own resistivity range.
    assert np.any(np.isclose(values, np.log10(1000.0)))
    assert np.any(np.isclose(values, np.log10(20000.0)))


def test_pattern_fill_without_a_pattern_stencil_stays_solid() -> None:
    """geology_fill="pattern" with no geology_patterns for a band must
    not change that band's rendering -- degrade to solid, not crash."""
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    fig = build_3d_map(
        data,
        VolumeMapOptions(
            mode="fence", geology=_NAMED_GEOLOGY_BANDS,
            geology_fill="pattern", geology_patterns=None,
        ),
    )
    solid = build_3d_map(
        data, VolumeMapOptions(mode="fence", geology=_NAMED_GEOLOGY_BANDS)
    )
    assert np.allclose(
        np.nan_to_num(np.asarray(fig.data[0].surfacecolor, dtype=float)),
        np.nan_to_num(np.asarray(solid.data[0].surfacecolor, dtype=float)),
    )


def test_pattern_fill_does_not_affect_block_or_surface_modes() -> None:
    """Block / iso-surface have no per-cell colour-axis override in
    Plotly -- geology_fill="pattern" must not attempt to touch them,
    and must not crash either."""
    data = MapData(
        sites=_VarSites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    patterns = {"Granodiorite": _checkerboard_stencil()}
    for mode in ("block", "surface"):
        fig = build_3d_map(
            data,
            VolumeMapOptions(
                mode=mode, geology=_NAMED_GEOLOGY_BANDS,
                geology_fill="pattern", geology_patterns=patterns,
            ),
        )
        assert fig.data
