"""Contracts for :mod:`pycsamt.ai.geology.topography`."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.ai.geology import (
    GeologyGrid,
    TopographicSurface,
    interpolate_topography,
    topography_from_sites,
)


def _fake_sites():
    return [
        SimpleNamespace(
            Head=SimpleNamespace(dataid="S01", lat=5.0, lon=-3.0, elev=100.0)
        ),
        SimpleNamespace(
            Head=SimpleNamespace(dataid="S02", lat=5.0, lon=-2.999, elev=120.0)
        ),
        SimpleNamespace(
            Head=SimpleNamespace(dataid="S03", lat=5.0, lon=-2.998, elev=110.0)
        ),
    ]


def test_surface_depth_masks_slope_summary_and_readonly():
    grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
    surface = TopographicSurface(grid, [120, 110, 100, 90], 120)
    np.testing.assert_array_equal(surface.surface_depth_m, [0, 10, 20, 30])
    assert surface.relief_m == 30
    assert surface.local_depth_m().shape == grid.shape
    np.testing.assert_array_equal(surface.air_mask(), ~surface.earth_mask())
    assert surface.summary()["air_cell_fraction"] >= 0
    assert not surface.elevation_m.flags.writeable
    assert np.all(surface.slope_degrees() > 0)
    assert len(surface.surface_hash) == 64


def test_2d_linear_nearest_and_cubic_interpolation():
    grid = GeologyGrid.regular_2d(nx=5, nz=3, dx_m=100, dz_m=50)
    x = np.array([0, 200, 500], dtype=float)
    elevation = np.array([100, 120, 110], dtype=float)
    linear = interpolate_topography(grid, x, elevation, interpolation_method="linear")
    nearest = interpolate_topography(grid, x, elevation, interpolation_method="nearest")
    assert linear.elevation_m.shape == nearest.elevation_m.shape == (5,)
    assert np.all(np.isfinite(linear.elevation_m))
    assert not np.array_equal(linear.elevation_m, nearest.elevation_m)

    cubic = interpolate_topography(
        grid,
        [0, 150, 350, 500],
        [100, 120, 90, 110],
        interpolation_method="cubic",
    )
    assert np.all(np.isfinite(cubic.elevation_m))


def test_3d_scattered_interpolation_fills_outside_convex_hull():
    grid = GeologyGrid.regular_3d(nx=5, ny=4, nz=3, dx_m=100, dy_m=100, dz_m=50)
    coordinates = np.array([[0, 0], [500, 0], [0, 400], [500, 400]], dtype=float)
    elevations = np.array([100, 110, 120, 90], dtype=float)
    surface = interpolate_topography(grid, coordinates, elevations)
    assert surface.elevation_m.shape == (4, 5)
    assert np.all(np.isfinite(surface.elevation_m))
    assert surface.reference_elevation_m == pytest.approx(np.max(surface.elevation_m))
    assert surface.earth_mask().shape == grid.shape


def test_sites_2d_extraction_alignment_and_provenance():
    grid = GeologyGrid.regular_2d(nx=5, nz=4, dx_m=50, dz_m=25)
    surface = topography_from_sites(
        _fake_sites(),
        grid,
        station_names=["S01", "S02", "S03"],
    )
    assert surface.source == "sites"
    assert surface.station_names == ("S01", "S02", "S03")
    np.testing.assert_array_equal(surface.sample_elevation_m, [100, 120, 110])
    assert surface.sample_coordinates_m.shape == (3, 1)
    assert surface.provenance()["vertical_datum"] == "metres above sea level"


def test_sites_projected_coordinates_support_3d_and_2d_chainage():
    sites = _fake_sites()
    coordinates = np.array([[0, 0], [100, 0], [200, 100]], dtype=float)
    grid2d = GeologyGrid.regular_2d(nx=6, nz=4, dx_m=50, dz_m=25)
    profile = topography_from_sites(sites, grid2d, coordinates_m=coordinates)
    expected_last = 100 + np.sqrt(2) * 100
    assert profile.sample_coordinates_m[-1, 0] == pytest.approx(expected_last)

    grid3d = GeologyGrid.regular_3d(nx=5, ny=5, nz=3, dx_m=50, dy_m=50, dz_m=25)
    surface = topography_from_sites(
        sites,
        grid3d,
        coordinates_m=coordinates,
        interpolation_method="nearest",
    )
    assert surface.elevation_m.shape == (5, 5)
    np.testing.assert_array_equal(surface.sample_coordinates_m, coordinates)


def test_topography_npz_roundtrip_preserves_samples_and_hash(tmp_path):
    grid = GeologyGrid.regular_2d(nx=4, nz=3, dx_m=100, dz_m=50)
    surface = interpolate_topography(
        grid,
        [0, 400],
        [100, 120],
        station_names=["A", "B"],
        source="survey",
    )
    restored = TopographicSurface.from_npz(surface.to_npz(tmp_path / "topo.npz"))
    assert restored.surface_hash == surface.surface_hash
    np.testing.assert_array_equal(
        restored.sample_coordinates_m, surface.sample_coordinates_m
    )
    assert restored.station_names == surface.station_names


def test_topography_validation_rejects_ambiguous_or_invalid_samples():
    grid2d = GeologyGrid.regular_2d(nx=4, nz=3, dx_m=100, dz_m=50)
    with pytest.raises(ValueError, match="unique"):
        interpolate_topography(grid2d, [0, 0], [100, 110])
    with pytest.raises(ValueError, match="four samples"):
        interpolate_topography(
            grid2d, [0, 1, 2], [1, 2, 1], interpolation_method="cubic"
        )
    with pytest.raises(ValueError, match="missing requested"):
        topography_from_sites(_fake_sites(), grid2d, station_names=["missing"])
    zeros = [
        SimpleNamespace(Head=SimpleNamespace(dataid="A", lat=5, lon=-3, elev=0)),
        SimpleNamespace(Head=SimpleNamespace(dataid="B", lat=5, lon=-2.9, elev=0)),
    ]
    with (
        pytest.warns(UserWarning),
        pytest.raises(ValueError, match="all zero"),
    ):
        topography_from_sites(zeros, grid2d)
    grid3d = GeologyGrid.regular_3d(nx=3, ny=3, nz=3, dx_m=100, dy_m=100, dz_m=50)
    with pytest.raises(ValueError, match="must be unique"):
        interpolate_topography(
            grid3d,
            [[0, 0], [0, 0], [1, 1]],
            [1, 2, 3],
            interpolation_method="nearest",
        )
    with pytest.raises(ValueError, match="requires projected"):
        topography_from_sites(_fake_sites(), grid3d)
