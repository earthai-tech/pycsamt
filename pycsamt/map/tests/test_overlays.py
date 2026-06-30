# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for map overlay helpers."""

from __future__ import annotations

import numpy as np

from pycsamt.map import (
    CRSConfig,
    build_basemap_layout,
    build_profile_line_overlay,
    build_station_label_overlay,
    build_topography_overlay,
    interpolate_overlay_grid,
    normalize_epsg,
    transform_xy,
)


def test_normalize_epsg_and_identity_transform() -> None:
    assert normalize_epsg(4326) == "EPSG:4326"
    x, y = transform_xy(
        np.array([1.0, 2.0]),
        np.array([2.0, 3.0]),
        crs=CRSConfig(source=4326, target=4326),
    )
    assert x.tolist() == [1.0, 2.0]
    assert y.tolist() == [2.0, 3.0]


def test_interpolate_overlay_grid_shape() -> None:
    xi, yi, zz = interpolate_overlay_grid(
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([10.0, 20.0, 30.0]),
        grid_size=5,
    )
    assert xi.shape == (5,)
    assert yi.shape == (5,)
    assert zz.shape == (5, 5)


def test_topography_supports_mesh_and_surface() -> None:
    mesh = build_topography_overlay(
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([10.0, 20.0, 30.0]),
    )
    surface = build_topography_overlay(
        np.array([0.0, 1.0]),
        np.array([0.0, 1.0]),
        np.ones((2, 2)),
    )
    assert mesh.type == "mesh3d"
    assert surface.type == "surface"


def test_station_and_profile_overlays_build_traces() -> None:
    labels = build_station_label_overlay(
        np.array([0.0, 1.0]),
        np.array([0.0, 1.0]),
        np.array(["S0", "S1"]),
    )
    profile = build_profile_line_overlay(
        np.array([0.0, 1.0]),
        np.array([0.0, 1.0]),
    )
    assert labels.mode == "text"
    assert profile.mode == "lines"


def test_basemap_layout_handles_finite_points() -> None:
    layout = build_basemap_layout(
        np.array([2.0, 2.2]),
        np.array([1.0, 1.2]),
        dark=True,
        bearing=20.0,
    )
    assert layout.style == "carto-darkmatter"
    assert layout.center["lat"] == 1.1
    assert layout.bearing == 20.0
