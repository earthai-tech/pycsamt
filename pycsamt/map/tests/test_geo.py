# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# ruff: noqa: E501
"""Tests for ESRI basemaps and the geo contour-image overlay."""

from __future__ import annotations

import numpy as np

from pycsamt.map import (
    basemap_style_and_layers,
    build_geo_contour_image,
)


def test_esri_style_resolves_to_raster_layer() -> None:
    style, layers = basemap_style_and_layers("esri-satellite")
    assert style == "white-bg"
    assert len(layers) == 1
    assert layers[0]["sourcetype"] == "raster"
    assert "arcgisonline" in layers[0]["source"][0]


def test_native_style_has_no_layers() -> None:
    style, layers = basemap_style_and_layers("open-street-map")
    assert style == "open-street-map"
    assert layers == []


def test_default_style_depends_on_theme() -> None:
    assert basemap_style_and_layers(None, dark=True)[0] == "carto-darkmatter"
    assert basemap_style_and_layers(None, dark=False)[0] == "open-street-map"


def test_geo_contour_image_builds_png_and_corners() -> None:
    rng = np.random.RandomState(0)
    lon = 1.0 + rng.rand(12) * 0.1
    lat = 6.0 + rng.rand(12) * 0.1
    val = 50.0 + rng.rand(12) * 500.0
    out = build_geo_contour_image(lon, lat, val, n_levels=8)
    assert out is not None
    assert out["image"].startswith("data:image/png;base64,")
    assert len(out["coordinates"]) == 4
    # corners are [lon, lat] and span the data extent
    lons = [c[0] for c in out["coordinates"]]
    lats = [c[1] for c in out["coordinates"]]
    assert min(lons) < float(lon.min()) + 1e-6
    assert max(lats) > float(lat.max()) - 1e-6


def test_geo_contour_image_needs_three_points() -> None:
    assert build_geo_contour_image([1.0], [6.0], [10.0]) is None


def test_geo_contour_image_interp_and_smoothing() -> None:
    rng = np.random.RandomState(1)
    lon = 1.0 + rng.rand(15) * 0.1
    lat = 6.0 + rng.rand(15) * 0.1
    val = 50.0 + rng.rand(15) * 500.0
    for interp in ("cubic", "linear", "nearest"):
        out = build_geo_contour_image(
            lon,
            lat,
            val,
            interp=interp,
            smooth_sigma=1.5,
            grid_res=120,
            n_levels=10,
        )
        assert out is not None
        assert out["image"].startswith("data:image/png;base64,")


def test_geo_contour_image_log_scale_skips_nonpositive() -> None:
    lon = np.array([1.0, 1.1, 1.2, 1.3])
    lat = np.array([6.0, 6.1, 6.2, 6.3])
    val = np.array([-1.0, 10.0, 100.0, 1000.0])
    out = build_geo_contour_image(lon, lat, val, log_scale=True)
    # 3 positive points remain -> still builds
    assert out is not None
