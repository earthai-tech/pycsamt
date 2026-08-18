from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import NavigationTrack


def test_navigation_track_basic_and_bbox():
    nav = NavigationTrack(
        sample_ids=("P1", "P2", "P3"),
        latitude=[5.0, 5.1, 5.2],
        longitude=[-3.0, -3.1, -3.2],
        terrain_elevation=[100.0, 101.0, 102.0],
        platform_elevation=[180.0, 181.0, 184.0],
    )
    assert nav.n_samples == 3
    assert nav.has_geographic_coordinates
    np.testing.assert_allclose(nav.clearance_values, [80.0, 80.0, 82.0])
    assert nav.bbox is not None
    assert nav.bbox.lat_min == pytest.approx(5.0)
    assert nav.bbox.lon_min == pytest.approx(-3.2)


def test_navigation_explicit_clearance_takes_precedence():
    nav = NavigationTrack(
        sample_ids=("P1", "P2"),
        terrain_elevation=[100.0, 100.0],
        platform_elevation=[200.0, 200.0],
        clearance=[60.0, 61.0],
    )
    np.testing.assert_allclose(nav.clearance_values, [60.0, 61.0])


def test_navigation_allows_missing_sample_coordinates():
    nav = NavigationTrack(
        sample_ids=("P1", "P2"),
        latitude=[5.0, np.nan],
        longitude=[-3.0, np.nan],
    )
    assert nav.bbox is not None
    assert nav.bbox.centre == pytest.approx((5.0, -3.0))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"sample_ids": ()},
        {"sample_ids": ("P1", "P1")},
        {"sample_ids": ("P1", "P2"), "latitude": [5.0]},
        {"sample_ids": ("P1",), "latitude": [5.0]},
        {"sample_ids": ("P1",), "easting": [100.0]},
        {
            "sample_ids": ("P1",),
            "latitude": [95.0],
            "longitude": [0.0],
        },
    ],
)
def test_navigation_rejects_invalid_alignment(kwargs):
    with pytest.raises(ValueError):
        NavigationTrack(**kwargs)


def test_navigation_index_lookup():
    nav = NavigationTrack(sample_ids=("P1", "P2"))
    assert nav.index_of("P2") == 1
    with pytest.raises(KeyError):
        nav.index_of("P3")
