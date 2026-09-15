from __future__ import annotations

import builtins
from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.models.mare2dem.geom.area_of_interest import (
    estimate_area_of_interest,
    survey_points,
)
from pycsamt.models.mare2dem.geom.triangle_regions import get_triangle_regions
from pycsamt.models.mare2dem.geom.utm import (
    _ll_to_utm_pure,
    _utm_to_ll_pure,
    lonlat_to_utm,
    utm_to_lonlat,
)


def _group(**values):
    defaults = {
        "receivers": np.empty((0, 3)),
        "transmitters": np.empty((0, 3)),
        "tx_electrodes": np.empty((0, 3)),
        "rx_electrodes": np.empty((0, 3)),
    }
    defaults.update(values)
    return SimpleNamespace(**defaults)


def _em(*, csem=None, mt=None, dc=None, data=None):
    return SimpleNamespace(csem=csem, mt=mt, dc=dc, data=data)


def test_survey_points_collects_all_supported_geometries():
    em = _em(
        csem=_group(
            receivers=np.array([[0, 1, 2], [0, 3, 4]]),
            transmitters=np.array([[0, 5, 6]]),
        ),
        mt=_group(receivers=np.array([[0, 7, 8]])),
        dc=_group(
            tx_electrodes=np.array([[0, 9, 10]]),
            rx_electrodes=np.array([[0, 11, 12]]),
        ),
    )
    np.testing.assert_array_equal(
        survey_points(em),
        [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12]],
    )
    assert survey_points(_em()).shape == (0, 2)


def test_area_of_interest_empty_single_and_square_variants():
    assert estimate_area_of_interest(_em()) == (None, None)

    single = _em(mt=_group(receivers=np.array([[0, 10, 20]])))
    ylim, zlim = estimate_area_of_interest(single)
    np.testing.assert_allclose(ylim, [-990, 1010])
    np.testing.assert_allclose(zlim, [-980, 1020])

    wider = _em(mt=_group(receivers=np.array([[0, 0, 0], [0, 10, 8]])))
    ylim, zlim = estimate_area_of_interest(wider)
    np.testing.assert_allclose(ylim, [-1, 11])
    np.testing.assert_allclose(zlim, [-2, 10])

    taller = _em(mt=_group(receivers=np.array([[0, 0, 0], [0, 8, 10]])))
    ylim, zlim = estimate_area_of_interest(taller)
    np.testing.assert_allclose(zlim, [-1, 11])
    np.testing.assert_allclose(ylim, [-2, 10])


def test_area_of_interest_elongated_mt_and_csem_range_limit():
    mt = _em(mt=_group(receivers=np.array([[0, 0, 2], [0, 100, 4]])))
    ylim, zlim = estimate_area_of_interest(mt)
    np.testing.assert_allclose(ylim, [-10, 110])
    np.testing.assert_allclose(zlim, [-8, 79])

    csem = _em(
        csem=_group(
            transmitters=np.array([[0, 0, 0]]),
            receivers=np.array([[0, 50, 1], [0, 100, 2]]),
        ),
        data=np.array([[1, 0, 1, 2], [101, 0, 1, 1]]),
    )
    ylim, zlim = estimate_area_of_interest(csem)
    np.testing.assert_allclose(ylim, [-10, 110])
    np.testing.assert_allclose(zlim, [-10, 102])

    invalid_links = _em(
        csem=csem.csem,
        data=np.array([[1, 0, 99, 99]]),
    )
    _, fallback_zlim = estimate_area_of_interest(invalid_links)
    np.testing.assert_allclose(fallback_zlim, [-10, 77])


def _square_mesh():
    points = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    triangles = np.array([[1, 2, 3], [1, 3, 4]])
    return points, triangles


def test_triangle_regions_connected_and_boundary_separated():
    points, triangles = _square_mesh()
    labels, mapping = get_triangle_regions(
        points, triangles, np.array([[1, 2]])
    )
    np.testing.assert_array_equal(labels, [1, 1])
    np.testing.assert_array_equal(mapping, [0])

    labels, mapping = get_triangle_regions(
        points, triangles, np.array([[1, 3]])
    )
    assert set(labels) == {1, 2}
    np.testing.assert_array_equal(mapping, [0, 0])


def test_triangle_regions_seed_mapping_duplicate_and_outside_seeds():
    points, triangles = _square_mesh()
    seeds = np.array([[0.75, 0.25], [0.25, 0.75], [2.0, 2.0]])
    labels, mapping = get_triangle_regions(
        points, triangles, np.array([[1, 3]]), seeds
    )
    assert set(labels) == {1, 2}
    # scipy's Delaunay simplex ordering can differ from the caller's triangle
    # ordering; one seeded component is therefore retained and the other is
    # discovered by the unassigned-triangle pass.
    assert sorted(mapping.tolist()) == [0, 1]

    # With no wall, the second seed lands in an already-filled component.
    labels, mapping = get_triangle_regions(
        points, triangles, np.array([[1, 2]]), seeds[:2]
    )
    np.testing.assert_array_equal(labels, [1, 1])
    np.testing.assert_array_equal(mapping, [1])


def test_triangle_regions_accepts_zero_based_connectivity():
    points, triangles = _square_mesh()
    labels, _ = get_triangle_regions(
        points, triangles - 1, np.array([[0, 2]])
    )
    assert set(labels) == {1, 2}


@pytest.mark.parametrize(
    ("lon", "lat", "zone", "south"),
    [(-70.0, 42.0, 19, False), (18.0, -34.0, 34, True)],
)
def test_pure_utm_roundtrip_north_and_south(lon, lat, zone, south):
    east, north = _ll_to_utm_pure(
        np.array([lat]), np.array([lon]), zone, south, "wgs84"
    )
    lat2, lon2 = _utm_to_ll_pure(east, north, zone, south, "wgs84")
    np.testing.assert_allclose(lon2, [lon], atol=2e-3)
    np.testing.assert_allclose(lat2, [lat], atol=2e-3)


def test_public_utm_normalizes_longitude_and_hemisphere_strings():
    east, north, zone, south = lonlat_to_utm(290.0, 42.0)
    assert zone == 19 and south is False
    lon, lat = utm_to_lonlat(east, north, zone, "North")
    np.testing.assert_allclose(lon, [-70], atol=2e-3)
    np.testing.assert_allclose(lat, [42], atol=2e-3)

    east, north, zone, south = lonlat_to_utm(18.0, -34.0)
    assert south is True
    lon, lat = utm_to_lonlat(east, north, zone, "South")
    np.testing.assert_allclose([lon[0], lat[0]], [18, -34], atol=2e-6)


def test_public_utm_falls_back_when_pyproj_is_unavailable(monkeypatch):
    original_import = builtins.__import__

    def without_pyproj(name, *args, **kwargs):
        if name == "pyproj":
            raise ImportError("disabled")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_pyproj)
    east, north, zone, south = lonlat_to_utm(-70.0, 42.0)
    lon, lat = utm_to_lonlat(east, north, zone, south)
    np.testing.assert_allclose(lon, [-70], atol=2e-3)
    np.testing.assert_allclose(lat, [42], atol=2e-3)
