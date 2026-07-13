# License: LGPL-3.0

from __future__ import annotations

import math

import pytest

import pycsamt.loc as locmod


@pytest.fixture(autouse=True)
def stub_ll_utm(monkeypatch):
    """
    Stub deterministic converters so tests don't depend on
    real geodesy. Inverse is exact.
    """

    def ll_to_utm_stub(
        *,
        reference_ellipsoid: int,
        lat: float,
        lon: float,
    ) -> tuple[str, float, float]:
        n = int(math.floor((lon + 180.0) / 6.0) + 1)
        n = min(60, max(1, n))
        hem = "N" if float(lat) >= 0.0 else "S"
        z = f"{n}{hem}"
        e = (float(lon) + 180.0) * 1000.0
        nrt = (float(lat) + 90.0) * 1000.0
        return z, e, nrt

    def utm_to_ll_stub(
        *,
        reference_ellipsoid: int,
        northing: float,
        easting: float,
        zone: str,
    ) -> tuple[float, float]:
        lat = float(northing) / 1000.0 - 90.0
        lon = float(easting) / 1000.0 - 180.0
        return lat, lon

    monkeypatch.setattr(locmod, "ll_to_utm", ll_to_utm_stub)
    monkeypatch.setattr(locmod, "utm_to_ll", utm_to_ll_stub)


# --------------------------- Location ---------------------------------


def test_location_init_latlon_sets_utm():
    loc = locmod.Location(latitude=10, longitude=5, elevation=12.5)
    assert loc.latitude == 10.0
    assert loc.longitude == 5.0
    z, e, n = loc.to_utm()
    assert isinstance(z, str) and z.endswith("N")
    assert e == pytest.approx((5.0 + 180.0) * 1000.0)
    assert n == pytest.approx((10.0 + 90.0) * 1000.0)


def test_location_init_utm_requires_zone():
    with pytest.raises(locmod.LocationError):
        locmod.Location(easting=1000, northing=2000)


def test_location_init_with_utm_and_zone_to_latlon():
    loc = locmod.Location(
        utm_zone="31N",
        easting=181000.0,
        northing=100000.0,
    )
    la, lo = loc.to_latlon()
    assert la == pytest.approx(10.0)
    assert lo == pytest.approx(1.0)


def test_property_validators_and_zone_normalization():
    loc = locmod.Location(latitude=0, longitude=0, utm_zone="31n")
    assert loc.utm_zone == "31N"
    with pytest.raises(locmod.LocationError):
        loc.latitude = 100.0
    with pytest.raises(locmod.LocationError):
        loc.longitude = -200.0


def test_to_utm_in_vectorized_basic():
    e, n = locmod.Location.to_utm_in(
        lats=[0.0, 10.0],
        lons=[0.0, 5.0],
    )
    assert e.shape == (2,)
    assert n.shape == (2,)
    assert e[0] == pytest.approx(180000.0)
    assert n[0] == pytest.approx(90000.0)
    assert e[1] == pytest.approx(185000.0)
    assert n[1] == pytest.approx(100000.0)


def test_to_latlon_in_vectorized_requires_zone():
    with pytest.raises(locmod.LocationError):
        locmod.Location.to_latlon_in(
            easts=[181000.0],
            norths=[100000.0],
        )


def test_to_latlon_in_vectorized_ok():
    la, lo = locmod.Location.to_latlon_in(
        easts=[181000.0, 186000.0],
        norths=[100000.0, 95000.0],
        utm_zone="31N",
    )
    assert la.tolist() == pytest.approx([10.0, 5.0])
    assert lo.tolist() == pytest.approx([1.0, 6.0])


def test_baseloc_batch_converters():
    e, n, z = locmod.BaseLoc.to_utm_batch(
        lats=[0.0, -5.0],
        lons=[0.0, 12.0],
    )
    assert list(z) == ["31N", "33S"]
    la, lo = locmod.BaseLoc.to_latlon_batch(
        easts=e,
        norths=n,
        zones=z,
    )
    assert la.tolist() == pytest.approx([0.0, -5.0])
    assert lo.tolist() == pytest.approx([0.0, 12.0])


def test_distance_and_nearly_equals():
    a = locmod.Location(latitude=10.0, longitude=10.0)
    b = locmod.Location(latitude=10.0001, longitude=10.0001)
    d = a.distance_to(b)
    assert d > 0.0
    assert a.nearly_equals(b, tol_m=20.0) is True
    c = locmod.Location(latitude=11.0, longitude=10.0)
    assert a.nearly_equals(c, tol_m=50.0) is False


def test_dict_roundtrip_and_equality_hash():
    a = locmod.Location(latitude=3.5, longitude=7.25, utm_zone="31N")
    # force both systems to exist
    a.to_utm()
    d = a.to_dict()
    b = locmod.Location.from_dict(d)
    assert a == b
    s = {a, b}
    assert len(s) == 1


# ----------------------------- UTMZone --------------------------------


def test_utmzone_basic_and_from_latlon():
    z1 = locmod.UTMZone("31n")
    assert z1.as_str() == "31N"
    z2 = locmod.UTMZone.from_latlon(10.0, 5.0)
    assert z2.as_str() == "31N"


def test_norm_zone_helper_variants():
    assert locmod._norm_zone("32s") == "32S"
    assert locmod._norm_zone("60N") == "60N"
    assert locmod._norm_zone(31) == "31N"
    assert locmod._norm_zone(-31) == "31S"
    with pytest.raises(locmod.LocationError):
        locmod._norm_zone("0N")
    with pytest.raises(locmod.LocationError):
        locmod._norm_zone("61N")


# ------------------------------ Bounds --------------------------------


def test_bounds_from_points_and_ops():
    la = [0.0, 1.0, -2.0]
    lo = [10.0, 12.0, 9.5]
    b = locmod.Bounds.from_points(la, lo)
    assert b.to_tuple() == (-2.0, 9.5, 1.0, 12.0)
    assert b.contains(0.0, 10.0) is True
    assert b.contains(5.0, 10.0) is False
    c = b.buffer_m(1000.0)
    # buffer expands both min and max
    assert c.min_lat < b.min_lat and c.max_lat > b.max_lat
    u = b.union(locmod.Bounds(-1.0, 8.0, 2.0, 9.0))
    assert u.to_tuple() == (-2.0, 8.0, 2.0, 12.0)
    i = b.intersection(locmod.Bounds(-1.0, 8.0, 2.0, 10.0))
    assert i is not None
    assert i.to_tuple() == (-1.0, 9.5, 1.0, 10.0)


def test_bounds_to_utm_rect_uses_zone_center_and_warns(caplog):
    b = locmod.Bounds(0.0, 0.0, 1.0, 10.0)
    z, e0, n0, e1, n1 = b.to_utm_rect()
    assert z in {"31N", "32N"}
    assert e1 > e0 and n1 > n0
    # explicit zone normalization path
    z2, *_ = b.to_utm_rect(zone="31n")
    assert z2 == "31N"


# ------------------------------ GeoPath -------------------------------


def test_geopath_append_extend_bbox_length():
    g = locmod.GeoPath()
    g.append(0.0, 0.0)
    g.extend([(0.0, 1.0), (1.0, 1.0)])
    assert len(g) == 3
    bb = g.bbox()
    assert bb.to_tuple() == (0.0, 0.0, 1.0, 1.0)
    # length > 0 and roughly consistent
    L = g.length_m()
    assert L > 0.0


def test_geopath_to_utm_single_zone_check():
    # lon=5 -> zone 31; lon=10 -> zone 32 (with stub)
    g = locmod.GeoPath([(0.0, 5.0), (0.0, 10.0)])
    zones, e, n = g.to_utm()
    assert zones == ["31N", "32N"]
    with pytest.raises(locmod.LocationError):
        g.to_utm(ensure_single_zone=True)


def test_geopath_simplify_reduces_points():
    pts = [
        (0.0, 0.0),
        (0.0, 0.001),
        (0.0, 0.002),
        (0.0, 0.010),
    ]
    g = locmod.GeoPath(pts)
    g2 = g.simplify(eps_m=100.0)
    # expect endpoints only after aggressive simplify
    assert len(g2) == 2
    # keep endpoints
    bb = g2.bbox()
    assert bb.to_tuple() == (0.0, 0.0, 0.0, 0.01)
