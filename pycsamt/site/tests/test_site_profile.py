import math

import pytest

from pycsamt.site.location import Coord, chainage_along
from pycsamt.site.profile import (
    Profile,
    infer_line_orientation,
)

# ----- tiny EDI-like mocks ---------------------------------------


class _Head:
    def __init__(self, lat=None, lon=None, long=None, dataid=None):
        # support both lon/long spellings
        self.lat = lat
        if long is not None:
            self.long = long
        if lon is not None:
            self.lon = lon
            # keep "long" in sync if present on the real class
            try:
                self.long = lon
            except Exception:
                pass
        self.dataid = dataid


class _EDI:
    def __init__(self, name, lat, lon):
        self.station = name
        self._head = _Head(lat=lat, lon=lon, dataid=name)

    def get_section(self, key):
        return self._head if key.lower() == "head" else None


# helper to generate a line of mock sites
def _make_line(prefix, n, lat0, lon0, d_lon_deg=0.0, d_lat_deg=0.0):
    edis = []
    for i in range(n):
        lat = lat0 + i * d_lat_deg
        lon = lon0 + i * d_lon_deg
        edis.append(_EDI(f"{prefix}{i + 1:02d}", lat, lon))
    return edis


# ----- tests ------------------------------------------------------
def test_infer_line_orientation_cardinals():
    def assert_axis(a, target, tol=5.0):
        # compare line orientation up to 180° ambiguity
        d = (a - target) % 180.0
        d = min(d, 180.0 - d)
        assert d <= tol

    # East-West line (constant lat, increasing lon) -> ~90°
    east_line = _make_line("E", 5, lat0=0.0, lon0=0.0, d_lon_deg=0.01)
    az_e = infer_line_orientation(east_line)
    assert az_e == pytest.approx(90.0, abs=3.0)

    # North-South line (constant lon, increasing lat) -> ~0°
    north_line = _make_line("N", 5, lat0=0.0, lon0=0.0, d_lat_deg=0.01)
    az_n = infer_line_orientation(north_line)
    # allow small numerical drift around 0/360
    assert min(az_n, 360.0 - az_n) == pytest.approx(0.0, abs=3.0)

    # ~45° line (NE)
    diag_line = _make_line("D", 5, 0.0, 0.0, d_lon_deg=0.01, d_lat_deg=0.01)
    az_d = infer_line_orientation(diag_line)
    assert_axis(az_d, 45.0, tol=5.0)


def test_profile_from_sites_chainage_and_stats():
    # 1 km spacing eastward near equator:
    # 1 deg lon ≈ 111 km => 0.009009... deg ≈ 1 km
    dlon_1km = 1.0 / 111.0
    sites = _make_line("S", 4, lat0=0.0, lon0=0.0, d_lon_deg=dlon_1km)
    prof = Profile.from_sites(sites)

    # chainages should be ≈ 0, 1k, 2k, 3k (meters)
    svals = [prof.chainages[s.station] for s in sites]
    for i, s in enumerate(svals):
        assert s == pytest.approx(1000.0 * i, rel=0.03)

    # spacing stats around 1k
    mean = prof.spacing_stats.get("spacing_mean")
    med = prof.spacing_stats.get("spacing_med")
    assert mean == pytest.approx(1000.0, rel=0.05)
    assert med == pytest.approx(1000.0, rel=0.05)

    # sort_sites returns original order here
    sorted_sites = prof.sort_sites(sites)
    assert [s.station for s in sorted_sites] == [s.station for s in sites]

    # slice in meters
    window = prof.slice(500.0, 2500.0)
    assert list(window.keys()) == ["S02", "S03"]

    # resample every 750 m from min..max (0..3000) → 0,750,1500,2250,3000
    grid = prof.resample(750.0)
    assert grid.size >= 5
    assert grid[0] == pytest.approx(0.0, abs=1e-6)
    assert grid[-1] == pytest.approx(3000.0, rel=1e-6)

    # summary looks coherent
    summary = prof.summary()
    assert summary["n_sites"] == 4.0
    assert summary["s_min"] == pytest.approx(0.0, abs=1.0)
    assert summary["s_max"] == pytest.approx(3000.0, abs=5.0)


def test_profile_handles_missing_coords_gracefully():
    good = _EDI("A01", 0.0, 0.0)
    bad1 = _EDI("A02", float("nan"), 0.0)
    bad2 = _EDI("A03", 0.0, float("nan"))
    prof = Profile.from_sites([good, bad1, bad2])

    # Only the good one should get a chainage
    assert set(prof.chainages.keys()) == {"A01"}
    assert math.isfinite(list(prof.chainages.values())[0])

    # Sorting should keep only finite
    sorted_sites = prof.sort_sites([good, bad1, bad2])
    assert [s.station for s in sorted_sites] == ["A01"]


def test_chainage_along_consistency_with_profile():
    # Along 90° (east), chainage should match lon delta * scale
    origin = Coord(0.0, 0.0, 0.0)
    p = (0.0, 1.0 / 111.0)  # ~1 km east
    s = chainage_along((origin.lat, origin.lon), 90.0, p)
    assert s == pytest.approx(1000.0, rel=0.03)
