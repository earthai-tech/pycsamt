# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Branch-level coverage for ``pycsamt.site.location``.

Complements ``test_site_location.py`` (real EDI + happy paths) with
targeted tests for the defensive branches a well-formed EDI never
triggers: locked/partial head objects, the ``Sites``-container path
through ``apply_topography``, the ``flat``/``utm`` distance/bearing
models, and the private matching helpers.
"""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from pycsamt.site.location import (
    Coord,
    apply_topography,
    bearing,
    chainage_along,
    distance,
    ensure_head_coords,
    parse_lat,
    parse_lon,
    project,
)
from pycsamt.site.location import (
    _get_station,
    _match_row,
    _parse_angle,
    _pick_col,
    _set_coords_from_row,
)

pd = pytest.importorskip("pandas")


# ---------------------------------------------------------------------------
# ensure_head_coords -- attribute-setting failures
# ---------------------------------------------------------------------------


class _SelectivelyLockedHead:
    """A head object whose named attributes reject assignment."""

    def __init__(self, locked=()):
        self._locked = set(locked)

    def __setattr__(self, name, value):
        if name != "_locked" and name in self._locked:
            raise AttributeError(f"{name} is locked")
        object.__setattr__(self, name, value)


class _StubEdiSetSection:
    def __init__(self, head):
        self._head = head
        self.set_section_calls = 0
        self.raise_on_set_section = False

    def get_section(self, name):
        return self._head if name == "head" else None

    def set_section(self, name, value):
        self.set_section_calls += 1
        if self.raise_on_set_section:
            raise RuntimeError("set_section failed")
        self._head = value


def test_ensure_head_coords_swallows_lat_lon_elev_setter_errors():
    h = _SelectivelyLockedHead(locked={"lat", "lon", "long", "elev"})
    ed = _StubEdiSetSection(h)
    # lon/long both locked -> wrote_lon stays False -> the fallback
    # Head-rebuild block runs and succeeds via set_section.
    out = ensure_head_coords(ed, lat=10.0, lon=20.0, elev=5.0)
    assert out.lon == pytest.approx(20.0)
    assert out.long == pytest.approx(20.0)
    assert ed.set_section_calls == 1


def test_ensure_head_coords_falls_back_to_edi_head_attr_when_set_section_fails():
    h = _SelectivelyLockedHead(locked={"lon", "long"})
    ed = _StubEdiSetSection(h)
    ed.raise_on_set_section = True
    # set_section raises -> except branch assigns ed.Head = nh directly.
    out = ensure_head_coords(ed, lat=1.0, lon=2.0, elev=3.0)
    assert out.lon == pytest.approx(2.0)
    assert ed.Head.lon == pytest.approx(2.0)


def test_ensure_head_coords_swallows_total_assignment_failure():
    class _TotallyLockedEdi(_StubEdiSetSection):
        def __setattr__(self, name, value):
            if name == "Head":
                raise AttributeError("no Head slot")
            object.__setattr__(self, name, value)

    h = _SelectivelyLockedHead(locked={"lon", "long"})
    ed = _TotallyLockedEdi(h)
    ed.raise_on_set_section = True
    # Both the set_section() and ed.Head = nh fallbacks fail silently;
    # ensure_head_coords must not raise.
    out = ensure_head_coords(ed, lat=1.0, lon=2.0, elev=3.0)
    assert out is h


# ---------------------------------------------------------------------------
# apply_topography -- Sites-like container (``._items``) path
# ---------------------------------------------------------------------------


class _FakeEdi:
    def __init__(self, station):
        self.station = station
        self.lat = 0.0
        self.lon = 0.0
        self.elev = 0.0

    def get_section(self, name):
        return self if name == "head" else None

    def set_section(self, name, value):
        pass


class _FakeSite:
    def __init__(self, edi):
        self.edi = edi


class _FakeSitesContainer:
    def __init__(self, items):
        self._items = items


def _topo_frame():
    return pd.DataFrame(
        {
            "station": ["S1"],
            "latitude": [12.5],
            "longitude": [45.5],
            "elevation": [321.0],
        }
    )


def test_apply_topography_sites_container_inplace():
    site = _FakeSite(_FakeEdi("S1"))
    container = _FakeSitesContainer([site])
    out = apply_topography(container, _topo_frame(), inplace=True)
    assert out is container
    assert site.edi.lat == pytest.approx(12.5)
    assert site.edi.lon == pytest.approx(45.5)
    assert site.edi.elev == pytest.approx(321.0)


def test_apply_topography_sites_container_not_inplace_deep_copies():
    site = _FakeSite(_FakeEdi("S1"))
    container = _FakeSitesContainer([site])
    out = apply_topography(container, _topo_frame(), inplace=False)
    assert out is not container
    assert out._items[0].edi.lat == pytest.approx(12.5)
    # original left untouched
    assert container._items[0].edi.lat == pytest.approx(0.0)


def test_apply_topography_sites_container_exception_is_swallowed():
    class _BadContainer:
        @property
        def _items(self):
            raise RuntimeError("no items for you")

    out = apply_topography(_BadContainer(), _topo_frame(), inplace=True)
    assert isinstance(out, _BadContainer)


def test_apply_topography_list_deepcopy_failure_is_swallowed():
    class _Undeepcopyable:
        def __deepcopy__(self, memo):
            raise RuntimeError("no copy")

        def get_section(self, name):
            return None

    item = _Undeepcopyable()
    out = apply_topography([item], _topo_frame(), inplace=False)
    # falls through to the single-EDI branch; must not raise
    assert out is not None


def test_apply_topography_no_matching_row_leaves_object_untouched():
    edi = _FakeEdi("UNKNOWN")
    out = apply_topography(edi, _topo_frame(), inplace=True)
    assert out.lat == 0.0
    assert out.lon == 0.0


# ---------------------------------------------------------------------------
# project() -- pyproj-missing fallback and multi-point input
# ---------------------------------------------------------------------------


def test_project_raises_without_pyproj_or_gdal(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def _blocked_import(name, *args, **kwargs):
        if name == "pyproj":
            raise ImportError("blocked for test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _blocked_import)
    monkeypatch.setattr("pycsamt.site.location.HAS_GDAL", False)
    with pytest.raises(RuntimeError, match="pyproj or GDAL"):
        project((0.0, 0.0), crs_from="EPSG:4326", crs_to="EPSG:3857")


def test_project_multi_point_sequence():
    pts = [(0.0, 0.0), (1.0, 1.0)]
    X, Y = project(pts, crs_from="EPSG:4326", crs_to="EPSG:4326")
    assert X.shape == (2,)
    assert Y.shape == (2,)
    np.testing.assert_allclose(X, [0.0, 1.0])
    np.testing.assert_allclose(Y, [0.0, 1.0])


# ---------------------------------------------------------------------------
# distance() / bearing() -- flat and utm modes
# ---------------------------------------------------------------------------


def test_distance_flat_mode_matches_geodetic_closely():
    a = Coord(35.0, 10.0, 0.0)
    b = Coord(35.01, 10.01, 0.0)
    d_flat = distance(a, b, mode="flat")
    d_geo = distance(a, b, mode="geodetic")
    assert d_flat == pytest.approx(d_geo, rel=0.05)


def test_distance_utm_mode_matches_geodetic_closely():
    a = Coord(35.0, 10.0, 0.0)
    b = Coord(35.01, 10.01, 0.0)
    d_utm = distance(a, b, mode="utm")
    d_geo = distance(a, b, mode="geodetic")
    assert d_utm == pytest.approx(d_geo, rel=0.02)


def test_distance_invalid_mode_raises():
    with pytest.raises(ValueError):
        distance((0.0, 0.0), (1.0, 1.0), mode="bogus")


def test_bearing_flat_mode():
    a = Coord(0.0, 0.0)
    east = Coord(0.0, 1.0)
    north = Coord(1.0, 0.0)
    assert bearing(a, east, mode="flat") == pytest.approx(90.0, abs=1e-6)
    assert bearing(a, north, mode="flat") == pytest.approx(0.0, abs=1e-6)


def test_bearing_utm_mode_roughly_matches_geodetic():
    a = Coord(35.0, 10.0)
    b = Coord(35.01, 10.01)
    b_utm = bearing(a, b, mode="utm")
    b_geo = bearing(a, b, mode="geodetic")
    assert b_utm == pytest.approx(b_geo, abs=5.0)


def test_bearing_invalid_mode_raises():
    with pytest.raises(ValueError):
        bearing((0.0, 0.0), (1.0, 1.0), mode="bogus")


def test_bearing_utm_explicit_crs():
    a = Coord(35.0, 10.0)
    b = Coord(35.01, 10.01)
    # UTM zone 32N (10E is on the boundary of 31N/32N -- use an explicit
    # EPSG rather than relying on auto-inference).
    b_utm = bearing(a, b, mode="utm", crs_to="EPSG:32632")
    assert 0.0 <= b_utm < 360.0


# ---------------------------------------------------------------------------
# chainage_along() -- Coord inputs and multi-point sequences
# ---------------------------------------------------------------------------


def test_chainage_along_accepts_coord_origin_and_points():
    origin = Coord(0.0, 0.0)
    pts = [Coord(0.0, 1.0), Coord(1.0, 0.0), (0.0, 0.5)]
    out = chainage_along(origin, 90.0, pts)
    assert isinstance(out, np.ndarray)
    assert out.shape == (3,)
    assert out[0] == pytest.approx(111_000.0, rel=0.02)
    assert out[1] == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# Private matching helpers
# ---------------------------------------------------------------------------


def test_match_row_returns_none_when_df_has_no_len():
    assert _match_row(42, "S1") is None


def test_match_row_returns_none_without_id_column():
    df = pd.DataFrame({"latitude": [1.0], "longitude": [2.0]})
    assert _match_row(df, "S1") is None


def test_match_row_returns_none_on_internal_error():
    class _WeirdColumn:
        def astype(self, *_a, **_k):
            raise RuntimeError("no astype here")

    df = pd.DataFrame({"station": [_WeirdColumn()]})
    assert _match_row(df, "S1") is None


def test_pick_col_matches_case_insensitively():
    df = pd.DataFrame({"Station": [1], "Lat": [2]})
    assert _pick_col(df, ("station",)) == "Station"
    assert _pick_col(df, ("missing",)) is None


def test_set_coords_from_row_returns_none_for_none_row():
    # No exception, no return value; purely defensive.
    assert _set_coords_from_row(object(), None, empty=0.0) is None


def test_get_station_falls_through_head_and_object_attrs():
    class _EmptyHead:
        pass

    class _Ed:
        def get_section(self, name):
            return _EmptyHead() if name == "head" else None

    assert _get_station(_Ed()) == ""

    class _EdWithName:
        name = "FromObjectAttr"

        def get_section(self, name):
            return None

    assert _get_station(_EdWithName()) == "FromObjectAttr"


# ---------------------------------------------------------------------------
# _parse_angle -- sign-without-hemisphere and non-DMS fallback branches
# ---------------------------------------------------------------------------


def test_parse_angle_negative_without_hemisphere_letter():
    assert parse_lat("-10 30 0") == pytest.approx(-10.5)


def test_parse_angle_fallback_non_dms_with_hemisphere_suffix():
    # "1.5e1S" does not match the DMS regex (exponent notation), so it
    # falls through to the plain-numeric-plus-trailing-letter branch.
    assert _parse_angle("1.5e1S", {"N": 1.0, "S": -1.0}) == pytest.approx(-15.0)


def test_parse_angle_fallback_non_dms_without_hemisphere():
    assert _parse_angle("1.5e1", {"N": 1.0, "S": -1.0}) == pytest.approx(15.0)


def test_parse_lon_rejects_garbage_via_fallback_branch():
    assert math.isnan(parse_lon("1.5eX"))
