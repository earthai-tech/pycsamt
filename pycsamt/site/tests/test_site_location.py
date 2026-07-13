from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.site.location import (
    Coord,
    apply_topography,
    bearing,
    chainage_along,
    distance,
    ensure_head_coords,
    parse_elev,
    parse_lat,
    parse_lon,
    project,
)


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


def _has_pyproj_or_gdal() -> bool:
    try:
        import pyproj  # noqa: F401

        return True
    except Exception:
        try:
            from osgeo import osr  # noqa: F401

            return True
        except Exception:
            return False


def test_parse_lat_lon_elev_variants():
    assert parse_lat("45N") == pytest.approx(45.0)
    assert parse_lat("45.5") == pytest.approx(45.5)
    assert parse_lat("45 30 0 S") == pytest.approx(-45.5)
    assert parse_lon("123W") == pytest.approx(-123.0)
    assert parse_lon("123 15 30 E") == pytest.approx(
        123 + 15 / 60 + 30 / 3600
    )

    assert parse_elev("1200") == 1200.0
    assert math.isnan(parse_elev("oops"))
    assert math.isnan(parse_lat("??"))
    assert math.isnan(parse_lon(None))


def test_ensure_head_coords_on_stub():
    class Stub:
        def __init__(self):
            self._sections = {}

        def set_section(self, name, sec):
            self._sections[name] = sec

        def get_section(self, name):
            return self._sections.get(name)

    ed = Stub()
    h = ensure_head_coords(
        ed,
        lat="10N",
        lon="20E",
        elev="300",
    )
    assert h.lat == pytest.approx(10.0)
    assert h.lon == pytest.approx(20.0)
    assert h.elev == pytest.approx(300.0)
    assert h.long == pytest.approx(20.0)


def test_ensure_head_coords_on_real_edi(simulated_edi: Path):
    ed = _load_edi(simulated_edi)
    h0 = ed.get_section("head")  # type: ignore
    ensure_head_coords(ed, lat=26.25, lon=11.5, elev=777)
    h1 = ed.get_section("head")  # type: ignore
    assert h1.lat == pytest.approx(26.25)
    assert h1.lon == pytest.approx(11.5)
    assert h1.elev == pytest.approx(777.0)
    # ensure we didn't drop other head fields
    assert getattr(h1, "dataid", None) == getattr(h0, "dataid", None)


@pytest.mark.parametrize(
    "use_list,inplace",
    [
        (False, True),
        (True, True),
        (True, False),
    ],
)
def test_apply_topography_updates_coords(
    tmp_path: Path,
    simulated_edi: Path,
    use_list: bool,
    inplace: bool,
):
    pd = pytest.importorskip("pandas")

    p = _dup_edi(tmp_path, simulated_edi, "T01")
    ed = _load_edi(p)

    sid = ed.get_section("head").dataid  # type: ignore

    df = pd.DataFrame(
        {
            "station": [sid, "OTHER"],
            "latitude": [35.125, 0.0],
            "longitude": [12.75, 0.0],
            "elevation": [1234.0, 0.0],
        }
    )

    obj = [ed] if use_list else ed
    out = apply_topography(obj, df, inplace=inplace)

    def _coords(e):
        h = e.get_section("head")  # type: ignore
        return float(h.lat), float(h.lon), float(h.elev)

    if use_list and not inplace:
        # new list returned; original unchanged
        la0, lo0, ev0 = _coords(ed)
        la1, lo1, ev1 = _coords(out[0])
        assert (la0, lo0, ev0) != (35.125, 12.75, 1234.0)
        assert (la1, lo1, ev1) == (35.125, 12.75, 1234.0)
    else:
        target = out[0] if use_list else out
        la, lo, ev = _coords(target)
        assert (la, lo, ev) == (35.125, 12.75, 1234.0)


def test_distance_bearing_chainage_simple():
    a = Coord(0.0, 0.0, 0.0)
    b_east = Coord(0.0, 1.0, 0.0)
    b_north = Coord(1.0, 0.0, 0.0)

    d_lat1 = distance(a, b_north)
    d_lon1 = distance(a, b_east)

    # ~111 km per degree (tolerate small ellipsoid diffs)
    assert d_lat1 == pytest.approx(111_000, rel=0.02)
    assert d_lon1 == pytest.approx(111_000, rel=0.02)

    # bearings: north=0, east=90
    assert bearing(a, b_north) == pytest.approx(0.0, abs=1e-6)
    assert bearing(a, b_east) == pytest.approx(90.0, abs=1e-6)

    # chainage along east-directed line from origin
    ch_e = chainage_along((0.0, 0.0), 90.0, (0.0, 1.0))
    ch_n = chainage_along((0.0, 0.0), 90.0, (1.0, 0.0))
    assert ch_e == pytest.approx(111_000, rel=0.02)
    assert ch_n == pytest.approx(0.0, abs=1e-6)


@pytest.mark.skipif(
    not _has_pyproj_or_gdal(),
    reason="requires pyproj or GDAL",
)
def test_project_identity_epsg4326():
    # lon/lat order (always_xy=True in implementation)
    x, y = project((12.5, 35.25), crs_from="EPSG:4326", crs_to="EPSG:4326")
    assert isinstance(x, np.ndarray) and isinstance(y, np.ndarray)
    assert x.size == 1 and y.size == 1
    assert x[0] == pytest.approx(12.5, abs=1e-9)
    assert y[0] == pytest.approx(35.25, abs=1e-9)
