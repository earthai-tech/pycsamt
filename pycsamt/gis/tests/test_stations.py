# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.gis.stations`."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pycsamt.gis.stations import (
    compare_crs,
    embed_station_coordinates,
    read_station_table,
    station_key,
    to_wgs84,
)

_ROWS = [("L1-0", 2830561.529, 792692.448, 373.721),
         ("L1-30", 2830543.800, 792716.693, 366.805),
         ("L1-60", 2830526.524, 792741.351, 343.662)]


def _write(tmp_path, text, name="t.csv"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_headerless_bom_and_empty_columns(tmp_path):
    text = "﻿" + "".join(
        f"{s},{n},{e},{z},,\n" for s, n, e, z in _ROWS)
    t = read_station_table(_write(tmp_path, text))
    assert t.columns.tolist() == ["station", "northing", "easting",
                                  "elev"]
    assert t.station.tolist() == ["L1-0", "L1-30", "L1-60"]
    np.testing.assert_allclose(t.northing, [r[1] for r in _ROWS])


def test_layout_inferred_when_columns_swapped(tmp_path):
    # easting before northing: inference uses magnitudes, not order
    text = "".join(f"{s} {e} {n} {z}\n" for s, n, e, z in _ROWS)
    t = read_station_table(_write(tmp_path, text, "t.txt"))
    np.testing.assert_allclose(t.northing, [r[1] for r in _ROWS])
    np.testing.assert_allclose(t.easting, [r[2] for r in _ROWS])


def test_header_detected_and_explicit_columns(tmp_path):
    text = "name,x,y,z\n" + "".join(
        f"{s},{e},{n},{z}\n" for s, n, e, z in _ROWS)
    p = _write(tmp_path, text)
    with pytest.raises(ValueError, match="station"):
        read_station_table(p)          # header has no 'station'
    t = read_station_table(p, columns=("station", "easting",
                                       "northing", "elev"))
    assert len(t) == 3 and t.easting.iloc[0] == pytest.approx(
        792692.448)


def test_ambiguous_layout_raises(tmp_path):
    text = "A,1,2,3\nB,4,5,6\n"
    with pytest.raises(ValueError, match="columns="):
        read_station_table(_write(tmp_path, text))


def test_station_key():
    assert station_key("l2-780.atimpd") == "L2-780"
    assert station_key(" S01 ") == "S01"


def test_to_wgs84_and_invalid_epsg():
    t = pd.DataFrame(_ROWS, columns=["station", "northing", "easting",
                                     "elev"])
    ll = to_wgs84(t, 4508)
    assert ll.lat.iloc[0] == pytest.approx(25.55334, abs=1e-5)
    assert ll.lon.iloc[0] == pytest.approx(113.91176, abs=1e-5)
    with pytest.raises(ValueError, match="not a valid"):
        to_wgs84(t, 15921)


def test_compare_crs_utm_vs_gauss_kruger():
    t = pd.DataFrame(_ROWS, columns=["station", "northing", "easting",
                                     "elev"])
    c = compare_crs(t, [32649, 4508])
    assert c.offset_m.iloc[0] == 0.0
    assert 1000.0 < c.offset_m.iloc[1] < 1300.0


class _Head(SimpleNamespace):
    pass


class _FakeEDI:
    def __init__(self, path, dataid):
        self.path = path
        self._head = _Head(dataid=dataid, lat=0.0, lon=0.0, elev=0.0)
        self._dm = SimpleNamespace(reflat=0.0, reflong=0.0, refelev=0.0)
        self._info = SimpleNamespace(info_text=["    old"])

    def get_section(self, key):
        return {"head": self._head, "definemeas": self._dm,
                "info": self._info}.get(key)


def _ll_table():
    t = pd.DataFrame(_ROWS, columns=["station", "northing", "easting",
                                     "elev"])
    return to_wgs84(t, 4508)


def test_embed_by_name_independent_of_order(tmp_path):
    from pathlib import Path

    edis = [_FakeEDI(Path("L1-60.edi"), "L1-60.atimpd"),
            _FakeEDI(Path("L1-0.edi"), "L1-0.atimpd"),
            _FakeEDI(Path("L1-30.edi"), "L1-30-1.atimpd")]
    t = _ll_table()
    rep = embed_station_coordinates(edis, t, note="CRS EPSG 4508")
    assert rep.matched.all()
    assert rep.station.tolist() == ["L1-60", "L1-0", "L1-30"]
    h = edis[1].get_section("head")
    assert h.lat == pytest.approx(t.lat.iloc[0])
    assert h.elev == pytest.approx(373.721)
    dm = edis[1].get_section("definemeas")
    assert dm.reflat == pytest.approx(t.lat.iloc[0])
    assert edis[1].get_section("info").info_text[0].strip() == \
        "CRS EPSG 4508"


def test_embed_unmatched_and_order_mode():
    from pathlib import Path

    t = _ll_table()
    edis = [_FakeEDI(Path("X-9.edi"), "X-9")]
    with pytest.raises(ValueError, match="no table row"):
        embed_station_coordinates(edis, t)
    rep = embed_station_coordinates(edis, t, require_all=False)
    assert not rep.matched.iloc[0]
    with pytest.raises(ValueError, match="stations but"):
        embed_station_coordinates(edis, t, match="order")


def test_coordinate_injector_default_and_invalid_epsg():
    from pycsamt.exceptions import ValidationError
    from pycsamt.stratagem import CoordinateInjector

    assert CoordinateInjector().epsg is None
    with pytest.raises(ValidationError, match="not a valid"):
        CoordinateInjector(epsg=15921)
