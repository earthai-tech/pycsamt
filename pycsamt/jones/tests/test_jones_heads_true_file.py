# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

from pycsamt.jones.heads import (
    Info,
    Head,
    Heads,
    Banner,
)
from pycsamt.jones.property import JSiteProperty
from pycsamt.jones.config import (
    RE_STATION,
    RE_DATATYPE_UNITS,
)


def _assert_head_valid(h: Head) -> None:
    assert h.station is not None
    assert RE_STATION.match(h.station)
    assert h.dtype is not None
    token = f"{h.dtype.kind}{h.dtype.comp}"
    assert RE_DATATYPE_UNITS.match(token)
    assert h.n is not None and h.n >= 0


def test_info_from_file_true(j_single_file: Path):
    info = Info.from_file(j_single_file)
    assert info.__has_read__() is True
    for k in info.items.keys():
        assert k == k.upper()
    # convenience props
    lat = info.latitude
    lon = info.longitude
    az = info.azimuth
    el = info.elevation
    assert lat is not None and -90.0 <= lat <= 90.0
    assert lon is not None and -180.0 <= lon < 180.0
    if az is not None:
        assert 0.0 <= az < 360.0 or -360.0 < az < 360.0
    if el is not None:
        assert isinstance(el, float)


def test_banner_from_file_true(j_single_file: Path):
    b = Banner.from_file(j_single_file)
    assert b.__has_read__() is True
    assert isinstance(b.software, str) or b.software is None
    if b.software:
        assert b.software.strip() != ""
    if b.station_hint:
        assert isinstance(b.station_hint, str)


def test_head_parse_station_with_azimuth_token():
    hdr = ["kb0001  -30", "RXY", "29"]
    h = Head.from_lines(hdr)
    assert h.station == "KB0001"
    _assert_head_valid(h)
    assert h.n == 29


def test_heads_from_file_true_properties(j_single_file: Path):
    heads = Heads.from_file(j_single_file)
    assert heads.__has_read__() is True
    assert heads.n in (0, 1)
    if heads.n == 1:
        _assert_head_valid(heads.head)
        # station exposure
        assert heads.station is not None
        assert RE_STATION.match(heads.station)
        # site props exposure
        if heads.latitude is not None:
            assert -90.0 <= heads.latitude <= 90.0
        if heads.longitude is not None:
            assert -180.0 <= heads.longitude < 180.0
        # azimuth may come from info or station hint
        if heads.azimuth is not None:
            assert -360.0 < heads.azimuth < 360.0
        # banner/software present or None
        if heads.software is not None:
            assert isinstance(heads.software, str)


def test_heads_write_roundtrip(j_single_file: Path):
    heads = Heads.from_file(j_single_file)
    out = heads.write()
    assert isinstance(out, list)
    assert all(isinstance(s, str) for s in out)
    # Look for a station/dtype/count triple in the output
    st, dt, ct = None, None, None
    for s in out:
        if st is None and RE_STATION.match(s):
            st = s
            continue
        if st is not None and dt is None and (
            RE_DATATYPE_UNITS.match(s)
        ):
            dt = s
            continue
        if st is not None and dt is not None and ct is None:
            if s.strip().isdigit():
                ct = s
                break
    assert st is not None and dt is not None and ct is not None


def test_info_write_roundtrip(j_single_file: Path):
    info = Info.from_file(j_single_file)
    lines = info.write()
    info2 = Info.from_lines(lines)
    assert info2.latitude == info.latitude
    assert info2.longitude == info.longitude
    assert info2.elevation == info.elevation
    assert info2.azimuth == info.azimuth


def test_property_from_file_true(j_single_file: Path):
    prop = JSiteProperty.from_file(j_single_file, strict=False)
    if prop.azimuth is not None:
        assert 0.0 <= prop.azimuth < 360.0 or (
            -360.0 < prop.azimuth < 360.0
        )
    if prop.latitude is not None:
        assert -90.0 <= prop.latitude <= 90.0
    if prop.longitude is not None:
        assert -180.0 <= prop.longitude < 180.0
