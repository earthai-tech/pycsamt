# -*- coding: utf-8 -*-
from pathlib import Path

from pycsamt.jones.heads import Info, Head, Heads

def _sample_lines():
    # comment + info block
    lines = [
        "# a",
        "# b",
        ">AZIMUTH = 45.0",
        ">LATITUDE = 10.0",
        ">LONGITUDE = -20.0",
        "",
        # first block head (ZXY SI) + 2 rows
        "S01",
        "ZXY SI",
        "2",
        " 1.0  1.0  0.0  0.1  1.0",
        " 2.0  0.5 -0.1  0.2  1.0",
        # second block head (RXY) + 1 row
        "S02",
        "RXY",
        "1",
        " 1.0  100  45  110  90  50  40  1  1",
    ]
    return lines


def test_info_from_file(tmp_path: Path):
    p = tmp_path / "j.j"
    p.write_text("\n".join(_sample_lines()), encoding="utf-8")

    info = Info.from_file(p)
    assert info.__has_read__() is True
    assert "AZIMUTH" in info.items
    assert info.items["AZIMUTH"].startswith("45.0")
    assert len(info.comments) == 2


def test_head_read_write_roundtrip():
    head_list = ["S01", "ZXY SI", "3"]
    h = Head().read(j_header_list=head_list)
    assert h.station == "S01"
    assert h.dtype.kind == "Z" and h.dtype.comp == "XY"
    assert h.dtype.units in ("ohms", "mV/km/nT")
    assert h.n == 3

    out = h.write()
    assert out[0] == "S01"
    assert out[1].startswith("ZXY")
    assert out[2] == "3"


def test_heads_from_lines_and_write():
    lines = _sample_lines()
    heads = Heads.from_lines(lines)

    assert heads.__has_read__() is True
    assert heads.info is not None
    # With the new API, Heads holds a single Head
    assert heads.n == 1
    assert heads.head.station == "S01"

    rendered = heads.write()
    joined = "\n".join(rendered)
    # should contain info header and the first head triple
    assert ">AZIMUTH" in joined
    assert "S01" in joined
