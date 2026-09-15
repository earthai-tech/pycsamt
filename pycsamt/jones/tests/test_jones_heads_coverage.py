from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.jones.heads import (
    Banner,
    Head,
    HeadMixin,
    Heads,
    Info,
    InfoMixin,
    _extract_first_head_list,
)

# ─────────────────────────────────────────────────────────────────────────
# Banner
# ─────────────────────────────────────────────────────────────────────────


def test_banner_constructor_with_top_lines_parses_immediately():
    b = Banner(top_lines=["#WRITTEN BY GEOTOOLS: KB0001 10/06/95 RAW RECS"])
    assert b.software == "GEOTOOLS"
    assert b.station_hint == "KB0001"


def test_banner_read_raises_when_top_lines_missing():
    with pytest.raises(ValueError):
        Banner().read()


def test_banner_write_new_false_uses_existing_software_and_date():
    b = Banner.from_lines(
        ["#WRITTEN BY GEOTOOLS: KB0001 10/06/95 RAW RECS"]
    )
    out = b.write(new=False, include_origin=True)
    assert out[0].startswith("#WRITTEN BY GEOTOOLS:")
    assert "10/06/95" in out[0]
    assert out[1].startswith("#FROM")


def test_banner_write_new_false_defaults_when_no_prior_parse():
    b = Banner()
    out = b.write(new=False)
    assert out[0].startswith("#WRITTEN BY PYCSAMT:")


def test_banner_date_parsed_valid_and_invalid_and_none():
    b1 = Banner.from_lines(["#WRITTEN BY X: S01 10/06/95 NOTE"])
    assert b1.date_parsed is not None

    b2 = Banner()
    assert b2.date_parsed is None

    b3 = Banner.from_lines(["#WRITTEN BY X: S01 not-a-date NOTE"])
    assert b3.date_parsed is None


# ─────────────────────────────────────────────────────────────────────────
# Info
# ─────────────────────────────────────────────────────────────────────────


def _info_lines():
    return [
        "# leading comment",
        ">AZIMUTH = 30.0",
        ">LATITUDE = 12.0",
        ">LONGITUDE = -15.0",
    ]


def test_info_from_lines_raises_when_none():
    with pytest.raises(ValueError):
        Info.from_lines(None)


def test_info_read_raises_when_none():
    with pytest.raises(ValueError):
        Info().read(None)


def test_info_write_with_explicit_lines_reads_first():
    info = Info()
    out = info.write(_info_lines())
    assert any("AZIMUTH" in ln for ln in out)


def test_info_accessors_and_dunders():
    info = Info.from_lines(_info_lines())
    assert info.get("AZIMUTH").startswith("30.0")
    assert info.get("NOPE", "default") == "default"
    assert "AZIMUTH" in info.keys()
    assert any(v.startswith("30.0") for v in info.values())
    im = info.items_map()
    assert im["AZIMUTH"].startswith("30.0")
    assert info.lat == info.latitude
    assert info.lon == info.longitude
    assert "azimuth" in info  # __contains__ is case-insensitive
    assert info["AZIMUTH"].startswith("30.0")
    assert str(info) == "Info(items=3)"
    assert "Info(items=3" in repr(info)


# ─────────────────────────────────────────────────────────────────────────
# Head
# ─────────────────────────────────────────────────────────────────────────


def test_head_from_file(tmp_path: Path):
    p = tmp_path / "h.j"
    p.write_text("S01\nZXY SI\n3\n", encoding="utf-8")
    h = Head.from_file(p)
    assert h.station == "S01"
    assert h.n == 3


def test_head_from_lines_raises_when_none():
    with pytest.raises(ValueError):
        Head.from_lines(None)


def test_head_from_lines_extracts_triple_from_longer_slice():
    lines = ["# comment", ">INFO = 1", "S01", "ZXY SI", "3"]
    h = Head.from_lines(lines)
    assert h.station == "S01" and h.n == 3


def test_head_read_raises_when_none():
    with pytest.raises(ValueError):
        Head().read(None)


def test_head_read_skips_leading_comment_info_blank_lines():
    lines = ["# c", ">AZIMUTH = 1", "", "S01", "ZXY SI", "3"]
    h = Head().read(lines)
    assert h.station == "S01"


def test_head_read_raises_when_no_station_line():
    with pytest.raises(ValueError):
        Head().read(["not a station!!", "ZXY SI", "3"])


def test_head_read_raises_when_missing_dtype_line():
    with pytest.raises(ValueError):
        Head().read(["S01"])


def test_head_read_raises_when_dtype_line_invalid():
    with pytest.raises(ValueError):
        Head().read(["S01", "@@@not-a-dtype@@@", "3"])


def test_head_read_raises_when_count_line_missing():
    with pytest.raises(ValueError):
        Head().read(["S01", "ZXY SI"])


def test_head_read_tolerates_blank_lines_between_fields():
    h = Head().read(["S01", "", "ZXY SI", "", "3"])
    assert h.station == "S01" and h.n == 3


def test_head_write_with_head_list_infos_reads_first():
    h = Head()
    out = h.write(["S01", "ZXY SI", "3"])
    assert out == ["S01", "ZXY SI", "3"]


def test_head_write_raises_when_not_fully_defined():
    with pytest.raises(ValueError):
        Head().write()


def test_head_properties_and_header_tuple():
    h = Head().read(["S01", "ZXY SI", "3"])
    assert h.units is not None
    assert h.tensor_hint is None or isinstance(h.tensor_hint, str)
    assert h.header == ("S01", "ZXY", 3)


def test_head_header_property_none_when_incomplete():
    h = Head()
    assert h.header is None
    h.station = "S01"
    assert h.header is None  # dtype still None


def test_head_str_and_repr():
    h = Head().read(["S01", "ZXY SI", "3"])
    assert str(h) == "Head(station='S01', n=3)"
    assert repr(h) == str(h)


# ─────────────────────────────────────────────────────────────────────────
# Heads
# ─────────────────────────────────────────────────────────────────────────


def test_heads_read_accepts_raw_text_string():
    text = "\n".join(
        [
            ">AZIMUTH = 5.0",
            "S01",
            "ZXY SI",
            "2",
        ]
    )
    heads = Heads().read(text)
    assert heads.station == "S01"


def test_heads_str_and_repr():
    heads = Heads()
    assert str(heads) == "Heads(n=0)"
    assert repr(heads) == str(heads)


# ─────────────────────────────────────────────────────────────────────────
# HeadMixin / InfoMixin
# ─────────────────────────────────────────────────────────────────────────


def test_head_mixin_from_file(tmp_path: Path):
    p = tmp_path / "h.j"
    p.write_text("S01\nZXY SI\n3\n", encoding="utf-8")
    h = HeadMixin.from_file(p)
    assert h.station == "S01"


def test_head_mixin_read_and_write_lazily_create_head():
    class Host(HeadMixin):
        pass

    host = Host()
    assert not hasattr(host, "head") or host.head is None
    out_head = host.read(["S01", "ZXY SI", "3"])
    assert out_head.station == "S01"
    lines = host.write()
    assert lines == ["S01", "ZXY SI", "3"]


def test_info_mixin_from_file(tmp_path: Path):
    p = tmp_path / "i.j"
    p.write_text(">AZIMUTH = 5.0\n", encoding="utf-8")
    info = InfoMixin.from_file(p)
    assert info.get("AZIMUTH").startswith("5.0")


def test_info_mixin_read_and_write_lazily_create_info():
    class Host(InfoMixin):
        pass

    host = Host()
    out_info = host.read([">AZIMUTH = 5.0"])
    assert out_info.get("AZIMUTH").startswith("5.0")
    lines = host.write()
    assert any("AZIMUTH" in ln for ln in lines)


# ─────────────────────────────────────────────────────────────────────────
# _extract_first_head_list: case A / case B / error branches
# ─────────────────────────────────────────────────────────────────────────


def test_extract_first_head_list_raises_when_no_station_at_all():
    with pytest.raises(ValueError):
        _extract_first_head_list(["# only a comment", ">INFO = 1"])


def test_extract_first_head_list_raises_when_missing_header_after_station():
    with pytest.raises(ValueError):
        _extract_first_head_list(["S01"])


def test_extract_first_head_list_raises_bad_header_after_station():
    with pytest.raises(ValueError):
        _extract_first_head_list(["S01", "not_valid_at_all!!"])


def test_extract_first_head_list_case_b_count_then_dtype():
    lines = ["S01", "5", "garbage row that is not a dtype", "ZXY SI"]
    out = _extract_first_head_list(lines)
    assert out == ["S01", "ZXY SI", "5"]


def test_extract_first_head_list_case_b_raises_when_no_dtype_found():
    lines = ["S01", "5", "garbage1", "garbage2"]
    with pytest.raises(ValueError):
        _extract_first_head_list(lines)


def test_extract_first_head_list_skips_leading_noise_before_station():
    lines = ["# c", ">INFO = 1", "", "S01", "ZXY SI", "3"]
    out = _extract_first_head_list(lines)
    assert out == ["S01", "ZXY SI", "3"]


def test_extract_first_head_list_tolerates_blank_between_dtype_and_count():
    lines = ["S01", "ZXY SI", "", "3"]
    out = _extract_first_head_list(lines)
    assert out == ["S01", "ZXY SI", "3"]
