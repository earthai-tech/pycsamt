"""Focused branch coverage for Jones helpers and validation."""

from __future__ import annotations

import io

import pytest

from pycsamt.exceptions import JError, JParseError
from pycsamt.jones.components import JComponentMixin
from pycsamt.jones.utils import (
    _fmt_err,
    _normalize_period,
    _normalize_units,
    _to_float,
    is_blank,
    is_comment,
    iter_info,
    iter_lines,
    iter_rows,
    parse_datatype_units,
    parse_info,
    parse_npoints,
    parse_row,
    parse_station,
    strip_nondata,
)
from pycsamt.jones.validation import IsJ, _scan_j_blocks, is_j_file
from pycsamt.jones.property import (
    JSiteProperty,
    _coerce_lat,
    _coerce_lon,
    _maybe_missing_float,
    _norm_azimuth,
    _parse_latitude,
    _parse_longitude,
    _to_float_safe,
)


class Bag(JComponentMixin):
    pass


class Writer:
    def __init__(self, value): self.value = value
    def write(self): return self.value


class BrokenWriter:
    def write(self): raise RuntimeError("bad")


def test_component_mixin_crud_accessors_and_headers():
    bag = Bag()
    assert bag.cget("missing", 3) == 3
    for name in ("banner", "info", "head", "heads", "blocks", "z", "tip", "res"):
        value = object()
        getattr(bag, f"set_{name}")(value)
        assert getattr(bag, f"get_{name}")() is value
        assert bag.chas(name.upper())
    assert set(bag.snapshot([" banner ", "absent"])) == {"banner", "absent"}
    bag.cdrop("RES")
    assert not bag.chas("res")

    bag.set_heads(Writer(["combined", 2]))
    assert bag.compose_headers() == ["combined", "2"]
    bag.set_banner(Writer("banner\nline2"))
    bag.set_head(Writer(["head"]))
    bag.set_info(BrokenWriter())
    assert bag.compose_headers(prefer_heads=False, join=True) == "banner\nline2\nhead"
    bag.cset("info", object())
    assert "head" in bag.compose_headers(prefer_heads=False)

    alternate = Bag()
    alternate.components = {"x": 1}
    assert alternate.cget("x") == 1
    alternate.components = []
    assert alternate._bag() == {}


def test_line_and_scalar_parsers(tmp_path):
    path = tmp_path / "lines.txt"
    path.write_text("a\nb\n")
    assert list(iter_lines(path)) == ["a", "b"]
    assert list(iter_lines(io.StringIO("a\nb\n"), keepends=True)) == ["a\n", "b\n"]
    assert list(iter_lines(["a\n", 2])) == ["a", "2"]
    assert is_comment("# hi") and is_blank("   ")
    assert parse_info(">AZIMUTH = 45") == ("azimuth", "45")
    assert parse_station("ab12") == "AB12"
    assert parse_npoints("3") == 3
    dtype = parse_datatype_units("RTE S.I.")
    assert dtype.kind == "R" and dtype.comp == "TE" and dtype.tensor_hint
    assert _normalize_units("field") == "mV/km/nT"
    assert _normalize_period(-2) == (0.5, 2, True)
    assert _normalize_period(4) == (4, 0.25, False)
    with pytest.raises(JParseError, match="positive"):
        _normalize_period(0)
    assert _to_float("1.5", "x", None) == 1.5
    with pytest.raises(JParseError, match="Bad float"):
        _to_float("bad", "x", 4)
    assert "..." in _fmt_err("bad", 2, "x" * 100)


@pytest.mark.parametrize(
    ("func", "value"),
    [(parse_info, "bad"), (parse_station, "###"),
     (parse_datatype_units, "bad type"), (parse_npoints, "x")],
)
def test_malformed_scalar_parsers(func, value):
    with pytest.raises(JParseError):
        func(value, lineno=9)


def test_row_parsing_and_iterators():
    rline = "-2 100 45 110 90 50 40 1 1"
    tline = "2 1 -2 .1 1"
    r = parse_row("R", rline)
    assert r.freq == 2 and r.flags["stored_as_freq"]
    assert not r.flags["rejected"]
    assert parse_row("Z", tline).period == 2
    with pytest.raises(JParseError, match="Malformed R/S"):
        parse_row("R", "bad")
    with pytest.raises(JParseError, match="Malformed TF"):
        parse_row("Z", "bad")
    with pytest.raises(JParseError, match="Unsupported"):
        parse_row("X", tline)
    assert list(strip_nondata(["# c", "data", "! c"])) == ["data", "! c"]
    assert list(iter_info(["", "# c", ">AZIMUTH=1", "SITE", ">LAT=2"])) == [("azimuth", "1")]
    assert len(list(iter_rows("Z", ["", "# c", tline]))) == 1


def test_deep_and_shallow_j_validation(tmp_path):
    valid = tmp_path / "valid.j"
    valid.write_text("# banner\n>AZIMUTH=0\nS01\nZXY SI\n1\n1 2 3 .1 1\n")
    assert is_j_file(valid)
    assert is_j_file(valid, deep=False)
    assert _scan_j_blocks(valid.read_text().splitlines()) == 1

    count_first = tmp_path / "count.txt"
    count_first.write_text("S01\n1\n\nRXY\n1 100 45 110 90 50 40 1 1\n")
    assert _scan_j_blocks(count_first.read_text().splitlines()) == 0
    wrong_ext = tmp_path / "valid.csv"
    wrong_ext.write_text(valid.read_text())
    with pytest.raises(JError, match="extension"):
        is_j_file(wrong_ext, deep=False)
    invalid = tmp_path / "invalid.j"
    invalid.write_text("S01\nZXY SI\n1\nnot a row\n")
    with pytest.raises(JError, match="no valid"):
        is_j_file(invalid)
    with pytest.raises(FileNotFoundError):
        is_j_file(tmp_path / "missing")
    with pytest.raises(Exception, match="NoneType"):
        is_j_file(None)


def test_site_property_mapping_and_numeric_helpers():
    prop = JSiteProperty.from_mapping(
        {"azimuth": 370, "latitude": 95, "longitude": 190,
         "elevation": 12, "custom": "value"}, verbose=1
    )
    assert prop.azimuth == 10
    assert prop.location == (90, -170)
    assert prop.azimuth_rad == pytest.approx(0.1745329)
    assert prop.asdict()["elevation"] == 12
    assert prop.extra
    assert _maybe_missing_float("bad") is None
    assert _maybe_missing_float("1.5") == 1.5
    with pytest.raises(ValueError, match="Invalid DMS"):
        _parse_latitude("bad", strict=False, verbose=0)
    with pytest.raises(ValueError, match="Invalid DMS"):
        _parse_longitude("bad", strict=False, verbose=0)
    assert _coerce_lat(100) == 90 and _coerce_lat(-100) == -90
    assert _coerce_lon(190) == -170
    assert _norm_azimuth(-10) == 350
    assert _to_float_safe("2") == 2
    with pytest.raises(ValueError, match="empty"):
        _to_float_safe("")
    empty = JSiteProperty()
    assert empty.location is None and empty.azimuth_rad is None
