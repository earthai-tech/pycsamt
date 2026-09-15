# test_seg_heads.py

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.exceptions import EdIDataError, FileHandlingError, HeaderError
from pycsamt.seg.heads import (
    Head,
    HeadMixin,
    Heads,
    Info,
    InfoMixin,
    _is_tag,
    _norm_key,
    _slice_section,
    _unquote,
)


@pytest.fixture()
def sample_edi(tmp_path: Path) -> Path:
    lines = [
        ">HEAD",
        "  DATAID=E1_2",
        "  LAT=26:23:00N",
        "  LONG=106:33:00W",
        "  ELEV=1234",
        "  STDVERS=SEG 1.0",
        "",
        ">INFO",
        "  PROJECT=DEMO",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        "  SECTID=E1_2",
        "  NFREQ=2",
        "",
        ">!****FREQUENCIES****!",
        ">FREQ  //2",
        "  1.000000E+02  2.000000E+02",
        "",
        ">ZXXR ROT=ZROT  //2",
        "  1.000000E+00  1.000000E+00",
        "",
        ">END",
    ]
    p = tmp_path / "demo.edi"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


def test_head_from_file_parses_coordinates(sample_edi: Path):
    h = Head.from_file(sample_edi)
    assert h.dataid == "E1_2"
    # 26°23'00"N  => 26 + 23/60 = 26.383333...
    assert h.lat == pytest.approx(26.3833, abs=1e-3)
    # 106°33'00"W => -(106 + 33/60) = -106.55
    assert h.long == pytest.approx(-106.55, abs=1e-3)
    assert h.elev == 1234.0

    out = h.write()
    assert out[0].strip() == ">HEAD"
    # quoted keys preserved
    assert any('DATAID="E1_2"' in ln for ln in out)
    # lat/long serialized (format may vary), but present
    assert any(ln.strip().startswith("LAT=") for ln in out)
    assert any(ln.strip().startswith("LONG=") for ln in out)


def test_info_from_file_routes_fields(sample_edi: Path):
    info = Info.from_file(sample_edi)
    assert info.Source.project == "DEMO"
    assert info.Processing.processedby == "pyCSAMT"
    sw = info.Processing.ProcessingSoftware.name
    assert sw == "pyCSAMT"

    out = info.write()
    assert out[0].strip() == ">INFO"
    assert any("PROJECT=DEMO" in ln for ln in out)
    # quoted value for processedby / processingsoftware
    assert any('PROCESSEDBY="pyCSAMT"' in ln for ln in out)
    assert any('PROCESSINGSOFTWARE="pyCSAMT"' in ln for ln in out)


def test_heads_aggregator_read_write(sample_edi: Path):
    hs = Heads.from_file(sample_edi)
    assert hs.head.dataid == "E1_2"
    assert hs.info.Source.project == "DEMO"

    out = hs.write()
    # Ensure >HEAD comes before >INFO
    text = "".join(out)
    assert text.index(">HEAD") < text.index(">INFO")


def test_head_mixin_instance_read():
    class Host(HeadMixin):
        def __init__(self):
            self.head = Head()

    host = Host()
    kv = ["  DATAID=H1", "  LAT=26:00:00N", "  LONG=10:00:00E"]
    head = host.read(kv)
    assert isinstance(head, Head)
    assert host.head.dataid == "H1"
    assert host.head.lat == pytest.approx(26.0, abs=1e-6)
    assert host.head.long == pytest.approx(10.0, abs=1e-6)

    out = host.write()
    assert out[0].strip() == ">HEAD"


def test_head_mixin_class_from_file(sample_edi: Path):
    # classmethod on mixin returns a Head
    class Host(HeadMixin):
        pass

    h = Host.from_file(sample_edi)
    assert isinstance(h, Head)
    assert h.dataid == "E1_2"


def test_info_mixin_instance_read():
    class Host(InfoMixin):
        def __init__(self):
            self.info = Info()

    host = Host()
    kv = ["  PROJECT=Mix", "  PROCESSEDBY=tool", "  RUNLIST=A,B"]
    info = host.read(kv)
    assert isinstance(info, Info)
    assert host.info.Source.project == "Mix"
    assert host.info.Processing.processedby == "tool"

    out = host.write()
    assert out[0].strip() == ">INFO"


def test_info_mixin_class_from_file(sample_edi: Path):
    # classmethod on mixin returns an Info
    class Host(InfoMixin):
        pass

    info = Host.from_file(sample_edi)
    assert isinstance(info, Info)
    assert info.Source.project == "DEMO"


# ─────────────────────────────────────────────────────────────────────────
# Module-private helpers
# ─────────────────────────────────────────────────────────────────────────


def test_norm_key_lon_alias():
    assert _norm_key("LON") == "long"
    assert _norm_key(" Lat ") == "lat"


def test_unquote_strips_matching_quotes_only():
    assert _unquote('"a b"') == "a b"
    assert _unquote("'a b'") == "a b"
    assert _unquote("no quotes") == "no quotes"
    assert _unquote('"mismatched\'') == '"mismatched\''


def test_is_tag_requires_leading_gt():
    assert _is_tag(">HEAD", ">HEAD") is True
    assert _is_tag("  >info", ">INFO") is True
    assert _is_tag("HEAD", ">HEAD") is False


def test_slice_section_raises_when_start_tag_missing():
    with pytest.raises(EdIDataError):
        _slice_section([">INFO", "x=1"], ">HEAD", after_tags=[">INFO"])


def test_slice_section_falls_back_to_scanning_for_any_tag():
    lines = [">HEAD", "DATAID=X", "MORE=1", ">END"]
    payload, i_start, i_stop = _slice_section(
        lines, ">HEAD", after_tags=[">INFO"],
    )
    assert payload == ["DATAID=X", "MORE=1"]
    assert i_start == 0
    assert i_stop == 3


def test_slice_section_scans_to_end_of_lines_when_no_tag_found_at_all():
    lines = [">HEAD", "DATAID=X", "MORE=1"]  # no closing tag anywhere
    payload, i_start, i_stop = _slice_section(
        lines, ">HEAD", after_tags=[">INFO"],
    )
    assert payload == ["DATAID=X", "MORE=1"]
    assert i_start == 0
    assert i_stop == 3  # == len(lines)


# ─────────────────────────────────────────────────────────────────────────
# Head: setters, from_file errors, read errors, write(explicit lines),
# compute_chainage, as_dict/update, __repr__
# ─────────────────────────────────────────────────────────────────────────


def test_head_kwargs_override_in_constructor():
    h = Head(dataid="X1", country="Testland")
    assert h.dataid == "X1"
    assert h.country == "Testland"


def test_head_lat_long_setters_none_and_dms_verbose(caplog):
    h = Head(verbose=1)
    h.lat = None
    assert h.lat is None
    h.long = None
    assert h.long is None
    h.lat = "26:23:00N"
    assert h.lat == pytest.approx(26.3833, abs=1e-3)
    h.long = "106:33:00W"
    assert h.long == pytest.approx(-106.55, abs=1e-3)


def test_head_lon_alias_getter_and_setter():
    h = Head()
    h.lon = "10:00:00E"
    assert h.lon == pytest.approx(10.0, abs=1e-6)
    assert h.long == h.lon


def test_head_elev_setter_empty_string_is_none():
    h = Head()
    h.elev = "   "
    assert h.elev is None
    h.elev = "100"
    assert h.elev == 100.0


def test_head_from_file_rejects_none_path():
    with pytest.raises(FileHandlingError):
        Head.from_file(None)


def test_head_read_raises_when_no_items():
    h = Head()
    h.edi_header = None
    with pytest.raises(HeaderError):
        h.read()


def test_head_read_handles_chainage_and_lon_and_skips_blank(sample_edi: Path):
    h = Head()
    h.read(
        [
            "",
            ">HEAD",
            "  DATAID=X",
            "  CHAINAGE=12.5",
            "  LON=5:00:00E",
        ]
    )
    assert h.chainage == 12.5
    assert h.long == pytest.approx(5.0, abs=1e-6)


def test_head_read_chainage_invalid_value_sets_none():
    h = Head()
    h.read(["  CHAINAGE=not-a-number"])
    assert h.chainage is None


def test_head_write_with_explicit_lines():
    h = Head()
    out = h.write(
        [
            "  DATAID=X1",
            "  LAT=1.0",
            "  BOGUS_NO_KV_LINE",
            ">SKIP",
        ]
    )
    text = "".join(out)
    assert 'DATAID="X1"' in text
    assert "LAT=1.0" in text


def test_head_write_chainage_invalid_is_skipped():
    h = Head(dataid="X1")
    h.chainage = "not-a-number"
    out = "".join(h.write())
    assert "CHAINAGE=" not in out


def test_head_write_decimal_to_dms_fallback_on_error(monkeypatch):
    import pycsamt.seg.heads as heads_mod

    h = Head(dataid="X1")
    h.lat = 10.0
    monkeypatch.setattr(
        heads_mod,
        "decimal_to_dms",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("boom")),
    )
    out = "".join(h.write())
    assert "LAT=10.0" in out


def test_head_as_dict_and_update():
    h = Head(dataid="X1")
    d = h.as_dict()
    assert d["dataid"] == "X1"
    assert "lat" in d and "long" in d and "elev" in d

    result = h.update(dataid="X2", country="Y")
    assert result is h
    assert h.dataid == "X2"
    assert h.country == "Y"


def test_head_compute_chainage_missing_coordinates_returns_nan():
    h = Head()
    ch = h.compute_chainage(origin=(0.0, 0.0), azimuth=0.0)
    assert ch != ch  # NaN
    assert h.chainage != h.chainage


def test_head_compute_chainage_computes_forward_distance():
    h = Head()
    h.lat = 1.0  # ~111 km north of origin
    h.long = 0.0
    ch = h.compute_chainage(origin=(0.0, 0.0), azimuth=0.0)
    assert ch == pytest.approx(111_000.0, rel=1e-3)
    assert h.chainage == pytest.approx(111_000.0, rel=1e-3)


def test_head_compute_chainage_set_attr_false_does_not_store():
    h = Head(dataid="X1")
    h.chainage = 5.0
    h.compute_chainage(origin=(0.0, 0.0), azimuth=0.0, set_attr=False)
    assert h.chainage == 5.0  # unchanged


def test_head_repr():
    h = Head(dataid="X1")
    assert "Head dataid='X1'" in repr(h)


# ─────────────────────────────────────────────────────────────────────────
# Info: kwargs override, from_file errors/fallback, read edge cases,
# write(explicit lines), as_dict/update, __repr__
# ─────────────────────────────────────────────────────────────────────────


def test_info_kwargs_override_in_constructor():
    info = Info(filter="F1")
    assert info.filter == "F1"


def test_info_from_file_rejects_none_path():
    with pytest.raises(FileHandlingError):
        Info.from_file(None)


def test_info_from_file_falls_back_to_scanning_when_no_explicit_stop(
    tmp_path: Path, monkeypatch,
):
    # IsEdi._assert_edi requires the file's last tag to be >END, which
    # means a genuinely *valid* EDI always gives Info.from_file's own
    # forward tag-scan something to stop on before end-of-file. The
    # "stop is None" scan-to-end-of-file fallback only fires for a file
    # with no closing tag at all, so bypass validation here to isolate
    # that parsing branch on its own.
    monkeypatch.setattr(
        "pycsamt.seg.heads.IsEdi._assert_edi", lambda *a, **k: True,
    )
    text = "\n".join(
        [
            ">HEAD",
            "  DATAID=X",
            "",
            ">INFO",
            "  PROJECT=DEMO",
            "  EXTRA=1",
        ]
    )
    p = tmp_path / "no_explicit_stop.edi"
    p.write_text(text, encoding="utf-8")
    info = Info.from_file(p)
    assert info.Source.project == "DEMO"


def test_info_read_skips_blank_and_tag_lines():
    info = Info()
    info.read(["", ">INFO", "  PROJECT=X"])
    assert info.Source.project == "X"


def test_info_read_invalid_maxinfo_falls_back_to_default():
    info = Info()
    info.read(["  MAXINFO=not-a-number"])
    assert info.maxinfo == 999


def test_info_write_with_explicit_lines_and_text_passthrough():
    info = Info()
    out = info.write(
        [
            "",
            ">SKIP_ME",
            "  PROJECT=demo",
            "  PROCESSEDBY=tool",
            "free text line",
        ]
    )
    text = "".join(out)
    assert "PROJECT=DEMO" in text
    assert 'PROCESSEDBY="tool"' in text
    assert "free text line" in text


def test_info_as_dict():
    info = Info()
    info.Source.project = "P1"
    d = info.as_dict()
    assert d["project"] == "P1"
    assert "processingsoftware" in d


def test_info_update_routes_to_nested_containers():
    info = Info()
    info.update(
        project="P1",
        processingsoftware="SW1",
        processedby="me",
        info_text=["line1"],
        filter="F1",
    )
    assert info.Source.project == "P1"
    assert info.Processing.ProcessingSoftware.name == "SW1"
    assert info.Processing.processedby == "me"
    assert info.info_text == ["line1"]
    assert info.filter == "F1"


def test_info_repr():
    info = Info()
    info.Source.project = "P1"
    text = repr(info)
    assert "Info project='P1'" in text


# ─────────────────────────────────────────────────────────────────────────
# Heads: from_file errors/fallback, read(), to_text, __repr__
# ─────────────────────────────────────────────────────────────────────────


def test_heads_from_file_rejects_none_path():
    with pytest.raises(FileHandlingError):
        Heads.from_file(None)


def test_heads_from_file_missing_info_creates_empty_info(tmp_path: Path):
    text = "\n".join(
        [
            ">HEAD",
            "  DATAID=X",
            "",
            ">=DEFINEMEAS",
            "  MAXCHAN=1",
            "",
            ">!****FREQUENCIES****!",
            ">FREQ  //1",
            "  1.000000E+02",
            "",
            ">END",
        ]
    )
    p = tmp_path / "no_info.edi"
    p.write_text(text, encoding="utf-8")
    hs = Heads.from_file(p)
    assert hs.head.dataid == "X"
    assert hs.info.Source.project is None


def test_heads_read_from_string_and_list(sample_edi: Path):
    text = sample_edi.read_text(encoding="utf-8")
    hs = Heads()
    result = hs.read(text)
    assert result is hs
    assert hs.head.dataid == "E1_2"
    assert hs.info.Source.project == "DEMO"

    hs2 = Heads()
    hs2.read(text.splitlines())
    assert hs2.head.dataid == "E1_2"


def test_heads_read_missing_info_leaves_empty(tmp_path: Path):
    lines = [">HEAD", "  DATAID=X", "", ">=DEFINEMEAS", "  MAXCHAN=1", ""]
    hs = Heads()
    hs.read(lines)
    assert hs.head.dataid == "X"
    assert hs.info.Source.project is None


def test_heads_to_text_and_repr(sample_edi: Path):
    hs = Heads.from_file(sample_edi)
    text = hs.to_text()
    assert ">HEAD" in text and ">INFO" in text
    assert "Heads dataid='E1_2'" in repr(hs)


# ─────────────────────────────────────────────────────────────────────────
# HeadMixin / InfoMixin: auto-create-on-missing-attribute branches
# ─────────────────────────────────────────────────────────────────────────


def test_head_mixin_auto_creates_head_when_missing():
    class Host(HeadMixin):
        pass

    host = Host()
    assert not hasattr(host, "head")
    out = host.write()
    assert out[0].strip() == ">HEAD"
    assert isinstance(host.head, Head)


def test_head_mixin_auto_creates_head_on_read_when_missing():
    class Host(HeadMixin):
        pass

    host = Host()
    head = host.read(["  DATAID=Z1"])
    assert isinstance(head, Head)
    assert host.head.dataid == "Z1"


def test_info_mixin_auto_creates_info_when_missing():
    class Host(InfoMixin):
        pass

    host = Host()
    assert not hasattr(host, "info")
    out = host.write()
    assert out[0].strip() == ">INFO"
    assert isinstance(host.info, Info)


def test_info_mixin_auto_creates_info_on_read_when_missing():
    class Host(InfoMixin):
        pass

    host = Host()
    info = host.read(["  PROJECT=Z1"])
    assert isinstance(info, Info)
    assert host.info.Source.project == "Z1"
