# test_seg_heads.py

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.seg.heads import (
    Head,
    HeadMixin,
    Heads,
    Info,
    InfoMixin,
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
