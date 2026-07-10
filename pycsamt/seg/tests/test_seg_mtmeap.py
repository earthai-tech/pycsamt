# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import textwrap
from pathlib import Path

from pycsamt.seg.mtemap import (
    MTEMAP,
    EMAPComponents,
    EMAPMixin,
)


def _write(tmp_path: Path, name: str, body: str) -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


def test_mt_from_file_parses_minimal_header(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="S01"

        >=MTSECT
          SECTID="S01"
          NFREQ=60
          HX=251.025
          HY=252.025
          HZ=253.025
          EX=254.025
          EY=255.025

        >FREQ
          1.0  2.0  3.0

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "mt.edi", edi)

    m = MTEMAP.from_file(str(p))
    assert m.sectid == "S01"
    assert m.nfreq == 60
    assert m.hx == "251.025"
    assert m.hy == "252.025"
    assert m.hz == "253.025"
    assert m.ex == "254.025"
    assert m.ey == "255.025"
    assert m.start_data_lines_num is not None


def test_mt_from_file_uses_dataid_as_fallback(tmp_path: Path):
    # SECTID looks numeric -> should fallback to DATAID
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="SITE_A"

        >=MTSECT
          SECTID=123.4
          NFREQ=2
          HX=1
          HY=2

        >FREQ
          1.0  2.0

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "mt_num_sectid.edi", edi)

    m = MTEMAP.from_file(str(p))
    assert m.sectid == "SITE_A"


def test_emap_from_file_parses_and_writes_header(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="EM01"

        >=EMAPSECT
          SECTID="EM01"
          NFREQ=3
          HX=11
          HY=12
          RX=61
          RY=62
          NDIPOLE=3
          TYPE=PROFILE

        >FREQ
          10  20  30

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "emap.edi", edi)

    m = MTEMAP.from_file(str(p))
    assert m.sectid == "EM01"
    assert m.nfreq == 3
    assert m.hx == "11"
    assert m.hy == "12"
    assert m.rx == "61"
    assert m.ry == "62"
    assert m.ndipole == 3
    assert m.type.upper() == "PROFILE"

    out = m.write()
    # First line must be EMAP header
    assert out[0].strip() == ">=EMAPSECT"
    # EMAP write omits EX/EY/HZ lines
    joined = "".join(out).upper()
    assert " EX=" not in joined
    assert " EY=" not in joined
    assert " HZ=" not in joined
    # Must include NDIPOLE/TYPE
    assert " NDIPOLE=" in joined
    assert " TYPE=" in joined


def test_emap_components_from_mtemap(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="S02"

        >=MTSECT
          SECTID="S02"
          NFREQ=1
          HX=101
          HY=102
          HZ=103
          EX=201
          EY=202
          RX=301
          RY=302

        >FREQ
          1

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "mt_comp.edi", edi)

    m = MTEMAP.from_file(str(p))
    comp = EMAPComponents.from_mtemap(m)
    d = comp.as_dict()

    assert d["hx"] == "101"
    assert d["hy"] == "102"
    assert d["hz"] == "103"
    assert d["ex"] == "201"
    assert d["ey"] == "202"
    assert d["rx"] == "301"
    assert d["ry"] == "302"


def test_emap_mixin_from_file_returns_mtemap(tmp_path: Path):
    # Enforce consistent API: mixin should call MTEMAP.from_file
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="MX"

        >=MTSECT
          SECTID="MX"
          NFREQ=1

        >FREQ
          1

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "mixin.edi", edi)

    class Host(EMAPMixin):
        pass

    obj = Host.from_file(str(p))
    assert isinstance(obj, MTEMAP)
    assert obj.sectid == "MX"
