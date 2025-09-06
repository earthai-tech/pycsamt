# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from pycsamt.seg.other import OtherSECT, OtherIO, OtherMixin
from pycsamt.exceptions import EdIDataError


def _write(tmp: Path, name: str, text: str) -> Path:
    p = tmp / name
    p.write_text(text, encoding="utf-8")
    return p


def test_othersect_from_file_parses_header_and_ids(
    tmp_path: Path,
):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="OTH_A"

        >=OTHERSECT
          SECTID=A1
          NCHAN=3
          NFREQ=2
          MAXBLKS=7
          TYPE=CUSTOM
            HX HY
            HZ

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "other_hdr.edi", edi)

    sect = OtherSECT.from_file(str(p))
    assert sect.sectid == "A1"
    assert sect.nchan == 3
    assert sect.nfreq == 2
    assert sect.maxblks == 7
    assert sect.type == "CUSTOM"
    # measurement IDs gathered from free lines
    assert sect.meas_ids == ["HX", "HY", "HZ"]
    # body pointer points past header
    assert isinstance(sect.start_data_lines_num, int)


def test_otherio_from_file_parses_numeric_and_raw_blocks(
    tmp_path: Path,
):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="OTH_B"

        >=OTHERSECT
          SECTID=B1
          NCHAN=2
          NFREQ=1
          MAXBLKS=4
          TYPE=FREE
            HX HY

        >COH ID=HX ROT=0 // 4
          0.1 0.2
          0.3 0.4

        >ANNO NOTE=TEXT // 2
          alpha beta
          gamma

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "other_blocks.edi", edi)

    sect = OtherSECT.from_file(str(p))
    io = OtherIO.from_file(
        str(p),
        start_line=sect.start_data_lines_num,
    )
    assert len(io.blocks) == 2

    b0, b1 = io.blocks
    # first, numeric block
    assert b0.keyword == ">COH"
    # options are typed best-effort
    assert b0.options["id"] == "HX"
    assert b0.options["rot"] == 0
    assert b0.nitems_hint == 4
    assert b0.values == [0.1, 0.2, 0.3, 0.4]
    assert b0.raw_lines == []

    # second, raw text block
    assert b1.keyword == ">ANNO"
    assert b1.options["note"] == "TEXT"
    assert b1.values == []
    assert [r.strip() for r in b1.raw_lines] == ["alpha beta", "gamma"]
    assert b1.nitems_hint == 2


def test_otherio_write_roundtrip(tmp_path: Path):
    # Build an IO with two blocks and serialize it
    io = OtherIO()
    # numeric block
    nb = io.blocks.append
    from pycsamt.seg.other import _OtherBlock

    b0 = _OtherBlock()
    b0.keyword = ">COH"
    b0.options = {"id": "HX", "rot": 1}
    b0.nitems_hint = 3
    b0.values = [1.0, 2.0, 3.0]
    nb(b0)

    # raw block
    b1 = _OtherBlock()
    b1.keyword = ">META"
    b1.options = {"note": "X"}
    b1.raw_lines = ["line-a", "line-b"]
    nb(b1)

    lines = io.write()
    text = "".join(lines)

    # headers keep key order (ID/ROT first), keys upper-cased
    assert ">COH ID=HX ROT=1 // 3" in text
    # numeric values formatted using Base.FLOAT_FMT
    assert "  1.000000E+00" in text
    # raw block uses raw_lines and counts them
    assert ">META NOTE=X // 2" in text
    assert "line-a" in text and "line-b" in text


def test_other_mixin_helpers(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="OTH_C"

        >=OTHERSECT
          SECTID=C1
          NCHAN=1
            HX

        >COH ID=HX // 2
          9 10

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "other_mixin.edi", edi)

    class Host(OtherMixin):
        pass

    hdr = Host.read_other_header(str(p))
    assert hdr.sectid == "C1"
    io = Host.read_other_blocks(str(p))
    assert len(io.blocks) == 1
    assert io.blocks[0].values == [9.0, 10.0]


def test_otherio_raises_when_no_blocks(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="EMPTY"

        >=OTHERSECT
          SECTID=E0
          NCHAN=0

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "other_empty.edi", edi)

    sect = OtherSECT.from_file(str(p))
    with pytest.raises(EdIDataError):
        OtherIO.from_file(
            str(p),
            start_line=sect.start_data_lines_num + 1,
        )
