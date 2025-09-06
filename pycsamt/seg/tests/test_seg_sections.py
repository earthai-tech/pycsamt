# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import textwrap
from pathlib import Path

from pycsamt.seg.sections import (
    SECT_REGISTRY,
    iter_sections,
)
from pycsamt.seg.mtemap import MTEMAP
from pycsamt.seg.spectra import SpectraSECT
from pycsamt.seg.time_series import TSect
from pycsamt.seg.other import OtherSECT


def _write(tmp: Path, name: str, text: str) -> Path:
    p = tmp / name
    p.write_text(text, encoding="utf-8")
    return p


def test_registry_contains_expected_tags():
    keys = set(SECT_REGISTRY.keys())
    assert ">=MTSECT" in keys
    assert ">=EMAPSECT" in keys
    assert ">=SPECTRASECT" in keys
    assert ">=TSERIESSECT" in keys
    assert ">=OTHERSECT" in keys

    assert SECT_REGISTRY[">=MTSECT"] is MTEMAP
    assert SECT_REGISTRY[">=EMAPSECT"] is MTEMAP
    assert SECT_REGISTRY[">=SPECTRASECT"] is SpectraSECT
    assert SECT_REGISTRY[">=TSERIESSECT"] is TSect
    assert SECT_REGISTRY[">=OTHERSECT"] is OtherSECT


def test_iter_sections_yields_in_order_and_types(
    tmp_path: Path,
):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="S"

        >=MTSECT
          SECTID=M1
          NFREQ=2
          HX=1
          HY=2
        >FREQ
          1 2

        >=SPECTRASECT
          SECTID=SP
          NCHAN=3
        >SPECTRA FREQ=1.0 // 0

        >=TSERIESSECT
          SECTID=TS
          NCHAN=1
        >TSERIES ID=HX // 0

        >=OTHERSECT
          SECTID=OT
          NCHAN=1
        >COH ID=HX // 0

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "sections.edi", edi)

    got = list(iter_sections(str(p)))
    tags = [t for (t, _, _) in got]
    assert tags == [
        ">=MTSECT",
        ">=SPECTRASECT",
        ">=TSERIESSECT",
        ">=OTHERSECT",
    ]

    # type checks
    assert isinstance(got[0][1], MTEMAP)
    assert isinstance(got[1][1], SpectraSECT)
    assert isinstance(got[2][1], TSect)
    assert isinstance(got[3][1], OtherSECT)

    # MTEMAP.read() should fill sectid
    assert getattr(got[0][1], "sectid") == "M1"

    # body_start points to a new block line
    with open(p, "r", encoding="utf-8") as f:
        lines = f.readlines()
    for _, _, start in got:
        assert lines[start].lstrip().startswith(">")


def test_iter_sections_include_filter(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="S"

        >=MTSECT
          SECTID=A
        >FREQ
          1

        >=SPECTRASECT
          SECTID=B
        >SPECTRA FREQ=1.0 // 0

        >=OTHERSECT
          SECTID=C
        >COH // 0

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "filter.edi", edi)

    only = list(
        iter_sections(
            str(p),
            include=[">=spectrasect", ">=othersect"],
        )
    )
    tags = [t for (t, _, _) in only]
    assert tags == [">=SPECTRASECT", ">=OTHERSECT"]


def test_iter_sections_ignores_unknown_blocks(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="X"

        >=FOOBAR
          SOME=THING

        >=OTHERSECT
          SECTID=OT
        >COH // 0

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "unknown.edi", edi)

    got = list(iter_sections(str(p)))
    # Only the known OTHERSECT should be yielded
    assert len(got) == 1
    assert got[0][0] == ">=OTHERSECT"
    assert isinstance(got[0][1], OtherSECT)


def test_iter_sections_empty_file_yields_nothing(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="EMPTY"
        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "empty.edi", edi)

    got = list(iter_sections(str(p)))
    assert got == []
