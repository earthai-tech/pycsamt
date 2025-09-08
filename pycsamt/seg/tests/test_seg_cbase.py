# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

from pycsamt.seg.cbase import CoreParser, CBBase
from pycsamt.seg.edi import EDIFile


def _mk_edi(tmp: Path, station: str, nf: int = 2) -> Path:
    fvals = np.geomspace(1.0, 10.0, nf)
    fstr = "  " + "  ".join(f"{v: .6E}" for v in fvals)
    lines: List[str] = [
        ">HEAD",
        f"  DATAID={station}",
        "  LAT=26:00:00N",
        "  LONG=010:00:00E",
        "  ELEV=1000",
        "",
        ">INFO",
        "  PROJECT=SIM",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        f"  SECTID={station}",
        f"  NFREQ={nf}",
        "",
        ">!****FREQUENCIES****!",
        f">FREQ  //{nf}",
        fstr,
        "",
        ">END",
    ]
    p = tmp / f"{station}.edi"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


def test_coreparser_parse_file_and_dir(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "A1", nf=3)
    f2 = _mk_edi(tmp_path, "A2", nf=4)

    pr = CoreParser(recursive=True, strict=False, on_dup="replace")
    items = pr.parse(tmp_path)

    # returns EDIFile objects
    assert isinstance(items, list)
    assert all(isinstance(x, EDIFile) for x in items)
    assert {ed.station for ed in items} == {"A1", "A2"}

    # parse a single file path too
    items2 = pr.parse([str(f1)])
    assert len(items2) == 1
    assert items2[0].station == "A1"

    # no unexpected hard errors collected
    assert pr.errors() == []


def test_coreparser_glob_and_errors(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "B1", nf=2)
    _mk_edi(tmp_path, "B2", nf=2)

    # include a bogus path to exercise error capture
    pr = CoreParser(recursive=False, strict=False, on_dup="replace")
    items = pr.parse([tmp_path / "*.edi", tmp_path / "NOPE.edi"])
    assert len(items) == 2
    assert len(pr.errors()) >= 1  # NOPE.edi recorded


def test_cbbase_add_iter_summary(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "C1", nf=2)
    f2 = _mk_edi(tmp_path, "C2", nf=5)

    ed1 = EDIFile(f1)
    ed2 = EDIFile(f2)

    col = CBBase(items=[ed1])
    assert len(col) == 1
    col.add(ed2)
    assert len(col) == 2

    st = col.stations()
    assert set(st) == {"C1", "C2"}

    sm = col.summary()
    by = {r["station"]: r for r in sm}
    assert by["C1"]["n_freq"] == 2
    assert by["C2"]["n_freq"] == 5

    # iteration yields EDIFile
    assert all(isinstance(x, EDIFile) for x in col)


def test_cbbase_repr_str(tmp_path: Path) -> None:
    f = _mk_edi(tmp_path, "D1", nf=3)
    ed = EDIFile(f)
    col = CBBase(items=[ed])
    s = str(col)
    r = repr(col)
    # light sanity: contains station and count
    assert "D1" in s
    assert "1" in r
