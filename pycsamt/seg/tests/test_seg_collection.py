# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

from pycsamt.seg.collection import EDICollection


def _mk_edi(tmp: Path, station: str, nf: int = 3) -> Path:
    fvals = np.geomspace(1.0, 100.0, nf)
    fstr = "  " + "  ".join(f"{v: .6E}" for v in fvals)
    lines: List[str] = [
        ">HEAD",
        f"  DATAID={station}",
        "  LAT=26:00:00N",
        "  LONG=010:00:00E",
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


def test_from_sources_and_add_from(tmp_path: Path) -> None:
    a = _mk_edi(tmp_path, "S1", nf=2)
    b = _mk_edi(tmp_path, "S2", nf=4)

    # from directory
    col = EDICollection.from_sources(tmp_path)
    assert len(col) == 2
    assert set(col.stations()) == {"S1", "S2"}

    # add from glob + bogus path (tolerated)
    col2 = EDICollection(verbose=0)
    col2.add_from([tmp_path / "*.edi", tmp_path / "NOPE.edi"])
    assert len(col2) == 2
    assert set(col2.stations()) == {"S1", "S2"}


def test_select_where_sort_and_stats(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "A", nf=2)
    _mk_edi(tmp_path, "B", nf=5)
    _mk_edi(tmp_path, "C", nf=3)

    col = EDICollection.from_sources(tmp_path)

    sel = col.select(["A", "C"])
    assert set(sel.stations()) == {"A", "C"}

    wh = col.where(lambda ed: ed.Z.n_freq >= 3)
    assert set(wh.stations()) == {"B", "C"}

    s1 = col.sort("station")
    assert s1.stations() == sorted(col.stations())

    s2 = col.sort("n_freq", reverse=True)
    nf = [r["n_freq"] for r in s2.summary()]
    assert nf == sorted(nf, reverse=True)

    stats = col.nf_stats()
    assert stats["min"] >= 2
    assert stats["max"] <= 5
    assert stats["mean"] >= 2.0


def test_merge_on_dup(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "M1", nf=2)
    _mk_edi(tmp_path, "M2", nf=3)

    c1 = EDICollection.from_sources(tmp_path)
    # make an overlapping collection
    c2 = EDICollection.from_sources(tmp_path).select(["M2"])

    kept = c1.merge(c2, on_dup="keep")
    assert len(kept) == len(c1)

    repl = c1.merge(c2, on_dup="replace")
    # still same count; but station M2 comes from c2
    assert len(repl) == len(c1)
    assert "M2" in repl.stations()


def test_paths_and_repr_str(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "P1", nf=2)
    _mk_edi(tmp_path, "P2", nf=2)

    col = EDICollection.from_sources(tmp_path)
    ps = col.paths
    assert all(p.endswith(".edi") for p in ps)

    s = str(col)
    r = repr(col)
    # light sanity: names appear
    assert "P1" in s and "P2" in s
    assert "stations" in r
