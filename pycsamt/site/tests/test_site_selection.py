# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np

from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Sites
from pycsamt.site import selection as sel
from pycsamt.site import edit as ed


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"),
                   encoding="utf-8")
    return dst


def _mk_edifiles(
    tmp_path: Path,
    simulated_edi: Path,
    *stems: str,
) -> List[EDIFile]:
    out: List[EDIFile] = []
    for st in stems:
        p = _dup_edi(tmp_path, simulated_edi, st)
        out.append(_load_edi(p))
    return out


def test_by_names_exact_wildcard_regex_callable(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2, e3 = _mk_edifiles(tmp_path, simulated_edi,
                              "A01", "A02", "B01")
    s = Sites([e1, e2, e3])

    # exact (case-insensitive)
    out = sel.by_names(s, "a01", case=False)
    assert isinstance(out, Sites)
    assert len(out) == 1
    assert out.by_index(0).name == "A01"

    # wildcard
    out = sel.by_names(s, "A*", case=False)
    assert len(out) == 2
    got = sorted([x.name for x in out])
    assert got == ["A01", "A02"]

    # regex
    import re as _re
    out = sel.by_names(s, _re.compile(r"^B\d+"))
    assert len(out) == 1 and out.by_index(0).name == "B01"

    # callable
    out = sel.by_names(s, lambda n: n.endswith("1"))
    assert len(out) == 2
    got = sorted([x.name for x in out])
    assert got == ["A01", "B01"]


def test_by_index_positive_and_negative(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2, e3 = _mk_edifiles(tmp_path, simulated_edi,
                              "I01", "I02", "I03")
    s = Sites([e1, e2, e3])

    out = sel.by_index(s, [0, -1, 10])
    assert isinstance(out, Sites)
    got = [x.name for x in out]
    assert got == ["I01", "I03"]

    out = sel.by_index(s, -2)
    assert len(out) == 1 and out.by_index(0).name == "I02"


def test_by_chainage_head_and_fallback_attr(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2 = _mk_edifiles(tmp_path, simulated_edi, "C01",
                          "C02")
    # set HEAD.chainage
    h1 = e1.get_section("head")  # type: ignore
    h2 = e2.get_section("head")  # type: ignore
    setattr(h1, "chainage", 100.0)
    setattr(h2, "chainage", 150.0)

    s = Sites([e1, e2])
    out = sel.by_chainage(s, 120.0, 200.0)
    assert len(out) == 1 and out.by_index(0).name == "C02"

    # fallback: remove head.chainage, set edi.chainage
    delattr(h2, "chainage")
    setattr(e2, "chainage", 180.0)
    out = sel.by_chainage(s, 170.0, 190.0)
    assert len(out) == 1 and out.by_index(0).name == "C02"


def test_by_freq_range_intersection(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2 = _mk_edifiles(tmp_path, simulated_edi, "F01",
                          "F02")
    s = Sites([e1, e2])

    # simulated EDI has 100 and 200 Hz
    out = sel.by_freq(s, 150.0, 250.0)
    assert len(out) == 2

    out = sel.by_freq(s, 1000.0, 2000.0)
    assert len(out) == 0


def test_by_bbox_with_coordinate_update(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2 = _mk_edifiles(tmp_path, simulated_edi, "B01",
                          "B02")
    # move stations apart
    e1 = ed.set_coords(e1, lat=0.0, lon=0.0, elev=0.0,
                       inplace=False)
    e2 = ed.set_coords(e2, lat=10.0, lon=10.0, elev=0.0,
                       inplace=False)

    s = Sites([e1, e2])
    out = sel.by_bbox(s, -1.0, -1.0, 1.0, 1.0)
    assert len(out) == 1 and out.by_index(0).name == "B01"

    out = sel.by_bbox(s, 9.0, 9.0, 11.0, 11.0)
    assert len(out) == 1 and out.by_index(0).name == "B02"


def test_by_predicate_lambda(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2, e3 = _mk_edifiles(tmp_path, simulated_edi,
                              "P01", "Q02", "P03")
    s = Sites([e1, e2, e3])

    out = sel.by_predicate(s, lambda edf: ed.station_name(
        edf
    ).startswith("P"))
    # Note: station_name is also in utils; we reuse bound
    # from selection's import through the Sites wrapper at
    # call time.
    assert len(out) == 2
    got = sorted([x.name for x in out])
    assert got == ["P01", "P03"]


def test_keep_finite_z_and_drop_empty(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2 = _mk_edifiles(tmp_path, simulated_edi, "Z01",
                          "Z02")
    # Make Z01 have finite z; Z02 becomes empty
    e1 = ed.fill_missing(e1, how="zero", components=("Z",),
                         inplace=False)
    # empty Z section object, no arrays
    e2.Z = type("Z", (), {})()  # type: ignore

    s = Sites([e1, e2])
    kept = sel.keep_finite_z(s)
    assert len(kept) == 1 and kept.by_index(0).name == "Z01"

    # drop_empty removes e2
    dropped = sel.drop_empty(s)
    assert len(dropped) == 1 and dropped.by_index(0).name == "Z01"


def test_mask_large_phase_err_threshold_and_no_err(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2, e3 = _mk_edifiles(tmp_path, simulated_edi,
                              "E01", "E02", "E03")
    # Allocate finite z and phase_err for E01/E02
    e1 = ed.fill_missing(e1, how="zero", components=("Z",),
                         inplace=False)
    e2 = ed.fill_missing(e2, how="zero", components=("Z",),
                         inplace=False)

    # Inject phase_err
    z1 = getattr(e1, "Z", None)
    z2 = getattr(e2, "Z", None)
    assert z1 is not None and z2 is not None

    f = np.asarray([100.0, 200.0], float)
    # small errors -> should be kept under thresh=10
    setattr(z1, "_phase_err", np.full((2, 2, 2), 5.0, float))
    # large errors -> should be dropped under thresh=10
    setattr(z2, "_phase_err", np.full((2, 2, 2), 50.0, float))

    # E03 has no phase_err -> should be kept
    s = Sites([e1, e2, e3])
    out = sel.mask_large_phase_err(s, thresh=10.0)
    got = sorted([x.name for x in out])
    assert got == ["E01", "E03"]
