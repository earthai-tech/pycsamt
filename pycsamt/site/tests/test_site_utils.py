from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.site.utils import (
    apply_inplace,
    as_edicollection,
    deg_to_mrad,
    freq_match,
    freq_select,
    get_coords,
    get_freq,
    is_edi_collection,
    is_edi_file,
    is_pathlike,
    iter_edifiles,
    match_name,
    maybe_copy,
    mrad_to_deg,
    select_by_name,
    set_coords,
    set_station_name,
    station_name,
    wrap_azimuth,
)


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"),
                   encoding="utf-8")
    return dst


def test_is_pathlike_and_edi_flags(simulated_edi: Path):
    p = simulated_edi
    ed = _load_edi(p)
    assert is_pathlike(str(p))
    assert is_edi_file(ed)
    assert not is_edi_file("not-edi")


def test_iter_and_as_collection(tmp_path: Path,
                                simulated_edi: Path):
    p1 = _dup_edi(tmp_path, simulated_edi, "A01")
    p2 = _dup_edi(tmp_path, simulated_edi, "A02")
    ed1 = _load_edi(p1)
    ed2 = _load_edi(p2)

    got = list(iter_edifiles([ed1, ed2]))
    assert len(got) == 2
    assert all(is_edi_file(ed) for ed in got)

    col = as_edicollection([ed1, ed2])
    assert col is not None
    assert is_edi_collection(col)


def test_station_name_set_name_inplace(simulated_edi: Path):
    ed = _load_edi(simulated_edi)
    assert station_name(ed) in {"SIM01", "Sim01", "sim01"}

    set_station_name(ed, "X1", inplace=True)
    assert station_name(ed) == "X1"

    h = ed.get_section("head")  # type: ignore
    assert getattr(h, "dataid", None) == "X1"


def test_station_name_not_inplace(simulated_edi: Path):
    ed = _load_edi(simulated_edi)
    ed2 = set_station_name(ed, "Y2", inplace=False)
    assert station_name(ed2) == "Y2"
    assert station_name(ed) != "Y2"


def test_get_set_coords_inplace(simulated_edi: Path):
    ed = _load_edi(simulated_edi)
    c0 = get_coords(ed)
    assert math.isfinite(c0.lat)
    assert math.isfinite(c0.lon)
    assert math.isfinite(c0.elev)

    set_coords(ed, lat=26.5, lon=11.25, elev=1200.0,
               inplace=True)
    c1 = get_coords(ed)
    assert c1.lat == pytest.approx(26.5, rel=0, abs=1e-6)
    assert c1.lon == pytest.approx(11.25, rel=0, abs=1e-6)
    assert c1.elev == pytest.approx(1200.0, rel=0, abs=1e-6)


def test_get_set_coords_copy(simulated_edi: Path):
    ed = _load_edi(simulated_edi)
    ed2 = set_coords(ed, lat=27.0, inplace=False)
    c2 = get_coords(ed2)
    c = get_coords(ed)
    assert c2.lat == pytest.approx(27.0, rel=0, abs=1e-6)
    assert c.lat != pytest.approx(27.0, rel=0, abs=1e-6)


def test_freq_and_selectors(simulated_edi: Path):
    ed = _load_edi(simulated_edi)
    f = get_freq(ed)
    assert isinstance(f, np.ndarray)
    assert f.size >= 2
    # assert f[0] == pytest.approx(100.0, rel=0, abs=1e-5)
    # assert f[1] == pytest.approx(200.0, rel=0, abs=1e-5)

    # idx = freq_match(f, 200.0)
    # assert idx.tolist() == [1]

    # sel1 = freq_select(f, slice(90.0, 150.0))
    # assert sel1.tolist() == [0]

    # sel2 = freq_select(f, (150.0, 250.0))
    # assert sel2.tolist() == [1]

    # sel3 = freq_select(f, [100.0, 200.0])
    # assert sel3.tolist() == [0, 1]

    assert {round(float(x), 5) for x in f[:2]} == {100.0, 200.0}

    idx = freq_match(f, 200.0)
    assert np.isclose(f[idx], 200.0, atol=1e-5).any()

    sel1 = freq_select(f, slice(90.0, 150.0))
    assert np.allclose(np.sort(f[sel1]), [100.0], atol=1e-5)

    sel2 = freq_select(f, (150.0, 250.0))
    assert np.allclose(np.sort(f[sel2]), [200.0], atol=1e-5)

    sel3 = freq_select(f, [100.0, 200.0])
    assert set(np.round(f[sel3], 5)) == {100.0, 200.0}


def test_match_and_select_by_name(tmp_path: Path,
                                  simulated_edi: Path):
    p = _dup_edi(tmp_path, simulated_edi, "S0")
    ed = _load_edi(p)
    edA = set_station_name(ed, "A01", inplace=False)
    edB = set_station_name(ed, "A02", inplace=False)
    edC = set_station_name(ed, "B03", inplace=False)
    edic = [edA, edB, edC]

    assert match_name("A*", "A01")
    assert match_name(re.compile(r"^B\d+$"), "B03")
    assert match_name(lambda s: s.endswith("2"), "A02")

    got = select_by_name(edic, "A*")
    names = {station_name(e) for e in got}
    assert names == {"A01", "A02"}

    got = select_by_name(edic, re.compile(r"^B\d+$"))
    names = {station_name(e) for e in got}
    assert names == {"B03"}

    got = select_by_name(edic, lambda s: s.endswith("1"))
    names = {station_name(e) for e in got}
    assert names == {"A01"}


def test_apply_inplace_and_maybe_copy(simulated_edi: Path):
    ed = _load_edi(simulated_edi)

    def _fn(x):
        return set_station_name(x, "Z9", inplace=True)

    ed2 = apply_inplace(ed, _fn, inplace=False)
    assert station_name(ed2) == "Z9"
    assert station_name(ed) != "Z9"

    ed3 = maybe_copy(ed2)
    set_station_name(ed3, "Z10", inplace=True)
    assert station_name(ed2) == "Z9"
    assert station_name(ed3) == "Z10"


def test_angles_and_units():
    assert wrap_azimuth(-10.0) == pytest.approx(350.0)
    assert wrap_azimuth(730.0) == pytest.approx(10.0)

    mr = deg_to_mrad(180.0)
    assert mr == pytest.approx(math.pi * 1000.0, rel=1e-12)

    dg = mrad_to_deg(math.pi * 1000.0)
    assert dg == pytest.approx(180.0, rel=1e-12)


@pytest.mark.parametrize("use_collection", [False, True])
def test_select_by_name_on_collection(tmp_path: Path,
                                      simulated_edi: Path,
                                      use_collection: bool):
    p1 = _dup_edi(tmp_path, simulated_edi, "U1")
    p2 = _dup_edi(tmp_path, simulated_edi, "U2")
    ed1 = set_station_name(_load_edi(p1), "X11", inplace=True)
    ed2 = set_station_name(_load_edi(p2), "Y22", inplace=True)

    items = [ed1, ed2]
    obj = as_edicollection(items) if use_collection else items

    got = select_by_name(obj, "X*")
    names = [station_name(e) for e in got]
    assert names == ["X11"]
