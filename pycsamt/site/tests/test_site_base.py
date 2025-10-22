# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile 
from pycsamt.site.base import Site, Sites


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"),
                   encoding="utf-8")
    return dst


def _mk_two_sites(
    tmp_path: Path,
    simulated_edi: Path,
    n1: str = "S01",
    n2: str = "S02",
) -> Tuple[Site, Site]:

    p1 = _dup_edi(tmp_path, simulated_edi, n1)
    p2 = _dup_edi(tmp_path, simulated_edi, n2)
    e1 = _load_edi(p1)
    e2 = _load_edi(p2)
    s1 = Site(e1)
    s2 = Site(e2)
    return s1, s2


def test_sites_container_index_get(tmp_path: Path,
                                   simulated_edi: Path) -> None:
    s1, s2 = _mk_two_sites(tmp_path, simulated_edi, "A01",
                           "A02")
    sites = Sites([s1.edi, s2.edi])
    assert len(sites) == 2
    # int index
    assert isinstance(sites[0], Site)
    # name lookup is case-insensitive
    assert sites["a01"].name.lower() == "a01"
    # by_index/get helpers
    assert sites.by_index(1).name == "A02"
    assert sites.get("A02").name == "A02"  # type: ignore
    assert sites.get("missing") is None
    # as_list returns EDIFile objects
    lst = sites.as_list()
    assert len(lst) == 2
    assert hasattr(lst[0], "get_section")


def test_site_coords_set_and_summary(tmp_path: Path,
                                     simulated_edi: Path) -> None:
    s1, _ = _mk_two_sites(tmp_path, simulated_edi, "C01",
                          "C02")
    # in-place update
    s1.set_coords(10.0, 20.0, 100.0, inplace=True)
    la, lo, ev = s1.coords
    assert (la, lo, ev) == (10.0, 20.0, 100.0)
    summ = s1.summary()
    assert summ["name"] == "C01"
    assert np.isfinite(summ["nfreq"])


def test_sites_closest(tmp_path: Path,
                       simulated_edi: Path) -> None:
    s1, s2 = _mk_two_sites(tmp_path, simulated_edi, "D01",
                           "D02")
    s1.set_coords(0.0, 0.0, 0.0, inplace=True)
    s2.set_coords(0.0, 2.0, 0.0, inplace=True)
    sites = Sites([s1.edi, s2.edi])

    a = sites.closest(0.0, 0.1)
    b = sites.closest(0.0, 1.9)
    assert a and a.name == "D01"
    assert b and b.name == "D02"

    # tolerance: ~1110 m for 0.01 degree at equator
    near = sites.closest(0.0, 0.01, tol=500.0)
    assert near is None


def test_sites_edit_all_rename_and_slice(tmp_path: Path,
                                         simulated_edi: Path) -> None:
    s1, s2 = _mk_two_sites(tmp_path, simulated_edi, "E01",
                           "E02")
    sites = Sites([s1.edi, s2.edi])

    def rnm(n: str) -> str:
        return f"X_{n}"

    out = sites.edit_all(rename=rnm, inplace=False)
    assert out["X_E01"].name == "X_E01"
    assert sites["E01"].name == "E01"

    # frequency slicing keeps alignment across arrays
    pre = Site(sites["E01"].edi).freq
    if pre is None or len(pre) < 2:
        pytest.skip("simulated EDI has <2 freq")
    sl = slice(1, None)
    out2 = sites.edit_all(freq_slice=sl, inplace=False)
    post = Site(out2["E01"].edi).freq
    assert post is not None
    assert len(post) == len(pre) - 1


def test_sites_edit_all_mask(tmp_path: Path,
                             simulated_edi: Path) -> None:
    s1, _ = _mk_two_sites(tmp_path, simulated_edi, "F01",
                          "F02")
    sites = Sites([s1.edi])

    df = Site(sites["F01"].edi).to_dataframe("z")
    if df.empty:
        pytest.skip("no Z data in simulated EDI")

    # mask out half of rows
    def mk_mask(frame):
        m = np.ones(len(frame), dtype=bool)
        m[: len(frame) // 2] = False
        return m

    out = sites.edit_all(mask=mk_mask, inplace=False)
    df2 = Site(out["F01"].edi).to_dataframe("z")
    if df2.empty:
        pytest.skip("no Z data after edit")

    # first half must be NaN now
    first_row = df2.iloc[0]
    assert np.all(np.isnan(first_row.values))


def test_sites_with_topography(tmp_path: Path,
                               simulated_edi: Path) -> None:
    pd = pytest.importorskip("pandas")
    s1, s2 = _mk_two_sites(tmp_path, simulated_edi, "T01",
                           "T02")
    sites = Sites([s1.edi, s2.edi])

    # grab station ids from HEAD.dataid
    def sid(site: Site) -> str:
        h = site.edi.get_section("head")  # type: ignore
        return str(getattr(h, "dataid"))

    df = pd.DataFrame(
        {
            "station": [sid(Site(s1.edi)), sid(Site(s2.edi))],
            "latitude": [35.125, 36.0],
            "longitude": [12.75, 13.25],
            "elevation": [1234.0, 50.0],
        }
    )

    out = sites.with_topography(df, inplace=False)
    h1 = out["T01"].edi.get_section("head")  # type: ignore
    h2 = out["T02"].edi.get_section("head")  # type: ignore
    got1 = (float(h1.lat), float(h1.lon), float(h1.elev))
    got2 = (float(h2.lat), float(h2.lon), float(h2.elev))
    assert got1 == (35.125, 12.75, 1234.0)
    assert got2 == (36.0, 13.25, 50.0)

    # original unchanged
    h1o = sites["T01"].edi.get_section("head")  # type: ignore
    assert (float(h1o.lat), float(h1o.lon)) != (35.125, 12.75)


def test_sites_to_profile_order(tmp_path: Path,
                                simulated_edi: Path) -> None:
    s1, s2 = _mk_two_sites(tmp_path, simulated_edi, "P01",
                           "P02")
    s1.set_coords(0.0, 0.0, 0.0, inplace=True)
    s2.set_coords(0.0, 1.0, 0.0, inplace=True)
    sites = Sites([s1.edi, s2.edi])

    prof = sites.to_profile((0.0, 0.0), 90.0)

    # Accept both return shapes: Profile or fallback dict
    if hasattr(prof, "chainages"):
        ch = prof.chainages  # type: ignore[attr-defined]
        # P02 should be farther east -> larger chainage
        assert ch["P02"] > ch["P01"]
    else:
        seq: List[Site] = prof["sites"]  # type: ignore[index]
        assert [s.name for s in seq] == ["P01", "P02"]


def test_sites_write(tmp_path: Path, simulated_edi: Path) -> None:
    s1, s2 = _mk_two_sites(tmp_path, simulated_edi, "W01",
                           "W02")
    sites = Sites([s1.edi, s2.edi])

    out = sites.write(tmp_path, template="{station}.edi",
                      exist_ok=True)
    assert len(out) == 2
    for p in out:
        assert p.exists()

