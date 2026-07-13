from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.jones.blocks import JBlocks
from pycsamt.jones.cbase import (
    JCBBase,
    JCoreParser,
    JParseMixin,
)
from pycsamt.jones.collection import JCollection
from pycsamt.jones.heads import Heads
from pycsamt.jones.j import JFile


class _Dummy(JParseMixin):
    """Lightweight concrete to access JParseMixin helpers."""

    pass


def test_as_path_and_is_j_path(j_single_file: Path, tmp_path: Path):
    d = _Dummy()
    # _as_path coerces and resolves
    out = d._as_path(str(j_single_file))
    assert isinstance(out, Path) and out.exists() and out.is_file()

    # _is_j_path checks extension and file-ness
    assert d._is_j_path(out) is True
    not_j = tmp_path / "foo.bin"
    not_j.write_text("x", encoding="utf-8")
    assert d._is_j_path(not_j) is False
    assert d._is_j_path(tmp_path) is False  # directory


def test_iter_paths_and_iter_j_files(project_root: Path, j_single_file: Path):
    d = _Dummy()
    # _iter_paths yields Path objects
    got = list(d._iter_paths([j_single_file, str(j_single_file)]))
    assert all(isinstance(p, Path) for p in got)

    # _iter_j_files accepts folder and glob patterns
    data_j = project_root / "data" / "j"
    paths = list(d._iter_j_files([data_j]))
    # assert any(p == j_single_file for p in paths)

    assert paths  # we found at least one candidate in the folder
    # and ensure they look like Jones candidates
    assert all(p.suffix.lower() in {".j", ".txt", ".dat"} for p in paths)

    # glob within folder
    patterns = [
        str(data_j / "*.j"),
        str(data_j / "*.txt"),
        str(data_j / "*.dat"),
    ]
    paths2 = list(d._iter_j_files(patterns))
    assert paths2  # we should find at least one candidate
    assert all(p.suffix.lower() in {".j", ".txt", ".dat"} for p in paths2)


def test_fast_station_scan(j_single_file: Path):
    d = _Dummy()
    sta = d._fast_station(j_single_file)
    assert isinstance(sta, str) and len(sta) >= 2


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_coreparser_parse_single(j_single_file: Path):
    pr = JCoreParser(
        recursive=True, strict=False, on_dup="replace", verbose=0
    )
    items = pr.parse([j_single_file])
    assert isinstance(items, list) and len(items) >= 1
    assert all(isinstance(it, JFile) for it in items)

    errs = pr.errors()
    assert isinstance(errs, list)


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_jcbbase_load_and_write(tmp_path: Path, j_single_file: Path):
    col = JCBBase.load([j_single_file], recursive=False, verbose=0)
    assert isinstance(col, JCBBase)
    assert len(col) >= 1

    # summary is a list of dicts with expected keys
    rows = col.summary()
    assert isinstance(rows, list) and rows
    for r in rows:
        assert "station" in r and "n_freq" in r

    # write produces files; keep it light (R only to be safe)
    out = col.write(tmp_path, pattern="{station}.j", datatype="R")
    assert isinstance(out, list) and out
    for p in out:
        assert Path(p).exists()


def test_collection_from_single_file(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    assert isinstance(col, JCollection)
    assert len(col) >= 1

    for it in list(col):
        assert isinstance(it, JFile)
        # headers parse
        hd = Heads.from_file(it.path, verbose=0)  # type: ignore[arg-type]
        assert hd.__has_read__() is True and hd.station is not None
        # blocks parse (some files may have only R)
        bl = JBlocks.from_file(it.path, verbose=0)  # type: ignore[arg-type]
        assert bl.__has_read__() is True
        assert bl.n >= 0


def test_collection_from_folder(project_root: Path):
    data_j = project_root / "data" / "j"
    if not data_j.exists():
        pytest.skip("Jones data folder not present")

    col = JCollection.from_sources([data_j], verbose=0)
    assert len(col) >= 1

    paths = [jf.path for jf in col if hasattr(jf, "path")]
    assert len(paths) == len(set(paths))
    for p in paths:
        assert Path(p).exists()


def test_collection_metadata_and_filtering(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    assert len(col) == 1
    it = list(col)[0]

    # metadata exposure
    site = getattr(it, "site", None) or getattr(it, "name", None)
    assert isinstance(site, str)

    # select by station
    sub = col.select([site])
    assert isinstance(sub, JCollection) and len(sub) >= 1

    # where with predicate
    sub2 = col.where(lambda jf: getattr(jf, "site", None) == site)
    assert isinstance(sub2, JCollection) and len(sub2) >= 1


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_collection_roundtrip_headers(tmp_path: Path, j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    it = list(col)[0]
    assert isinstance(it, JFile)

    # compose header-only lines and re-parse heads
    hdr = it.compose_headers()
    assert isinstance(hdr, list) and all(isinstance(s, str) for s in hdr)
    hd2 = Heads.from_lines(hdr, verbose=0)
    assert hd2.__has_read__() is True

    # Write a minimal J (R-only) and read back via collection
    out = it.write(
        savepath=tmp_path,
        new_jfn="col_roundtrip.j",
        datatype="R",
        overwrite=True,
    )
    outp = Path(out)
    assert outp.exists()

    reread = JCollection.from_sources([outp], verbose=0)
    assert len(reread) >= 1


def test_collection_repr_and_str(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    s = str(col)
    r = repr(col)
    assert isinstance(s, str) and isinstance(r, str)
    assert "JCollection" in s or "JCollection" in r
