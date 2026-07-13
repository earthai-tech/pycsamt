from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from pycsamt.jones.blocks import JBlocks
from pycsamt.jones.heads import Heads
from pycsamt.jones.j import JFile


def _read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def _has_token(text: str, token: str) -> bool:
    return re.search(rf"\b{re.escape(token)}\b", text) is not None


def _assert_blocks_have_kinds(col: JBlocks, required: set[str]) -> None:
    kinds = {b.head.dtype.kind for b in col.blocks}  # type: ignore
    missing = required.difference(kinds)
    assert not missing, f"Missing kinds: {missing}"


def _assert_blocks_have_comps(col: JBlocks, tokens: set[str]) -> None:
    # tokens like "ZXY", "RXY", ...
    seen = {
        f"{b.head.dtype.kind}{b.head.dtype.comp}"  # type: ignore
        for b in col.blocks
    }
    missing = tokens.difference(seen)
    assert not missing, f"Missing components: {missing}"


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_roundtrip_R_only(tmp_path: Path, j_single_file: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    out = jf.write(
        savepath=tmp_path,
        new_jfn="roundtrip_R.j",
        datatype="R",
        overwrite=True,
    )
    outp = Path(out)
    assert outp.exists()
    txt = _read_text(outp)
    # banner + R
    assert txt.startswith("#WRITTEN BY")
    assert "RXY" in txt or "RYX" in txt

    # structural parse
    col = JBlocks.from_file(outp, verbose=0)
    assert col.__has_read__() is True
    assert col.n >= 1
    _assert_blocks_have_kinds(col, {"R"})

    # heads/info preserved
    hd = Heads.from_file(outp, verbose=0)
    assert hd.__has_read__() is True
    assert hd.station is not None
    # frequency/period sanity
    assert jf.freq is not None and jf.freq.size >= 1
    # periods computed from freq should match writer periods
    p_written = []
    for b in col.blocks:
        # first column is period
        p_written.extend([float(r.period) for r in b.rows])  # type: ignore
        break
    if p_written:
        pr = np.asarray(p_written, float)
        # if original jf.freq finite, compare one value
        fi = float(jf.freq[0])  # type: ignore
        if math.isfinite(fi) and fi > 0:
            assert math.isclose(1.0 / fi, pr[0], rel_tol=1e-6, abs_tol=1e-9)


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_write_ZR_selection(tmp_path: Path, j_single_file: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    # ensure Z can be derived from R
    assert jf.Z is not None or jf.Res is not None

    out = jf.write(
        savepath=tmp_path,
        new_jfn="zr.j",
        datatype="ZR",
        overwrite=True,
    )
    outp = Path(out)
    assert outp.exists()
    col = JBlocks.from_file(outp, verbose=0)
    assert col.__has_read__() is True
    assert col.n >= 2
    _assert_blocks_have_kinds(col, {"Z", "R"})
    # check that at least one off-diagonal appears
    txt = _read_text(outp)
    assert _has_token(txt, "ZXY") or _has_token(txt, "ZYX")


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_write_Z_only(tmp_path: Path, j_single_file: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    # For files without TF, Z is derived from R; should exist.
    assert jf.Z is not None

    out = jf.write(
        savepath=tmp_path,
        new_jfn="z_only.j",
        datatype="Z",
        overwrite=True,
    )
    outp = Path(out)
    assert outp.exists()

    col = JBlocks.from_file(outp, verbose=0)
    assert col.__has_read__() is True
    assert col.n >= 1
    _assert_blocks_have_kinds(col, {"Z"})
    _assert_blocks_have_comps(col, {"ZXX", "ZXY", "ZYX", "ZYY"})


def test_write_no_overwrite_suffix(tmp_path: Path, j_single_file: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    base = tmp_path / "dup.j"
    out1 = jf.write(
        savepath=tmp_path, new_jfn=base.name, datatype="R", overwrite=False
    )
    out2 = jf.write(
        savepath=tmp_path, new_jfn=base.name, datatype="R", overwrite=False
    )
    p1 = Path(out1)
    p2 = Path(out2)
    assert p1.exists() and p2.exists()
    assert p1 != p2
    assert p2.stem.startswith(p1.stem) or p2.stem.endswith("_1")


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_writer_headers_and_station(tmp_path: Path, j_single_file: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    out = jf.write(
        savepath=tmp_path,
        new_jfn="hdr_check.j",
        datatype="R",
        overwrite=True,
    )
    outp = Path(out)
    txt = _read_text(outp)
    assert txt.startswith("#WRITTEN BY")
    hd = Heads.from_file(outp, verbose=0)
    assert hd.__has_read__() is True
    # station carried through and uppercase
    assert isinstance(hd.station, str) and hd.station == hd.station.upper()
