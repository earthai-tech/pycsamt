from __future__ import annotations

import importlib.util as _imp
import math
from pathlib import Path

import pytest

from pycsamt.jones.blocks import (
    JBlocks,
    RBlock,
    TFBlock,
)
from pycsamt.jones.config import (
    RE_DATATYPE_UNITS,
    RE_NPOINTS,
    RE_STATION,
)

# --------------------------- helpers ---------------------------


def _r_block_lines() -> list[str]:
    # Head triple (RXY) + 3 rows (9 floats each)
    return [
        "S01",
        "RXY",
        "3",
        # p    rho   pha  rhomax rhomin phamax phamin wrho wpha
        "-1.0  100   45   110    90     50     40     1    1",
        " 2.0  -5    30   35     25     40     20     1    1",
        " 3.0   50   10   -999   40     20     -999   -1   1",
    ]


def _tf_block_lines() -> list[str]:
    # Head triple (ZXY SI) + 2 rows (5 floats each)
    return [
        "S02",
        "ZXY SI",
        "2",
        # p      real   imag  err   w
        "-1.0e+1  1.0   -2.0  0.1   1.0",
        " 5.0     -999  3.0   -999  -1.0",
    ]


def _join(*chunks: list[str]) -> list[str]:
    out: list[str] = []
    for c in chunks:
        out.extend(c)
    return out


# ---------------------------- tests ----------------------------


def test_rblock_read_normalize_numpy_qa():
    lines = _r_block_lines()
    blk = RBlock.from_lines(lines)
    assert blk.__has_read__() is True
    assert blk.head is not None
    assert blk.head.station == "S01"
    assert blk.nrows == 3

    # normalization checks
    a = blk.to_numpy()
    # period sign: -1.0 Hz -> 1.0 s
    assert math.isclose(a["period"][0], 1.0, rel_tol=0, abs_tol=1e-9)
    # -999 -> NaN
    assert math.isnan(a["rhomax"][2])
    assert math.isnan(a["phamin"][2])
    # rejection rules
    # RBlock test
    assert bool(a["rej"][1]) is True  # rho < 0
    assert bool(a["rej"][2]) is True  # wrho < 0
    assert bool(a["rej"][0]) is False  # not rejected

    # QA summary shape/content
    qa = blk.qa_summary()
    assert qa["kind"] == "R/S"
    assert qa["nrows"] == 3
    assert qa["rejected"] == 2
    assert qa["kept"] == 1
    assert qa["period_min"] > 0.0
    assert qa["period_max"] >= qa["period_min"]


def test_tfblock_read_normalize_numpy_qa():
    lines = _tf_block_lines()
    blk = TFBlock.from_lines(lines)
    assert blk.__has_read__() is True
    assert blk.head is not None
    assert blk.head.station == "S02"
    assert blk.nrows == 2

    a = blk.to_numpy()
    # period sign: -10 Hz -> 0.1 s
    assert math.isclose(a["period"][0], 0.1, rel_tol=0, abs_tol=1e-9)
    # -999 -> NaN
    assert math.isnan(a["real"][1])
    assert math.isnan(a["error"][1])
    # rejection by weight < 0

    # TFBlock test
    assert bool(a["rej"][0]) is False
    assert bool(a["rej"][1]) is True

    qa = blk.qa_summary()
    assert qa["kind"] == "TF"
    assert qa["nrows"] == 2
    assert qa["rejected"] == 1
    assert qa["kept"] == 1
    assert qa["period_min"] > 0.0
    assert qa["period_max"] >= qa["period_min"]
    assert isinstance(qa["error_rms"], float)


def test_blocks_from_lines_multiple_blocks():
    lines = _join(_r_block_lines(), _tf_block_lines())
    col = JBlocks.from_lines(lines)
    assert col.__has_read__() is True
    assert col.n == 2

    # write() must reproduce two header triples
    out = col.write()
    txt = "\n".join(out)
    assert "S01" in txt

    # basic header sanity in output
    st_seen = sum(1 for s in out if RE_STATION.match(s))
    dt_seen = sum(1 for s in out if RE_DATATYPE_UNITS.match(s))
    ct_seen = sum(1 for s in out if RE_NPOINTS.match(s))
    assert st_seen >= 2 and dt_seen >= 2 and ct_seen >= 2


def test_blocks_from_file_true(j_single_file: Path):
    col = JBlocks.from_file(j_single_file)
    assert col.__has_read__() is True
    assert col.n >= 1

    # all blocks should yield positive periods after normalize
    for b in col.blocks:
        a = b.to_numpy()
        assert (a["period"] > 0).all()

        # kind token must be one of J families
        k = b.head.dtype.kind
        assert k in {"R", "S", "Z", "Q", "C", "T"}

    # QA list shape matches number of blocks
    qas = col.qa_summary()
    assert isinstance(qas, list) and len(qas) == col.n


@pytest.mark.skipif(
    not _imp.find_spec("pandas"),
    reason="pandas not installed",
)
def test_to_dataframe_optional_dep():
    blk = RBlock.from_lines(_r_block_lines())
    df = blk.to_dataframe()
    need = {
        "period",
        "rho",
        "pha",
        "rhomax",
        "rhomin",
        "phamax",
        "phamin",
        "wrho",
        "wpha",
        "rej",
    }
    assert set(df.columns) >= need


def test_write_roundtrip_counts():
    # R-block write length = 3 header + n rows
    r = RBlock.from_lines(_r_block_lines())
    out_r = r.write()
    assert len(out_r) == 3 + r.nrows

    # TF-block write length = 3 header + n rows
    t = TFBlock.from_lines(_tf_block_lines())
    out_t = t.write()
    assert len(out_t) == 3 + t.nrows
