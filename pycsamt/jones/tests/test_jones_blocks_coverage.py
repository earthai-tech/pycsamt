from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pycsamt.jones.blocks import (
    JBlock,
    JBlocks,
    RBlock,
    TFBlock,
    _block_factory,
    _extract_all_blocks,
    _extract_first_head_and_body,
    _index_of_subsequence,
    _safe_div,
)
from pycsamt.jones.heads import Head


def _r_block_lines() -> list[str]:
    return [
        "S01",
        "RXY",
        "3",
        "-1.0  100   45   110    90     50     40     1    1",
        " 2.0  -5    30   35     25     40     20     1    1",
        " 3.0   50   10   -999   40     20     -999   -1   1",
    ]


def _tf_block_lines() -> list[str]:
    return [
        "S02",
        "ZXY SI",
        "2",
        "-1.0e+1  1.0   -2.0  0.1   1.0",
        " 5.0     -999  3.0   -999  -1.0",
    ]


# ─────────────────────────────────────────────────────────────────────────
# JBlock: abstract-base behaviors
# ─────────────────────────────────────────────────────────────────────────


def test_jblock_from_file(tmp_path: Path):
    p = tmp_path / "b.j"
    p.write_text("\n".join(_r_block_lines()), encoding="utf-8")
    blk = JBlock.from_file(p)
    assert isinstance(blk, RBlock)
    assert blk.head.station == "S01"


def test_jblock_read_and_write_raise_not_implemented():
    blk = JBlock()
    with pytest.raises(NotImplementedError):
        blk.read((None, []))
    with pytest.raises(NotImplementedError):
        blk.write()


def test_jblock_to_numpy_raises_not_implemented():
    blk = JBlock()
    with pytest.raises(NotImplementedError):
        blk.to_numpy()


def test_jblock_normalize_returns_self():
    blk = JBlock()
    assert blk.normalize() is blk


def test_jblock_qa_summary_default_empty_dict():
    blk = JBlock()
    assert blk.qa_summary() == {}


def test_jblock_properties_with_no_head():
    blk = JBlock()
    assert blk.station is None
    assert blk.kind is None
    assert blk.comp is None
    assert blk.units is None
    assert blk.has_units is False
    assert blk.columns == ()
    assert blk.shape == (0, 0)


def test_jblock_str_with_no_head():
    blk = JBlock()
    s = str(blk)
    assert "JBlock" in s and "station=None" in s


# ─────────────────────────────────────────────────────────────────────────
# RBlock: error paths and columns
# ─────────────────────────────────────────────────────────────────────────


def test_rblock_read_raises_when_data_none():
    with pytest.raises(ValueError):
        RBlock().read(None)


def test_rblock_read_skips_blank_and_stops_on_garbage():
    lines = [
        "S01",
        "RXY",
        "3",
        "",
        "-1.0  100   45   110    90     50     40     1    1",
        "not a valid data row at all",
    ]
    head, body = _extract_first_head_and_body(lines)
    blk = RBlock(head=head).read((head, body))
    assert blk.nrows == 1


def test_rblock_write_raises_when_missing_head():
    with pytest.raises(ValueError):
        RBlock().write()


def test_rblock_columns_property():
    blk = RBlock.from_lines(_r_block_lines())
    assert blk.columns == (
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
    )


def test_rblock_periods_property():
    blk = RBlock.from_lines(_r_block_lines())
    assert len(blk.periods) == 3


# ─────────────────────────────────────────────────────────────────────────
# TFBlock: error paths and columns
# ─────────────────────────────────────────────────────────────────────────


def test_tfblock_read_raises_when_data_none():
    with pytest.raises(ValueError):
        TFBlock().read(None)


def test_tfblock_read_skips_blank_and_stops_on_garbage():
    lines = [
        "S02",
        "ZXY SI",
        "2",
        "",
        "-1.0e+1  1.0   -2.0  0.1   1.0",
        "not a valid data row at all",
    ]
    head, body = _extract_first_head_and_body(lines)
    blk = TFBlock(head=head).read((head, body))
    assert blk.nrows == 1


def test_tfblock_write_raises_when_missing_head():
    with pytest.raises(ValueError):
        TFBlock().write()


def test_tfblock_columns_property():
    blk = TFBlock.from_lines(_tf_block_lines())
    assert blk.columns == ("period", "real", "imag", "error", "weight", "rej")


# ─────────────────────────────────────────────────────────────────────────
# JBlocks: read/write/exports/selection/period_range
# ─────────────────────────────────────────────────────────────────────────


def test_jblocks_read_raises_when_lines_none():
    with pytest.raises(ValueError):
        JBlocks().read(None)


def test_jblocks_read_populates_blocks():
    col = JBlocks()
    col.read(_r_block_lines())
    assert col.n == 1
    assert col.__has_read__() is True


def test_jblocks_to_numpy_list():
    col = JBlocks.from_lines(_r_block_lines() + _tf_block_lines())
    out = col.to_numpy()
    assert isinstance(out, list) and len(out) == 2


def test_jblocks_to_dataframe_empty_and_nonempty():
    pd = pytest.importorskip("pandas")
    empty_col = JBlocks()
    df_empty = empty_col.to_dataframe()
    assert isinstance(df_empty, pd.DataFrame)
    assert df_empty.empty

    col = JBlocks.from_lines(_r_block_lines())
    df = col.to_dataframe()
    assert not df.empty


def test_jblocks_stations_property():
    col = JBlocks.from_lines(_r_block_lines())
    assert col.stations == ["S01"]


def test_jblocks_station_property_empty_and_nonempty():
    assert JBlocks().station is None
    col = JBlocks.from_lines(_r_block_lines())
    assert col.station == "S01"


def test_jblocks_kinds_property():
    col = JBlocks.from_lines(_r_block_lines() + _tf_block_lines())
    assert set(col.kinds) == {"R", "Z"}


def test_jblocks_select_by_kind_and_comp():
    col = JBlocks.from_lines(_r_block_lines() + _tf_block_lines())
    by_kind = col.select(kind="Z")
    assert len(by_kind) == 1 and by_kind[0].kind == "Z"
    by_comp = col.select(comp="XY")
    assert len(by_comp) == 2
    by_neither_match = col.select(kind="Z", comp="XX")
    assert by_neither_match == []


def test_jblocks_period_range_empty_and_nonempty():
    assert JBlocks().period_range() is None
    col = JBlocks.from_lines(_r_block_lines() + _tf_block_lines())
    rng = col.period_range()
    assert rng is not None
    lo, hi = rng
    assert lo <= hi


def test_jblocks_period_range_skips_blocks_with_zero_rows():
    empty_head = Head().read(["S01", "RXY", "0"])
    empty_block = RBlock(head=empty_head, rows=[])
    col = JBlocks(blocks=[empty_block])
    assert col.period_range() is None


# ─────────────────────────────────────────────────────────────────────────
# Module-level helpers
# ─────────────────────────────────────────────────────────────────────────


def test_safe_div_normal_and_zero_division():
    assert _safe_div(1.0, 2.0) == 0.5
    assert math.isnan(_safe_div(1.0, 0.0))


def test_index_of_subsequence_empty_needle_and_no_match_and_match():
    assert _index_of_subsequence(["a", "b"], []) == -1
    assert _index_of_subsequence(["a", "b"], ["x"]) == -1
    assert _index_of_subsequence(["a", "b", "c"], ["b", "c"]) == 1


def test_block_factory_raises_for_unsupported_kind():
    from pycsamt.jones.config import KIND_COMPLEX_TF, KIND_RHO_PHI

    class _FakeDType:
        kind = "!"

    class _FakeHead:
        dtype = _FakeDType()

    assert "!" not in KIND_RHO_PHI and "!" not in KIND_COMPLEX_TF
    with pytest.raises(ValueError):
        _block_factory(_FakeHead(), verbose=0)


def test_extract_first_head_and_body_fallback_when_no_exact_match():
    lines = ["S01", "RXY", "3", "-1.0 1 1 1 1 1 1 1 1"]
    head, body = _extract_first_head_and_body(lines)
    assert head.station == "S01"
    assert body == ["-1.0 1 1 1 1 1 1 1 1"]


def test_extract_all_blocks_returns_empty_when_no_station():
    assert _extract_all_blocks(["# just a comment", ">INFO = 1"]) == []


def test_extract_all_blocks_skips_leading_noise_before_station():
    lines = ["# c", ">INFO = 1", ""] + _r_block_lines()
    blocks = _extract_all_blocks(lines)
    assert len(blocks) == 1
    assert blocks[0].station == "S01"


def test_extract_all_blocks_stops_when_malformed_header_after_dtype():
    # dtype found but the next non-blank line is neither a count nor
    # a valid dtype -> "malformed header" break, no block emitted.
    lines = ["S01", "RXY", "not_a_count_or_dtype"]
    blocks = _extract_all_blocks(lines)
    assert blocks == []


def test_extract_all_blocks_no_station_verbose_logs_warning():
    assert _extract_all_blocks(["# just a comment"], verbose=1) == []


def test_extract_all_blocks_no_count_verbose_logs_warning():
    lines = ["S01", "RXY", "not_a_count_or_dtype"]
    assert _extract_all_blocks(lines, verbose=1) == []


def test_extract_all_blocks_skips_blank_lines_while_scanning_for_dtype():
    lines = ["S01", "", "", "RXY", "3"] + _r_block_lines()[3:]
    blocks = _extract_all_blocks(lines)
    assert len(blocks) == 1 and blocks[0].nrows == 3


def test_extract_all_blocks_skips_blank_lines_while_scanning_for_count():
    lines = ["S01", "RXY", "", "", "3"] + _r_block_lines()[3:]
    blocks = _extract_all_blocks(lines)
    assert len(blocks) == 1 and blocks[0].nrows == 3


def test_extract_all_blocks_multi_block_same_station():
    lines = _r_block_lines() + ["ZXY SI", "2"] + _tf_block_lines()[3:]
    blocks = _extract_all_blocks(lines)
    assert len(blocks) == 2
    assert all(b.station == "S01" for b in blocks)
