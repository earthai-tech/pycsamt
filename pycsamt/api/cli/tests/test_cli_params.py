from __future__ import annotations

import sys
import types

import click
import pytest

from pycsamt.api.cli.params import (
    EDIDir,
    EDIPath,
    FreqRange,
    PipeStepList,
    StationList,
)


# ─────────────────────────────────────────────────────────────────────────
# EDIPath
# ─────────────────────────────────────────────────────────────────────────


def test_edi_path_accepts_valid_edi_file(tmp_path):
    edi = tmp_path / "station.edi"
    edi.write_text("dummy edi content", encoding="utf-8")
    result = EDIPath().convert(str(edi), None, None)
    assert result == edi


def test_edi_path_rejects_wrong_extension(tmp_path):
    bad = tmp_path / "station.txt"
    bad.write_text("x", encoding="utf-8")
    with pytest.raises(click.BadParameter):
        EDIPath().convert(str(bad), None, None)


def test_edi_path_rejects_missing_file(tmp_path):
    missing = tmp_path / "missing.edi"
    with pytest.raises(click.BadParameter):
        EDIPath().convert(str(missing), None, None)


# ─────────────────────────────────────────────────────────────────────────
# EDIDir
# ─────────────────────────────────────────────────────────────────────────


def test_edi_dir_accepts_directory_with_edi_files(tmp_path):
    (tmp_path / "a.edi").write_text("x", encoding="utf-8")
    result = EDIDir().convert(str(tmp_path), None, None)
    assert result == tmp_path


def test_edi_dir_rejects_directory_without_edi_files(tmp_path):
    (tmp_path / "readme.txt").write_text("x", encoding="utf-8")
    with pytest.raises(click.BadParameter):
        EDIDir().convert(str(tmp_path), None, None)


def test_edi_dir_rejects_missing_directory(tmp_path):
    missing = tmp_path / "nope"
    with pytest.raises(click.BadParameter):
        EDIDir().convert(str(missing), None, None)


# ─────────────────────────────────────────────────────────────────────────
# FreqRange
# ─────────────────────────────────────────────────────────────────────────


def test_freq_range_parses_valid_string():
    assert FreqRange().convert("0.1:1000", None, None) == (0.1, 1000.0)


def test_freq_range_passthrough_when_already_tuple():
    assert FreqRange().convert((1.0, 2.0), None, None) == (1.0, 2.0)


def test_freq_range_rejects_malformed_string():
    with pytest.raises(click.BadParameter):
        FreqRange().convert("not-a-range", None, None)


def test_freq_range_rejects_non_numeric_parts():
    with pytest.raises(click.BadParameter):
        FreqRange().convert("a:b", None, None)


def test_freq_range_rejects_inverted_bounds():
    with pytest.raises(click.BadParameter):
        FreqRange().convert("1000:0.1", None, None)


def test_freq_range_rejects_equal_bounds():
    with pytest.raises(click.BadParameter):
        FreqRange().convert("5:5", None, None)


# ─────────────────────────────────────────────────────────────────────────
# StationList
# ─────────────────────────────────────────────────────────────────────────


def test_station_list_parses_comma_separated():
    assert StationList().convert("S01,S02,S03", None, None) == [
        "S01", "S02", "S03",
    ]


def test_station_list_strips_whitespace_and_drops_empties():
    assert StationList().convert(" S01 , , S02 ", None, None) == [
        "S01", "S02",
    ]


def test_station_list_passthrough_when_already_list():
    assert StationList().convert(["S01"], None, None) == ["S01"]


def test_station_list_rejects_empty_result():
    with pytest.raises(click.BadParameter):
        StationList().convert("   ,  ,", None, None)


# ─────────────────────────────────────────────────────────────────────────
# PipeStepList
# ─────────────────────────────────────────────────────────────────────────


def test_pipe_step_list_passthrough_when_already_list():
    assert PipeStepList().convert(["NR001"], None, None) == ["NR001"]


def test_pipe_step_list_rejects_empty_string():
    with pytest.raises(click.BadParameter):
        PipeStepList().convert(" , ,", None, None)


def test_pipe_step_list_resolves_known_codes_and_aliases():
    result = PipeStepList().convert("NR001,notch_powerline", None, None)
    assert result == ["NR001", "NR001"]


def test_pipe_step_list_rejects_unknown_step():
    with pytest.raises(click.BadParameter):
        PipeStepList().convert("bogus_step_xyz", None, None)


def test_pipe_step_list_falls_back_when_pipeline_not_importable(monkeypatch):
    fake_module = types.ModuleType("pycsamt.pipeline")
    monkeypatch.setitem(sys.modules, "pycsamt.pipeline", fake_module)
    result = PipeStepList().convert("anything,goes", None, None)
    assert result == ["anything", "goes"]
