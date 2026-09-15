# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import AvgDataError, AvgFileError, StationError
from pycsamt.zonge.utils import (
    _block_to_dict,
    _dict_to_lines,
    _find_col,
    _first_present,
    _get_weight_bool,
    _next_block,
    _norm_comp,
    _parse_kind1,
    _parse_kind2,
    _to_complex,
    _to_float,
    _to_num,
    _to_numeric_percent,
    chunk_by_frequency,
    classify_avg_format,
    detect_stn_header,
    extract_core_columns,
    find_and_rename_column,
    load_avg,
    number_stations,
    read_stn,
    round_dipole_length,
    split_by_station,
    to_numeric_if_possible,
    to_xarray,
    write_avg,
)


# ─────────────────────────────────────────────────────────────────────────
# _to_float
# ─────────────────────────────────────────────────────────────────────────


def test_to_float_passthrough_for_numeric_types():
    assert _to_float(5) == 5.0
    assert _to_float(5.5) == 5.5


def test_to_float_leading_dot():
    assert _to_float(".5") == 0.5


def test_to_float_trailing_dot():
    assert _to_float("5.") == 5.0


def test_to_float_unparseable_returns_nan():
    assert np.isnan(_to_float("abc"))


def test_to_float_sentinel_returns_nan():
    assert np.isnan(_to_float("*"))


# ─────────────────────────────────────────────────────────────────────────
# classify_avg_format: fallback heuristics
# ─────────────────────────────────────────────────────────────────────────


def test_classify_avg_format_dot_notation_keyword_fallback():
    lines = ["$Survey.Type=CSAMT", "1,2,3"]
    assert classify_avg_format(lines) == 2


def test_classify_avg_format_comma_heavy_data_fallback():
    lines = ["some text", "1,2,3,4,5"]
    assert classify_avg_format(lines) == 2


def test_classify_avg_format_raises_when_no_indicators():
    lines = ["just some text", "no clear format here"]
    with pytest.raises(AvgFileError):
        classify_avg_format(lines)


def test_classify_avg_format_skips_blank_lines():
    lines = ["", "   ", "$Survey.Type=CSAMT"]
    assert classify_avg_format(lines) == 2


# ─────────────────────────────────────────────────────────────────────────
# _parse_kind1
# ─────────────────────────────────────────────────────────────────────────


def test_parse_kind1_raises_when_no_header_found():
    with pytest.raises(AvgFileError, match="Header row not found"):
        _parse_kind1(["just some text", "no header here"])


def test_parse_kind1_raises_when_no_data_rows():
    lines = [
        "skp Station Freq Comp Amps Emag Ephz Hmag Hphz "
        "Resistivity Phase %Emag sEphz %Hmag sHphz %Rho sPhz",
        "",  # blank immediately stops row collection
    ]
    with pytest.raises(AvgDataError, match="No data rows"):
        _parse_kind1(lines)


def test_parse_kind1_stops_at_blank_line_and_skips_comment_rows():
    header = (
        "skp Station Freq Comp Amps Emag Ephz Hmag Hphz "
        "Resistivity Phase %Emag sEphz %Hmag sHphz %Rho sPhz"
    )
    row = (
        " 2   150.0   8192 ExHy  5.00  3.1e+2  1371.6  9.2e-2  "
        "1953.2  2.7e+2  -581.6   13.8   84.7    9.8   73.6   14.7  136.0"
    )
    lines = [header, r"\ a comment row to skip", row, "", "extra ignored"]
    df = _parse_kind1(lines)
    assert len(df) == 1


def test_parse_kind1_appends_generic_columns_for_surplus_fields():
    header = "skp Station Freq Comp Amps"
    row1 = "2 150.0 8192 ExHy 5.00 99.9"  # one extra field beyond header
    row2 = "2 200.0 4096 ExHy 6.00 88.8"
    df = _parse_kind1([header, row1, row2])
    assert "extra_1" in df.columns


# ─────────────────────────────────────────────────────────────────────────
# _next_block
# ─────────────────────────────────────────────────────────────────────────


def test_next_block_returns_none_when_no_header_found():
    start, end = _next_block(["no header here", "still nothing"], 0)
    assert start is None


def test_next_block_advances_past_non_header_lines():
    lines = ["not a header", "Z.mwgt,Z.pwgt,Freq", "1,1,1", ""]
    start, end = _next_block(lines, 0)
    assert start == 1


def test_next_block_ends_on_dollar_line_without_blank():
    lines = ["Z.mwgt,Z.pwgt,Freq", "1,1,1", "$Rx.Stn=5"]
    start, end = _next_block(lines, 0)
    assert start == 0
    assert end == 2  # stops before the $ line, no trailing blank needed


# ─────────────────────────────────────────────────────────────────────────
# _parse_kind2 edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_parse_kind2_raises_when_no_blocks_found():
    with pytest.raises(AvgDataError, match="Data block"):
        _parse_kind2(["$Survey.Type=CSAMT", "no table here at all"])


def test_parse_kind2_skips_standalone_comment_lines():
    lines = [
        "$Survey.Type=CSAMT",
        r"\ a standalone comment between meta and table",
        "$Rx.Stn=10",
        "Z.mwgt,Z.pwgt,Freq",
        "1,1,1",
        "",
    ]
    df, meta = _parse_kind2(lines)
    assert len(df) == 1
    assert meta["Survey.Type"] == "CSAMT"


def test_parse_kind2_station_fallback_when_to_float_raises(monkeypatch):
    import pycsamt.zonge.utils as u

    real_to_float = u._to_float

    def _boom(val):
        if val == "BOOMVAL":
            raise RuntimeError("simulated failure")
        return real_to_float(val)

    monkeypatch.setattr(u, "_to_float", _boom)
    lines = [
        "$Rx.Stn=BOOMVAL",
        "Z.mwgt,Z.pwgt,Freq",
        "1,1,1",
        "",
    ]
    df, meta = u._parse_kind2(lines)
    assert df["station"].iloc[0] == "BOOMVAL"


# ─────────────────────────────────────────────────────────────────────────
# to_numeric_if_possible
# ─────────────────────────────────────────────────────────────────────────


def test_to_numeric_if_possible_converts_clean_numeric_strings():
    out = to_numeric_if_possible(pd.Series(["1", "2", "3"]))
    assert pd.api.types.is_numeric_dtype(out)


def test_to_numeric_if_possible_keeps_original_when_conversion_creates_new_na():
    out = to_numeric_if_possible(pd.Series(["1", "not_a_number", "3"]))
    assert list(out) == ["1", "not_a_number", "3"]


def test_to_numeric_if_possible_returns_input_on_exception():
    class _Bad:
        def __iter__(self):
            raise RuntimeError("cannot iterate")

    bad = _Bad()
    assert to_numeric_if_possible(bad) is bad


# ─────────────────────────────────────────────────────────────────────────
# split_by_station
# ─────────────────────────────────────────────────────────────────────────


def test_split_by_station_raises_when_column_missing():
    with pytest.raises(AvgDataError, match="'station' column missing"):
        split_by_station(pd.DataFrame({"freq": [1.0]}))


def test_split_by_station_groups_and_uses_int_keys_for_integral_floats():
    df = pd.DataFrame({"station": [25.0, 25.0, 75.0], "freq": [1.0, 2.0, 1.0]})
    out = split_by_station(df)
    assert set(out.keys()) == {25, 75}
    assert len(out[25]) == 2


def test_split_by_station_coerces_non_numeric_station_column():
    df = pd.DataFrame({"station": ["25", "25", "75"], "freq": [1.0, 2.0, 1.0]})
    out = split_by_station(df)
    assert set(out.keys()) == {25, 75}


def test_split_by_station_handles_nan_station_group():
    df = pd.DataFrame({"station": [25.0, np.nan], "freq": [1.0, 2.0]})
    out = split_by_station(df)
    assert len(out) == 2
    nan_keys = [k for k in out if isinstance(k, float) and np.isnan(k)]
    assert len(nan_keys) == 1


# ─────────────────────────────────────────────────────────────────────────
# to_xarray edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_to_xarray_raises_when_no_coordinate_columns_found():
    pytest.importorskip("xarray")
    df = pd.DataFrame({"rho": [1.0, 2.0]})
    with pytest.raises(AvgDataError, match="No coordinate columns"):
        to_xarray(df, coords=())


def test_to_xarray_default_data_vars_and_raises_when_none_found():
    pytest.importorskip("xarray")
    df = pd.DataFrame({"station": [1.0], "freq": [1.0]})
    with pytest.raises(AvgDataError, match="No data variables"):
        to_xarray(df)


def test_to_xarray_default_data_vars_selected_when_numeric_present():
    pytest.importorskip("xarray")
    df = pd.DataFrame(
        {"station": [1.0, 2.0], "freq": [1.0, 1.0], "rho": [10.0, 20.0]}
    )
    ds = to_xarray(df)
    assert "rho" in ds.data_vars


# ─────────────────────────────────────────────────────────────────────────
# write_avg edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_write_avg_drops_fully_empty_extra_columns(tmp_path: Path):
    core = pd.DataFrame(
        {
            "station": [1, 1],
            "freq": [1.0, 2.0],
            "rho": [10.0, 20.0],
            "coh": [np.nan, np.nan],
        }
    )
    out = write_avg(core, None, {}, tmp_path / "out.avg", stamp=False)
    txt = out.read_text(encoding="utf-8")
    assert "Choer" not in txt


def test_write_avg_skips_nan_station_group(tmp_path: Path):
    core = pd.DataFrame(
        {"station": [1.0, np.nan], "freq": [1.0, 2.0], "rho": [10.0, 20.0]}
    )
    out = write_avg(core, None, {}, tmp_path / "out2.avg", stamp=False)
    txt = out.read_text(encoding="utf-8")
    assert "$Rx.Stn=1" in txt


def test_write_avg_single_block_when_no_station_column(tmp_path: Path):
    core = pd.DataFrame({"freq": [1.0, 2.0], "rho": [10.0, 20.0]})
    out = write_avg(core, None, {}, tmp_path / "out3.avg", stamp=False)
    txt = out.read_text(encoding="utf-8")
    assert "$Rx.Stn" not in txt
    assert "ARes.mag" in txt or "rho" in txt.lower()


def test_write_avg_stamps_rx_gdpstn_and_rx_length(tmp_path: Path):
    core = pd.DataFrame(
        {
            "station": [1, 1],
            "freq": [1.0, 2.0],
            "rho": [10.0, 20.0],
            "rx_gdpstn": [99, 99],
            "rx_length": [50, 50],
        }
    )
    out = write_avg(core, None, {}, tmp_path / "out4.avg", stamp=False)
    txt = out.read_text(encoding="utf-8")
    assert "$Rx.GdpStn=99" in txt
    assert "$Rx.Length=50" in txt


# ─────────────────────────────────────────────────────────────────────────
# load_avg edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_load_avg_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        load_avg("no_such_file_at_all.avg")


def test_load_avg_computes_utm_when_lat_lon_present(tmp_path: Path, monkeypatch):
    import pycsamt.zonge.utils as u

    monkeypatch.setattr(
        u,
        "to_utm",
        lambda lat, lon, utm_zone=None: (
            np.array([500000.0]),
            np.array([4000000.0]),
            "12T",
        ),
    )
    text = "\n".join(
        [
            "$Survey.Type=CSAMT",
            "Z.mwgt,Z.pwgt,Freq,latitude,longitude",
            "1,1,1,40.0,-110.0",
            "",
        ]
    )
    p = tmp_path / "ll.avg"
    p.write_text(text, encoding="utf-8")
    df, meta = u.load_avg(p)
    assert "easting" in df.columns
    assert "northing" in df.columns
    assert np.allclose(df["easting"], 500000.0)


def test_load_avg_utm_conversion_failure_is_caught_and_logged(
    tmp_path: Path, monkeypatch
):
    import pycsamt.zonge.utils as u

    def _boom(lat, lon, utm_zone=None):
        raise RuntimeError("simulated projection failure")

    monkeypatch.setattr(u, "to_utm", _boom)
    text = "\n".join(
        [
            "$Survey.Type=CSAMT",
            "Z.mwgt,Z.pwgt,Freq,latitude,longitude",
            "1,1,1,40.0,-110.0",
            "",
        ]
    )
    p = tmp_path / "ll_fail.avg"
    p.write_text(text, encoding="utf-8")
    df, meta = u.load_avg(p)
    assert "easting" not in df.columns  # conversion failed, no crash


# ─────────────────────────────────────────────────────────────────────────
# read_stn edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_read_stn_raises_when_file_missing():
    with pytest.raises(StationError, match="Cannot read STN file"):
        read_stn("no_such_stn_file.stn")


def test_read_stn_raises_when_empty_or_comment_only(tmp_path: Path):
    p = tmp_path / "empty.stn"
    p.write_text("! just a comment\n; another\n", encoding="utf-8")
    with pytest.raises(StationError, match="Empty or comment-only"):
        read_stn(p)


def test_read_stn_splits_header_and_data_on_same_line(tmp_path: Path):
    p = tmp_path / "inline.stn"
    p.write_text(
        "station,easting,northing,elevation,roll 100,500000,4000000,10,0.5\n"
        "200,500100,4000010,11,0.6\n",
        encoding="utf-8",
    )
    df = read_stn(p)
    assert len(df) == 2
    assert "roll" in [c.lower() for c in df.columns]


def test_read_stn_no_header_assumes_four_legacy_columns(tmp_path: Path):
    p = tmp_path / "legacy.stn"
    p.write_text("100 500000 4000000 10\n200 500100 4000010 11\n", encoding="utf-8")
    df = read_stn(p)
    assert list(df.columns) == ["station", "easting", "northing", "elevation"]


def test_read_stn_raises_when_header_found_but_no_data(tmp_path: Path):
    p = tmp_path / "headeronly.stn"
    p.write_text("station,easting,northing,elevation\n", encoding="utf-8")
    with pytest.raises(StationError, match="No data rows"):
        read_stn(p)


def test_read_stn_raises_when_all_rows_empty(tmp_path: Path):
    p = tmp_path / "allnan.stn"
    p.write_text("station,easting,northing,elevation\n,,,\n", encoding="utf-8")
    with pytest.raises(StationError, match="empty DataFrame"):
        read_stn(p)


def test_read_stn_whitespace_header_and_data(tmp_path: Path):
    p = tmp_path / "ws.stn"
    p.write_text(
        "station easting northing elevation\n100 500000 4000000 10\n",
        encoding="utf-8",
    )
    df = read_stn(p)
    assert len(df) == 1


# ─────────────────────────────────────────────────────────────────────────
# detect_stn_header
# ─────────────────────────────────────────────────────────────────────────


def test_detect_stn_header_no_match_returns_zero_score():
    score, matches = detect_stn_header(["! comment", "random text no labels"])
    assert score == 0
    assert matches == []


def test_detect_stn_header_comma_header_matches_labels():
    score, matches = detect_stn_header(["station,easting,northing,elevation"])
    assert score == 4
    labels = [m[0] for m in matches]
    assert labels == ["station", "easting", "northing", "elevation"]


def test_detect_stn_header_whitespace_header_matches_labels():
    score, matches = detect_stn_header(["dot e n h heading pitch roll"])
    assert score == 7


def test_detect_stn_header_skips_comments():
    score, matches = detect_stn_header(
        ["! ignore this station easting", "station,easting,northing"]
    )
    assert score == 3


def test_detect_stn_header_header_and_data_same_line():
    score, matches = detect_stn_header(["station,easting,northing 100,500000,4000000"])
    assert score >= 2


def test_detect_stn_header_picks_best_scoring_line():
    lines = [
        "station",  # score 1
        "station,easting,northing,elevation",  # score 4 (best)
    ]
    score, matches = detect_stn_header(lines)
    assert score == 4


def test_detect_stn_header_explicit_splitter():
    score, matches = detect_stn_header(["station;easting;northing"], splitter=";")
    assert score == 3


def test_detect_stn_header_label_equals_value_stripped():
    # The embedded "label=value" segment is stripped down to its label
    # before tokenizing, so "line=1" still contributes the "line" token.
    score, matches = detect_stn_header(["station easting line=1"])
    assert score == 3
    labels = [m[0] for m in matches]
    assert "line" in labels


# ─────────────────────────────────────────────────────────────────────────
# round_dipole_length
# ─────────────────────────────────────────────────────────────────────────


def test_round_dipole_length_low_remainder_rounds_down():
    assert round_dipole_length(52.0) == 50.0


def test_round_dipole_length_mid_remainder_rounds_up_to_next_five():
    assert round_dipole_length(54.0) == 55.0


def test_round_dipole_length_high_remainder_rounds_up_ten():
    assert round_dipole_length(58.0) == 60.0


# ─────────────────────────────────────────────────────────────────────────
# extract_core_columns
# ─────────────────────────────────────────────────────────────────────────


def test_extract_core_columns_default_keep():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [1.0],
            "rho": [10.0],
            "unrelated_extra": [99],
        }
    )
    core, extra = extract_core_columns(df)
    assert "rho" in core.columns
    assert "unrelated_extra" in extra.columns


def test_extract_core_columns_custom_keep_inserts_station():
    df = pd.DataFrame({"station": [1], "rho": [10.0], "other": [5]})
    core, extra = extract_core_columns(df, keep=["rho"])
    assert "station" in core.columns
    assert "rho" in core.columns
    assert "other" in extra.columns


def test_extract_core_columns_uses_first_column_when_no_station():
    df = pd.DataFrame({"freq": [1.0], "rho": [10.0]})
    core, extra = extract_core_columns(df, keep=["rho"])
    assert "freq" in core.columns  # fallback: first column stands in


# ─────────────────────────────────────────────────────────────────────────
# _block_to_dict / _dict_to_lines
# ─────────────────────────────────────────────────────────────────────────


def test_block_to_dict_parses_key_value_lines_and_skips_others():
    out = _block_to_dict(["Key1=Val1", "no equals here", "KEY2 = Val2"])
    assert out == {"key1": "Val1", "key2": "Val2"}


def test_dict_to_lines_from_dict_skips_none_values():
    lines = _dict_to_lines({"a": 1, "b": None, "c": "x"})
    assert lines == ["a=1\n", "c=x\n"]


def test_dict_to_lines_from_json_string():
    lines = _dict_to_lines('{"a": 1}')
    assert lines == ["a=1\n"]


def test_dict_to_lines_from_arbitrary_mapping_like():
    lines = _dict_to_lines([("a", 1), ("b", 2)])
    assert lines == ["a=1\n", "b=2\n"]


# ─────────────────────────────────────────────────────────────────────────
# number_stations / chunk_by_frequency
# ─────────────────────────────────────────────────────────────────────────


def test_number_stations_raises_for_invalid_counts():
    with pytest.raises(ValueError):
        number_stations(0, 1)
    with pytest.raises(ValueError):
        number_stations(1, 0)


def test_chunk_by_frequency_basic():
    out = chunk_by_frequency([0, 1, 2, 3, 4], 2)
    assert [list(c) for c in out] == [[0, 1], [2, 3], [4]]


def test_chunk_by_frequency_drop_remainder():
    out = chunk_by_frequency([0, 1, 2, 3, 4], 2, drop_remainder=True)
    assert [list(c) for c in out] == [[0, 1], [2, 3]]


def test_chunk_by_frequency_raises_for_invalid_n_per_chunk():
    with pytest.raises(ValueError):
        chunk_by_frequency([1, 2, 3], 0)


# ─────────────────────────────────────────────────────────────────────────
# _find_col / _to_num / _norm_comp / _first_present / _to_numeric_percent /
# _to_complex
# ─────────────────────────────────────────────────────────────────────────


def test_find_col_case_insensitive_match_and_none():
    df = pd.DataFrame({"Rho": [1]})
    assert _find_col(df, ["rho"]) == "Rho"
    assert _find_col(df, ["nope"]) is None


def test_to_num_none_and_sentinels_and_valid_and_invalid():
    assert np.isnan(_to_num(None))
    assert np.isnan(_to_num("*"))
    assert _to_num("3.5") == 3.5
    assert np.isnan(_to_num("garbage"))


def test_norm_comp_none_and_blank_return_default():
    assert _norm_comp(None) == "ExHy"
    assert _norm_comp("   ") == "ExHy"


def test_norm_comp_zxx_family_capitalized():
    assert _norm_comp("zxy") == "Zxy"
    assert _norm_comp("ZYX") == "Zyx"


def test_norm_comp_fallback_capitalizes_first_letter_only():
    assert _norm_comp("weird_label") == "Weird_label"


def test_first_present_matches_and_none():
    df = pd.DataFrame({"Freq": [1.0]})
    assert _first_present(df, ["freq", "frequency"]) == "Freq"
    assert _first_present(df, ["nope"]) is None


def test_to_numeric_percent_handles_blank_and_asterisk():
    out = _to_numeric_percent(pd.Series(["5", "", "*", "10.5"]))
    assert out.iloc[0] == 5.0
    assert pd.isna(out.iloc[1])
    assert pd.isna(out.iloc[2])
    assert out.iloc[3] == 10.5


def test_to_complex_none_and_sentinel_and_valid_and_invalid():
    assert np.isnan(_to_complex(None))
    assert np.isnan(_to_complex("*"))
    assert _to_complex("1+2j") == complex(1, 2)
    assert np.isnan(_to_complex("not-complex-at-all!"))


def test_norm_comp_exhy_exhx_eyhx_eyhy_family():
    assert _norm_comp("exhy") == "ExHy"
    assert _norm_comp("EX-HY") == "ExHy"
    assert _norm_comp("exhx") == "ExHx"
    assert _norm_comp("eyhx") == "EyHx"
    assert _norm_comp("eyhy") == "EyHy"


# ─────────────────────────────────────────────────────────────────────────
# find_and_rename_column
# ─────────────────────────────────────────────────────────────────────────


def test_find_and_rename_column_renames_known_alias():
    df = pd.DataFrame({"Resistivity": [1.0]})
    out = find_and_rename_column(df, "rho")
    assert "rho" in out.columns


def test_find_and_rename_column_noop_when_no_alias_present():
    df = pd.DataFrame({"unrelated": [1.0]})
    out = find_and_rename_column(df, "rho")
    assert list(out.columns) == ["unrelated"]


def test_find_and_rename_column_noop_when_already_canonical():
    df = pd.DataFrame({"rho": [1.0]})
    out = find_and_rename_column(df, "rho")
    assert list(out.columns) == ["rho"]


# ─────────────────────────────────────────────────────────────────────────
# _get_weight_bool
# ─────────────────────────────────────────────────────────────────────────


def test_get_weight_bool_default_when_column_missing():
    df = pd.DataFrame({"other": [1.0]})
    out = _get_weight_bool(df, "z_mwgt")
    assert out is True  # default value 1 > 0


def test_get_weight_bool_from_column_with_nan_fill():
    df = pd.DataFrame({"z_mwgt": [1.0, np.nan, 0.0]})
    out = _get_weight_bool(df, "z_mwgt")
    assert list(out) == [True, True, False]


# ─────────────────────────────────────────────────────────────────────────
# to_numeric_if_possible: second exception branch
# ─────────────────────────────────────────────────────────────────────────


def test_to_numeric_if_possible_returns_input_when_isna_check_fails(monkeypatch):
    import pycsamt.zonge.utils as u

    original_isna = u.pd.isna
    calls = {"n": 0}

    def _flaky_isna(x):
        calls["n"] += 1
        if calls["n"] == 2:
            raise RuntimeError("simulated isna failure")
        return original_isna(x)

    monkeypatch.setattr(u.pd, "isna", _flaky_isna)
    values = pd.Series(["1", "2"])
    assert u.to_numeric_if_possible(values) is values


# ─────────────────────────────────────────────────────────────────────────
# to_xarray: duplicate averaging + attrs
# ─────────────────────────────────────────────────────────────────────────


def test_to_xarray_averages_duplicate_coordinate_rows():
    pytest.importorskip("xarray")
    df = pd.DataFrame(
        {
            "station": [1.0, 1.0],
            "freq": [1.0, 1.0],
            "comp": ["ExHy", "ExHy"],
            "rho": [10.0, 20.0],
        }
    )
    ds = to_xarray(df)
    val = ds["rho"].sel(station=1.0, freq=1.0, comp="ExHy").item()
    assert np.isclose(val, 15.0)


def test_to_xarray_merges_attrs_and_drops_blocks_key():
    pytest.importorskip("xarray")
    df = pd.DataFrame(
        {"station": [1.0], "freq": [1.0], "comp": ["ExHy"], "rho": [10.0]}
    )
    ds = to_xarray(df, attrs={"blocks": [1, 2], "Survey.Type": "CSAMT"})
    assert "blocks" not in ds.attrs
    assert ds.attrs.get("Survey.Type") == "CSAMT"


# ─────────────────────────────────────────────────────────────────────────
# write_avg: banner lines / stamp / comp stamping
# ─────────────────────────────────────────────────────────────────────────


def test_write_avg_with_banner_lines_and_stamp(tmp_path: Path):
    core = pd.DataFrame(
        {
            "station": [1],
            "freq": [1.0],
            "rho": [10.0],
            "comp": ["ExHy"],
        }
    )
    out = write_avg(
        core,
        None,
        {},
        tmp_path / "banner.avg",
        stamp=True,
        banner_lines=[r"\ AMTAVG banner line"],
    )
    txt = out.read_text(encoding="utf-8")
    assert "AMTAVG banner line" in txt
    assert "$Written=" in txt
    assert "$Rx.Cmp=ExHy" in txt


def test_write_avg_header_spaces_option(tmp_path: Path):
    core = pd.DataFrame({"station": [1], "freq": [1.0], "rho": [10.0]})
    out = write_avg(
        core, None, {"Survey.Type": "CSAMT"}, tmp_path / "spaced.avg",
        stamp=False, header_spaces=True,
    )
    txt = out.read_text(encoding="utf-8")
    assert "$Survey.Type = CSAMT" in txt


# ─────────────────────────────────────────────────────────────────────────
# read_stn: comma-only-in-data delimiter detection
# ─────────────────────────────────────────────────────────────────────────


def test_read_stn_detects_comma_from_first_data_line(tmp_path: Path):
    p = tmp_path / "mixed.stn"
    p.write_text(
        "station easting northing elevation\n100,500000,4000000,10\n",
        encoding="utf-8",
    )
    df = read_stn(p)
    assert len(df) == 1


def test_read_stn_wraps_parse_exception(tmp_path: Path, monkeypatch):
    import pycsamt.zonge.utils as u

    def _boom(*a, **k):
        raise ValueError("simulated tokenizer failure")

    monkeypatch.setattr(u.pd, "read_csv", _boom)
    p = tmp_path / "any.stn"
    p.write_text("station,easting,northing,elevation\n1,2,3,4\n", encoding="utf-8")
    with pytest.raises(StationError, match="Failed to parse STN data"):
        u.read_stn(p)


# ─────────────────────────────────────────────────────────────────────────
# detect_stn_header: empty-tokens continue
# ─────────────────────────────────────────────────────────────────────────


def test_detect_stn_header_skips_lines_with_no_letter_tokens():
    score, matches = detect_stn_header(["12345", "station,easting"])
    assert score == 2


# ─────────────────────────────────────────────────────────────────────────
# number_stations: happy path
# ─────────────────────────────────────────────────────────────────────────


def test_number_stations_happy_path():
    ids, expanded = number_stations(2, 3, prefix="S")
    assert ids == ["S00", "S01"]
    assert expanded == ["S00", "S00", "S00", "S01", "S01", "S01"]
