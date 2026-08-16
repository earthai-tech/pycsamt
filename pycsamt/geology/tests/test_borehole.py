# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.geology.borehole — Interval and Borehole."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.geology.borehole import Borehole, Interval

# ─────────────────────────────────────────────────────────────────────────────
# Interval
# ─────────────────────────────────────────────────────────────────────────────


def test_interval_thickness():
    iv = Interval(top=10.0, bottom=25.0, lithology="Clay")
    assert iv.thickness == 15.0


def test_interval_contains_top_inclusive_bottom_exclusive():
    iv = Interval(top=10.0, bottom=20.0, lithology="Clay")
    assert iv.contains(10.0) is True
    assert iv.contains(19.999) is True
    assert iv.contains(20.0) is False
    assert iv.contains(9.999) is False


def test_interval_invalid_raises():
    with pytest.raises(ValueError, match="must be > top"):
        Interval(top=20.0, bottom=20.0, lithology="Clay")
    with pytest.raises(ValueError, match="must be > top"):
        Interval(top=20.0, bottom=10.0, lithology="Clay")


def test_interval_resistivity_default_none():
    iv = Interval(top=0.0, bottom=1.0, lithology="Sand")
    assert iv.resistivity is None


# ─────────────────────────────────────────────────────────────────────────────
# Borehole — construction / accessors
# ─────────────────────────────────────────────────────────────────────────────


def _intervals():
    return [
        Interval(top=20.0, bottom=40.0, lithology="Sand", resistivity=150.0),
        Interval(top=0.0, bottom=20.0, lithology="Clay", resistivity=10.0),
        Interval(top=40.0, bottom=80.0, lithology="Granite"),  # no TRES
    ]


def test_borehole_sorts_intervals_by_top():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    tops = [iv.top for iv in bh.intervals]
    assert tops == sorted(tops)
    assert bh.intervals[0].lithology == "Clay"


def test_interval_at_depth_found_and_missing():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    assert bh.interval_at_depth(10.0).lithology == "Clay"
    assert bh.interval_at_depth(30.0).lithology == "Sand"
    assert bh.interval_at_depth(500.0) is None


def test_tres_at_depth():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    assert bh.tres_at_depth(10.0) == 10.0
    assert bh.tres_at_depth(60.0) is None  # Granite has no resistivity
    assert bh.tres_at_depth(500.0) is None  # outside any interval


def test_lithology_at_depth():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    assert bh.lithology_at_depth(30.0) == "Sand"
    assert bh.lithology_at_depth(500.0) is None


def test_tres_column():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    z = np.array([10.0, 30.0, 60.0, 500.0])
    col = bh.tres_column(z)
    np.testing.assert_allclose(col[:2], [10.0, 150.0])
    assert np.isnan(col[2])  # Granite: no TRES
    assert np.isnan(col[3])  # outside any interval


def test_min_max_depth():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    assert bh.min_depth == 0.0
    assert bh.max_depth == 80.0


def test_min_max_depth_empty_borehole():
    bh = Borehole("Empty", x=0.0, intervals=[])
    assert bh.min_depth == 0.0
    assert bh.max_depth == 0.0


def test_collar_elevation_default_and_custom():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    assert bh.collar_elevation == 0.0
    bh2 = Borehole("Bo2", x=0.0, intervals=[], collar_elevation=350.5)
    assert bh2.collar_elevation == 350.5


# ─────────────────────────────────────────────────────────────────────────────
# Borehole.from_csv
# ─────────────────────────────────────────────────────────────────────────────


def test_from_csv_basic(tmp_path):
    p = tmp_path / "Bo17.csv"
    p.write_text(
        "top,bottom,lithology,resistivity\n"
        "0,20,Clay,10.0\n"
        "20,40,Sand,150.0\n"
        "40,80,Granite,5000.0\n"
    )
    bh = Borehole.from_csv(p, x=1050.0)
    assert bh.name == "Bo17"  # defaults to file stem
    assert bh.x == 1050.0
    assert len(bh.intervals) == 3
    assert bh.tres_at_depth(10.0) == 10.0


def test_from_csv_explicit_name_and_collar(tmp_path):
    p = tmp_path / "raw.csv"
    p.write_text("top,bottom,lithology,resistivity\n0,10,Clay,10.0\n")
    bh = Borehole.from_csv(p, name="MyWell", x=5.0, collar_elevation=200.0)
    assert bh.name == "MyWell"
    assert bh.collar_elevation == 200.0


def test_from_csv_case_insensitive_headers(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("TOP,BOTTOM,LITHOLOGY,RESISTIVITY\n0,10,Clay,25.0\n")
    bh = Borehole.from_csv(p)
    assert bh.tres_at_depth(5.0) == 25.0


def test_from_csv_custom_column_names(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("z0,z1,unit,rho\n0,10,Clay,25.0\n")
    bh = Borehole.from_csv(
        p,
        top_col="z0",
        bottom_col="z1",
        lithology_col="unit",
        resistivity_col="rho",
    )
    assert bh.intervals[0].lithology == "Clay"
    assert bh.tres_at_depth(5.0) == 25.0


def test_from_csv_custom_delimiter(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("top;bottom;lithology;resistivity\n0;10;Clay;25.0\n")
    bh = Borehole.from_csv(p, delimiter=";")
    assert len(bh.intervals) == 1
    assert bh.tres_at_depth(5.0) == 25.0


def test_from_csv_missing_resistivity_column(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("top,bottom,lithology\n0,10,Clay\n")
    bh = Borehole.from_csv(p)
    assert bh.tres_at_depth(5.0) is None


@pytest.mark.parametrize("raw", ["", "nan", "None", "NA"])
def test_from_csv_blank_resistivity_sentinels(tmp_path, raw):
    p = tmp_path / "b.csv"
    p.write_text(f"top,bottom,lithology,resistivity\n0,10,Clay,{raw}\n")
    bh = Borehole.from_csv(p)
    assert bh.tres_at_depth(5.0) is None


def test_from_csv_unparseable_resistivity_becomes_none(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("top,bottom,lithology,resistivity\n0,10,Clay,not_a_number\n")
    bh = Borehole.from_csv(p)
    assert bh.tres_at_depth(5.0) is None


def test_from_csv_no_header_raises(tmp_path):
    p = tmp_path / "empty.csv"
    p.write_text("")
    with pytest.raises(ValueError, match="no header"):
        Borehole.from_csv(p)


# ─────────────────────────────────────────────────────────────────────────────
# Borehole.from_las
# ─────────────────────────────────────────────────────────────────────────────

_LAS_BASIC = """~VERSION INFORMATION
 VERS.                  2.0:   CWLS log ASCII Standard
 WRAP.                  NO:    ONE LINE PER DEPTH STEP
~WELL INFORMATION
 STRT.M          0.0:   START DEPTH
 STOP.M          4.0:   STOP DEPTH
 STEP.M          1.0:   STEP
 COMP.           pycsamt:  COMPANY
 WELL.           TestWell:  WELL
~CURVE INFORMATION
 DEPT.M                 :  DEPTH
 RESD.OHMM              :  RESISTIVITY
 LITH.                  :  LITHOLOGY CODE
~A  DEPT         RESD         LITH
   0.0000      100.00000       1
   1.0000      100.00000       1
   2.0000      -9999.25        1
   3.0000      500.00000       2
   4.0000      500.00000       2
"""


def test_from_las_basic_parses_well_name_and_groups(tmp_path):
    p = tmp_path / "well1.las"
    p.write_text(_LAS_BASIC)
    bh = Borehole.from_las(p, x=250.0)
    assert bh.name == "TestWell"
    assert bh.x == 250.0
    # two lithology groups: code 1 (depths 0-2) and code 2 (depths 3-4)
    assert len(bh.intervals) == 2
    assert bh.intervals[0].lithology == "1"
    assert bh.intervals[1].lithology == "2"


def test_from_las_null_value_masked_to_nan(tmp_path):
    p = tmp_path / "well1.las"
    p.write_text(_LAS_BASIC)
    bh = Borehole.from_las(p)
    # First interval spans the null sample (depth=2) — mean excludes it.
    assert bh.intervals[0].resistivity == pytest.approx(100.0)


def test_from_las_null_header_line_overrides_default(tmp_path):
    las = _LAS_BASIC.replace(
        " COMP.           pycsamt:  COMPANY\n",
        " COMP.           pycsamt:  COMPANY\n" " NULL.           -999:  NULL VALUE\n",
    ).replace("-9999.25", "-999.00")
    p = tmp_path / "well2.las"
    p.write_text(las)
    bh = Borehole.from_las(p)
    assert bh.intervals[0].resistivity == pytest.approx(100.0)


def test_from_las_no_lithology_curve_single_group(tmp_path):
    p = tmp_path / "well3.las"
    p.write_text(_LAS_BASIC)
    bh = Borehole.from_las(p, lithology_curve=None)
    # No lithology curve -> everything grouped into a single interval
    assert len(bh.intervals) == 1
    assert bh.intervals[0].lithology == "0"


def test_from_las_missing_resistivity_curve_all_nan(tmp_path):
    las = _LAS_BASIC.replace(" RESD.OHMM              :  RESISTIVITY\n", "")
    p = tmp_path / "well4.las"
    p.write_text(las)
    bh = Borehole.from_las(p, resistivity_curve="RESD")
    assert all(iv.resistivity is None for iv in bh.intervals)


def test_from_las_single_sample_returns_empty(tmp_path):
    las = """~WELL INFORMATION
 WELL.           OneSample:  WELL
~CURVE INFORMATION
 DEPT.M                 :  DEPTH
 RESD.OHMM              :  RESISTIVITY
~A  DEPT         RESD
   0.0000      100.00000
"""
    p = tmp_path / "well5.las"
    p.write_text(las)
    bh = Borehole.from_las(p)
    assert bh.intervals == []
    assert bh.name == "OneSample"


def test_from_las_missing_depth_curve_raises(tmp_path):
    las = """~CURVE INFORMATION
 RESD.OHMM              :  RESISTIVITY
~A  RESD
   100.00000
   100.00000
"""
    p = tmp_path / "bad.las"
    p.write_text(las)
    with pytest.raises(ValueError, match="Depth curve"):
        Borehole.from_las(p)


_LAS_EDGE = """~VERSION INFORMATION
 VERS.                  2.0:   CWLS log ASCII Standard
# a comment line, should be skipped
~WELL INFORMATION
 STRT.M          0.0:   START DEPTH
 NULL.           abc:  NULL VALUE
 WELL.           EdgeWell:  WELL
~CURVE INFORMATION
 DEPT.M                 :  DEPTH
no-dot-line-should-be-skipped
.
 RESD.OHMM              :  RESISTIVITY
~A  DEPT         RESD

   0.0000      100.00000
   1.0000      abc
"""


def test_from_las_edge_cases_comment_bad_null_bad_curve_lines(tmp_path):
    """Exercises: '#' comment skip, unparsable NULL value (except/pass),
    curve-section lines without a '.', an empty curve mnemonic, a blank
    line inside the data section, and an unparsable data value."""
    p = tmp_path / "edge.las"
    p.write_text(_LAS_EDGE)
    bh = Borehole.from_las(p, lithology_curve=None)
    assert bh.name == "EdgeWell"
    assert len(bh.intervals) == 1
    # NULL header failed to parse -> falls back to the default null_value,
    # and the unparsable "abc" data value is replaced by that same
    # sentinel, so it is masked to nan and excluded from the mean.
    assert bh.intervals[0].resistivity == pytest.approx(100.0)


def test_from_las_step_override(tmp_path):
    p = tmp_path / "well6.las"
    p.write_text(_LAS_BASIC)
    bh_default = Borehole.from_las(p)
    bh_custom = Borehole.from_las(p, step=5.0)
    # a larger step widens the bottom of the last interval
    assert bh_custom.max_depth > bh_default.max_depth


# ─────────────────────────────────────────────────────────────────────────────
# Serialisation
# ─────────────────────────────────────────────────────────────────────────────


def test_to_dataframe():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    df = bh.to_dataframe()
    assert list(df.columns) == [
        "top",
        "bottom",
        "thickness",
        "lithology",
        "resistivity",
    ]
    assert len(df) == 3


def test_to_dict():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals(), collar_elevation=12.0)
    d = bh.to_dict()
    assert d["name"] == "Bo1"
    assert d["x"] == 100.0
    assert d["collar_elevation"] == 12.0
    assert len(d["intervals"]) == 3
    assert d["intervals"][0]["lithology"] == "Clay"


def test_to_dataframe_without_pandas_raises():
    import sys
    from unittest import mock

    bh = Borehole("Bo1", x=0.0, intervals=_intervals())
    with mock.patch.dict(sys.modules, {"pandas": None}):
        with pytest.raises(ImportError, match="pandas is required"):
            bh.to_dataframe()


def test_repr():
    bh = Borehole("Bo1", x=100.0, intervals=_intervals())
    r = repr(bh)
    assert "Bo1" in r
    assert "3 intervals" in r
    assert "80.0" in r
