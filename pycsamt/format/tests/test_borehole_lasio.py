"""Tests for loss-explicit PCBH LAS subset adapters."""

from __future__ import annotations

from pycsamt.format.borehole import (
    Collar,
    borehole_from_las,
    write_las_subset,
)

LAS_TEXT = """~VERSION INFORMATION
 VERS. 2.0: CWLS LOG ASCII STANDARD
 WRAP. NO: ONE LINE PER DEPTH STEP
~WELL INFORMATION
 STRT.M 0: START
 STOP.M 2: STOP
 STEP.M 1: STEP
 NULL. -999.25: NULL
 WELL. TEST-1: WELL
~CURVE INFORMATION
 DEPT.M: MEASURED DEPTH
 RESD.OHMM: DEEP RESISTIVITY
 LITH.CODE: LITHOLOGY CODE
 GR.API: GAMMA RAY
~A DEPT RESD LITH GR
0 10 1 40
1 -999.25 1 45
2 100 2 50
"""


def test_las_import_preserves_curves_units_nulls_and_metadata(tmp_path):
    source = tmp_path / "source.las"
    source.write_text(LAS_TEXT, encoding="utf-8")
    document, report = borehole_from_las(
        source,
        collar=Collar(500.0, 600.0, 70.0),
        crs_horizontal="LOCAL:test",
    )

    hole = document.boreholes[0]
    curves = hole.extensions["pcbh:continuous_curves"]
    resd = next(item for item in curves if item["mnemonic"] == "RESD")
    assert resd["unit"] == "OHMM"
    assert resd["samples"][1] == [1.0, None]
    assert hole.extensions["pcbh:las_metadata"]["null_value"] == -999.25
    assert len(hole.interval_logs["lithology"]) == 2
    assert report.rows_accepted == 3


def test_las_subset_round_trip_and_loss_report(tmp_path):
    source = tmp_path / "source.las"
    source.write_text(LAS_TEXT, encoding="utf-8")
    document, _ = borehole_from_las(
        source,
        collar=Collar(0.0, 0.0, 0.0),
        crs_horizontal="LOCAL:test",
    )
    output, loss = write_las_subset(
        document.boreholes[0], tmp_path / "roundtrip.las"
    )
    restored, _ = borehole_from_las(
        output,
        collar=Collar(0.0, 0.0, 0.0),
        crs_horizontal="LOCAL:test",
    )

    names = {
        item["mnemonic"]
        for item in restored.boreholes[0].extensions[
            "pcbh:continuous_curves"
        ]
    }
    assert names == {"RESD", "LITH", "GR"}
    assert set(loss.curves_written) == {"DEPT", "RESD", "LITH", "GR"}
    assert loss.losses


def test_las_rejects_broken_width_and_sample_limit(tmp_path):
    broken = tmp_path / "broken.las"
    broken.write_text(LAS_TEXT + "3 20\n", encoding="utf-8")
    try:
        borehole_from_las(
            broken,
            collar=Collar(0.0, 0.0, 0.0),
            crs_horizontal="LOCAL:test",
        )
    except ValueError as error:
        assert "width" in str(error)
    else:
        raise AssertionError("broken LAS width was accepted")

    source = tmp_path / "limited.las"
    source.write_text(LAS_TEXT, encoding="utf-8")
    try:
        borehole_from_las(
            source,
            collar=Collar(0.0, 0.0, 0.0),
            crs_horizontal="LOCAL:test",
            max_samples=2,
        )
    except ValueError as error:
        assert "sample limit" in str(error)
    else:
        raise AssertionError("LAS sample limit was ignored")


def test_las_feet_depth_is_converted_and_reported(tmp_path):
    source = tmp_path / "feet.las"
    source.write_text(
        LAS_TEXT.replace("DEPT.M", "DEPT.FT").replace("STEP.M", "STEP.FT"),
        encoding="utf-8",
    )
    document, report = borehole_from_las(
        source,
        collar=Collar(0.0, 0.0, 0.0),
        crs_horizontal="LOCAL:test",
    )

    assert document.boreholes[0].total_depth_md == 3 * 0.3048
    assert report.unit_conversions == ["LAS depth FT converted to PCBH metres"]
