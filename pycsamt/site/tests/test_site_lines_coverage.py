"""Branch coverage for survey-line helpers."""

from pycsamt.site.lines import (
    apply_line_renames,
    detect_lines_from_filenames,
    detect_lines_from_station_ids,
    pick_representative_stations,
    resolve_display_name,
)


def test_detect_station_lines_numeric_text_empty_and_trivial():
    assert detect_lines_from_station_ids(
        ["18-001", "18_002", "west-1", "west-2", "-orphan"]
    ) == {
        "L18": ["18-001", "18_002"],
        "WEST": ["west-1", "west-2"],
        "unassigned": ["-orphan"],
    }
    assert detect_lines_from_station_ids(["A", "B"]) == {
        "auto": ["A", "B"]
    }


def test_detect_station_lines_merges_small_groups():
    groups = detect_lines_from_station_ids(
        ["A-1", "A-2", "B-1"], min_stations_per_line=2
    )
    assert groups == {"A": ["A-1", "A-2"], "unassigned": ["B-1"]}


def test_detect_filename_lines_by_stem_and_parent(tmp_path):
    paths = [
        str(tmp_path / "line-a" / "18-001.edi"),
        str(tmp_path / "line-a" / "18_002.edi"),
        str(tmp_path / "line-b" / "west-1.edi"),
        str(tmp_path / "line-b" / "-.edi"),
    ]
    assert set(detect_lines_from_filenames(paths)) == {
        "L18",
        "WEST",
        "unassigned",
    }
    assert set(detect_lines_from_filenames(paths, use_parent=True)) == {
        "line-a",
        "line-b",
    }
    assert detect_lines_from_filenames(["station.edi"], use_parent=True) == {
        "unassigned": ["station.edi"]
    }


def test_rename_resolve_and_representative_selection():
    original = [{"ID": "S1", "Line": "L1"}, {"ID": "S2"}]
    renamed = apply_line_renames(original, {"L1": "South", "": ""})
    assert renamed == [{"ID": "S1", "Line": "South"}, {"ID": "S2"}]
    assert renamed is not original
    assert resolve_display_name("L1", {"L1": "South"}) == "South"
    assert resolve_display_name("L2", None) == "L2"
    assert pick_representative_stations(["S3", "S1", "S2"], 5) == [
        "S1",
        "S2",
        "S3",
    ]
    assert pick_representative_stations(["S4", "S1", "S3", "S2"], 2) == [
        "S1",
        "S4",
    ]
    assert pick_representative_stations(["S2", "S1"], 1) == ["S1"]
