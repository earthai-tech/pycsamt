"""Edge-case coverage for dependency-light site helper contracts."""

from types import SimpleNamespace

import numpy as np
import pandas as pd

import pycsamt.site.location as location
import pycsamt.site.recompute as recompute
import pycsamt.site.selection as selection
import pycsamt.site.utils as utils
from pycsamt.site.recompute import EDIRecomputeRecord


def test_location_dataframe_and_coordinate_helpers():
    df = pd.DataFrame(
        [{"Station": " s01 ", "LAT": "5:30:00N", "LON": "3:15:00W"}]
    )
    assert location._lower_cols(df) == ["Station", "LAT", "LON"]
    assert location._lower_cols([" a "]) == ["a"]
    assert location._lower_cols(1) == []
    assert location._pick_col(df, ("station",)) == "Station"
    assert location._pick_col(df, ("missing",)) is None
    assert location._match_row(None, "S01") is None
    assert location._match_row(df.iloc[:0], "S01") is None
    assert location._match_row(df, "missing") is None
    assert location._match_row(df, "S01")["Station"] == " s01 "

    assert location._infer_utm_epsg(5, -3) == 32630
    assert location._infer_utm_epsg(-5, -3) == 32730
    dx, dy = location._flat_offsets_m(0, 0, 1, 1)
    assert dx > 0 and dy > 0
    assert location._parse_angle(None, {"S": -1}) != location._parse_angle(
        None, {"S": -1}
    )
    assert location._parse_angle("5:30:00S", {"S": -1}) == -5.5
    assert location._parse_angle("3.5W", {"W": -1}) == -3.5
    assert np.isnan(location._parse_angle("bad", {}))


def test_location_station_lookup_fallbacks():
    assert location._get_station(SimpleNamespace(name="fallback")) == "fallback"
    assert location._get_station(SimpleNamespace()) == ""


def test_selection_numeric_helpers(monkeypatch):
    assert selection._in_box(1, 2, 0, 0, 2, 3)
    assert not selection._in_box(np.nan, 2, 0, 0, 2, 3)
    assert selection._name_matches("abc", "ABC", case=False)

    monkeypatch.setattr(
        selection,
        "match_name",
        lambda *args: (_ for _ in ()).throw(ValueError()),
    )
    assert selection._name_matches("abc", "ABC", case=False)
    assert not selection._name_matches("abc", "ABC", case=True)

    assert selection._any_finite_z(SimpleNamespace(Z=None)) is False
    assert selection._any_finite_z(SimpleNamespace(Z=SimpleNamespace())) is False
    rho = SimpleNamespace(_resistivity=np.array([1.0, np.nan]))
    assert bool(selection._any_finite_z(SimpleNamespace(Z=rho)))
    z = SimpleNamespace(_z=np.array([np.nan + 1j]))
    assert bool(selection._any_finite_z(SimpleNamespace(Z=z)))
    assert np.isnan(selection._max_phase_err(SimpleNamespace(Z=None)))
    assert np.isnan(selection._max_phase_err(SimpleNamespace(Z=SimpleNamespace())))
    assert selection._max_phase_err(
        SimpleNamespace(Z=SimpleNamespace(_phase_err=[1, 4]))
    ) == 4


def test_recompute_naming_path_and_manifest_helpers(tmp_path, monkeypatch):
    assert recompute._render_name("{line}-{station}-{missing}", {
        "line": "L1", "station": "S1"
    }) == "L1-S1-"
    assert recompute._is_pathlike(tmp_path)
    assert not recompute._is_seq_of_pathlike("one.edi")
    assert recompute._is_seq_of_pathlike(["one.edi", tmp_path / "two.edi"])
    assert not recompute._is_seq_of_pathlike(1)
    assert not recompute._is_seq_of_pathlike([])
    assert recompute._source_path(tmp_path / "a.edi").name == "a.edi"
    assert recompute._source_path(SimpleNamespace()) is None
    assert recompute._source_stem(tmp_path / "a.edi") == "a"
    assert recompute._safe_station(tmp_path / "a.edi") == "a"

    records = [
        EDIRecomputeRecord(
            source=None,
            output=None,
            line=None,
            station="S1",
            status="skipped",
            message="none",
        )
    ]
    manifest = recompute._write_manifest(records, tmp_path / "out" / "m.csv")
    assert "station" in manifest.read_text()
    monkeypatch.setattr(recompute, "progress_enabled", lambda enabled: False)
    group = recompute._LineGroup(name="L1", sources=["a.edi"])
    assert list(recompute._iter_progress_groups([group], enabled=False, total=1))


def test_utils_matching_frequency_and_section_fallbacks(monkeypatch):
    assert utils.freq_select([1, 2, 3], 2).tolist() == [1]
    assert utils.freq_select([1, 2, 3], object()).size == 0
    assert utils.match_name(r"^A\d+$", "A12")
    assert not utils.match_name("[", "anything")

    broken = SimpleNamespace(get_section=lambda name: (_ for _ in ()).throw(KeyError()))
    assert utils._get_head(broken) is None
    assert utils._get_definemeas(broken) is None

    target = SimpleNamespace(get_section=lambda name: None)
    head = utils._ensure_head(target)
    assert head is target.Head
