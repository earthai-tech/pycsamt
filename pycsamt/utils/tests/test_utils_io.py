from __future__ import annotations

import numpy as np
import pytest

from pycsamt.exceptions import PycsamtError, ValidationError
from pycsamt.utils.io import parse_stn_profile, stn_separation

# ------------------------------ stn_separation -----------------------------


def test_stn_separation_normal_1d_arrays():
    easting = [0, 10, 20, 30]
    northing = [0, 0, 0, 0]
    seps, mean_sep = stn_separation(easting, northing)
    assert np.allclose(seps, [10, 10, 10])
    assert mean_sep == pytest.approx(10.0)


def test_stn_separation_interpolate_true_matches_point_count():
    easting = [0, 10, 20, 30]
    northing = [0, 0, 0, 0]
    seps, mean_sep = stn_separation(easting, northing, interpolate=True)
    assert seps.shape[0] == len(easting)
    assert mean_sep == pytest.approx(10.0)


def test_stn_separation_interpolate_false_length_n_minus_1():
    easting = [0, 10, 20, 30]
    northing = [0, 0, 0, 0]
    seps, _ = stn_separation(easting, northing, interpolate=False)
    assert seps.shape[0] == len(easting) - 1


def test_stn_separation_mismatched_shapes_raises_validation_error():
    with pytest.raises(ValidationError, match="shapes differ"):
        stn_separation([0, 1, 2], [0, 1])


def test_stn_separation_fewer_than_two_points_returns_empty():
    seps, mean_sep = stn_separation([0], [0])
    assert seps.size == 0
    assert mean_sep == 0.0

    seps2, mean_sep2 = stn_separation([], [])
    assert seps2.size == 0
    assert mean_sep2 == 0.0


def test_stn_separation_non_numeric_input_raises_validation_error():
    with pytest.raises(ValidationError, match="Invalid coordinate arrays"):
        stn_separation(["a", "b"], [0, 1])


# ----------------------------- parse_stn_profile ----------------------------


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content, encoding="utf8")
    return str(p)


def test_parse_stn_profile_normal(tmp_path):
    content = (
        "dot east north elev\n"
        "0 100.0 200.0 10.0\n"
        "1 110.0 210.0 11.0\n"
        "2 120.0 220.0 12.0\n"
    )
    path = _write(tmp_path, "normal.stn", content)
    result = parse_stn_profile(path)
    assert np.allclose(result["position"], [0, 1, 2])
    assert np.allclose(result["easting"], [100, 110, 120])
    assert np.allclose(result["northing"], [200, 210, 220])
    assert np.allclose(result["elevation"], [10, 11, 12])


def test_parse_stn_profile_comment_lines_skipped(tmp_path):
    content = (
        "> header comment\n"
        "! another comment\n"
        "dot east north elev\n"
        "0 100.0 200.0 10.0\n"
        "1 110.0 210.0 11.0\n"
    )
    path = _write(tmp_path, "comments.stn", content)
    result = parse_stn_profile(path)
    assert result["position"].tolist() == [0.0, 1.0]


def test_parse_stn_profile_custom_delimiter(tmp_path):
    content = (
        "dot,east,north,elev\n"
        "0,100.0,200.0,10.0\n"
        "1,110.0,210.0,11.0\n"
    )
    path = _write(tmp_path, "delim.stn", content)
    result = parse_stn_profile(path, delimiter=",")
    assert result["easting"].tolist() == [100.0, 110.0]


def test_parse_stn_profile_mismatched_field_count_raises(tmp_path):
    content = (
        "dot east north elev\n"
        "0 100.0 200.0 10.0\n"
        "1 110.0 210.0\n"  # missing one field
    )
    path = _write(tmp_path, "mismatch.stn", content)
    with pytest.raises(PycsamtError, match="fields but expected"):
        parse_stn_profile(path)


def test_parse_stn_profile_missing_required_column_raises(tmp_path):
    content = "dot east elev\n0 100.0 10.0\n1 110.0 11.0\n"
    path = _write(tmp_path, "missingcol.stn", content)
    with pytest.raises(PycsamtError, match="Missing required column"):
        parse_stn_profile(path)


def test_parse_stn_profile_file_not_found_raises(tmp_path):
    missing_path = str(tmp_path / "does_not_exist.stn")
    with pytest.raises(PycsamtError, match="File not found"):
        parse_stn_profile(missing_path)


def test_parse_stn_profile_no_data_after_stripping_comments_raises(tmp_path):
    content = "> only a comment\n! another comment\n"
    path = _write(tmp_path, "onlycomments.stn", content)
    with pytest.raises(PycsamtError, match="No data found"):
        parse_stn_profile(path)


def test_parse_stn_profile_alternate_column_name_aliases(tmp_path):
    content = (
        "station e n h\n"
        "0 100.0 200.0 10.0\n"
        "1 110.0 210.0 11.0\n"
    )
    path = _write(tmp_path, "aliases.stn", content)
    result = parse_stn_profile(path)
    assert result["position"].tolist() == [0.0, 1.0]
    assert result["easting"].tolist() == [100.0, 110.0]
    assert result["northing"].tolist() == [200.0, 210.0]
    assert result["elevation"].tolist() == [10.0, 11.0]


def test_parse_stn_profile_extra_column_preserved_under_raw_name(tmp_path):
    content = (
        'dot east north elev "resistivity"\n'
        "0 100.0 200.0 10.0 55.5\n"
        "1 110.0 210.0 11.0 66.6\n"
    )
    path = _write(tmp_path, "extra.stn", content)
    result = parse_stn_profile(path)
    assert result["resistivity"].tolist() == [55.5, 66.6]
