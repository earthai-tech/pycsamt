"""Deep native-file recognition and validation contracts."""

import json

import pytest

from pycsamt.models.occam1d import (
    Occam1DFileType,
    Occam1DValidationStatus,
    inspect_occam1d_file,
    scan_occam1d_directory,
    validate_occam1d_file,
)
from pycsamt.models.occam1d.tests.test_occam1d_results import _make_run
from pycsamt.models.occam1d.validation import (
    detect_file_type,
    is_data_file,
    is_iter_file,
    is_log_file,
    is_model_file,
    is_response_file,
    is_startup_file,
)


def test_all_native_run_files_are_classified_with_metadata(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    expected = {
        "Occam1DData": Occam1DFileType.DATA,
        "Occam1DModel": Occam1DFileType.MODEL,
        "Startup": Occam1DFileType.STARTUP,
        "ITER_01": Occam1DFileType.ITER,
        "RESP_02.resp": Occam1DFileType.RESPONSE,
        "occam1d.log": Occam1DFileType.LOG,
    }
    for name, kind in expected.items():
        report = inspect_occam1d_file(run / name)
        assert report.valid
        assert report.file_type is kind
        assert report.lines_examined > 0
    assert inspect_occam1d_file(run / "Startup").iteration == 0
    assert inspect_occam1d_file(run / "ITER_01").iteration == 1


def test_compatibility_predicates_delegate_to_unambiguous_detection(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    assert is_data_file(run / "Occam1DData")
    assert is_model_file(run / "Occam1DModel")
    assert is_startup_file(run / "Startup")
    assert is_iter_file(run / "ITER_02")
    assert is_response_file(run / "RESP_02.resp")
    assert is_log_file(run / "occam1d.log")
    assert detect_file_type(run / "missing") is Occam1DFileType.UNKNOWN


def test_missing_directory_and_binary_states_are_explicit(tmp_path):
    missing = inspect_occam1d_file(tmp_path / "missing")
    assert missing.status is Occam1DValidationStatus.MISSING
    directory = inspect_occam1d_file(tmp_path)
    assert directory.status is Occam1DValidationStatus.NOT_FILE
    binary = tmp_path / "binary.dat"
    binary.write_bytes(b"EMData_1.1\x00content")
    report = inspect_occam1d_file(binary)
    assert report.status is Occam1DValidationStatus.INVALID


def test_ambiguous_signatures_are_not_silently_prioritized(tmp_path):
    path = tmp_path / "ambiguous.txt"
    path.write_text(
        "Format: EMData_1.1\n"
        "# Frequencies: 1\n"
        "# Data: 0\n"
        "Format: Resistivity1DMod_1.0\n"
        "# Layers: 2\n",
        encoding="utf8",
    )
    report = inspect_occam1d_file(path)
    assert report.ambiguous
    assert set(report.matches) == {
        Occam1DFileType.DATA,
        Occam1DFileType.MODEL,
    }
    assert detect_file_type(path) is Occam1DFileType.UNKNOWN


def test_shallow_signature_can_be_rejected_by_full_parser(tmp_path):
    path = tmp_path / "truncated.data"
    path.write_text(
        "Format: EMData_1.1\n# Frequencies: 2\n# Data: 0\n",
        encoding="utf8",
    )
    shallow = validate_occam1d_file(path)
    deep = validate_occam1d_file(path, expected="data", deep=True)
    assert shallow.valid
    assert deep.status is Occam1DValidationStatus.INVALID
    assert any("Full parser rejected" in item for item in deep.reasons)


def test_expected_type_mismatch_returns_invalid_report(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    report = validate_occam1d_file(
        run / "Occam1DModel",
        expected=Occam1DFileType.DATA,
    )
    assert report.status is Occam1DValidationStatus.INVALID
    assert report.file_type is Occam1DFileType.MODEL


def test_report_is_json_safe(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    values = inspect_occam1d_file(run / "ITER_02").to_dict()
    assert json.loads(json.dumps(values))["iteration"] == 2
    assert values["file_type"] == "iter"


def test_bounded_inspection_reports_truncation(tmp_path):
    path = tmp_path / "large.txt"
    path.write_text("unrelated\n" * 20, encoding="utf8")
    report = inspect_occam1d_file(path, max_lines=3)
    assert report.truncated
    assert report.lines_examined == 3
    assert report.status is Occam1DValidationStatus.UNKNOWN


def test_directory_scan_is_sorted_and_can_include_unknown(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    (run / "notes.txt").write_text("not native", encoding="utf8")
    native = scan_occam1d_directory(run)
    all_reports = scan_occam1d_directory(run, include_unknown=True)
    assert len(all_reports) == len(native) + 1
    assert [item.path.name for item in all_reports] == sorted(
        (item.path.name for item in all_reports),
        key=str.lower,
    )


def test_public_resource_limits_are_validated(tmp_path):
    path = tmp_path / "candidate"
    path.write_text("text", encoding="utf8")
    with pytest.raises(ValueError, match="max_lines"):
        inspect_occam1d_file(path, max_lines=0)
    with pytest.raises(TypeError, match="max_bytes"):
        inspect_occam1d_file(path, max_bytes=True)
    with pytest.raises(ValueError, match="expected"):
        validate_occam1d_file(path, expected="mesh")
