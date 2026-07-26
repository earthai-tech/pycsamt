from __future__ import annotations

from verify_pytest_junit import passing_report


def test_accepts_complete_passing_report(tmp_path):
    report = tmp_path / "pass.xml"
    report.write_text(
        '<testsuites><testsuite tests="12" failures="0" errors="0" '
        'skipped="2" /></testsuites>',
        encoding="utf-8",
    )
    assert passing_report(report)


def test_rejects_failure_error_empty_and_malformed_reports(tmp_path):
    for name, xml in {
        "failure": '<testsuite tests="2" failures="1" errors="0" />',
        "error": '<testsuite tests="2" failures="0" errors="1" />',
        "empty": '<testsuite tests="0" failures="0" errors="0" />',
        "malformed": "<testsuite",
    }.items():
        report = tmp_path / f"{name}.xml"
        report.write_text(xml, encoding="utf-8")
        assert not passing_report(report)
