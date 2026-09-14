from __future__ import annotations

import numpy as np

from pycsamt.airborne import (
    AirborneEMDataset,
    AirborneEMLine,
    AirborneEMRecord,
    NavigationTrack,
    assess_airborne_qc,
    inspect_airborne,
)
from pycsamt.airborne import qc as _qc
from pycsamt.airborne.mobilemt import build_mobilemt_emtf
from pycsamt.airborne.ztem import (
    ZTEMReferenceStation,
    build_ztem_emtf,
)
from pycsamt.emtf import EMTF, StatisticalEstimate, TransferFunction
from pycsamt.metadata import SiteMeta


def _dataset(*, reference=False, sparse=True):
    nav = NavigationTrack(
        sample_ids=("P1", "P2", "P3"),
        latitude=(5.0, 5.1, 5.2),
        longitude=(-3.0, -3.1, -3.2),
        terrain_elevation=(100.0, 101.0, 102.0),
        platform_elevation=(180.0, 181.0, 182.0),
    )
    from pycsamt.airborne import AirborneEMLine

    line = AirborneEMLine(
        line_id="L1",
        navigation=nav,
        attrs={"technology": "ZTEM"},
    )
    ref = None
    if reference:
        ref = ZTEMReferenceStation(
            station_id="BASE01",
            site=SiteMeta(site_id="BASE01"),
        )
    ids = ("P1", "P3") if sparse else ("P1", "P2", "P3")
    for sample_id in ids:
        doc = build_ztem_emtf(
            np.ones((2, 2), dtype=complex),
            frequency=[30.0, 90.0],
            reference_station=ref,
        )
        line.add_emtf(sample_id, doc)
    dataset = AirborneEMDataset(
        name="qc_test",
        attrs={"technology": "ZTEM"},
    )
    dataset.add_line(line)
    return dataset


def test_inspection_inventory():
    dataset = _dataset(reference=True, sparse=True)
    info = inspect_airborne(dataset)
    assert info.object_type == "dataset"
    assert info.technologies == ("ztem",)
    assert info.n_lines == 1
    assert info.n_samples == 3
    assert info.n_records == 2
    assert "tipper" in info.transfer_function_names
    assert info.bbox is not None


def test_qc_reports_structural_completeness_without_failing_sparse_data():
    report = assess_airborne_qc(_dataset(reference=True, sparse=True))
    assert report.technologies == ("ztem",)
    assert report.metrics["record_coverage_fraction"] == 2 / 3
    assert report.metrics["emtf_coverage_fraction"] == 2 / 3
    assert report.metrics["finite_response_fraction"] == 1.0
    assert report.metrics["reference_metadata_fraction"] == 1.0
    assert report.line_metrics["L1"]["location_fraction"] == 1.0
    assert report.line_metrics["L1"]["clearance_fraction"] == 1.0
    assert any(issue.code == "missing_em_records" for issue in report.issues)
    assert report.status == "pass"


def test_qc_flags_missing_required_reference_as_warning():
    report = assess_airborne_qc(_dataset(reference=False, sparse=False))
    assert report.metrics["reference_metadata_fraction"] == 0.0
    assert any(
        issue.code == "missing_reference_station"
        for issue in report.warnings
    )
    assert report.status == "warning"


def test_qc_covariance_metrics_are_descriptive():
    report = assess_airborne_qc(_dataset(reference=True, sparse=False))
    assert report.metrics["variance_tf_fraction"] == 0.0
    assert report.metrics["inverse_signal_covariance_tf_fraction"] == 0.0
    assert report.metrics["residual_covariance_tf_fraction"] == 0.0
    assert report.status == "pass"


def test_dataset_convenience_inspect_and_qc():
    dataset = _dataset(reference=True, sparse=False)
    assert dataset.inspect().technologies == ("ztem",)
    assert dataset.qc().status == "pass"


# ─────────────────────────────────────────────────────────────────────────
# inspect_airborne for line / record / bare EMTF
# ─────────────────────────────────────────────────────────────────────────


def test_inspect_airborne_line_record_and_emtf():
    dataset = _dataset(reference=True, sparse=True)
    line = next(dataset.iter_lines())
    line_info = inspect_airborne(line)
    assert line_info.object_type == "line"
    assert line_info.n_lines == 1
    assert line_info.n_records == 2
    assert line_info.attrs["line_id"] == "L1"

    record = line.get_record("P1")
    record_info = inspect_airborne(record)
    assert record_info.object_type == "record"
    assert record_info.n_records == 1
    assert record_info.attrs["sample_id"] == "P1"

    doc = record.emtf
    doc_info = inspect_airborne(doc)
    assert doc_info.object_type == "emtf"
    assert "tipper" in doc_info.transfer_function_names


def test_inspect_airborne_rejects_unsupported_type():
    import pytest

    with pytest.raises(TypeError):
        inspect_airborne(object())


def test_assess_airborne_qc_rejects_non_dataset():
    import pytest

    with pytest.raises(TypeError):
        assess_airborne_qc(object())


def test_airborne_qc_issue_rejects_empty_code_and_bad_severity():
    from pycsamt.airborne.qc import AirborneQCIssue

    import pytest

    with pytest.raises(ValueError):
        AirborneQCIssue(code="   ", severity="info", message="x")
    with pytest.raises(ValueError):
        AirborneQCIssue(code="foo", severity="bogus", message="x")


# ─────────────────────────────────────────────────────────────────────────
# AirborneQCReport.status / .errors for a genuinely invalid frequency axis
# ─────────────────────────────────────────────────────────────────────────


def test_qc_report_error_status_from_invalid_frequency_axis():
    dataset = _dataset(reference=True, sparse=True)
    line = next(dataset.iter_lines())
    record = line.get_record("P1")
    # Bypass EMTF.validate() to simulate an internally inconsistent
    # document (never producible through the public constructor).
    record.emtf.periods = np.array([-1.0, 90.0 ** -1])

    report = assess_airborne_qc(dataset)
    assert report.status == "error"
    assert any(issue.code == "invalid_frequency_axis" for issue in report.errors)
    assert report.errors and all(
        issue.severity == "error" for issue in report.errors
    )


# ─────────────────────────────────────────────────────────────────────────
# Private helper functions, exercised directly for edge-case branches
# ─────────────────────────────────────────────────────────────────────────


def test_finite_complex_fraction_helper():
    assert _qc._finite_complex_fraction(np.array([])) == (0, 0)
    finite, total = _qc._finite_complex_fraction(
        np.array([1.0 + 1.0j, np.nan + 1.0j, 2.0 + np.nan * 1j])
    )
    assert (finite, total) == (1, 3)
    finite, total = _qc._finite_complex_fraction(np.array([1.0, np.nan]))
    assert (finite, total) == (1, 2)


def test_navigation_location_count_projected_only():
    nav = NavigationTrack(
        sample_ids=("A", "B"),
        easting=(500000.0, 500100.0),
        northing=(4000000.0, 4000100.0),
    )
    line = AirborneEMLine(line_id="LP", navigation=nav)
    assert _qc._navigation_location_count(line) == 2


def test_navigation_location_count_none_when_no_coordinates():
    nav = NavigationTrack(sample_ids=("A", "B"))
    line = AirborneEMLine(line_id="LN", navigation=nav)
    assert _qc._navigation_location_count(line) == 0


def test_clearance_count_helper():
    nav_none = NavigationTrack(sample_ids=("A", "B"))
    line_none = AirborneEMLine(line_id="LC1", navigation=nav_none)
    assert _qc._clearance_count(line_none) == 0

    nav_explicit = NavigationTrack(
        sample_ids=("A", "B"), clearance=(30.0, float("nan")),
    )
    line_explicit = AirborneEMLine(line_id="LC2", navigation=nav_explicit)
    assert _qc._clearance_count(line_explicit) == 1


def test_reference_present_default_true_for_unhandled_technology():
    doc = EMTF()
    assert _qc._reference_present(doc, "afmag") is True


# ─────────────────────────────────────────────────────────────────────────
# assess_airborne_qc: missing navigation coordinates / mixed line technology
# ─────────────────────────────────────────────────────────────────────────


def test_qc_flags_missing_navigation_coordinates():
    nav = NavigationTrack(sample_ids=("P1",))
    line = AirborneEMLine(line_id="LMISS", navigation=nav)
    line.add_emtf(
        "P1",
        build_ztem_emtf(np.ones((2, 2), dtype=complex), frequency=[30.0, 90.0]),
    )
    dataset = AirborneEMDataset(name="no_coords_test")
    dataset.add_line(line)
    report = assess_airborne_qc(dataset)
    assert any(
        issue.code == "missing_navigation_coordinates"
        for issue in report.warnings
    )


def test_qc_flags_mixed_line_technology():
    nav = NavigationTrack(
        sample_ids=("A", "B"), latitude=(5.0, 5.1), longitude=(-3.0, -3.1),
    )
    line = AirborneEMLine(line_id="LMIX", navigation=nav)
    line.add_emtf(
        "A",
        build_ztem_emtf(np.ones((2, 2), dtype=complex), frequency=[30.0, 90.0]),
    )
    line.add_emtf(
        "B",
        build_mobilemt_emtf(
            np.ones((1, 3, 2), dtype=complex), frequency=[100.0],
        ),
    )
    dataset = AirborneEMDataset(name="mixed_line_test")
    dataset.add_line(line)
    report = assess_airborne_qc(dataset)
    assert any(
        issue.code == "mixed_line_technology" for issue in report.warnings
    )


def test_qc_flags_record_without_emtf():
    nav = NavigationTrack(
        sample_ids=("A", "B"), latitude=(5.0, 5.1), longitude=(-3.0, -3.1),
    )
    line = AirborneEMLine(line_id="LNOEMTF", navigation=nav)
    line.add_record(AirborneEMRecord(sample_id="A", emtf=None))
    line.add_emtf(
        "B",
        build_ztem_emtf(np.ones((2, 2), dtype=complex), frequency=[30.0, 90.0]),
    )
    dataset = AirborneEMDataset(name="record_without_emtf_test")
    dataset.add_line(line)
    report = assess_airborne_qc(dataset)
    assert any(
        issue.code == "record_without_emtf" and issue.sample_id == "A"
        for issue in report.issues
    )


def test_qc_flags_mixed_dataset_technology():
    ztem_nav = NavigationTrack(
        sample_ids=("A",), latitude=(5.0,), longitude=(-3.0,),
    )
    ztem_line = AirborneEMLine(line_id="LZ", navigation=ztem_nav)
    ztem_line.add_emtf(
        "A",
        build_ztem_emtf(np.ones((2, 2), dtype=complex), frequency=[30.0, 90.0]),
    )
    mmt_nav = NavigationTrack(
        sample_ids=("B",), latitude=(6.0,), longitude=(-4.0,),
    )
    mmt_line = AirborneEMLine(line_id="LM", navigation=mmt_nav)
    mmt_line.add_emtf(
        "B",
        build_mobilemt_emtf(
            np.ones((1, 3, 2), dtype=complex), frequency=[100.0],
        ),
    )
    dataset = AirborneEMDataset(name="mixed_dataset_test")
    dataset.add_line(ztem_line)
    dataset.add_line(mmt_line)
    report = assess_airborne_qc(dataset)
    assert len(report.technologies) > 1
    assert any(
        issue.code == "mixed_dataset_technology" for issue in report.issues
    )


def test_qc_skips_reference_check_for_technology_not_requiring_it():
    from pycsamt.airborne.afmag import build_original_afmag_emtf

    nav = NavigationTrack(
        sample_ids=("A",), latitude=(5.0,), longitude=(-3.0,),
    )
    line = AirborneEMLine(line_id="LAFM", navigation=nav)
    line.add_emtf(
        "A",
        build_original_afmag_emtf([1.0, 2.0], frequency=[150.0, 510.0]),
    )
    dataset = AirborneEMDataset(name="afmag_no_reference_required_test")
    dataset.add_line(line)
    report = assess_airborne_qc(dataset)
    assert "afmag" in report.technologies
    assert not any(
        issue.code == "missing_reference_station" for issue in report.issues
    )
    assert report.metrics["reference_metadata_fraction"] == 1.0


def test_qc_derived_transfer_function_and_estimate_coverage():
    doc = EMTF(periods=[1.0 / 30.0, 1.0 / 90.0])
    doc.add_transfer_function(
        TransferFunction(
            name="T",
            data=np.ones((2, 1, 2), dtype=complex),
            input_channels=("Hx", "Hy"),
            output_channels=("Hz",),
            periods=[1.0 / 30.0, 1.0 / 90.0],
        )
    )
    tipper_tf = doc.get_transfer_function("tipper")
    tipper_tf.add_estimate(
        StatisticalEstimate(name="VAR", data=np.ones(2), kind="variance")
    )
    tipper_tf.add_estimate(
        StatisticalEstimate(
            name="INVSIGCOV", data=np.ones(2), kind="inverse_signal_covariance",
        )
    )
    tipper_tf.add_estimate(
        StatisticalEstimate(
            name="RESIDCOV", data=np.ones(2), kind="residual_covariance",
        )
    )
    # Derived transfer function (auto-tagged intention="derived").
    doc.add_transfer_function(
        TransferFunction(
            name="TIPMAG", data=np.ones(2), periods=[1.0 / 30.0, 1.0 / 90.0],
        )
    )

    nav = NavigationTrack(
        sample_ids=("A",), latitude=(5.0,), longitude=(-3.0,),
    )
    line = AirborneEMLine(line_id="LDER", navigation=nav)
    line.add_emtf("A", doc)
    dataset = AirborneEMDataset(name="derived_tf_test")
    dataset.add_line(line)

    report = assess_airborne_qc(dataset)
    assert report.metrics["n_primary_transfer_functions"] == 1
    assert report.metrics["n_derived_transfer_functions"] == 1
    assert report.metrics["variance_tf_fraction"] == 1.0
    assert report.metrics["inverse_signal_covariance_tf_fraction"] == 1.0
    assert report.metrics["residual_covariance_tf_fraction"] == 1.0
