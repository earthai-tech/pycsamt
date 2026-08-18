from __future__ import annotations

import numpy as np

from pycsamt.airborne import (
    AirborneEMDataset,
    NavigationTrack,
    assess_airborne_qc,
    inspect_airborne,
)
from pycsamt.airborne.ztem import (
    ZTEMReferenceStation,
    build_ztem_emtf,
)
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
