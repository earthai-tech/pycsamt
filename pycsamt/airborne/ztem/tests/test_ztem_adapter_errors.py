from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import AirborneEMLine, NavigationTrack
from pycsamt.airborne.ztem import (
    ZTEMReferenceStation,
    ZTEMSystemSpec,
    ZTEMValidationError,
    build_ztem_dataset,
    build_ztem_emtf,
    build_ztem_line,
    validate_ztem_transfer_function,
)
from pycsamt.airborne.ztem.adapter import (
    _normalize_tipper,
    _processing_for_reference,
    _reference_mapping,
    _sample_axis_tipper,
    _xml_notes_mapping,
)
from pycsamt.emtf import TransferFunction
from pycsamt.metadata import (
    LocationMeta,
    ProcessingMeta,
    RemoteReferenceMeta,
    SiteMeta,
)


def _tipper(nf: int = 2) -> np.ndarray:
    base = np.arange(nf * 2, dtype=float).reshape(nf, 2)
    return base + 1j * (base + 0.25)


def _nav(n: int = 2) -> NavigationTrack:
    return NavigationTrack(
        sample_ids=tuple(f"S{i:03d}" for i in range(n)),
        latitude=[10.0 + 0.001 * i for i in range(n)],
        longitude=[20.0 + 0.001 * i for i in range(n)],
    )


# ─────────────────────────────────────────────────────────────────────────
# _normalize_tipper
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_tipper_rejects_bad_shape():
    with pytest.raises(ZTEMValidationError):
        _normalize_tipper(np.ones((2, 3)))


def test_normalize_tipper_rejects_non_numeric():
    with pytest.raises(ZTEMValidationError):
        _normalize_tipper(np.array(["a", "b"]))


# ─────────────────────────────────────────────────────────────────────────
# _reference_mapping / _xml_notes_mapping
# ─────────────────────────────────────────────────────────────────────────


def test_reference_mapping_includes_attrs_when_present():
    ref = ZTEMReferenceStation(station_id="BASE01", attrs={"note": "fixed"})
    out = _reference_mapping(ref)
    assert out["attrs"] == {"note": "fixed"}


def test_xml_notes_mapping_includes_full_reference_location():
    ref = ZTEMReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(
                latitude=5.0, longitude=-3.0, elevation=120.0, datum="WGS84",
            ),
        ),
    )
    notes = _xml_notes_mapping(ZTEMSystemSpec(), ref)
    ztem = notes["ZTEM"]
    assert ztem["ReferenceLatitude"] == 5.0
    assert ztem["ReferenceLongitude"] == -3.0
    assert ztem["ReferenceElevation"] == 120.0
    assert ztem["ReferenceDatum"] == "WGS84"


def test_xml_notes_mapping_reference_station_without_id_or_site():
    ref = ZTEMReferenceStation()
    notes = _xml_notes_mapping(ZTEMSystemSpec(), ref)
    ztem = notes["ZTEM"]
    assert "ReferenceStationId" not in ztem
    assert "ReferenceLatitude" not in ztem


def test_xml_notes_mapping_reference_location_with_no_optional_fields():
    ref = ZTEMReferenceStation(
        site=SiteMeta(site_id="BASE02", location=LocationMeta(datum=None)),
    )
    notes = _xml_notes_mapping(ZTEMSystemSpec(), ref)
    ztem = notes["ZTEM"]
    assert "ReferenceLatitude" not in ztem
    assert "ReferenceLongitude" not in ztem
    assert "ReferenceElevation" not in ztem
    assert "ReferenceDatum" not in ztem


# ─────────────────────────────────────────────────────────────────────────
# _processing_for_reference
# ─────────────────────────────────────────────────────────────────────────


def test_processing_for_reference_rejects_bad_processing_type():
    ref = ZTEMReferenceStation(station_id="BASE01")
    with pytest.raises(TypeError):
        _processing_for_reference(ref, "not-a-processing-meta")


def test_processing_for_reference_merges_when_no_existing_remote():
    ref = ZTEMReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(processed_by="Alice", run_list=["run1"])
    merged = _processing_for_reference(ref, processing)
    assert merged.processed_by == "Alice"
    assert merged.run_list == ["run1"]
    assert merged.remote_reference.site == "BASE01"


def test_processing_for_reference_passthrough_when_site_matches():
    ref = ZTEMReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(
        remote_reference=RemoteReferenceMeta(
            reference_type="fixed_ground_horizontal_magnetic", site="BASE01",
        ),
    )
    merged = _processing_for_reference(ref, processing)
    assert merged is processing


def test_processing_for_reference_skips_conflict_check_when_existing_site_missing():
    ref = ZTEMReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(
        remote_reference=RemoteReferenceMeta(
            reference_type="fixed_ground_horizontal_magnetic", site=None,
        ),
    )
    merged = _processing_for_reference(ref, processing)
    assert merged is processing


# ─────────────────────────────────────────────────────────────────────────
# validate_ztem_transfer_function
# ─────────────────────────────────────────────────────────────────────────


def test_validate_ztem_tf_rejects_non_tf():
    with pytest.raises(TypeError):
        validate_ztem_transfer_function("not-a-tf")


def test_validate_ztem_tf_rejects_wrong_name():
    tf = TransferFunction(
        name="Z",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    with pytest.raises(ZTEMValidationError):
        validate_ztem_transfer_function(tf)


def test_validate_ztem_tf_rejects_wrong_channels():
    from pycsamt.airborne.ztem.constants import ZTEM_TIPPER_TAG

    bad_input = TransferFunction(
        name=ZTEM_TIPPER_TAG,
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Hz",),
        periods=[1.0],
    )
    with pytest.raises(ZTEMValidationError):
        validate_ztem_transfer_function(bad_input)

    bad_output = TransferFunction(
        name=ZTEM_TIPPER_TAG,
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ez",),
        periods=[1.0],
    )
    with pytest.raises(ZTEMValidationError):
        validate_ztem_transfer_function(bad_output)


def test_validate_ztem_tf_rejects_mutated_bad_shape():
    from pycsamt.airborne.ztem.constants import ZTEM_TIPPER_TAG

    tf = TransferFunction(
        name=ZTEM_TIPPER_TAG,
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0],
    )
    tf.data = np.ones((1, 1, 3), dtype=complex)
    with pytest.raises(ZTEMValidationError):
        validate_ztem_transfer_function(tf)


# ─────────────────────────────────────────────────────────────────────────
# build_ztem_emtf
# ─────────────────────────────────────────────────────────────────────────


def test_build_ztem_emtf_rejects_frequency_count_mismatch():
    with pytest.raises(ZTEMValidationError):
        build_ztem_emtf(_tipper(2), frequency=[10.0, 20.0, 30.0])


def test_build_ztem_emtf_rejects_bad_system_spec_type():
    with pytest.raises(TypeError):
        build_ztem_emtf(_tipper(1), frequency=[10.0], system_spec="bad")


def test_build_ztem_emtf_rejects_bad_reference_station_type():
    with pytest.raises(TypeError):
        build_ztem_emtf(
            _tipper(1), frequency=[10.0], reference_station="bad",
        )


def test_build_ztem_emtf_rejects_bad_site_type():
    with pytest.raises(TypeError):
        build_ztem_emtf(_tipper(1), frequency=[10.0], site="bad")


def test_build_ztem_emtf_rejects_bad_orientation_type():
    with pytest.raises(TypeError):
        build_ztem_emtf(_tipper(1), frequency=[10.0], orientation="bad")


# ─────────────────────────────────────────────────────────────────────────
# _sample_axis_tipper
# ─────────────────────────────────────────────────────────────────────────


def test_sample_axis_tipper_single_sample_vector_promotion():
    arr = _sample_axis_tipper(np.array([1.0 + 1.0j, 2.0 + 2.0j]), n_samples=1)
    assert arr.shape == (1, 1, 1, 2)


def test_sample_axis_tipper_single_sample_stacked_promotion():
    arr = _sample_axis_tipper(_tipper(3), n_samples=1)
    assert arr.shape == (1, 3, 1, 2)


def test_sample_axis_tipper_single_sample_canonical_3d_promotion():
    arr = _sample_axis_tipper(_tipper(3)[:, None, :], n_samples=1)
    assert arr.shape == (1, 3, 1, 2)


def test_sample_axis_tipper_batched_component_axis_wrong_size():
    bad = np.ones((2, 3, 3))
    with pytest.raises(ZTEMValidationError):
        _sample_axis_tipper(bad, n_samples=2)


def test_sample_axis_tipper_rejects_bad_shape_for_n_samples():
    with pytest.raises(ZTEMValidationError):
        _sample_axis_tipper(np.ones((2, 2)), n_samples=2)


def test_sample_axis_tipper_rejects_bad_matrix_tail_shape():
    bad = np.ones((2, 3, 2, 2), dtype=complex)
    with pytest.raises(ZTEMValidationError):
        _sample_axis_tipper(bad, n_samples=2)


def test_sample_axis_tipper_rejects_non_numeric():
    bad = np.array([[["a", "b"]], [["c", "d"]]])
    with pytest.raises(ZTEMValidationError):
        _sample_axis_tipper(bad, n_samples=2)


# ─────────────────────────────────────────────────────────────────────────
# build_ztem_line
# ─────────────────────────────────────────────────────────────────────────


def test_build_ztem_line_rejects_non_navigation():
    with pytest.raises(TypeError):
        build_ztem_line(
            "L1", "not-a-navigation", _tipper(1), frequency=[10.0],
        )


# ─────────────────────────────────────────────────────────────────────────
# build_ztem_dataset
# ─────────────────────────────────────────────────────────────────────────


def test_build_ztem_dataset_rejects_bad_survey_type():
    line = build_ztem_line("L1", _nav(1), _tipper(1), frequency=[10.0])
    with pytest.raises(TypeError):
        build_ztem_dataset("SURVEY", [line], survey="not-a-survey")


def test_build_ztem_dataset_rejects_bad_system_spec_type():
    line = build_ztem_line("L1", _nav(1), _tipper(1), frequency=[10.0])
    with pytest.raises(TypeError):
        build_ztem_dataset("SURVEY", [line], system_spec="bad")


def test_build_ztem_dataset_rejects_non_line_items():
    with pytest.raises(TypeError):
        build_ztem_dataset("SURVEY", ["not-a-line"])


def test_build_ztem_dataset_rejects_mismatched_technology_tag():
    nav = _nav(1)
    line = AirborneEMLine(
        line_id="L1", navigation=nav, attrs={"technology": "MobileMT"},
    )
    with pytest.raises(ZTEMValidationError):
        build_ztem_dataset("SURVEY", [line])
