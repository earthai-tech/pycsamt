from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import NavigationTrack
from pycsamt.airborne.afmag import (
    AFMAG_AP_TAG,
    AFMAG_TENSOR_INPUT_CHANNELS,
    AFMAG_TENSOR_OUTPUT_CHANNELS,
    AFMAG_TENSOR_TAG,
    AFMAG_TILT_TAG,
    AFMAGReferenceStation,
    AFMAGValidationError,
    AirMtSystemSpec,
    OriginalAFMAGSystemSpec,
    build_airmt_dataset,
    build_airmt_emtf,
    build_airmt_line,
    build_original_afmag_dataset,
    build_original_afmag_emtf,
    build_original_afmag_line,
    compute_airmt_amplification_parameter,
    register_afmag_datatypes,
)
from pycsamt.emtf import get_emtf_datatype
from pycsamt.metadata import LocationMeta, OrientationMeta, SiteMeta


def _tensor(nf: int = 3) -> np.ndarray:
    data = np.zeros((nf, 3, 2), dtype=complex)
    for index in range(nf):
        scale = index + 1.0
        data[index] = np.array(
            [
                [1.0 + 0.1j * scale, 0.2 + 0.05j],
                [0.1 - 0.03j, 1.2 + 0.08j * scale],
                [0.3 + 0.02j, -0.15 + 0.04j],
            ]
        )
    return data


def test_afmag_datatypes_reuse_ti_and_register_derived_products():
    register_afmag_datatypes()
    ti = get_emtf_datatype("TI")
    assert ti is not None
    assert ti.tag == AFMAG_TENSOR_TAG
    ap = get_emtf_datatype("AP")
    assert ap is not None
    assert ap.tag == AFMAG_AP_TAG
    assert ap.is_derived
    tilt = get_emtf_datatype("AFTILT")
    assert tilt is not None
    assert tilt.tag == AFMAG_TILT_TAG
    assert tilt.is_derived


def test_original_afmag_system_spec_is_descriptive():
    spec = OriginalAFMAGSystemSpec()
    assert spec.typical_frequencies_hz == (150.0, 510.0)
    assert spec.coil_count == 2
    assert spec.coil_tilt_deg == 45.0
    assert spec.coil_separation_deg == 45.0
    assert spec.digital_recording is False


def test_airmt_system_spec_frequency_mask():
    spec = AirMtSystemSpec()
    mask = spec.practical_frequency_mask([10.0, 20.0, 800.0, 900.0])
    assert np.array_equal(mask, [False, True, True, False])
    instrument = spec.to_instrument_meta(serial="AIRMT-TEST")
    assert instrument.system == "AirMt / tensor AFMAG"
    assert instrument.electric_sensor is None


def test_original_afmag_tilt_angle_defaults_to_degrees():
    doc = build_original_afmag_emtf(
        [2.0, -1.5],
        frequency=[150.0, 510.0],
    )
    tf = doc.get_transfer_function(AFMAG_TILT_TAG)
    assert tf is not None
    assert tf.shape == (2, 1, 1)
    assert tf.units == "deg"
    assert not tf.input_channels
    assert not tf.output_channels
    assert tf.attrs["response_kind"] == "tilt_angle"
    assert doc.subtype == "afmag_original"


def test_original_afmag_comparator_deflection_does_not_invent_units():
    doc = build_original_afmag_emtf(
        [10.0, 12.0],
        frequency=[150.0, 510.0],
        response_kind="comparator_deflection",
    )
    tf = doc.get_transfer_function(AFMAG_TILT_TAG)
    assert tf.units is None
    assert tf.attrs["calibrated_angle"] is False


def test_original_afmag_tilt_must_be_real():
    with pytest.raises(AFMAGValidationError):
        build_original_afmag_emtf(
            [1.0 + 2.0j],
            frequency=[150.0],
        )


def test_airmt_tensor_reuses_interstation_magnetic_tf():
    freq = [24.0, 75.0, 300.0]
    doc = build_airmt_emtf(_tensor(), frequency=freq)
    tf = doc.get_transfer_function(AFMAG_TENSOR_TAG)
    assert tf is not None
    assert tf.shape == (3, 3, 2)
    assert tf.input_channels == AFMAG_TENSOR_INPUT_CHANNELS
    assert tf.output_channels == AFMAG_TENSOR_OUTPUT_CHANNELS
    assert np.allclose(tf.frequency, freq)
    assert doc.subtype == "afmag_airmt"


def test_airmt_amplification_parameter_identity_is_one():
    tensor = np.array(
        [
            [1.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, 1.0 + 0.0j],
            [0.0 + 0.0j, 0.0 + 0.0j],
        ]
    )
    ap = compute_airmt_amplification_parameter(tensor)
    assert ap.shape == ()
    assert np.allclose(ap, 1.0 + 0.0j)


def test_airmt_amplification_parameter_is_added_as_scalar_derived_product():
    doc = build_airmt_emtf(_tensor(2), frequency=[24.0, 75.0])
    ap = doc.get_transfer_function(AFMAG_AP_TAG)
    assert ap is not None
    assert ap.shape == (2, 1, 1)
    assert not ap.input_channels
    assert not ap.output_channels
    expected = compute_airmt_amplification_parameter(_tensor(2))
    assert np.allclose(ap.data[:, 0, 0], expected)


def test_airmt_amplification_parameter_zero_projection_policy():
    zero = np.zeros((3, 2), dtype=complex)
    result = compute_airmt_amplification_parameter(zero)
    assert np.isnan(result.real)
    assert np.isnan(result.imag)
    with pytest.raises(AFMAGValidationError):
        compute_airmt_amplification_parameter(zero, zero_policy="raise")


def test_airmt_covariance_shapes_are_preserved():
    nf = 2
    doc = build_airmt_emtf(
        _tensor(nf),
        frequency=[24.0, 75.0],
        variance=np.ones((nf, 3, 2)),
        inverse_signal_covariance=np.tile(
            np.eye(2), (nf, 1, 1)
        ).astype(complex),
        residual_covariance=np.tile(
            np.eye(3), (nf, 1, 1)
        ).astype(complex),
    )
    tf = doc.get_transfer_function(AFMAG_TENSOR_TAG)
    assert tf.get_estimate("VAR").shape == (nf, 3, 2)
    assert tf.get_estimate("INVSIGCOV").shape == (nf, 2, 2)
    assert tf.get_estimate("RESIDCOV").shape == (nf, 3, 3)


def test_afmag_reference_station_maps_to_processing_metadata():
    ref = AFMAGReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(latitude=5.2, longitude=-3.1),
        ),
    )
    doc = build_airmt_emtf(
        _tensor(2),
        frequency=[24.0, 75.0],
        reference_station=ref,
    )
    remote = doc.processing.remote_reference
    assert remote.reference_type == "fixed_ground_magnetic"
    assert remote.site == "BASE01"
    assert remote.extra["measured_channels"] == ["Hx", "Hy", "Hz"]
    assert remote.extra["transfer_input_channels"] == ["Hx", "Hy"]


def test_airmt_does_not_infer_orientation():
    doc = build_airmt_emtf(_tensor(2), frequency=[24.0, 75.0])
    assert doc.orientation is None


def test_airmt_rotation_preserves_recomputed_amplification_parameter():
    source = build_airmt_emtf(
        _tensor(3),
        frequency=[24.0, 75.0, 300.0],
        orientation=OrientationMeta(
            mode="orthogonal",
            angle_to_geographic_north=0.0,
        ),
    )
    original_ap = source.get_transfer_function(AFMAG_AP_TAG).data[:, 0, 0]
    rotated = source.rotate(37.0, derived_policy="drop")
    assert rotated.get_transfer_function(AFMAG_AP_TAG) is None
    rotated_tf = rotated.get_transfer_function(AFMAG_TENSOR_TAG)
    recomputed = compute_airmt_amplification_parameter(rotated_tf.data)
    assert np.allclose(recomputed, original_ap, atol=1.0e-12, rtol=1.0e-12)


def test_build_airmt_line_shared_grid_and_sparse_records():
    nav = NavigationTrack(sample_ids=("P001", "P002", "P003"))
    data = np.stack([_tensor(2) + index for index in range(3)])
    line = build_airmt_line(
        "L001",
        nav,
        data,
        frequency=[24.0, 75.0],
        record_mask=[True, False, True],
    )
    assert line.n_samples == 3
    assert line.n_records == 2
    assert line.missing_sample_ids == ("P002",)
    assert line.get_record("P001").emtf.product_id == "L001.P001"
    assert line.get_record("P003").attrs["navigation_index"] == 2


def test_build_airmt_line_supports_per_sample_frequency_grids():
    nav = NavigationTrack(sample_ids=("A", "B"))
    data = np.stack([_tensor(2), _tensor(2) + 10.0])
    freq = np.array([[24.0, 75.0], [30.0, 90.0]])
    line = build_airmt_line("L002", nav, data, frequency=freq)
    assert np.allclose(line.get_record("A").emtf.frequency, freq[0])
    assert np.allclose(line.get_record("B").emtf.frequency, freq[1])


def test_build_original_afmag_line_preserves_missing_samples():
    nav = NavigationTrack(sample_ids=("A", "B", "C"))
    tilt = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    line = build_original_afmag_line(
        "H001",
        nav,
        tilt,
        frequency=[150.0, 510.0],
        record_mask=[True, False, True],
    )
    assert line.n_records == 2
    assert line.missing_sample_ids == ("B",)
    tf = line.get_record("C").emtf.get_transfer_function(AFMAG_TILT_TAG)
    assert np.allclose(tf.data[:, 0, 0], [5.0, 6.0])


def test_afmag_datasets_reuse_common_airborne_containers():
    nav = NavigationTrack(sample_ids=("A",))
    airmt_line = build_airmt_line(
        "T001",
        nav,
        _tensor(2),
        frequency=[24.0, 75.0],
    )
    airmt = build_airmt_dataset(
        "AIRMT_SURVEY",
        [airmt_line],
        instrument_serial="A-TEST",
    )
    assert airmt.method == "AEM"
    assert airmt.instrument.system == "AirMt / tensor AFMAG"
    assert airmt.attrs["technology"] == "AirMt"

    original_line = build_original_afmag_line(
        "H001",
        nav,
        [1.0, 2.0],
        frequency=[150.0, 510.0],
    )
    original = build_original_afmag_dataset(
        "AFMAG_SURVEY",
        [original_line],
    )
    assert original.instrument.system == "AFMAG (original comparator)"
    assert original.attrs["technology"] == "AFMAG"


def test_airmt_amplification_parameter_accepts_patent_3_by_3_form():
    tensor = np.eye(3, dtype=complex)
    ap = compute_airmt_amplification_parameter(tensor)
    assert np.allclose(ap, 1.0 + 0.0j)
