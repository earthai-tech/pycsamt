from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import NavigationTrack
from pycsamt.airborne.ztem import (
    ZTEM_INPUT_CHANNELS,
    ZTEM_OUTPUT_CHANNELS,
    ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ,
    ZTEM_TIPPER_TAG,
    ZTEMReferenceStation,
    ZTEMSystemSpec,
    ZTEMValidationError,
    build_ztem_dataset,
    build_ztem_emtf,
    build_ztem_line,
    validate_ztem_transfer_function,
)
from pycsamt.emtf import get_emtf_datatype
from pycsamt.metadata import LocationMeta, OrientationMeta, SiteMeta


def _tipper(nf: int = 3) -> np.ndarray:
    base = np.arange(nf * 2, dtype=float).reshape(nf, 2)
    return base + 1j * (base + 0.25)


def test_ztem_reuses_standard_tipper_datatype():
    definition = get_emtf_datatype("T")
    assert definition is not None
    assert definition.tag == ZTEM_TIPPER_TAG == "tipper"
    assert definition.input_kind == "H"
    assert definition.output_kind == "H"


def test_ztem_system_spec_is_descriptive_not_grid_validation():
    spec = ZTEMSystemSpec()
    assert spec.practical_frequency_range_hz == (
        ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ
    )
    mask = spec.practical_frequency_mask([10.0, 30.0, 720.0, 900.0])
    assert np.array_equal(mask, [False, True, True, False])


def test_ztem_instrument_metadata_has_no_electric_sensor():
    instrument = ZTEMSystemSpec().to_instrument_meta(serial="ZT-TEST")
    assert instrument.system == "ZTEM"
    assert instrument.serial == "ZT-TEST"
    assert instrument.magnetic_sensor.sensor_type == "induction_coil"
    assert instrument.electric_sensor is None


def test_build_ztem_emtf_accepts_nf_by_two_tipper():
    freq = np.array([30.0, 90.0, 360.0])
    doc = build_ztem_emtf(_tipper(), frequency=freq)
    tf = doc.tipper_tf
    assert tf is not None
    assert tf.shape == (3, 1, 2)
    assert tf.input_channels == ZTEM_INPUT_CHANNELS
    assert tf.output_channels == ZTEM_OUTPUT_CHANNELS
    assert tf.units == "[]"
    assert np.allclose(tf.frequency, freq)
    assert doc.subtype == "ztem"
    validate_ztem_transfer_function(tf)


def test_build_ztem_emtf_accepts_single_tipper_vector():
    doc = build_ztem_emtf(
        np.array([1.0 + 2.0j, 3.0 + 4.0j]),
        frequency=[100.0],
    )
    assert doc.tipper_tf.shape == (1, 1, 2)


def test_ztem_requires_exactly_one_period_axis():
    with pytest.raises(ZTEMValidationError):
        build_ztem_emtf(_tipper())
    with pytest.raises(ZTEMValidationError):
        build_ztem_emtf(
            _tipper(),
            frequency=[30.0, 90.0, 360.0],
            periods=[1 / 30, 1 / 90, 1 / 360],
        )


def test_ztem_covariance_shapes_are_preserved():
    nf = 3
    var = np.ones((nf, 1, 2))
    inv = np.tile(np.eye(2), (nf, 1, 1)).astype(complex)
    res = np.ones((nf, 1, 1), dtype=complex) * 2.0
    doc = build_ztem_emtf(
        _tipper(nf),
        frequency=[30.0, 90.0, 360.0],
        variance=var,
        inverse_signal_covariance=inv,
        residual_covariance=res,
    )
    tf = doc.tipper_tf
    assert tf.get_estimate("VAR").shape == (nf, 1, 2)
    assert tf.get_estimate("INVSIGCOV").shape == (nf, 2, 2)
    assert tf.get_estimate("RESIDCOV").shape == (nf, 1, 1)


def test_ztem_missing_and_physical_zero_remain_distinct():
    data = _tipper(2)
    data[0, 0] = 0.0 + 0.0j
    data[1, 1] = np.nan + 1j * np.nan
    doc = build_ztem_emtf(data, frequency=[30.0, 90.0])
    tf = doc.tipper_tf
    assert tf.data[0, 0, 0] == 0.0
    assert np.isnan(tf.data[1, 0, 1].real)
    assert np.isnan(tf.data[1, 0, 1].imag)


def test_ztem_reference_station_maps_to_remote_reference_metadata():
    ref = ZTEMReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(latitude=5.2, longitude=-3.1),
        ),
    )
    doc = build_ztem_emtf(
        _tipper(2),
        frequency=[30.0, 90.0],
        reference_station=ref,
    )
    remote = doc.processing.remote_reference
    assert remote.reference_type == "fixed_ground_horizontal_magnetic"
    assert remote.site == "BASE01"
    assert remote.extra["channels"] == ["Hx", "Hy"]
    assert doc.metadata["notes"]["ZTEM"]["ReferenceStationId"] == "BASE01"


def test_ztem_does_not_infer_orientation():
    doc = build_ztem_emtf(_tipper(2), frequency=[30.0, 90.0])
    assert doc.orientation is None


def test_ztem_can_carry_explicit_orthogonal_orientation():
    orientation = OrientationMeta(
        mode="orthogonal",
        angle_to_geographic_north=0.0,
    )
    doc = build_ztem_emtf(
        _tipper(2),
        frequency=[30.0, 90.0],
        orientation=orientation,
    )
    assert doc.orientation is orientation


def test_build_ztem_line_shared_grid_and_sparse_records():
    nav = NavigationTrack(
        sample_ids=("P001", "P002", "P003"),
        latitude=(5.0, 5.001, 5.002),
        longitude=(-3.0, -3.001, -3.002),
    )
    data = np.stack([_tipper(2) + index for index in range(3)])
    line = build_ztem_line(
        "L001",
        nav,
        data,
        frequency=[30.0, 90.0],
        record_mask=[True, False, True],
    )
    assert line.n_samples == 3
    assert line.n_records == 2
    assert line.missing_sample_ids == ("P002",)
    assert line.get_record("P001").emtf.product_id == "L001.P001"
    assert line.get_record("P003").attrs["navigation_index"] == 2


def test_build_ztem_line_supports_per_sample_frequency_grids():
    nav = NavigationTrack(sample_ids=("A", "B"))
    data = np.stack([_tipper(2), _tipper(2) + 10])
    freq = np.array([[30.0, 90.0], [25.0, 75.0]])
    line = build_ztem_line("L002", nav, data, frequency=freq)
    assert np.allclose(line.get_record("A").emtf.frequency, freq[0])
    assert np.allclose(line.get_record("B").emtf.frequency, freq[1])


def test_build_ztem_line_propagates_covariance():
    nav = NavigationTrack(sample_ids=("A", "B"))
    nf = 2
    data = np.stack([_tipper(nf), _tipper(nf) + 10])
    var = np.ones((2, nf, 1, 2))
    inv = np.tile(np.eye(2), (2, nf, 1, 1)).astype(complex)
    res = np.ones((2, nf, 1, 1), dtype=complex)
    line = build_ztem_line(
        "L003",
        nav,
        data,
        frequency=[30.0, 90.0],
        variance=var,
        inverse_signal_covariance=inv,
        residual_covariance=res,
    )
    tf = line.get_record("B").emtf.tipper_tf
    assert tf.get_estimate("VAR").shape == (nf, 1, 2)
    assert tf.get_estimate("INVSIGCOV").shape == (nf, 2, 2)
    assert tf.get_estimate("RESIDCOV").shape == (nf, 1, 1)


def test_ztem_dataset_reuses_common_airborne_metadata():
    nav = NavigationTrack(sample_ids=("A", "B"))
    line = build_ztem_line(
        "L010",
        nav,
        np.stack([_tipper(2), _tipper(2) + 2]),
        frequency=[30.0, 90.0],
    )
    dataset = build_ztem_dataset(
        "ZTEM_SURVEY",
        [line],
        instrument_serial="ZT-TEST",
    )
    assert dataset.method == "AEM"
    assert dataset.n_lines == 1
    assert dataset.instrument.system == "ZTEM"
    assert dataset.instrument.serial == "ZT-TEST"
    assert dataset.instrument.electric_sensor is None
    assert dataset.attrs["technology"] == "ZTEM"


def test_build_ztem_line_single_sample_accepts_unbatched_tipper():
    nav = NavigationTrack(sample_ids=("A",))
    line = build_ztem_line(
        "L011",
        nav,
        _tipper(2),
        frequency=[30.0, 90.0],
    )
    assert line.n_records == 1
    assert line.get_record("A").emtf.tipper_tf.shape == (2, 1, 2)


def test_ztem_conflicting_remote_reference_metadata_is_rejected():
    from pycsamt.metadata import ProcessingMeta, RemoteReferenceMeta

    ref = ZTEMReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(
        remote_reference=RemoteReferenceMeta(
            reference_type="fixed_ground_horizontal_magnetic",
            site="OTHER_BASE",
        )
    )
    with pytest.raises(ZTEMValidationError):
        build_ztem_emtf(
            _tipper(2),
            frequency=[30.0, 90.0],
            reference_station=ref,
            processing=processing,
        )
