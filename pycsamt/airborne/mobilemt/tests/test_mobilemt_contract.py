from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import NavigationTrack
from pycsamt.airborne.mobilemt import (
    MOBILEMT_ADMITTANCE_TAG,
    MOBILEMT_APPARENT_CONDUCTIVITY_FIELD,
    MOBILEMT_INPUT_CHANNELS,
    MOBILEMT_OUTPUT_CHANNELS,
    MobileMTReferenceStation,
    MobileMTSystemSpec,
    MobileMTValidationError,
    build_mobilemt_dataset,
    build_mobilemt_emtf,
    build_mobilemt_line,
    build_mobilemt_record,
    validate_mobilemt_transfer_function,
)
from pycsamt.emtf import get_emtf_datatype
from pycsamt.metadata import LocationMeta, SiteMeta


def _admittance(nf: int = 3) -> np.ndarray:
    base = np.arange(nf * 6, dtype=float).reshape(nf, 3, 2)
    return base + 1j * (base + 0.5)


def test_mobilemt_datatype_is_registered_explicitly():
    definition = get_emtf_datatype("Y")
    assert definition is not None
    assert definition.tag == MOBILEMT_ADMITTANCE_TAG
    assert definition.input_kind == "E"
    assert definition.output_kind == "H"


def test_system_spec_defaults_are_descriptive_not_data_limits():
    spec = MobileMTSystemSpec()
    assert spec.nominal_frequency_range_hz == (19.0, 26000.0)
    assert spec.nominal_max_frequency_windows == 30
    assert spec.nominal_sampling_rate_hz == 73728.0
    mask = spec.nominal_frequency_mask([10.0, 19.0, 26000.0, 30000.0])
    assert mask.tolist() == [False, True, True, False]


def test_reference_station_reuses_site_metadata():
    site = SiteMeta(
        site_id="BASE01",
        location=LocationMeta(latitude=5.2, longitude=-3.1),
    )
    ref = MobileMTReferenceStation(site=site)
    assert ref.preferred_id == "BASE01"
    assert ref.electric_channels == ("Ex", "Ey")


def test_build_mobilemt_emtf_constructs_3x2_tf():
    freq = np.array([20.0, 100.0, 1000.0])
    doc = build_mobilemt_emtf(
        _admittance(),
        frequency=freq,
        units="[nT]/[mV/km]",
        product_id="L1.P001",
    )
    tf = doc.get_transfer_function(MOBILEMT_ADMITTANCE_TAG)
    assert tf is not None
    assert tf.shape == (3, 3, 2)
    assert tf.input_channels == MOBILEMT_INPUT_CHANNELS
    assert tf.output_channels == MOBILEMT_OUTPUT_CHANNELS
    assert np.allclose(tf.frequency, freq)
    assert doc.subtype == "mobilemt"
    validate_mobilemt_transfer_function(tf)


def test_build_mobilemt_emtf_accepts_one_frequency_matrix():
    doc = build_mobilemt_emtf(
        _admittance(1)[0],
        frequency=[100.0],
    )
    tf = doc.get_transfer_function("Y")
    assert tf is not None
    assert tf.shape == (1, 3, 2)


def test_mobilemt_emtf_requires_exactly_one_period_axis():
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_emtf(_admittance())
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_emtf(
            _admittance(),
            frequency=[1.0, 2.0, 3.0],
            periods=[1.0, 0.5, 1.0 / 3.0],
        )


def test_bad_admittance_shape_is_rejected():
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_emtf(
            np.ones((3, 2, 2), dtype=complex),
            frequency=[1.0, 2.0, 3.0],
        )


def test_mobilemt_covariance_shapes_are_preserved():
    nf = 3
    var = np.ones((nf, 3, 2))
    inv = np.zeros((nf, 2, 2), dtype=complex)
    res = np.zeros((nf, 3, 3), dtype=complex)
    for i in range(nf):
        inv[i] = np.eye(2) * (i + 1)
        res[i] = np.eye(3) * (i + 2)
    doc = build_mobilemt_emtf(
        _admittance(nf),
        frequency=[20.0, 100.0, 1000.0],
        variance=var,
        inverse_signal_covariance=inv,
        residual_covariance=res,
    )
    tf = doc.get_transfer_function("Y")
    assert tf.get_estimate("VAR").shape == (nf, 3, 2)
    assert tf.get_estimate("INVSIGCOV").shape == (nf, 2, 2)
    assert tf.get_estimate("RESIDCOV").shape == (nf, 3, 3)


def test_missing_and_physical_zero_remain_distinct():
    data = _admittance(2)
    data[0, 0, 0] = 0.0 + 0.0j
    data[1, 2, 1] = np.nan + 1j * np.nan
    doc = build_mobilemt_emtf(data, frequency=[20.0, 40.0])
    tf = doc.get_transfer_function("Y")
    assert tf.data[0, 0, 0] == 0.0
    assert np.isnan(tf.data[1, 2, 1].real)
    assert np.isnan(tf.data[1, 2, 1].imag)


def test_apparent_conductivity_is_preserved_not_recomputed():
    sigma = np.array([12.0, 8.0, 5.0])
    record = build_mobilemt_record(
        "P001",
        _admittance(),
        frequency=[20.0, 100.0, 1000.0],
        apparent_conductivity=sigma,
    )
    assert np.array_equal(
        record.fields[MOBILEMT_APPARENT_CONDUCTIVITY_FIELD],
        sigma,
    )
    assert record.attrs["apparent_conductivity_units"] == "mS/m"


def test_build_mobilemt_line_shared_grid_and_sparse_records():
    nav = NavigationTrack(
        sample_ids=("P001", "P002", "P003"),
        latitude=(5.0, 5.001, 5.002),
        longitude=(-3.0, -3.001, -3.002),
    )
    data = np.stack([_admittance(2) + i for i in range(3)])
    line = build_mobilemt_line(
        "L001",
        nav,
        data,
        frequency=[20.0, 100.0],
        record_mask=[True, False, True],
    )
    assert line.n_samples == 3
    assert line.n_records == 2
    assert line.missing_sample_ids == ("P002",)
    assert line.get_record("P001").emtf.product_id == "L001.P001"
    assert line.get_record("P003").attrs["navigation_index"] == 2


def test_build_mobilemt_line_supports_per_sample_frequency_grids():
    nav = NavigationTrack(sample_ids=("A", "B"))
    data = np.stack([_admittance(2), _admittance(2) + 10])
    freq = np.array([[20.0, 100.0], [25.0, 125.0]])
    line = build_mobilemt_line(
        "L002",
        nav,
        data,
        frequency=freq,
    )
    assert np.allclose(line.get_record("A").emtf.frequency, freq[0])
    assert np.allclose(line.get_record("B").emtf.frequency, freq[1])


def test_build_mobilemt_line_propagates_covariance_and_conductivity():
    nav = NavigationTrack(sample_ids=("A", "B"))
    nf = 2
    data = np.stack([_admittance(nf), _admittance(nf) + 10])
    var = np.ones((2, nf, 3, 2))
    inv = np.tile(np.eye(2), (2, nf, 1, 1)).astype(complex)
    res = np.tile(np.eye(3), (2, nf, 1, 1)).astype(complex)
    sigma = np.array([[1.0, 2.0], [3.0, 4.0]])
    line = build_mobilemt_line(
        "L003",
        nav,
        data,
        frequency=[20.0, 100.0],
        variance=var,
        inverse_signal_covariance=inv,
        residual_covariance=res,
        apparent_conductivity=sigma,
    )
    rec = line.get_record("B")
    tf = rec.emtf.get_transfer_function("Y")
    assert tf.get_estimate("VAR").shape == (nf, 3, 2)
    assert tf.get_estimate("INVSIGCOV").shape == (nf, 2, 2)
    assert tf.get_estimate("RESIDCOV").shape == (nf, 3, 3)
    assert np.array_equal(
        rec.fields[MOBILEMT_APPARENT_CONDUCTIVITY_FIELD],
        sigma[1],
    )


def test_mobilemt_dataset_reuses_common_airborne_and_instrument_metadata():
    nav = NavigationTrack(sample_ids=("A", "B"))
    data = np.stack([_admittance(2), _admittance(2) + 2])
    line = build_mobilemt_line(
        "L010",
        nav,
        data,
        frequency=[20.0, 100.0],
    )
    dataset = build_mobilemt_dataset(
        "MOBILE_SURVEY",
        [line],
        instrument_serial="MMT-TEST",
    )
    assert dataset.method == "AEM"
    assert dataset.n_lines == 1
    assert dataset.instrument.system == "MobileMT"
    assert dataset.instrument.serial == "MMT-TEST"
    assert dataset.instrument.magnetic_sensor.sensor_type == "induction_coil"
    assert dataset.instrument.electric_sensor.sensor_type == "electrode"
    assert dataset.attrs["technology"] == "MobileMT"
