from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import AirborneEMLine, NavigationTrack
from pycsamt.airborne.mobilemt import (
    MobileMTReferenceStation,
    MobileMTSystemSpec,
    MobileMTValidationError,
    build_mobilemt_dataset,
    build_mobilemt_emtf,
    build_mobilemt_line,
    build_mobilemt_record,
    validate_mobilemt_transfer_function,
)
from pycsamt.airborne.mobilemt.adapter import (
    _as_frequency,
    _normalize_estimate,
    _reference_mapping,
    _sample_axis_array,
    _xml_notes_mapping,
)
from pycsamt.emtf import TransferFunction
from pycsamt.metadata import LocationMeta, SiteMeta


def _admittance(nf: int = 3) -> np.ndarray:
    base = np.arange(nf * 6, dtype=float).reshape(nf, 3, 2)
    return base + 1j * (base + 0.5)


# ─────────────────────────────────────────────────────────────────────────
# _as_frequency
# ─────────────────────────────────────────────────────────────────────────


def test_as_frequency_scalar_promotion():
    arr = _as_frequency(100.0)
    assert arr.shape == (1,)


def test_as_frequency_rejects_bad_ndim_empty_and_values():
    with pytest.raises(MobileMTValidationError):
        _as_frequency([[1.0, 2.0]])
    with pytest.raises(MobileMTValidationError):
        _as_frequency([])
    with pytest.raises(MobileMTValidationError):
        _as_frequency([1.0, -1.0])
    with pytest.raises(MobileMTValidationError):
        _as_frequency([1.0, float("nan")])


# ─────────────────────────────────────────────────────────────────────────
# _normalize_estimate
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_estimate_promotes_2d():
    arr = _normalize_estimate(
        np.ones((3, 2)), n_frequency=1, tail=(3, 2), name="VAR",
    )
    assert arr.shape == (1, 3, 2)


def test_normalize_estimate_rejects_bad_shape_and_dtype():
    with pytest.raises(MobileMTValidationError):
        _normalize_estimate(
            np.ones((2, 3, 2)), n_frequency=1, tail=(3, 2), name="VAR",
        )
    with pytest.raises(MobileMTValidationError):
        _normalize_estimate(
            np.array([[["a", "b"]]]), n_frequency=1, tail=(1, 2), name="VAR",
        )


# ─────────────────────────────────────────────────────────────────────────
# _reference_mapping / _xml_notes_mapping
# ─────────────────────────────────────────────────────────────────────────


def test_reference_mapping_includes_attrs_when_present():
    ref = MobileMTReferenceStation(
        station_id="BASE01", attrs={"note": "fixed"},
    )
    out = _reference_mapping(ref)
    assert out["attrs"] == {"note": "fixed"}


def test_xml_notes_mapping_includes_full_location_metadata():
    ref = MobileMTReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(
                latitude=5.0, longitude=-3.0, elevation=120.0, datum="WGS84",
            ),
        ),
    )
    notes = _xml_notes_mapping(MobileMTSystemSpec(), ref)
    mobile = notes["MobileMT"]
    assert mobile["ReferenceLatitude"] == 5.0
    assert mobile["ReferenceLongitude"] == -3.0
    assert mobile["ReferenceElevation"] == 120.0
    assert mobile["ReferenceDatum"] == "WGS84"


def test_xml_notes_mapping_reference_station_without_id_or_site():
    ref = MobileMTReferenceStation()
    notes = _xml_notes_mapping(MobileMTSystemSpec(), ref)
    mobile = notes["MobileMT"]
    assert "ReferenceStationId" not in mobile
    assert "ReferenceLatitude" not in mobile


def test_xml_notes_mapping_reference_location_with_no_optional_fields():
    ref = MobileMTReferenceStation(
        site=SiteMeta(site_id="BASE02", location=LocationMeta(datum=None)),
    )
    notes = _xml_notes_mapping(MobileMTSystemSpec(), ref)
    mobile = notes["MobileMT"]
    assert "ReferenceLatitude" not in mobile
    assert "ReferenceLongitude" not in mobile
    assert "ReferenceElevation" not in mobile
    assert "ReferenceDatum" not in mobile


# ─────────────────────────────────────────────────────────────────────────
# validate_mobilemt_transfer_function
# ─────────────────────────────────────────────────────────────────────────


def test_validate_transfer_function_rejects_non_tf():
    with pytest.raises(TypeError):
        validate_mobilemt_transfer_function("not-a-tf")


def test_validate_transfer_function_rejects_wrong_name():
    tf = TransferFunction(
        name="T",
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0],
    )
    with pytest.raises(MobileMTValidationError):
        validate_mobilemt_transfer_function(tf)


def test_validate_transfer_function_rejects_wrong_channels_and_shape():
    from pycsamt.airborne.mobilemt.constants import MOBILEMT_ADMITTANCE_TAG

    bad_input = TransferFunction(
        name=MOBILEMT_ADMITTANCE_TAG,
        data=np.ones((1, 3, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hx", "Hy", "Hz"),
        periods=[1.0],
    )
    with pytest.raises(MobileMTValidationError):
        validate_mobilemt_transfer_function(bad_input)

    bad_output = TransferFunction(
        name=MOBILEMT_ADMITTANCE_TAG,
        data=np.ones((1, 3, 2), dtype=complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Ex", "Ey", "Ez"),
        periods=[1.0],
    )
    with pytest.raises(MobileMTValidationError):
        validate_mobilemt_transfer_function(bad_output)


def test_validate_transfer_function_rejects_mutated_bad_data_shape():
    from pycsamt.airborne.mobilemt.constants import MOBILEMT_ADMITTANCE_TAG

    tf = TransferFunction(
        name=MOBILEMT_ADMITTANCE_TAG,
        data=np.ones((1, 3, 2), dtype=complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Hx", "Hy", "Hz"),
        periods=[1.0],
    )
    # Bypass TransferFunction.validate() to simulate a data array whose
    # matrix shape no longer matches its own channel counts (never
    # producible through the public constructor).
    tf.data = np.ones((1, 3, 3), dtype=complex)
    with pytest.raises(MobileMTValidationError):
        validate_mobilemt_transfer_function(tf)


# ─────────────────────────────────────────────────────────────────────────
# build_mobilemt_emtf
# ─────────────────────────────────────────────────────────────────────────


def test_build_mobilemt_emtf_from_periods():
    doc = build_mobilemt_emtf(_admittance(2), periods=[0.1, 0.05])
    assert np.allclose(doc.frequency, [10.0, 20.0])


def test_build_mobilemt_emtf_rejects_frequency_count_mismatch():
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_emtf(_admittance(2), frequency=[10.0, 20.0, 30.0])


def test_build_mobilemt_emtf_rejects_non_numeric_admittance():
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_emtf(
            np.array([[["a", "b"], ["c", "d"], ["e", "f"]]]),
            frequency=[10.0],
        )


def test_build_mobilemt_emtf_rejects_bad_system_spec_type():
    with pytest.raises(TypeError):
        build_mobilemt_emtf(
            _admittance(1), frequency=[10.0], system_spec="not-a-spec",
        )


def test_build_mobilemt_emtf_rejects_bad_reference_station_type():
    with pytest.raises(TypeError):
        build_mobilemt_emtf(
            _admittance(1),
            frequency=[10.0],
            reference_station="not-a-reference-station",
        )


# ─────────────────────────────────────────────────────────────────────────
# build_mobilemt_record: apparent_conductivity validation
# ─────────────────────────────────────────────────────────────────────────


def test_build_mobilemt_record_apparent_conductivity_scalar_promotion():
    record = build_mobilemt_record(
        "S1", _admittance(1), frequency=[10.0], apparent_conductivity=5.0,
    )
    from pycsamt.airborne.mobilemt.constants import (
        MOBILEMT_APPARENT_CONDUCTIVITY_FIELD,
    )

    assert record.fields[MOBILEMT_APPARENT_CONDUCTIVITY_FIELD].shape == (1,)


def test_build_mobilemt_record_apparent_conductivity_size_mismatch():
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_record(
            "S1",
            _admittance(2),
            frequency=[10.0, 20.0],
            apparent_conductivity=[1.0, 2.0, 3.0],
        )


def test_build_mobilemt_record_apparent_conductivity_rejects_infinite():
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_record(
            "S1",
            _admittance(1),
            frequency=[10.0],
            apparent_conductivity=[float("inf")],
        )


# ─────────────────────────────────────────────────────────────────────────
# _sample_axis_array
# ─────────────────────────────────────────────────────────────────────────


def test_sample_axis_array_promotes_single_sample():
    arr = _sample_axis_array(
        np.ones((3, 2)), name="admittance", n_samples=1, tail_ndim=2,
    )
    assert arr.shape == (1, 3, 2)


def test_sample_axis_array_rejects_bad_shape():
    with pytest.raises(MobileMTValidationError):
        _sample_axis_array(
            np.ones((3, 2)), name="admittance", n_samples=2, tail_ndim=2,
        )


# ─────────────────────────────────────────────────────────────────────────
# build_mobilemt_line
# ─────────────────────────────────────────────────────────────────────────


def _nav(n: int = 2) -> NavigationTrack:
    return NavigationTrack(
        sample_ids=tuple(f"S{i:03d}" for i in range(n)),
        latitude=[10.0 + 0.001 * i for i in range(n)],
        longitude=[20.0 + 0.001 * i for i in range(n)],
    )


def test_build_mobilemt_line_rejects_non_navigation():
    with pytest.raises(TypeError):
        build_mobilemt_line(
            "L1", "not-a-navigation-track", _admittance(1), frequency=[10.0],
        )


def test_build_mobilemt_line_rejects_bad_admittance_tail_shape():
    nav = _nav(1)
    bad = np.ones((1, 3, 2, 2))  # wrong tail shape
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_line(
            "L1", nav, bad, frequency=[10.0, 20.0, 30.0],
        )


def test_build_mobilemt_line_rejects_shared_frequency_length_mismatch():
    nav = _nav(1)
    data = np.stack([_admittance(3)])
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_line("L1", nav, data, frequency=[10.0, 20.0])


def test_build_mobilemt_line_rejects_bad_frequency_shape():
    nav = _nav(2)
    data = np.stack([_admittance(2), _admittance(2)])
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_line(
            "L1", nav, data, frequency=np.ones((2, 2, 2)),
        )


def test_build_mobilemt_line_rejects_bad_record_mask_shape():
    nav = _nav(2)
    data = np.stack([_admittance(2), _admittance(2)])
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_line(
            "L1", nav, data, frequency=[10.0, 20.0],
            record_mask=[True, False, True],
        )


def test_build_mobilemt_line_apparent_conductivity_single_sample_promotion():
    nav = _nav(1)
    data = np.stack([_admittance(2)])
    line = build_mobilemt_line(
        "L1", nav, data, frequency=[10.0, 20.0],
        apparent_conductivity=[1.0, 2.0],
    )
    assert line.n_records == 1


def test_build_mobilemt_line_apparent_conductivity_shape_mismatch():
    nav = _nav(2)
    data = np.stack([_admittance(2), _admittance(2)])
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_line(
            "L1", nav, data, frequency=[10.0, 20.0],
            apparent_conductivity=np.ones((2, 3)),
        )


def test_build_mobilemt_line_optional_sample_shape_mismatch():
    nav = _nav(2)
    data = np.stack([_admittance(2), _admittance(2)])
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_line(
            "L1", nav, data, frequency=[10.0, 20.0],
            variance=np.ones((2, 3, 3, 2)),
        )


def test_build_mobilemt_line_per_sample_frequency_rows():
    nav = _nav(2)
    data = np.stack([_admittance(2), _admittance(2)])
    freq_rows = np.array([[10.0, 20.0], [11.0, 21.0]])
    line = build_mobilemt_line("L1", nav, data, frequency=freq_rows)
    assert line.n_records == 2


# ─────────────────────────────────────────────────────────────────────────
# build_mobilemt_dataset
# ─────────────────────────────────────────────────────────────────────────


def test_build_mobilemt_dataset_rejects_bad_survey_type():
    line = build_mobilemt_line(
        "L1", _nav(1), np.stack([_admittance(1)]), frequency=[10.0],
    )
    with pytest.raises(TypeError):
        build_mobilemt_dataset("SURVEY", [line], survey="not-a-survey")


def test_build_mobilemt_dataset_rejects_bad_system_spec_type():
    line = build_mobilemt_line(
        "L1", _nav(1), np.stack([_admittance(1)]), frequency=[10.0],
    )
    with pytest.raises(TypeError):
        build_mobilemt_dataset("SURVEY", [line], system_spec="not-a-spec")


def test_build_mobilemt_dataset_rejects_non_line_items():
    with pytest.raises(TypeError):
        build_mobilemt_dataset("SURVEY", ["not-a-line"])


def test_build_mobilemt_dataset_rejects_mismatched_technology_tag():
    nav = _nav(1)
    line = AirborneEMLine(
        line_id="L1", navigation=nav, attrs={"technology": "ZTEM"},
    )
    with pytest.raises(MobileMTValidationError):
        build_mobilemt_dataset("SURVEY", [line])


def test_build_mobilemt_dataset_accepts_mapping_of_lines():
    line = build_mobilemt_line(
        "L1", _nav(1), np.stack([_admittance(1)]), frequency=[10.0],
    )
    dataset = build_mobilemt_dataset("SURVEY", {"L1": line})
    assert dataset.n_lines == 1
