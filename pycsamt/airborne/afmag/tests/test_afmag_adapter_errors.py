from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import NavigationTrack
from pycsamt.airborne.afmag import (
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
    validate_airmt_transfer_function,
    validate_original_afmag_tilt,
)
from pycsamt.airborne.afmag.adapter import (
    _airmt_notes,
    _line_frequency_rows,
    _normalize_scalar_variance,
    _normalize_tilt,
)
from pycsamt.emtf import TransferFunction
from pycsamt.metadata import LocationMeta, OrientationMeta, SiteMeta


def _tensor(nf: int = 2) -> np.ndarray:
    base = np.arange(nf * 6, dtype=float).reshape(nf, 3, 2)
    return base + 1j * (base + 0.5)


def _nav(n: int = 2) -> NavigationTrack:
    return NavigationTrack(
        sample_ids=tuple(f"S{i:03d}" for i in range(n)),
        latitude=[10.0 + 0.001 * i for i in range(n)],
        longitude=[20.0 + 0.001 * i for i in range(n)],
    )


# ─────────────────────────────────────────────────────────────────────────
# _normalize_tilt
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_tilt_scalar_promotion():
    arr = _normalize_tilt(5.0)
    assert arr.shape == (1,)


def test_normalize_tilt_rejects_bad_ndim():
    with pytest.raises(AFMAGValidationError):
        _normalize_tilt([[1.0, 2.0]])


def test_normalize_tilt_rejects_non_numeric():
    with pytest.raises(AFMAGValidationError):
        _normalize_tilt(np.array(["a", "b"]))


# ─────────────────────────────────────────────────────────────────────────
# _normalize_scalar_variance
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_scalar_variance_scalar_promotion():
    arr = _normalize_scalar_variance(5.0, n_frequency=1)
    assert arr.shape == (1, 1, 1)


def test_normalize_scalar_variance_accepts_1x1_matrix():
    arr = _normalize_scalar_variance(np.ones((1, 1)), n_frequency=1)
    assert arr.shape == (1, 1, 1)


def test_normalize_scalar_variance_rejects_bad_shape():
    with pytest.raises(AFMAGValidationError):
        _normalize_scalar_variance(np.ones((3,)), n_frequency=2)


def test_normalize_scalar_variance_rejects_complex_with_imag():
    with pytest.raises(AFMAGValidationError):
        _normalize_scalar_variance(
            np.array([1.0 + 1.0j]), n_frequency=1,
        )


# ─────────────────────────────────────────────────────────────────────────
# compute_airmt_amplification_parameter
# ─────────────────────────────────────────────────────────────────────────


def test_amplification_parameter_rejects_bad_shape():
    with pytest.raises(AFMAGValidationError):
        compute_airmt_amplification_parameter(np.ones((2, 4)))


def test_amplification_parameter_rejects_non_numeric():
    with pytest.raises(AFMAGValidationError):
        compute_airmt_amplification_parameter(
            np.array([["a", "b"], ["c", "d"], ["e", "f"]])
        )


def test_amplification_parameter_rejects_bad_zero_policy():
    with pytest.raises(ValueError):
        compute_airmt_amplification_parameter(
            np.ones((3, 2)), zero_policy="bogus",
        )


# ─────────────────────────────────────────────────────────────────────────
# validate_airmt_transfer_function
# ─────────────────────────────────────────────────────────────────────────


def test_validate_airmt_tf_rejects_non_tf():
    with pytest.raises(TypeError):
        validate_airmt_transfer_function("not-a-tf")


def test_validate_airmt_tf_rejects_wrong_name():
    tf = TransferFunction(
        name="T",
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0],
    )
    with pytest.raises(AFMAGValidationError):
        validate_airmt_transfer_function(tf)


def test_validate_airmt_tf_rejects_wrong_channels():
    from pycsamt.airborne.afmag.constants import AFMAG_TENSOR_TAG

    bad_input = TransferFunction(
        name=AFMAG_TENSOR_TAG,
        data=np.ones((1, 3, 2), dtype=complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Hx", "Hy", "Hz"),
        periods=[1.0],
    )
    with pytest.raises(AFMAGValidationError):
        validate_airmt_transfer_function(bad_input)

    bad_output = TransferFunction(
        name=AFMAG_TENSOR_TAG,
        data=np.ones((1, 3, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey", "Ez"),
        periods=[1.0],
    )
    with pytest.raises(AFMAGValidationError):
        validate_airmt_transfer_function(bad_output)


def test_validate_airmt_tf_rejects_mutated_bad_shape():
    from pycsamt.airborne.afmag.constants import AFMAG_TENSOR_TAG

    tf = TransferFunction(
        name=AFMAG_TENSOR_TAG,
        data=np.ones((1, 3, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hx", "Hy", "Hz"),
        periods=[1.0],
    )
    tf.data = np.ones((1, 3, 3), dtype=complex)
    with pytest.raises(AFMAGValidationError):
        validate_airmt_transfer_function(tf)


# ─────────────────────────────────────────────────────────────────────────
# validate_original_afmag_tilt
# ─────────────────────────────────────────────────────────────────────────


def test_validate_original_afmag_tilt_rejects_non_tf():
    with pytest.raises(TypeError):
        validate_original_afmag_tilt("not-a-tf")


def test_validate_original_afmag_tilt_rejects_wrong_name():
    tf = TransferFunction(name="T", data=np.ones(2), periods=[1.0, 2.0])
    with pytest.raises(AFMAGValidationError):
        validate_original_afmag_tilt(tf)


def test_validate_original_afmag_tilt_rejects_nonscalar_channels():
    from pycsamt.airborne.afmag.constants import AFMAG_TILT_TAG

    tf = TransferFunction(
        name=AFMAG_TILT_TAG,
        data=np.ones((1, 1, 2)),
        input_channels=("Hx", "Hy"),
        periods=[1.0],
    )
    with pytest.raises(AFMAGValidationError):
        validate_original_afmag_tilt(tf)


def test_validate_original_afmag_tilt_rejects_mutated_bad_shape():
    from pycsamt.airborne.afmag.constants import AFMAG_TILT_TAG

    tf = TransferFunction(name=AFMAG_TILT_TAG, data=np.ones(2), periods=[1.0, 2.0])
    tf.data = np.ones((2, 1, 2))
    with pytest.raises(AFMAGValidationError):
        validate_original_afmag_tilt(tf)


def test_validate_original_afmag_tilt_rejects_mutated_complex():
    from pycsamt.airborne.afmag.constants import AFMAG_TILT_TAG

    tf = TransferFunction(name=AFMAG_TILT_TAG, data=np.ones(2), periods=[1.0, 2.0])
    tf.data = (np.ones((2, 1, 1)) + 1.0j).astype(complex)
    with pytest.raises(AFMAGValidationError):
        validate_original_afmag_tilt(tf)


# ─────────────────────────────────────────────────────────────────────────
# _airmt_notes: full reference-location metadata
# ─────────────────────────────────────────────────────────────────────────


def test_airmt_notes_includes_full_reference_location():
    ref = AFMAGReferenceStation(
        station_id="REF01",
        site=SiteMeta(
            site_id="REF01",
            location=LocationMeta(
                latitude=1.0, longitude=2.0, elevation=300.0, datum="WGS84",
            ),
        ),
    )
    notes = _airmt_notes(AirMtSystemSpec(), ref)
    afmag = notes["AFMAG"]
    assert afmag["ReferenceLatitude"] == 1.0
    assert afmag["ReferenceLongitude"] == 2.0
    assert afmag["ReferenceElevation"] == 300.0
    assert afmag["ReferenceDatum"] == "WGS84"


# ─────────────────────────────────────────────────────────────────────────
# build_airmt_emtf: type checks
# ─────────────────────────────────────────────────────────────────────────


def test_build_airmt_emtf_rejects_bad_system_spec():
    with pytest.raises(TypeError):
        build_airmt_emtf(_tensor(1), frequency=[10.0], system_spec="bad")


def test_build_airmt_emtf_rejects_bad_reference_station():
    with pytest.raises(TypeError):
        build_airmt_emtf(
            _tensor(1), frequency=[10.0], reference_station="bad",
        )


def test_build_airmt_emtf_rejects_bad_site():
    with pytest.raises(TypeError):
        build_airmt_emtf(_tensor(1), frequency=[10.0], site="bad")


def test_build_airmt_emtf_rejects_bad_orientation():
    with pytest.raises(TypeError):
        build_airmt_emtf(_tensor(1), frequency=[10.0], orientation="bad")


# ─────────────────────────────────────────────────────────────────────────
# build_original_afmag_emtf: type checks and value errors
# ─────────────────────────────────────────────────────────────────────────


def test_build_original_afmag_emtf_rejects_frequency_count_mismatch():
    with pytest.raises(AFMAGValidationError):
        build_original_afmag_emtf([1.0, 2.0, 3.0], frequency=[10.0, 20.0])


def test_build_original_afmag_emtf_rejects_bad_response_kind():
    with pytest.raises(ValueError):
        build_original_afmag_emtf(
            [1.0], frequency=[10.0], response_kind="bogus",
        )


def test_build_original_afmag_emtf_rejects_bad_system_spec():
    with pytest.raises(TypeError):
        build_original_afmag_emtf(
            [1.0], frequency=[10.0], system_spec="bad",
        )


def test_build_original_afmag_emtf_rejects_bad_site():
    with pytest.raises(TypeError):
        build_original_afmag_emtf([1.0], frequency=[10.0], site="bad")


def test_build_original_afmag_emtf_rejects_bad_orientation():
    with pytest.raises(TypeError):
        build_original_afmag_emtf(
            [1.0], frequency=[10.0], orientation="bad",
        )


# ─────────────────────────────────────────────────────────────────────────
# _line_frequency_rows
# ─────────────────────────────────────────────────────────────────────────


def test_line_frequency_rows_shared_length_mismatch():
    with pytest.raises(AFMAGValidationError):
        _line_frequency_rows([10.0], n_samples=2, n_frequency=2)


def test_line_frequency_rows_per_sample_grid_must_be_finite_positive():
    grid = np.array([[10.0, 20.0], [-1.0, 20.0]])
    with pytest.raises(AFMAGValidationError):
        _line_frequency_rows(grid, n_samples=2, n_frequency=2)


def test_line_frequency_rows_rejects_bad_shape():
    with pytest.raises(AFMAGValidationError):
        _line_frequency_rows(np.ones((2, 2, 2)), n_samples=2, n_frequency=2)


def test_line_frequency_rows_accepts_valid_per_sample_grid():
    grid = np.array([[10.0, 20.0], [11.0, 21.0]])
    common, rows = _line_frequency_rows(grid, n_samples=2, n_frequency=2)
    assert common is None
    assert rows.shape == (2, 2)


# ─────────────────────────────────────────────────────────────────────────
# build_airmt_line
# ─────────────────────────────────────────────────────────────────────────


def test_build_airmt_line_rejects_non_navigation():
    with pytest.raises(TypeError):
        build_airmt_line(
            "L1", "not-a-navigation", _tensor(1), frequency=[10.0],
        )


def test_build_airmt_line_single_sample_2d_promotion():
    nav = _nav(1)
    line = build_airmt_line("L1", nav, np.ones((3, 2), dtype=complex), frequency=[10.0])
    assert line.n_records == 1


def test_build_airmt_line_rejects_bad_ndim_shape():
    nav = _nav(2)
    with pytest.raises(AFMAGValidationError):
        build_airmt_line("L1", nav, np.ones((3, 2), dtype=complex), frequency=[10.0])


def test_build_airmt_line_rejects_bad_tail_shape():
    nav = _nav(2)
    bad = np.ones((2, 2, 3, 3), dtype=complex)
    with pytest.raises(AFMAGValidationError):
        build_airmt_line("L1", nav, bad, frequency=[10.0, 20.0])


def test_build_airmt_line_accepts_per_sample_variance():
    nav = _nav(2)
    tensor = np.stack([_tensor(2), _tensor(2)])
    line = build_airmt_line(
        "L1", nav, tensor, frequency=[10.0, 20.0],
        variance=np.ones((2, 2, 3, 2)),
    )
    assert line.n_records == 2


def test_build_airmt_line_rejects_non_numeric_tensor():
    nav = _nav(1)
    bad = np.array([[["a", "b"], ["c", "d"], ["e", "f"]]])
    with pytest.raises(AFMAGValidationError):
        build_airmt_line("L1", nav, bad, frequency=[10.0])


# ─────────────────────────────────────────────────────────────────────────
# build_original_afmag_line
# ─────────────────────────────────────────────────────────────────────────


def test_build_original_afmag_line_rejects_non_navigation():
    with pytest.raises(TypeError):
        build_original_afmag_line(
            "L1", "not-a-navigation", [1.0], frequency=[10.0],
        )


def test_build_original_afmag_line_scalar_single_sample_promotion():
    nav = _nav(1)
    line = build_original_afmag_line(
        "L1", nav, np.array(5.0), frequency=[10.0],
    )
    assert line.n_records == 1


def test_build_original_afmag_line_rejects_bad_shape():
    nav = _nav(2)
    with pytest.raises(AFMAGValidationError):
        build_original_afmag_line(
            "L1", nav, np.ones((2, 2, 2)), frequency=[10.0, 20.0],
        )


def test_build_original_afmag_line_rejects_non_numeric():
    nav = _nav(1)
    with pytest.raises(AFMAGValidationError):
        build_original_afmag_line(
            "L1", nav, np.array([["a", "b"]]), frequency=[10.0, 20.0],
        )


def test_build_original_afmag_line_variance_single_sample_promotion():
    nav = _nav(1)
    line = build_original_afmag_line(
        "L1", nav, [[1.0, 2.0]], frequency=[10.0, 20.0],
        variance=[0.1, 0.2],
    )
    assert line.n_records == 1


def test_build_original_afmag_line_variance_shape_mismatch():
    nav = _nav(2)
    with pytest.raises(AFMAGValidationError):
        build_original_afmag_line(
            "L1", nav, [[1.0, 2.0], [3.0, 4.0]], frequency=[10.0, 20.0],
            variance=np.ones((2, 3)),
        )


# ─────────────────────────────────────────────────────────────────────────
# build_airmt_dataset / build_original_afmag_dataset
# ─────────────────────────────────────────────────────────────────────────


def test_build_airmt_dataset_rejects_bad_system_spec():
    line = build_airmt_line("L1", _nav(1), np.ones((3, 2), dtype=complex), frequency=[10.0])
    with pytest.raises(TypeError):
        build_airmt_dataset("SURVEY", [line], system_spec="bad")


def test_build_original_afmag_dataset_rejects_bad_system_spec():
    line = build_original_afmag_line("L1", _nav(1), [5.0], frequency=[10.0])
    with pytest.raises(TypeError):
        build_original_afmag_dataset("SURVEY", [line], system_spec="bad")
