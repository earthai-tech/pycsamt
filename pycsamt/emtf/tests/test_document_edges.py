from __future__ import annotations

import numpy as np
import pytest

from pycsamt.core.base import TFBundle
from pycsamt.emtf import EMTF, TransferFunction
from pycsamt.emtf.document import (
    _as_frequency,
    _legacy_scalar_or_tensor,
    _tipper_to_matrix,
)


def _matrix_tf(
    name: str = "impedance", *, periods=(1.0,), shape=(1, 2, 2)
) -> TransferFunction:
    return TransferFunction(
        name=name,
        data=np.ones(shape, dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey") if shape[1] == 2 else ("E",),
        periods=periods,
    )


def test_frequency_normalization_and_validation():
    assert _as_frequency(None) is None
    np.testing.assert_allclose(_as_frequency(10.0), [10.0])
    np.testing.assert_allclose(_as_frequency([]), [])
    with pytest.raises(ValueError, match="1-D"):
        _as_frequency([[1.0]])
    for value in ([0.0], [np.nan], [np.inf]):
        with pytest.raises(ValueError, match="finite positive"):
            _as_frequency(value)


@pytest.mark.parametrize(
    ("value", "shape"),
    [
        (np.array([1.0, 2.0]), (1, 1, 2)),
        (np.ones((3, 2)), (3, 1, 2)),
        (np.ones((3, 1, 2)), (3, 1, 2)),
    ],
)
def test_legacy_tipper_shape_normalization(value, shape):
    assert _tipper_to_matrix(value).shape == shape


def test_legacy_tipper_rejects_ambiguous_shapes():
    with pytest.raises(ValueError, match="legacy tipper"):
        _tipper_to_matrix(np.ones((2, 2, 2)))


@pytest.mark.parametrize(
    ("value", "shape"),
    [
        (3.0, (1, 1, 1)),
        (np.array([1.0, 2.0]), (2, 1, 1)),
        (np.eye(2), (1, 2, 2)),
        (np.ones((3, 2, 2)), (3, 2, 2)),
    ],
)
def test_legacy_derived_data_normalization(value, shape):
    tf = _legacy_scalar_or_tensor("apparent_resistivity", value, periods=None)
    assert tf.shape == shape


def test_legacy_derived_data_rejects_unsupported_matrix():
    with pytest.raises(ValueError, match="legacy impedance_phase"):
        _legacy_scalar_or_tensor(
            "impedance_phase", np.ones((2, 3)), periods=None
        )


def test_empty_document_accessors_and_period_counts():
    doc = EMTF()
    assert doc.frequency is None
    assert doc.n_periods == 0
    assert doc.get_transfer_function("") is None
    assert doc.z is None and doc.z_err is None and doc.Z is None
    assert doc.tipper is None and doc.tipper_err is None and doc.Tip is None
    assert doc.rho is None and doc.phase is None
    assert doc.is_empty()

    tf = _matrix_tf(periods=(1.0, 2.0), shape=(2, 2, 2))
    tf.periods = None
    doc.transfer_functions["impedance"] = tf
    assert doc.n_periods == 2


def test_add_transfer_function_validates_type_key_duplicate_and_periods():
    doc = EMTF()
    with pytest.raises(TypeError, match="TransferFunction"):
        doc.add_transfer_function(object())

    tf = _matrix_tf()
    with pytest.raises(ValueError, match="key must be non-empty"):
        doc.add_transfer_function(tf, key=" ")
    doc.add_transfer_function(tf)
    with pytest.raises(ValueError, match="already exists"):
        doc.add_transfer_function(tf)

    adopted = EMTF()
    adopted.add_transfer_function(_matrix_tf(periods=(2.0,)))
    np.testing.assert_allclose(adopted.periods, [2.0])


def test_transfer_function_lookup_accepts_registered_codes_and_aliases():
    tf = _matrix_tf(name="Z")
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    assert doc.get_transfer_function("Z") is tf
    assert doc.get_transfer_function("IMPEDANCE") is tf
    assert doc.get_transfer_function("not-a-type") is None


def test_stored_scalar_and_matrix_derived_accessors():
    scalar = TransferFunction(
        name="apparent_resistivity", data=[10.0, 20.0], periods=[1.0, 2.0]
    )
    matrix = _matrix_tf(name="impedance_phase", periods=(1.0, 2.0), shape=(2, 2, 2))
    doc = EMTF(periods=[1.0, 2.0])
    doc.add_transfer_function(scalar).add_transfer_function(matrix)
    np.testing.assert_allclose(doc.rho, [10.0, 20.0])
    assert doc.phase.shape == (2, 2, 2)

    bundle = doc.to_bundle()
    np.testing.assert_allclose(bundle.rho, [10.0, 20.0])
    assert bundle.phase.shape == (2, 2, 2)


def test_z_compatibility_rejects_non_impedance_matrix_shape():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex",),
        periods=[1.0],
    )
    assert EMTF(periods=[1.0], transfer_functions={"impedance": tf}).Z is None


def test_write_rejects_unknown_format(tmp_path):
    with pytest.raises(ValueError, match="unsupported EMTF output format"):
        EMTF().write(tmp_path / "out.dat", format="binary")


def test_from_bundle_validates_type_and_accepts_missing_frequency():
    with pytest.raises(TypeError, match="TFBundle"):
        EMTF.from_bundle(object())
    doc = EMTF.from_bundle(TFBundle())
    assert doc.periods is None
    assert doc.is_empty()
