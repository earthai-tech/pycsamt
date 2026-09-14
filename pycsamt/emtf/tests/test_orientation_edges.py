from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.emtf import EMTF, StatisticalEstimate, TransferFunction
from pycsamt.emtf.orientation import (
    EMTFRotationError,
    EMTFRotationWarning,
    RotationMatrices,
    _angle,
    _angle_vector,
    _channel_pair_indices,
    _handle_derived,
    _handle_variance_without_full_covariance,
    _layout_channels,
    _masked_bilinear,
    _require_horizontal_pair,
    _source_angles_for_document,
    _source_override_for_tf,
    _variance_from_covariance,
    rotate_covariance,
    rotate_emtf,
    rotate_transfer_function,
)
from pycsamt.metadata import ChannelMeta, OrientationMeta, SiteLayout


def _tf(*, n_input: int = 2) -> TransferFunction:
    inputs = tuple(f"H{x}" for x in range(n_input))
    return TransferFunction(
        name="impedance",
        data=np.ones((1, 2, n_input), dtype=complex),
        input_channels=inputs,
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )


def _matrices() -> RotationMatrices:
    q = np.array([[[0.0, 1.0], [-1.0, 0.0]]])
    return RotationMatrices(
        input_matrix=q,
        output_matrix=q,
        source_mode="orthogonal",
        target_mode="orthogonal",
        source_angles=np.array([0.0]),
        target_angle=90.0,
    )


@pytest.mark.parametrize("value", ["bad", None, np.nan, np.inf])
def test_angle_rejects_nonfinite_values(value):
    with pytest.raises(EMTFRotationError, match="finite angle"):
        _angle(value, name="angle")


def test_angle_vector_rejects_bad_length_and_nonfinite_values():
    with pytest.raises(EMTFRotationError, match="scalar or length 2"):
        _angle_vector([1, 2, 3], 2, name="angles")
    with pytest.raises(EMTFRotationError, match="finite angles"):
        _angle_vector([1, np.nan], 2, name="angles")


def test_masked_bilinear_preserves_real_dtype_and_localizes_nan():
    data = np.array([[1.0, np.nan], [2.0, 3.0]])
    out = _masked_bilinear(np.eye(2), data, np.eye(2))
    assert out.dtype.kind == "f"
    assert out[0, 0] == 1.0
    assert np.isnan(out[0, 1])
    assert out[1, 1] == 3.0


def test_channel_pair_discovery_rejects_unpaired_horizontal_channel():
    assert _channel_pair_indices(("Hz",)) == []
    assert _channel_pair_indices(("Ex", "Ey")) == [(0, 1)]
    assert _channel_pair_indices(("Hx", "Hy", "Hz")) == [(0, 1)]
    with pytest.raises(EMTFRotationError, match="complete x/y"):
        _channel_pair_indices(("Hx", "Hz", "Aux"))


def test_layout_channel_and_horizontal_pair_validation():
    layout = SiteLayout(input_channels=[ChannelMeta("Hx", "magnetic")])
    with pytest.raises(EMTFRotationError, match="no input channel"):
        _layout_channels(layout, ("Hy",), role="input")
    with pytest.raises(EMTFRotationError, match="exactly two"):
        _require_horizontal_pair([], label="input")
    with pytest.raises(EMTFRotationError, match="no orientation"):
        _require_horizontal_pair(
            [ChannelMeta("Hx", "magnetic"), ChannelMeta("Hy", "magnetic", orientation=90)],
            label="input",
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"source_mode": "bad", "source_angles": 0}, "source_mode"),
        ({"source_mode": "orthogonal", "target_mode": "bad", "source_angles": 0}, "target_mode"),
        ({"source_mode": "orthogonal", "source_angles": None}, "source angle metadata"),
        ({"source_mode": "orthogonal", "source_angles": 0, "target_angle": None}, "requires target_angle"),
    ],
)
def test_transfer_rotation_rejects_invalid_frame_configuration(kwargs, message):
    with pytest.raises(EMTFRotationError, match=message):
        rotate_transfer_function(_tf(), **kwargs)


def test_transfer_rotation_requires_two_inputs_and_transfer_function_type():
    with pytest.raises(TypeError, match="TransferFunction"):
        rotate_transfer_function(object(), source_mode="orthogonal", source_angles=0)
    with pytest.raises(EMTFRotationError, match="exactly two input"):
        rotate_transfer_function(_tf(n_input=1), source_mode="orthogonal", source_angles=0)


def test_site_layout_modes_require_layout_and_identity_copy():
    tf = _tf()
    with pytest.raises(EMTFRotationError, match="SiteLayout metadata"):
        rotate_transfer_function(tf, source_mode="sitelayout", target_angle=0)

    copied = rotate_transfer_function(
        tf,
        source_mode="sitelayout",
        target_mode="sitelayout",
        target_angle=None,
        site_layout=SiteLayout(),
    )
    np.testing.assert_array_equal(copied.data, tf.data)
    assert copied is not tf


def test_covariance_rotation_validates_side_and_shapes():
    with pytest.raises(ValueError, match="side"):
        rotate_covariance(np.eye(2)[None], np.eye(2)[None], side="middle")
    with pytest.raises(EMTFRotationError, match="shape"):
        rotate_covariance(np.ones((2, 2)), np.eye(2)[None], side="input")
    with pytest.raises(EMTFRotationError, match="shape mismatch"):
        rotate_covariance(np.eye(2)[None], np.repeat(np.eye(2)[None], 2, axis=0), side="input")
    out = rotate_covariance(np.eye(2, dtype=float)[None], np.eye(2)[None], side="input")
    assert out.dtype.kind == "f"


def test_covariance_variance_warns_for_imaginary_and_negative_products():
    inverse = np.array([[[1 + 1j, 0], [0, -1 + 0j]]])
    residual = np.array([[[1 + 1j, 0], [0, 1 + 0j]]])
    with pytest.warns(EMTFRotationWarning) as caught:
        variance = _variance_from_covariance(inverse, residual)
    assert len(caught) == 2
    assert np.any(variance < 0)

    tiny = _variance_from_covariance(
        np.array([[[-1.0e-13, 0], [0, 1.0]]]),
        np.eye(2)[None],
    )
    assert tiny[0, 0, 0] == 0.0


def test_variance_policy_validation_and_shape_errors():
    matrices = _matrices()
    with pytest.raises(ValueError, match="variance_policy"):
        _handle_variance_without_full_covariance(np.ones((1, 2, 2)), matrices, policy="bad")
    with pytest.raises(EMTFRotationError, match="cannot be rotated exactly"):
        _handle_variance_without_full_covariance(np.ones((1, 2, 2)), matrices, policy="raise")
    with pytest.raises(EMTFRotationError, match="variance shape"):
        _handle_variance_without_full_covariance(np.ones((1, 1, 1)), matrices, policy="independent")
    with pytest.raises(EMTFRotationError, match="variance shape"):
        _handle_variance_without_full_covariance(np.ones((1, 1, 1)), matrices, policy="fcu")


@pytest.mark.parametrize(
    ("kind", "shape", "message"),
    [
        ("inverse_signal_covariance", (1, 1, 1), "INVSIGCOV shape"),
        ("residual_covariance", (1, 1, 1), "RESIDCOV shape"),
    ],
)
def test_transfer_rotation_rejects_covariance_shape_mismatch(kind, shape, message):
    tf = _tf()
    tf.add_estimate(StatisticalEstimate(name=kind, kind=kind, data=np.ones(shape)))
    with pytest.raises(EMTFRotationError, match=message):
        rotate_transfer_function(tf, source_mode="orthogonal", source_angles=0, target_angle=20)


def test_transfer_rotation_rejects_unknown_estimate_policy():
    tf = _tf()
    with pytest.raises(ValueError, match="unsupported_estimates"):
        rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0,
            target_angle=20,
            unsupported_estimates="maybe",
        )


def test_real_input_values_rotate_without_inventing_imaginary_components():
    tf = TransferFunction(
        name="impedance",
        data=np.array([[[1.0, 2.0], [3.0, 4.0]]]),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    out = rotate_transfer_function(
        tf, source_mode="orthogonal", source_angles=0, target_angle=30
    )
    np.testing.assert_allclose(out.data.imag, 0.0)


def test_source_angle_resolution_errors_and_case_insensitive_override():
    tf = _tf()
    doc = EMTF(periods=[1], transfer_functions={"impedance": tf})
    with pytest.raises(EMTFRotationError, match="ambiguous"):
        _source_angles_for_document(doc, tf, None, use_legacy_edi_rotation=False)
    with pytest.raises(EMTFRotationError, match="no orientation metadata"):
        _source_angles_for_document(doc, tf, None, use_legacy_edi_rotation=True)

    doc.orientation = OrientationMeta(mode="orthogonal")
    with pytest.raises(EMTFRotationError, match="no angle metadata"):
        _source_angles_for_document(doc, tf, None, use_legacy_edi_rotation=False)

    doc.orientation = OrientationMeta(extra={})
    with pytest.raises(EMTFRotationError, match="no historical EDI"):
        _source_angles_for_document(doc, tf, None, use_legacy_edi_rotation=True)
    assert _source_override_for_tf({"IMPEDANCE": 17}, tf) == 17
    assert _source_override_for_tf({"other": 17}, tf) is None

    tipper = TransferFunction(
        name="tipper",
        data=np.ones((1, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0],
    )
    doc.orientation = OrientationMeta(extra={"edi_zrot": [12.0]})
    with pytest.warns(match="EDI_ZROT"):
        mode, angles = _source_angles_for_document(
            doc, tipper, None, use_legacy_edi_rotation=True
        )
    assert mode == "orthogonal" and angles == [12.0]


def test_derived_and_document_policy_validation():
    tf = _tf()
    with pytest.raises(ValueError, match="derived_policy"):
        _handle_derived(tf, policy="bad")
    with pytest.raises(EMTFRotationError, match="recomputed"):
        _handle_derived(tf, policy="raise")
    with pytest.raises(TypeError, match="document must be an EMTF"):
        rotate_emtf(SimpleNamespace())
    with pytest.raises(ValueError, match="target must"):
        rotate_emtf(EMTF(), target="diagonal")
