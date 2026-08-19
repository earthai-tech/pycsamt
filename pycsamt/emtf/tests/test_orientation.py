# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 7 format-neutral transfer-function rotation tests."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtf import (
    ApproximateVarianceRotationWarning,
    DerivedDataRotationWarning,
    EMTF,
    EMTFRotationError,
    LegacyRotationAssumptionWarning,
    StatisticalEstimate,
    TransferFunction,
    UnsupportedEstimateRotationWarning,
    horizontal_inverse_rotation_matrix,
    horizontal_rotation_matrix,
    rotate_transfer_function,
)
from pycsamt.metadata import (
    ChannelMeta,
    OrientationMeta,
    SiteLayout,
)


def _layout(*, nonorthogonal: bool = True) -> SiteLayout:
    input_angles = (-10.0, 75.0) if nonorthogonal else (0.0, 90.0)
    output_angles = (20.0, 112.0) if nonorthogonal else (0.0, 90.0)
    return SiteLayout(
        input_channels=[
            ChannelMeta("Hx", "magnetic", orientation=input_angles[0]),
            ChannelMeta("Hy", "magnetic", orientation=input_angles[1]),
        ],
        output_channels=[
            ChannelMeta("Hz", "magnetic", orientation=0.0, tilt=90.0),
            ChannelMeta("Ex", "electric", orientation=output_angles[0]),
            ChannelMeta("Ey", "electric", orientation=output_angles[1]),
        ],
    )


def _fcu_impedance_tf() -> TransferFunction:
    z = np.array(
        [
            [
                [1.0 + 0.5j, 2.0 - 0.2j],
                [-3.0 + 0.7j, 4.0 + 1.1j],
            ]
        ]
    )
    tf = TransferFunction(
        name="impedance",
        data=z,
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="INVSIGCOV",
            kind="inverse_signal_covariance",
            data=np.array(
                [
                    [
                        [2.0 + 0.0j, 0.3 + 0.1j],
                        [0.3 - 0.1j, 1.5 + 0.0j],
                    ]
                ]
            ),
        )
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="RESIDCOV",
            kind="residual_covariance",
            data=np.array(
                [
                    [
                        [4.0 + 0.0j, 0.4 - 0.2j],
                        [0.4 + 0.2j, 3.0 + 0.0j],
                    ]
                ]
            ),
        )
    )
    return tf


def _fcu_document() -> EMTF:
    return EMTF(
        periods=[1.0],
        transfer_functions={"impedance": _fcu_impedance_tf()},
        site_layout=_layout(),
        orientation=OrientationMeta(mode="sitelayout"),
    )


def test_general_horizontal_rotation_and_inverse_are_consistent():
    q = horizontal_rotation_matrix(-10.0, 75.0, 13.0)
    qinv = horizontal_inverse_rotation_matrix(-10.0, 75.0, 13.0)
    np.testing.assert_allclose(q @ qinv, np.eye(2), atol=1.0e-14)
    np.testing.assert_allclose(qinv @ q, np.eye(2), atol=1.0e-14)


def test_parallel_horizontal_channels_are_rejected():
    with pytest.raises(EMTFRotationError, match="parallel"):
        horizontal_inverse_rotation_matrix(10.0, 190.0, 0.0)


def test_identity_rotation_is_exact_and_keeps_estimates():
    tf = _fcu_impedance_tf()
    out = rotate_transfer_function(
        tf,
        source_mode="orthogonal",
        source_angles=20.0,
        target_mode="orthogonal",
        target_angle=20.0,
    )
    np.testing.assert_array_equal(out.data, tf.data)
    assert set(out.estimates) == set(tf.estimates)
    for key in tf.estimates:
        np.testing.assert_array_equal(
            out.estimates[key].data,
            tf.estimates[key].data,
        )


def test_orthogonal_90_degree_impedance_rotation():
    tf = TransferFunction(
        name="impedance",
        data=np.array([[[1.0, 2.0], [3.0, 4.0]]], dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    out = rotate_transfer_function(
        tf,
        source_mode="orthogonal",
        source_angles=0.0,
        target_mode="orthogonal",
        target_angle=90.0,
    )
    expected = np.array([[[4.0, -3.0], [-2.0, 1.0]]], dtype=complex)
    np.testing.assert_allclose(out.data, expected, atol=1.0e-14)


def test_fcu_site_layout_reference_impedance_and_covariance():
    """Match FCU v4.1 rotate_data for a non-orthogonal site layout."""
    doc = _fcu_document()
    out = doc.rotate(0.0)

    expected_z = np.array(
        [
            [
                [
                    2.0672192481485361 + 0.09672526991822134j,
                    -0.15951222533783982 - 0.6470891938195381j,
                ],
                [
                    -1.5996961731763872 + 0.9609793686924775j,
                    4.976373178597356 + 0.7275697909457486j,
                ],
            ]
        ]
    )
    expected_s = np.array(
        [
            [
                [2.027294299892441 + 0.0j, 0.02862495604151122 + 0.10038198339177593j],
                [0.028624956041511274 - 0.10038198339177591j, 1.4468019181314853 + 0.0j],
            ]
        ]
    )
    expected_n = np.array(
        [
            [
                [3.6714671959715126 + 0.0j, 0.5408454890921977 - 0.19987816509193782j],
                [0.540845489092198 + 0.19987816509193784j, 3.3006131709420594 + 0.0j],
            ]
        ]
    )
    expected_var = np.array(
        [[[7.4431445186351315, 5.31188578148841],
          [6.691314267600752, 4.775333466729015]]]
    )

    np.testing.assert_allclose(out.z, expected_z, rtol=5.0e-7, atol=5.0e-7)
    np.testing.assert_allclose(
        out.impedance.get_estimate("INVSIGCOV").data,
        expected_s,
        rtol=5.0e-7,
        atol=5.0e-7,
    )
    np.testing.assert_allclose(
        out.impedance.get_estimate("RESIDCOV").data,
        expected_n,
        rtol=5.0e-7,
        atol=5.0e-7,
    )
    np.testing.assert_allclose(
        out.impedance.get_estimate("VAR").data,
        expected_var,
        rtol=5.0e-7,
        atol=5.0e-7,
    )
    assert out.orientation.mode == "orthogonal"
    assert out.orientation.angle_to_geographic_north == pytest.approx(0.0)


def test_fcu_site_layout_roundtrip_is_reversible_with_full_covariance():
    doc = _fcu_document()
    north = doc.rotate(0.0)
    back = north.rotate(target="sitelayout")

    np.testing.assert_allclose(back.z, doc.z, rtol=2.0e-13, atol=2.0e-13)
    for key in ("INVSIGCOV", "RESIDCOV"):
        np.testing.assert_allclose(
            back.impedance.get_estimate(key).data,
            doc.impedance.get_estimate(key).data,
            rtol=2.0e-13,
            atol=2.0e-13,
        )
    assert back.orientation.mode == "sitelayout"
    assert back.orientation.angle_to_geographic_north is None


def test_fcu_tipper_reference_and_variance():
    tf = TransferFunction(
        name="tipper",
        data=np.array([[[0.2 + 0.1j, -0.3 + 0.05j]]]),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="INVSIGCOV",
            kind="inverse_signal_covariance",
            data=np.array(
                [[[2.0 + 0j, 0.3 + 0.1j], [0.3 - 0.1j, 1.5 + 0j]]]
            ),
        )
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="RESIDCOV",
            kind="residual_covariance",
            data=np.array([[[2.5 + 0j]]]),
        )
    )
    out = rotate_transfer_function(
        tf,
        source_mode="sitelayout",
        target_mode="orthogonal",
        target_angle=0.0,
        site_layout=_layout(),
    )
    expected = np.array(
        [[[0.14162965549711895 + 0.10567712568196848j,
           -0.3485323973769504 + 0.02344771215613209j]]]
    )
    expected_var = np.array([[[5.068235749731103, 3.617004795328713]]])
    np.testing.assert_allclose(out.data, expected, rtol=5.0e-7, atol=5.0e-7)
    np.testing.assert_allclose(
        out.get_estimate("VAR").data,
        expected_var,
        rtol=5.0e-7,
        atol=5.0e-7,
    )


def test_site_layout_metadata_is_not_rotated_or_mutated():
    doc = _fcu_document()
    before = [
        channel.orientation
        for channel in doc.site_layout.input_channels
        + doc.site_layout.output_channels
    ]
    out = doc.rotate(37.0)
    after_source = [
        channel.orientation
        for channel in doc.site_layout.input_channels
        + doc.site_layout.output_channels
    ]
    after_copy = [
        channel.orientation
        for channel in out.site_layout.input_channels
        + out.site_layout.output_channels
    ]
    assert before == after_source == after_copy
    assert out.site_layout is not doc.site_layout


def test_variance_without_covariance_is_dropped_by_default():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="VAR",
            kind="variance",
            data=np.ones((1, 2, 2)),
        )
    )
    with pytest.warns(ApproximateVarianceRotationWarning, match="dropping"):
        out = rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0.0,
            target_angle=30.0,
        )
    assert out.get_estimate("VAR") is None


def test_independent_variance_policy_is_nonnegative():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="VAR",
            kind="variance",
            data=np.array([[[1.0, 2.0], [3.0, 4.0]]]),
        )
    )
    with pytest.warns(ApproximateVarianceRotationWarning, match="independent"):
        out = rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0.0,
            target_angle=31.0,
            variance_policy="independent",
        )
    assert np.all(out.get_estimate("VAR").data >= 0.0)


def test_fcu_variance_policy_matches_direct_matrix_rotation():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    variance = np.array([[[1.0, 2.0], [3.0, 4.0]]])
    tf.add_estimate(
        StatisticalEstimate(name="VAR", kind="variance", data=variance)
    )
    q = horizontal_rotation_matrix(0.0, 90.0, 30.0)
    expected = np.array([q @ variance[0] @ q.T])
    with pytest.warns(ApproximateVarianceRotationWarning, match="FCU"):
        out = rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0.0,
            target_angle=30.0,
            variance_policy="fcu",
        )
    np.testing.assert_allclose(out.get_estimate("VAR").data, expected)


def test_custom_three_output_admittance_rotates_xy_and_leaves_z_output():
    tf = TransferFunction(
        name="mobilemt_admittance",
        data=np.array(
            [[[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]],
            dtype=complex,
        ),
        input_channels=("Ex", "Ey"),
        output_channels=("Hx", "Hy", "Hz"),
        periods=[1.0],
    )
    layout = SiteLayout(
        input_channels=[
            ChannelMeta("Ex", "electric", orientation=0.0),
            ChannelMeta("Ey", "electric", orientation=90.0),
        ],
        output_channels=[
            ChannelMeta("Hx", "magnetic", orientation=0.0),
            ChannelMeta("Hy", "magnetic", orientation=90.0),
            ChannelMeta("Hz", "magnetic", orientation=0.0, tilt=90.0),
        ],
    )
    out = rotate_transfer_function(
        tf,
        source_mode="sitelayout",
        target_angle=90.0,
        site_layout=layout,
    )
    q = horizontal_rotation_matrix(0.0, 90.0, 90.0)
    expected = np.eye(3)
    expected[:2, :2] = q
    expected = expected @ tf.data[0] @ q.T
    np.testing.assert_allclose(out.data[0], expected, atol=1.0e-14)


def test_missing_component_does_not_contaminate_identity_dependencies():
    data = np.array([[[1.0 + 0j, np.nan + 1j * np.nan],
                      [2.0 + 0j, 3.0 + 0j]]])
    tf = TransferFunction(
        name="impedance",
        data=data,
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    out = rotate_transfer_function(
        tf,
        source_mode="orthogonal",
        source_angles=0.0,
        target_angle=0.0,
    )
    assert out.data[0, 0, 0] == pytest.approx(1.0 + 0j)
    assert np.isnan(out.data[0, 0, 1].real)
    assert out.data[0, 1, 1] == pytest.approx(3.0 + 0j)


def test_tilted_horizontal_input_is_rejected_for_site_layout_rotation():
    layout = _layout()
    layout.input_channels[0].tilt = 5.0
    tf = _fcu_impedance_tf()
    with pytest.raises(EMTFRotationError, match="tilted input"):
        rotate_transfer_function(
            tf,
            source_mode="sitelayout",
            target_angle=0.0,
            site_layout=layout,
        )


def test_frequency_dependent_explicit_source_angles_normalize_principal_axis():
    data = np.array(
        [
            [[1.0, 2.0], [3.0, 4.0]],
            [[1.0, 2.0], [3.0, 4.0]],
        ],
        dtype=complex,
    )
    tf = TransferFunction(
        name="impedance",
        data=data,
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0, 10.0],
    )
    doc = EMTF(
        periods=[1.0, 10.0],
        transfer_functions={"impedance": tf},
        orientation=None,
    )
    out = doc.rotate(
        0.0,
        source_angles={"impedance": [0.0, 45.0]},
    )
    np.testing.assert_allclose(out.z[0], data[0], atol=1.0e-14)
    q = horizontal_rotation_matrix(45.0, 135.0, 0.0)
    np.testing.assert_allclose(out.z[1], q @ data[1] @ q.T)
    assert out.orientation.mode == "orthogonal"
    assert out.orientation.angle_to_geographic_north == pytest.approx(0.0)


def test_legacy_edi_rotation_requires_opt_in_and_is_normalized():
    tf = TransferFunction(
        name="impedance",
        data=np.repeat(np.eye(2, dtype=complex)[None, :, :], 2, axis=0),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0, 2.0],
    )
    doc = EMTF(
        periods=[1.0, 2.0],
        transfer_functions={"impedance": tf},
        orientation=OrientationMeta(
            extra={"edi_zrot": [10.0, 25.0], "edi_had_zrot": True}
        ),
    )
    with pytest.raises(EMTFRotationError, match="ambiguous"):
        doc.rotate(0.0)
    with pytest.warns(LegacyRotationAssumptionWarning):
        out = doc.rotate(0.0, use_legacy_edi_rotation=True)
    assert out.orientation.extra["edi_zrot"] == [0.0, 0.0]
    assert out.orientation.extra["edi_rotation_metadata_normalized"] is True


def test_unsupported_estimate_policy_drop_keep_and_raise():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="COH",
            kind="coherence",
            data=np.ones((1, 2, 2)),
        )
    )
    with pytest.warns(UnsupportedEstimateRotationWarning, match="dropping"):
        dropped = rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0.0,
            target_angle=10.0,
        )
    assert dropped.get_estimate("COH") is None

    with pytest.warns(UnsupportedEstimateRotationWarning, match="keeping"):
        kept = rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0.0,
            target_angle=10.0,
            unsupported_estimates="keep",
        )
    assert kept.get_estimate("COH").attrs["stale_after_rotation"] is True

    with pytest.raises(EMTFRotationError, match="COH"):
        rotate_transfer_function(
            tf,
            source_mode="orthogonal",
            source_angles=0.0,
            target_angle=10.0,
            unsupported_estimates="raise",
        )


def test_derived_products_are_dropped_by_default_and_can_be_marked_stale():
    periods = [1.0]
    ztf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=periods,
    )
    rho = TransferFunction(
        name="apparent_resistivity",
        data=np.ones((1, 2, 2)),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=periods,
    )
    doc = EMTF(
        periods=periods,
        transfer_functions={"impedance": ztf, "apparent_resistivity": rho},
        orientation=OrientationMeta(mode="orthogonal", angle_to_geographic_north=0.0),
    )
    with pytest.warns(DerivedDataRotationWarning, match="dropping"):
        out = doc.rotate(10.0)
    assert out.get_transfer_function("apparent_resistivity") is None

    with pytest.warns(DerivedDataRotationWarning, match="stale"):
        kept = doc.rotate(10.0, derived_policy="keep")
    assert kept.get_transfer_function("apparent_resistivity").attrs[
        "stale_after_rotation"
    ] is True


def test_inplace_rotation_retains_exact_site_layout_object():
    doc = _fcu_document()
    layout = doc.site_layout
    returned = doc.rotate(12.0, inplace=True)
    assert returned is doc
    assert doc.site_layout is layout
    assert doc.orientation.mode == "orthogonal"
    assert doc.orientation.angle_to_geographic_north == pytest.approx(12.0)


def test_rotated_xml_roundtrip_preserves_data_covariance_and_orientation():
    from pycsamt.emtf.tests.test_xml_reader import FULL_XML

    doc = EMTF.from_xml(FULL_XML)
    rotated = doc.rotate(23.0)
    reread = EMTF.from_xml(rotated.to_xml())

    assert reread.orientation.mode == "orthogonal"
    assert reread.orientation.angle_to_geographic_north == pytest.approx(23.0)
    np.testing.assert_allclose(
        reread.z,
        rotated.z,
        equal_nan=True,
        rtol=2.0e-14,
        atol=2.0e-14,
    )
    for name in ("INVSIGCOV", "RESIDCOV", "VAR"):
        np.testing.assert_allclose(
            reread.impedance.get_estimate(name).data,
            rotated.impedance.get_estimate(name).data,
            equal_nan=True,
            rtol=2.0e-14,
            atol=2.0e-14,
        )


def test_rotated_edi_document_updates_rotation_vectors(tmp_path):
    from pycsamt.emtf.tests.test_edi_conversion import _edi_path

    doc = EMTF.from_edi(_edi_path(tmp_path))
    with pytest.warns(ApproximateVarianceRotationWarning):
        rotated = doc.rotate(30.0, variance_policy="independent")
    edi = rotated.to_edi(on_loss="ignore")

    np.testing.assert_allclose(edi.Z.rotation_angle, [30.0, 30.0])
    np.testing.assert_allclose(edi.Tip.rotation_angle, [30.0, 30.0])
    np.testing.assert_allclose(edi.Z.z, rotated.z)
    np.testing.assert_allclose(edi.Tip.tipper, rotated.tipper)
