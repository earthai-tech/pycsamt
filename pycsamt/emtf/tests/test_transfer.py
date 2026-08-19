from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtf import StatisticalEstimate, TransferFunction


def _periods(n=3):
    return np.geomspace(0.01, 1.0, n)


def test_impedance_transfer_function_matrix_contract():
    z = np.ones((3, 2, 2), dtype=complex) * (1.0 + 2.0j)
    tf = TransferFunction(
        name="Z",
        data=z,
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=_periods(),
    )
    assert tf.name == "impedance"
    assert tf.shape == (3, 2, 2)
    assert tf.n_input == 2
    assert tf.n_output == 2
    assert np.allclose(tf.frequency, 1.0 / _periods())
    assert tf.units == "[mV/km]/[nT]"


def test_tipper_and_arbitrary_matrix_contracts():
    tip = TransferFunction(
        name="tipper",
        data=np.zeros((4, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=np.arange(1.0, 5.0),
    )
    assert tip.shape == (4, 1, 2)

    # Deliberately not registered: demonstrates a future MobileMT-like 3x2
    # admittance can live in the core without changing the class.
    y = TransferFunction(
        name="admittance",
        data=np.zeros((5, 3, 2), dtype=complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Hx", "Hy", "Hz"),
    )
    assert y.shape == (5, 3, 2)


def test_scalar_derived_data_are_promoted_to_1x1_matrices():
    rho = TransferFunction(
        name="apparent_resistivity",
        data=np.array([10.0, 20.0, 30.0]),
        periods=_periods(),
    )
    assert rho.shape == (3, 1, 1)
    assert not np.iscomplexobj(rho.data)


def test_registered_real_type_rejects_nonzero_imaginary_data():
    with pytest.raises(TypeError):
        TransferFunction(
            name="apparent_resistivity",
            data=np.array([10.0 + 1.0j]),
            periods=[1.0],
        )


def test_invalid_matrix_shape_and_period_grid_are_rejected():
    with pytest.raises(ValueError):
        TransferFunction(
            name="impedance",
            data=np.zeros((3, 2, 1)),
            input_channels=("Hx", "Hy"),
            output_channels=("Ex", "Ey"),
        )

    with pytest.raises(ValueError):
        TransferFunction(
            name="impedance",
            data=np.zeros((3, 2, 2)),
            input_channels=("Hx", "Hy"),
            output_channels=("Ex", "Ey"),
            periods=[1.0, 2.0],
        )


def test_statistical_estimate_attachment_and_lookup():
    tf = TransferFunction(
        name="Z",
        data=np.zeros((3, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
    )
    var = StatisticalEstimate(
        name="VAR",
        kind="variance",
        data=np.ones((3, 2, 2)),
    )
    tf.add_estimate(var)
    assert tf.get_estimate("variance") is var
    assert tf.get_estimate("VAR") is var

    with pytest.raises(ValueError):
        tf.add_estimate(var)

    with pytest.raises(ValueError):
        tf.add_estimate(
            StatisticalEstimate(
                name="RESIDCOV",
                kind="residual_covariance",
                data=np.ones((2, 2, 2)),
            )
        )


def test_copy_detaches_arrays_and_estimates():
    tf = TransferFunction(
        name="tipper",
        data=np.ones((2, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0, 2.0],
        estimates={
            "variance": StatisticalEstimate(
                name="VAR",
                kind="variance",
                data=np.ones((2, 1, 2)),
            )
        },
    )
    cp = tf.copy()
    cp.data[0, 0, 0] = 99
    cp.estimates["variance"].data[0, 0, 0] = 77
    assert tf.data[0, 0, 0] == 1
    assert tf.estimates["variance"].data[0, 0, 0] == 1
