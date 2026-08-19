"""Scientific and numerical contracts for native Occam1D forward modelling."""

import json

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DData,
    Occam1DForwardModel,
    Occam1DForwardResponse,
    Occam1DModel,
)


def _uniform_model(resistivity=100.0):
    return Occam1DModel(
        [0.0, 10.0, 100.0, 1000.0],
        [resistivity] * 4,
    )


def test_uniform_earth_has_constant_resistivity_and_45_degree_phase():
    frequencies = np.logspace(-4, 4, 33)
    response = Occam1DForwardModel(_uniform_model()).predict(frequencies)
    np.testing.assert_allclose(response.resistivity, 100.0, rtol=2e-14)
    np.testing.assert_allclose(response.phase, 45.0, atol=2e-13)


def test_uniform_earth_is_independent_of_discretization():
    frequencies = np.logspace(-3, 3, 15)
    coarse = Occam1DForwardModel(
        Occam1DModel([0.0, 500.0], [30.0, 30.0])
    ).predict(frequencies)
    fine = Occam1DForwardModel(
        Occam1DModel(
            [0.0, 1.0, 5.0, 25.0, 100.0, 500.0],
            [30.0] * 6,
        )
    ).predict(frequencies)
    np.testing.assert_allclose(coarse.impedance, fine.impedance, rtol=2e-14)


def test_two_layer_limits_approach_surface_and_basement_resistivity():
    model = Occam1DModel([0.0, 100.0], [10.0, 1000.0])
    response = Occam1DForwardModel(model).predict([1e-8, 1e8])
    assert response.resistivity[0] == pytest.approx(1000.0, rel=0.02)
    assert response.resistivity[1] == pytest.approx(10.0, rel=1e-10)


def test_tensor_is_antisymmetric_and_determinant_recovers_scalar():
    forward = Occam1DForwardModel(_uniform_model())
    tensor = forward.tensor_impedance([0.1, 1.0, 10.0])
    np.testing.assert_allclose(tensor[:, 0, 0], 0.0)
    np.testing.assert_allclose(tensor[:, 1, 1], 0.0)
    np.testing.assert_allclose(tensor[:, 1, 0], -tensor[:, 0, 1])
    determinant = np.sqrt(-tensor[:, 0, 1] * tensor[:, 1, 0])
    np.testing.assert_allclose(determinant, tensor[:, 0, 1])


@pytest.mark.parametrize("mode", ["xy", "yx", "te", "tm", "determinant"])
def test_all_modes_use_consistent_1d_physical_response(mode):
    response = Occam1DForwardModel(_uniform_model()).predict([1.0], mode=mode)
    assert response.resistivity[0] == pytest.approx(100.0)
    assert response.phase[0] == pytest.approx(45.0)


def test_log10_override_does_not_mutate_bound_model():
    model = _uniform_model(100.0)
    forward = Occam1DForwardModel(model)
    response = forward.predict([1.0], log10_resistivity=[2.5] * 4)
    assert response.resistivity[0] == pytest.approx(10**2.5)
    np.testing.assert_allclose(model.resistivity, 100.0)
    np.testing.assert_allclose(forward.predict([1.0]).resistivity, 100.0)


def test_frequency_order_is_preserved_and_response_is_immutable():
    frequency = np.array([100.0, 0.1, 10.0])
    response = Occam1DForwardModel(_uniform_model()).predict(frequency)
    np.testing.assert_array_equal(response.frequency, frequency)
    with pytest.raises(ValueError):
        response.resistivity[0] = 1.0
    np.testing.assert_allclose(response.period, 1.0 / frequency)


def test_solver_vector_matches_native_row_order_and_missing_pattern():
    data = Occam1DData(
        [100.0, 10.0, 1.0],
        [100.0, np.nan, 100.0],
        [45.0, 45.0, np.nan],
        [0.05, np.nan, 0.05],
        [2.0, 2.0, np.nan],
        mode="xy",
    )
    values = Occam1DForwardModel(_uniform_model()).solver_vector(data)
    np.testing.assert_allclose(values, [2.0, 45.0, 45.0, 2.0])
    assert values.shape == (data.n_data,)


def test_response_records_and_diagnostics_are_json_safe():
    response = Occam1DForwardModel(_uniform_model()).predict([1.0, 10.0])
    assert isinstance(response, Occam1DForwardResponse)
    assert len(response.to_records()) == 2
    values = json.loads(json.dumps(response.diagnostics()))
    assert values["n_frequencies"] == 2
    assert values["metadata"]["earth"] == "isotropic_1d"


@pytest.mark.parametrize(
    "frequency, message",
    [([], "non-empty"), ([0.0], "strictly positive"), ([np.nan], "finite")],
)
def test_invalid_frequency_is_rejected(frequency, message):
    with pytest.raises(ValueError, match=message):
        Occam1DForwardModel(_uniform_model()).predict(frequency)


def test_model_and_parameter_contracts_are_explicit():
    forward = Occam1DForwardModel()
    with pytest.raises(RuntimeError, match="No layer model"):
        forward.predict([1.0])
    with pytest.raises(ValueError, match="not both"):
        Occam1DForwardModel(_uniform_model()).predict(
            [1.0],
            resistivity=[10.0] * 4,
            log10_resistivity=[1.0] * 4,
        )
    with pytest.raises(ValueError, match="shape"):
        Occam1DForwardModel(_uniform_model()).predict(
            [1.0],
            resistivity=[10.0],
        )


def test_invalid_permeability_and_mode_are_rejected():
    with pytest.raises(ValueError, match="permeability"):
        Occam1DForwardModel(permeability=0.0)
    forward = Occam1DForwardModel(_uniform_model())
    with pytest.raises(ValueError, match="mode"):
        forward.predict([1.0], mode="tipper")
