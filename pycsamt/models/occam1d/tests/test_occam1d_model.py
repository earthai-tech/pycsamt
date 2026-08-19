"""Geometry, parameter, and native-I/O tests for Occam1D models."""

import numpy as np
import pytest

from pycsamt.models.occam1d import Occam1DModel


def _model(**overrides):
    values = {
        "depth": [0.0, 10.0, 30.0],
        "resistivity": [100.0, 50.0, 200.0],
    }
    values.update(overrides)
    return Occam1DModel(**values)


def test_model_geometry_properties_and_diagnostics():
    model = _model()

    assert model.n_layers == 3
    assert model.n_finite_layers == 2
    assert model.depth_max == 30.0
    assert model.is_parameterized
    assert model.n_parameters == 3
    assert model.resistivity_bounds == (50.0, 200.0)
    np.testing.assert_allclose(model.finite_thickness, [10.0, 20.0])
    np.testing.assert_allclose(model.conductance, [0.1, 0.4])
    assert np.isinf(model.thickness[-1])
    assert len(model.to_records()) == 3
    assert model.diagnostics()["depth_max_m"] == 30.0
    assert "layers=3" in model.summary()


def test_constructor_copies_arrays_and_defaults_controls():
    depth = np.array([0.0, 10.0])
    model = Occam1DModel(depth, [100.0, 200.0])
    depth[1] = 20.0

    assert model.depth[1] == 10.0
    np.testing.assert_array_equal(model.penalty, [0.0, 0.0])
    np.testing.assert_array_equal(model.preference, [0.0, 0.0])
    np.testing.assert_array_equal(model.weight, [0.0, 0.0])


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"depth": [1.0, 10.0, 30.0]}, "start at zero"),
        ({"depth": [0.0, 10.0, 10.0]}, "increase strictly"),
        ({"resistivity": [100.0, 0.0, 20.0]}, "positive"),
        ({"resistivity": [100.0, np.inf, 20.0]}, "infinite"),
        ({"penalty": [0.0, -1.0, 0.0]}, "non-negative"),
        ({"weight": [0.0, np.nan, 0.0]}, "finite"),
        ({"preference": [0.0, 1.0]}, "exactly 3"),
    ],
)
def test_invalid_model_states_are_rejected(overrides, message):
    with pytest.raises(ValueError, match=message):
        _model(**overrides)


def test_geometric_build_hits_all_boundary_conditions():
    model = Occam1DModel.build(8, 5.0, 2000.0, resistivity=80.0)

    assert model.finite_thickness[0] == pytest.approx(5.0)
    assert model.depth_max == pytest.approx(2000.0)
    assert np.all(np.diff(model.finite_thickness) > 0)
    ratios = model.finite_thickness[1:] / model.finite_thickness[:-1]
    np.testing.assert_allclose(ratios, ratios[0], rtol=1e-10)


def test_uniform_build_and_two_layer_geometry():
    uniform = Occam1DModel.build(4, 10.0, 30.0)
    np.testing.assert_allclose(uniform.finite_thickness, 10.0)

    two = Occam1DModel.build(2, 10.0, 10.0)
    np.testing.assert_allclose(two.depth, [0.0, 10.0])
    with pytest.raises(ValueError, match="two-layer"):
        Occam1DModel.build(2, 10.0, 20.0)


def test_filled_and_replacement_models_are_independent():
    model = _model(resistivity=[100.0, np.nan, 200.0])

    assert not model.is_parameterized
    assert model.n_parameters == 2
    assert np.isnan(model.conductance[1])
    filled = model.filled(75.0)

    assert filled.is_parameterized
    assert filled.resistivity[1] == 75.0
    assert np.isnan(model.resistivity[1])
    replacement = model.with_resistivity([1.0, 2.0, 3.0])
    np.testing.assert_allclose(replacement.resistivity, [1.0, 2.0, 3.0])


def test_unset_model_cannot_be_written(tmp_path):
    model = _model(resistivity=[100.0, np.nan, 200.0])

    with pytest.raises(ValueError, match="filled"):
        model.write(tmp_path / "Model")


def test_native_roundtrip_preserves_control_columns_and_binding(tmp_path):
    model = _model(
        penalty=[0.0, 1.0, 2.0],
        preference=[2.0, 2.1, 2.2],
        weight=[0.0, 0.5, 1.0],
    )

    path = model.write(tmp_path / "Model")
    restored = Occam1DModel.read(path)

    assert model.is_bound
    assert restored.path == path.resolve()
    np.testing.assert_allclose(restored.depth, model.depth)
    np.testing.assert_allclose(restored.resistivity, model.resistivity)
    np.testing.assert_allclose(restored.penalty, model.penalty)
    np.testing.assert_allclose(restored.preference, model.preference)
    np.testing.assert_allclose(restored.weight, model.weight)


def test_reader_accepts_legacy_rows_and_fortran_exponents(tmp_path):
    path = tmp_path / "Model"
    path.write_text(
        "Format: Resistivity1DMod_1.0\n"
        " # LAYERS : 2\n"
        "0.0D+00 1.0D+02\n"
        "1.0D+01 2.0D+02 1 2\n",
        encoding="utf8",
    )

    model = Occam1DModel.read(path)

    np.testing.assert_allclose(model.depth, [0.0, 10.0])
    np.testing.assert_allclose(model.penalty, [0.0, 1.0])
    np.testing.assert_allclose(model.preference, [0.0, 2.0])
    np.testing.assert_allclose(model.weight, [0.0, 0.0])


@pytest.mark.parametrize(
    "body, message",
    [
        ("# Layers: nope\n", "Invalid # Layers"),
        ("# Layers: 1\n0 100\n", "at least two"),
        ("# Layers: 2\n0 100\n", "Declared # Layers"),
        ("# Layers: 2\n0\n10 100\n", "two to five"),
        ("# Layers: 2\n0 100 0 0 0 9\n10 100\n", "two to five"),
        ("# Layers: 2\n0 100\n10 bad\n", "Invalid layer"),
    ],
)
def test_reader_rejects_corrupt_native_models(tmp_path, body, message):
    path = tmp_path / "Model"
    path.write_text(
        f"Format: Resistivity1DMod_1.0\n{body}", encoding="utf8"
    )

    with pytest.raises(ValueError, match=message):
        Occam1DModel.read(path)
