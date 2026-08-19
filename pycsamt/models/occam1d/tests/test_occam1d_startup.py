"""Control, parameter, and native-format tests for OCCAMITER_FLEX."""

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DConfig,
    Occam1DModel,
    Occam1DStartup,
)


def _startup(**overrides):
    values = {
        "model_file": "Model",
        "data_file": "Data",
        "parameters": [2.0, 2.5, 3.0],
    }
    values.update(overrides)
    return Occam1DStartup(**values)


def _native(parameter_block="2.0 2.5 3.0\n", **overrides):
    values = {
        "format": "OCCAMITER_FLEX",
        "model": "Model",
        "data": "Data",
        "iterations": "30",
        "target": "1.0",
        "roughness": "1",
        "steps": "8",
        "debug": "1",
        "iteration": "0",
        "lagrange": "5.0",
        "count": "3",
    }
    values.update(overrides)
    return (
        f"Format: {values['format']}\n"
        f"Model File: {values['model']}\n"
        f"Data File: {values['data']}\n"
        f"Iterations to run: {values['iterations']}\n"
        f"Target Misfit: {values['target']}\n"
        f"Roughness Type: {values['roughness']}\n"
        f"Stepsize Cut Count: {values['steps']}\n"
        f"Debug Level: {values['debug']}\n"
        f"Iteration: {values['iteration']}\n"
        f"Lagrange Value: {values['lagrange']}\n"
        f"Param Count: {values['count']}\n"
        f"{parameter_block}"
    )


def test_startup_properties_physical_values_and_diagnostics():
    startup = _startup()

    assert startup.n_parameters == 3
    assert startup.is_initial
    assert not startup.is_iteration
    np.testing.assert_allclose(
        startup.physical_resistivity, [100.0, 10**2.5, 1000.0]
    )
    assert startup.resistivity_bounds == (100.0, 1000.0)
    assert startup.diagnostics()["target_misfit"] == 1.0
    assert "kind='startup'" in startup.summary()


def test_constructor_copies_parameter_vector():
    parameters = np.array([2.0, 2.0])
    startup = _startup(parameters=parameters)
    parameters[0] = 3.0

    assert startup.parameters[0] == 2.0


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"parameters": [2.0]}, "at least two"),
        ({"parameters": [2.0, np.nan]}, "finite log10"),
        ({"iteration": -1}, "iteration"),
        ({"target_misfit": 0.0}, "target_misfit"),
        ({"max_iterations": 0}, "max_iterations"),
        ({"roughness_type": 3}, "roughness_type"),
        ({"stepsize_cut_count": -1}, "stepsize_cut_count"),
        ({"debug_level": -1}, "debug_level"),
        ({"lagrange_start": np.inf}, "lagrange_start"),
        ({"model_file": "  "}, "model_file"),
        ({"data_file": "Data\nbad"}, "newline"),
    ],
)
def test_invalid_startup_controls_are_rejected(overrides, message):
    with pytest.raises((TypeError, ValueError), match=message):
        _startup(**overrides)


def test_extra_fields_are_validated_and_read_only_view_is_exposed():
    startup = _startup(extra_fields={"Description": "test run"})

    assert startup.extra_fields_view["Description"] == "test run"
    with pytest.raises(TypeError):
        startup.extra_fields_view["Description"] = "changed"
    with pytest.raises(ValueError, match="override"):
        _startup(extra_fields={"Target Misfit": "2"})
    with pytest.raises(ValueError, match="field name"):
        _startup(extra_fields={"Bad:Key": "x"})


def test_from_model_fills_unset_values_and_records_provenance():
    model = Occam1DModel(
        [0.0, 10.0, 30.0], [100.0, np.nan, 1000.0]
    )
    config = Occam1DConfig(
        n_layers=3,
        first_thickness=10.0,
        depth_max=30.0,
        starting_resistivity=50.0,
    )

    startup = Occam1DStartup.from_model(model, config)

    np.testing.assert_allclose(startup.physical_resistivity, [100, 50, 1000])
    assert startup.metadata["filled_resistivity_count"] == 1
    applied = startup.apply_to_model(model)
    np.testing.assert_allclose(applied.resistivity, [100, 50, 1000])
    with pytest.raises(ValueError, match="layer count"):
        startup.apply_to_model(Occam1DModel([0, 10], [1, 1]))
    with pytest.raises(TypeError, match="Occam1DModel"):
        startup.apply_to_model(object())


def test_physical_resistivity_overflow_is_reported():
    startup = _startup(parameters=[2.0, 1000.0])

    with pytest.raises(ValueError, match="exceed"):
        startup.physical_resistivity


def test_canonical_write_and_extra_field_roundtrip(tmp_path):
    startup = _startup(
        iteration=4,
        extra_fields={
            "Description": "solver iteration",
            "Misfit Value": "9.5D-01",
        },
    )

    path = startup.write(tmp_path / "ITER_04")
    text = path.read_text(encoding="utf8")
    restored = Occam1DStartup.read(path)

    assert "Iterations to run:" in text
    assert "Model Values:" not in text
    assert text.index("Param Count:") < text.index("2 2.5 3")
    assert startup.is_bound
    assert restored.is_bound
    assert restored.iteration == 4
    assert restored.is_iteration
    assert restored.extra_fields["Description"] == "solver iteration"
    assert restored.extra_fields["Misfit Value"] == "9.5D-01"
    np.testing.assert_allclose(restored.parameters, startup.parameters)


def test_reader_accepts_aliases_legacy_marker_and_fortran_exponents(tmp_path):
    text = _native(
        parameter_block="Model Values:\n2.0D+00 2.5D+00 3.0D+00\n",
        iterations="3.0D+01",
        target="1.0D+00",
    ).replace("Iterations to run:", "Max Iterations:")
    path = tmp_path / "Startup"
    path.write_text(text, encoding="utf8")

    startup = Occam1DStartup.read(path)

    assert startup.max_iterations == 30
    np.testing.assert_allclose(startup.parameters, [2.0, 2.5, 3.0])


@pytest.mark.parametrize(
    "text, message",
    [
        (_native(format="wrong"), "Expected format"),
        (_native(count="4"), "Declared Param Count"),
        (_native(parameter_block="2 bad 3\n"), "Invalid model parameter"),
        (_native(iteration="1.5"), "must be an integer"),
        (_native(count="1", parameter_block="2\n"), "at least two"),
        (
            _native().replace(
                "Target Misfit: 1.0\n",
                "Target Misfit: 1.0\nTarget Misfit: 2.0\n",
            ),
            "Duplicate native field",
        ),
        (_native().replace("Data File: Data\n", ""), "Missing required"),
    ],
)
def test_reader_rejects_corrupt_native_files(tmp_path, text, message):
    path = tmp_path / "Startup"
    path.write_text(text, encoding="utf8")

    with pytest.raises(ValueError, match=message):
        Occam1DStartup.read(path)
