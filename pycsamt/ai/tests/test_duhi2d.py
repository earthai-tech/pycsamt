from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.api.property import PyCSAMTObject
from pycsamt.compat.sklearn import InvalidParameterError
from pycsamt.ai.inversion.duhi2d import (
    DUHIPreparation,
    DUHIInverter2D,
    apply_observation_reliability,
    map_ai_grid_to_occam,
)
from pycsamt.ai.inversion.schema import (
    DUHI_INIT_SCHEMA,
    DUHI_PREPARE_SCHEMA,
    OBSERVATION_RELIABILITY_SCHEMA,
)
from pycsamt.models.occam2d import OccamConfig, OccamPrejudice


class _Writable:
    def __init__(self, **values):
        self.__dict__.update(values)

    def write(self, path):
        from pathlib import Path

        path = Path(path)
        path.write_text("written", encoding="ascii")
        self.path = path
        return path


def _builder(tmp_path):
    model = _Writable(
        layers=[
            {"n_cols": 3, "n_merge": 1, "params": np.array([1, 1, 1])},
            {"n_cols": 2, "n_merge": 1, "params": np.array([1, 2])},
        ],
        n_params=5,
        prejudice_file="none",
    )
    return SimpleNamespace(
        is_ready=True,
        workdir=tmp_path,
        config=OccamConfig(
            data_file="OccamData.dat",
            model_file="OccamModel",
            startup_file="Startup",
        ),
        data=_Writable(
            data_blocks=np.array(
                [[1, 1, 1, 2.0, 0.1], [1, 1, 2, 45.0, 1.0]], dtype=float
            )
        ),
        mesh=SimpleNamespace(
            x_widths=np.ones(3),
            z_widths=np.ones(2),
            n_airlayers=0,
        ),
        model=model,
        startup=_Writable(param_values=np.zeros(5)),
    )


def test_reliability_inflates_only_uncertain_errors():
    actual = apply_observation_reliability([2.0, 2.0], [1.0, 0.25])
    np.testing.assert_allclose(actual, [2.0, 4.0])
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        apply_observation_reliability([1.0], [1.1])


def test_ai_grid_mapping_matches_layer_major_parameter_count(tmp_path):
    builder = _builder(tmp_path)
    mapped = map_ai_grid_to_occam(
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        builder.model,
        builder.mesh,
    )
    assert mapped.shape == (5,)
    np.testing.assert_allclose(mapped[:3], [1.0, 1.5, 2.0])


def test_duhi_prepare_writes_native_inputs(tmp_path):
    builder = _builder(tmp_path)
    inverter = DUHIInverter2D(lambda_ai=2.0, sigma_ai_floor=0.1)
    result = inverter.prepare(
        builder,
        ai_mean=np.array([[1.0, 2.0], [3.0, 4.0]]),
        ai_std=np.full((2, 2), 0.2),
        observation_reliability=np.array([1.0, 0.25]),
    )

    assert isinstance(inverter, PyCSAMTObject)
    assert isinstance(result, DUHIPreparation)
    assert inverter.is_prepared
    assert inverter.preparation is result
    assert inverter.prejudice.n_prejudiced == 5
    assert "prepared" in inverter.summary()
    assert result.n_data == 2
    assert result.n_params == 5
    assert result.ai_initialized
    assert set(result.files) == {"data", "model", "startup", "prejudice"}
    np.testing.assert_allclose(builder.data.data_blocks[:, 4], [0.1, 2.0])
    assert builder.model.prejudice_file == "DUHIPrejudice"
    np.testing.assert_allclose(builder.startup.param_values[:3], [1.0, 1.5, 2.0])
    prejudice = OccamPrejudice.read(result.prejudice_file)
    np.testing.assert_allclose(prejudice.weights, 2.0 / np.sqrt(0.05))

    mean_parameters = inverter.ai_mean_parameters
    mean_parameters[:] = -999.0
    assert np.all(inverter.ai_mean_parameters != -999.0)
    assert inverter.ai_std_parameters.shape == (5,)
    assert inverter.prejudice_weights.shape == (5,)

    with pytest.raises(RuntimeError, match="one project"):
        inverter.prepare(
            builder,
            ai_mean=np.ones((2, 2)),
            ai_std=np.ones((2, 2)),
            observation_reliability=np.ones(2),
        )


def test_duhi_unprepared_state_and_control_validation():
    inverter = DUHIInverter2D()
    assert not inverter.is_prepared
    assert "unprepared" in str(inverter)
    with pytest.raises(RuntimeError, match=r"prepare\(\)"):
        _ = inverter.preparation

    with pytest.raises(ValueError, match="lambda_ai"):
        DUHIInverter2D(lambda_ai=-1.0)
    with pytest.raises(ValueError, match="sigma_ai_floor"):
        DUHIInverter2D(sigma_ai_floor=0.0)
    with pytest.raises(ValueError, match="reliability_floor"):
        DUHIInverter2D(reliability_floor=1.1)
    with pytest.raises(ValueError, match="local file"):
        DUHIInverter2D(prejudice_filename="nested/prejudice")
    with pytest.raises(TypeError, match="callable"):
        DUHIInverter2D(grid_mapper=1)


def test_duhi_custom_mapper_can_preserve_startup(tmp_path):
    builder = _builder(tmp_path)
    original = builder.startup.param_values.copy()
    calls = []

    def mapper(grid, model, mesh, x_coordinates, z_coordinates):
        calls.append((grid.shape, model.n_params))
        return np.linspace(float(grid.min()), float(grid.max()), model.n_params)

    inverter = DUHIInverter2D(grid_mapper=mapper)
    result = inverter.prepare(
        builder,
        ai_mean=np.array([[1.0, 2.0], [3.0, 4.0]]),
        ai_std=np.full((2, 2), 0.5),
        observation_reliability=np.ones(2),
        ai_initialize=False,
    )

    assert calls == [((2, 2), 5), ((2, 2), 5)]
    assert not result.ai_initialized
    np.testing.assert_allclose(builder.startup.param_values, original)


def test_duhi_rejects_incompatible_run_inputs(tmp_path):
    builder = _builder(tmp_path)
    inverter = DUHIInverter2D()
    with pytest.raises(ValueError, match="identical shape"):
        inverter.prepare(
            builder,
            ai_mean=np.ones((2, 2)),
            ai_std=np.ones((3, 2)),
            observation_reliability=np.ones(2),
        )

    builder = _builder(tmp_path)
    with pytest.raises(ValueError, match="per Occam datum"):
        DUHIInverter2D().prepare(
            builder,
            ai_mean=np.ones((2, 2)),
            ai_std=np.ones((2, 2)),
            observation_reliability=np.ones(3),
        )


def test_duhi_uses_compatibility_parameter_schemas(tmp_path):
    assert (
        apply_observation_reliability._skl_parameter_constraints
        is OBSERVATION_RELIABILITY_SCHEMA
    )
    assert (
        DUHIInverter2D.__init__._skl_parameter_constraints
        is DUHI_INIT_SCHEMA
    )
    assert (
        DUHIInverter2D.prepare._skl_parameter_constraints
        is DUHI_PREPARE_SCHEMA
    )
    with pytest.raises(InvalidParameterError, match="lambda_ai"):
        DUHIInverter2D(lambda_ai="strong")
    with pytest.raises(InvalidParameterError, match="ai_initialize"):
        DUHIInverter2D().prepare(
            _builder(tmp_path),
            ai_mean=np.ones((2, 2)),
            ai_std=np.ones((2, 2)),
            observation_reliability=np.ones(2),
            ai_initialize="yes",
        )
