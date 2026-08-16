import numpy as np
import pytest

from pycsamt.compat.sklearn import InvalidParameterError
from pycsamt.models.occam2d import OccamConfig, OccamPrejudice
from pycsamt.models.occam2d.base import OccamBase
from pycsamt.models.occam2d.schema import (
    OCCAM_PREJUDICE_DENSE_SCHEMA,
    OCCAM_PREJUDICE_INIT_SCHEMA,
)


def test_prejudice_sparse_roundtrip(tmp_path):
    config = OccamConfig()
    prejudice = OccamPrejudice.from_dense(
        [1.5, 2.0, 2.5],
        [4.0, 0.0, 1.0],
        config=config,
    )
    path = prejudice.write(tmp_path / "DUHIPrejudice")

    assert isinstance(prejudice, OccamBase)
    assert prejudice.config is config
    assert prejudice.path == path
    text = path.read_text(encoding="ascii")
    assert text.startswith("Format:           OCCAM2MTPREJ_2.0")
    assert "Param Count:      2" in text
    assert "1 6 4" in text
    assert "3 2.5 1" in text
    restored = OccamPrejudice.read(path, config=config)
    assert restored.path == path
    assert restored.config is config
    np.testing.assert_array_equal(restored.parameter_indices, [1, 3])
    np.testing.assert_allclose(restored.target_values, [1.5, 2.5])
    np.testing.assert_allclose(restored.weights, [4.0, 1.0])


def test_prejudice_validates_solver_parameter_count():
    prejudice = OccamPrejudice([1, 4], [2.0, 3.0], [1.0, 1.0])
    with pytest.raises(ValueError, match="beyond"):
        prejudice.validate_parameter_count(3)


def test_prejudice_empty_container_and_native_values():
    empty = OccamPrejudice()
    assert empty.n_prejudiced == 0
    assert empty.native_values.size == 0

    prejudice = OccamPrejudice(
        (index for index in [1, 2]),
        (target for target in [1.5, 2.5]),
        (weight for weight in [2.0, 0.0]),
    )
    np.testing.assert_allclose(prejudice.native_values, [3.0, 2.5])
    assert prejudice.validate_parameter_count(2) is prejudice


@pytest.mark.parametrize(
    ("indices", "targets", "weights", "message"),
    [
        ([0], [1.0], [1.0], "one-based"),
        ([1, 1], [1.0, 2.0], [1.0, 1.0], "unique"),
        ([1], [np.nan], [1.0], "finite"),
        ([1], [1.0], [-1.0], "non-negative"),
    ],
)
def test_prejudice_rejects_invalid_records(
    indices, targets, weights, message
):
    with pytest.raises(ValueError, match=message):
        OccamPrejudice(indices, targets, weights)


def test_prejudice_uses_compatibility_parameter_schemas(tmp_path):
    assert (
        OccamPrejudice.__init__._skl_parameter_constraints
        is OCCAM_PREJUDICE_INIT_SCHEMA
    )
    assert (
        OccamPrejudice.from_dense.__func__._skl_parameter_constraints
        is OCCAM_PREJUDICE_DENSE_SCHEMA
    )
    with pytest.raises(InvalidParameterError, match="parameter_indices"):
        OccamPrejudice(1, [2.0], [1.0])
    with pytest.raises(InvalidParameterError, match="n_params"):
        OccamPrejudice().validate_parameter_count(0)
    with pytest.raises(InvalidParameterError, match="path"):
        OccamPrejudice.read(1)
    with pytest.raises(InvalidParameterError, match="path"):
        OccamPrejudice().write(1)
