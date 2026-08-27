from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.compat.sklearn import InvalidParameterError
from pycsamt.ai.inversion.mapping2d import (
    build_occam_parameter_adjacency,
    map_ai_grid_to_occam,
)
from pycsamt.ai.inversion.schema import GRID_MAPPING_SCHEMA


def _model(layers):
    return SimpleNamespace(
        layers=layers,
        n_params=sum(int(layer["n_cols"]) for layer in layers),
    )


def _mesh(x_widths, z_widths, n_airlayers=0):
    return SimpleNamespace(
        x_widths=np.asarray(x_widths, dtype=float),
        z_widths=np.asarray(z_widths, dtype=float),
        n_airlayers=n_airlayers,
    )


def test_geometry_mapper_area_weights_merged_parameter():
    mesh = _mesh([1.0, 3.0], [2.0, 4.0])
    model = _model(
        [
            {
                "n_merge": 2,
                "n_cols": 1,
                "params": np.array([2]),
            }
        ]
    )
    x = np.array([0.5, 2.5])
    z = np.array([1.0, 4.0])
    grid = z[:, None] * 10.0 + x[None, :]

    mapped = map_ai_grid_to_occam(grid, model, mesh, x, z)

    np.testing.assert_allclose(mapped, [32.0])


def test_geometry_mapper_excludes_air_and_respects_layers():
    mesh = _mesh([1.0, 3.0], [5.0, 2.0, 4.0], n_airlayers=1)
    model = _model(
        [
            {
                "n_merge": 1,
                "n_cols": 1,
                "params": np.array([2]),
            },
            {
                "n_merge": 1,
                "n_cols": 1,
                "params": np.array([2]),
            },
        ]
    )
    x = np.array([0.5, 2.5])
    z = np.array([1.0, 4.0])
    grid = z[:, None] * 10.0 + x[None, :]

    mapped = map_ai_grid_to_occam(grid, model, mesh, x, z)

    np.testing.assert_allclose(mapped, [12.0, 42.0])


def test_geometry_mapper_respects_horizontal_spans():
    mesh = _mesh([1.0, 2.0, 3.0, 4.0], [1.0])
    model = _model(
        [
            {
                "n_merge": 1,
                "n_cols": 3,
                "params": np.array([1, 2, 1]),
            }
        ]
    )
    x = np.array([0.5, 2.0, 4.5, 8.0])
    grid = x.reshape(1, -1)

    mapped = map_ai_grid_to_occam(grid, model, mesh, x, [0.5])

    np.testing.assert_allclose(mapped, [0.5, 3.5, 8.0])


def test_geometry_mapper_infers_uniform_source_centres():
    mesh = _mesh([1.0, 1.0], [1.0, 1.0])
    model = _model(
        [
            {
                "n_merge": 1,
                "n_cols": 1,
                "params": np.array([2]),
            },
            {
                "n_merge": 1,
                "n_cols": 1,
                "params": np.array([2]),
            },
        ]
    )
    grid = np.array([[1.0, 3.0], [5.0, 7.0]])

    mapped = map_ai_grid_to_occam(grid, model, mesh)

    np.testing.assert_allclose(mapped, [2.0, 6.0])


@pytest.mark.parametrize(
    ("model", "message"),
    [
        (
            _model(
                [
                    {
                        "n_merge": 1,
                        "n_cols": 1,
                        "params": np.array([1]),
                    }
                ]
            ),
            "span the mesh width",
        ),
        (
            _model(
                [
                    {
                        "n_merge": 1,
                        "n_cols": 1,
                        "params": np.array([2]),
                    }
                ]
            ),
            "do not span",
        ),
    ],
)
def test_geometry_mapper_rejects_inconsistent_model(model, message):
    mesh = _mesh([1.0, 1.0], [1.0, 1.0])
    with pytest.raises(ValueError, match=message):
        map_ai_grid_to_occam(np.ones((2, 2)), model, mesh)


def test_geometry_mapper_rejects_invalid_coordinates():
    mesh = _mesh([1.0, 1.0], [1.0])
    model = _model(
        [
            {
                "n_merge": 1,
                "n_cols": 1,
                "params": np.array([2]),
            }
        ]
    )
    with pytest.raises(ValueError, match="strictly increasing"):
        map_ai_grid_to_occam(
            np.ones((1, 2)),
            model,
            mesh,
            x_coordinates=[1.0, 1.0],
        )


def test_geometry_mapper_uses_compatibility_parameter_schema():
    assert (
        map_ai_grid_to_occam._skl_parameter_constraints
        is GRID_MAPPING_SCHEMA
    )
    with pytest.raises(InvalidParameterError, match="grid"):
        map_ai_grid_to_occam(1.0, object(), object())


def test_parameter_adjacency_single_layer_is_consecutive_only():
    model = _model(
        [{"n_merge": 1, "n_cols": 3, "params": np.array([1, 2, 1])}]
    )
    assert build_occam_parameter_adjacency(model) == [(1, 2), (2, 3)]


def test_parameter_adjacency_aligned_layers_links_directly_below():
    model = _model(
        [
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
        ]
    )
    assert build_occam_parameter_adjacency(model) == [
        (1, 2), (1, 3), (2, 4), (3, 4),
    ]


def test_parameter_adjacency_handles_span_mismatch_between_layers():
    model = _model(
        [
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
            {"n_merge": 1, "n_cols": 1, "params": np.array([4])},
        ]
    )
    assert build_occam_parameter_adjacency(model) == [(1, 2), (1, 3), (2, 3)]


def test_parameter_adjacency_pairs_are_unique_and_sorted():
    model = _model(
        [
            {"n_merge": 1, "n_cols": 3, "params": np.array([2, 2, 2])},
            {"n_merge": 1, "n_cols": 3, "params": np.array([2, 2, 2])},
        ]
    )
    pairs = build_occam_parameter_adjacency(model)
    assert len(pairs) == len(set(pairs))
    assert pairs == sorted(pairs)
    assert all(i < j for i, j in pairs)


def test_parameter_adjacency_rejects_empty_model():
    with pytest.raises(ValueError, match="no parameter layers"):
        build_occam_parameter_adjacency(_model([]))


def test_parameter_adjacency_rejects_missing_layers_attribute():
    with pytest.raises(TypeError, match="OccamModel interface"):
        build_occam_parameter_adjacency(object())
