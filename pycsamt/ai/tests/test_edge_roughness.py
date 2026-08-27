from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.ai.inversion.edge_roughness import (
    compute_roughness_exceptions,
    roughness_exception_weight,
)


def _model(layers):
    return SimpleNamespace(
        layers=layers,
        n_params=sum(int(layer["n_cols"]) for layer in layers),
    )


# -----------------------------------------------------------------------
# roughness_exception_weight
# -----------------------------------------------------------------------


def test_weight_is_one_when_both_endpoints_uncertain():
    w = roughness_exception_weight(
        5.0, 5.0, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    assert w == pytest.approx(1.0, abs=1e-3)


def test_weight_approaches_floor_when_both_confident():
    w = roughness_exception_weight(
        1e-6, 1e-6, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    assert w == pytest.approx(0.3, abs=1e-3)


def test_weight_uses_weaker_endpoint():
    confident_pair = roughness_exception_weight(
        1e-6, 1e-6, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    mixed_pair = roughness_exception_weight(
        1e-6, 5.0, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    assert mixed_pair > confident_pair
    assert mixed_pair == pytest.approx(1.0, abs=1e-2)


def test_weight_is_monotonic_in_uncertainty():
    stds = [0.001, 0.01, 0.1, 1.0, 10.0]
    weights = [
        roughness_exception_weight(
            s, s, sigma_ai_floor=0.1, roughness_weight_floor=0.3
        )
        for s in stds
    ]
    assert weights == sorted(weights)


def test_weight_bounds_respect_floor():
    for s in (0.0, 1e-9, 0.05, 0.1, 1.0, 100.0):
        w = roughness_exception_weight(
            s, s, sigma_ai_floor=0.1, roughness_weight_floor=0.3
        )
        assert 0.3 - 1e-9 <= w <= 1.0 + 1e-9


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (dict(std_i=-1.0, std_j=0.0, sigma_ai_floor=0.1, roughness_weight_floor=0.3), "std_i"),
        (dict(std_i=0.0, std_j=float("nan"), sigma_ai_floor=0.1, roughness_weight_floor=0.3), "std_j"),
        (dict(std_i=0.0, std_j=0.0, sigma_ai_floor=0.0, roughness_weight_floor=0.3), "sigma_ai_floor"),
        (dict(std_i=0.0, std_j=0.0, sigma_ai_floor=0.1, roughness_weight_floor=1.5), "roughness_weight_floor"),
    ],
)
def test_weight_rejects_invalid_inputs(kwargs, match):
    with pytest.raises(ValueError, match=match):
        roughness_exception_weight(**kwargs)


# -----------------------------------------------------------------------
# compute_roughness_exceptions
# -----------------------------------------------------------------------


def _two_layer_model():
    return _model(
        [
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
        ]
    )


def test_compute_exceptions_confident_region_produces_relaxed_pairs():
    model = _two_layer_model()
    std = np.array([1e-6, 1e-6, 1e-6, 1e-6])
    exceptions = compute_roughness_exceptions(
        model, std, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    assert len(exceptions) == 4  # all four adjacent pairs relaxed
    for _, _, expen in exceptions:
        assert expen == pytest.approx(0.3, abs=1e-3)


def test_compute_exceptions_uncertain_region_produces_no_exceptions():
    model = _two_layer_model()
    std = np.full(4, 10.0)
    exceptions = compute_roughness_exceptions(
        model, std, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    assert exceptions == []


def test_compute_exceptions_is_selective_by_region():
    model = _two_layer_model()
    # Parameters 1 and 3 (top-left, bottom-left) confident; 2 and 4 uncertain.
    std = np.array([1e-6, 10.0, 1e-6, 10.0])
    exceptions = compute_roughness_exceptions(
        model, std, sigma_ai_floor=0.1, roughness_weight_floor=0.3
    )
    pairs = {(i, j) for i, j, _ in exceptions}
    assert pairs == {(1, 3)}  # only the confident-confident pair relaxed


def test_compute_exceptions_indices_are_valid_bricks():
    model = _two_layer_model()
    std = np.full(4, 1e-6)
    exceptions = compute_roughness_exceptions(model, std, sigma_ai_floor=0.1)
    for i, j, _ in exceptions:
        assert 1 <= i <= model.n_params
        assert 1 <= j <= model.n_params


def test_compute_exceptions_rejects_length_mismatch():
    model = _two_layer_model()
    with pytest.raises(ValueError, match="n_params"):
        compute_roughness_exceptions(model, np.zeros(3))


def test_compute_exceptions_rejects_negative_std():
    model = _two_layer_model()
    with pytest.raises(ValueError, match="non-negative"):
        compute_roughness_exceptions(model, np.array([-1.0, 0.0, 0.0, 0.0]))
