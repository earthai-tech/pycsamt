from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.ai.inversion.gradient_prejudice import (
    boost_prejudice_weights,
    gradient_boost_multiplier,
)


def _model(layers):
    return SimpleNamespace(
        layers=layers,
        n_params=sum(int(layer["n_cols"]) for layer in layers),
    )


def _two_layer_model():
    return _model(
        [
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
            {"n_merge": 1, "n_cols": 2, "params": np.array([2, 2])},
        ]
    )


# -----------------------------------------------------------------------
# gradient_boost_multiplier
# -----------------------------------------------------------------------


def test_multiplier_is_one_for_flat_neighbourhood():
    m = gradient_boost_multiplier(
        0.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0
    )
    assert m == pytest.approx(1.0)


def test_multiplier_grows_with_gradient():
    grads = [0.0, 0.1, 0.3, 1.0, 5.0]
    mult = [
        gradient_boost_multiplier(
            g, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0
        )
        for g in grads
    ]
    assert mult == sorted(mult)
    assert mult[0] == pytest.approx(1.0)


def test_multiplier_saturates_at_cap():
    near_cap = gradient_boost_multiplier(
        100.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0
    )
    far_past_cap = gradient_boost_multiplier(
        1_000_000.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0
    )
    assert near_cap == pytest.approx(far_past_cap)
    assert near_cap == pytest.approx(1.0 + 2.0 * 5.0)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (dict(local_gradient=-1.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0), "local_gradient"),
        (dict(local_gradient=float("nan"), gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0), "local_gradient"),
        (dict(local_gradient=0.0, gradient_scale=0.0, boost_gain=2.0, boost_cap=5.0), "gradient_scale"),
        (dict(local_gradient=0.0, gradient_scale=0.3, boost_gain=-1.0, boost_cap=5.0), "boost_gain"),
        (dict(local_gradient=0.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=-1.0), "boost_cap"),
    ],
)
def test_multiplier_rejects_invalid_inputs(kwargs, match):
    with pytest.raises(ValueError, match=match):
        gradient_boost_multiplier(**kwargs)


# -----------------------------------------------------------------------
# boost_prejudice_weights
# -----------------------------------------------------------------------


def test_boost_is_identity_for_uniform_ai_mean():
    model = _two_layer_model()
    mean = np.full(4, 2.0)
    weights = np.array([1.0, 1.0, 1.0, 1.0])
    boosted = boost_prejudice_weights(model, mean, weights)
    assert boosted == pytest.approx(weights)


def test_boost_amplifies_only_near_a_sharp_boundary():
    # A single row of six adjacent parameters with one sharp jump
    # between params 3 and 4 (indices 2, 3); everything else is flat.
    model = _model(
        [{"n_merge": 1, "n_cols": 6, "params": np.ones(6, dtype=int)}]
    )
    mean = np.array([1.0, 1.0, 1.0, 4.0, 4.0, 4.0])
    weights = np.ones(6)
    boosted = boost_prejudice_weights(
        model, mean, weights, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0
    )
    assert np.all(boosted >= weights)
    # Only the two params straddling the boundary are boosted.
    near_boundary = [2, 3]
    away_from_boundary = [0, 1, 4, 5]
    for idx in near_boundary:
        assert boosted[idx] > 1.0 + 1e-6
    for idx in away_from_boundary:
        assert boosted[idx] == pytest.approx(weights[idx])


def test_boost_never_decreases_weights():
    rng = np.random.default_rng(0)
    model = _two_layer_model()
    mean = rng.normal(size=4)
    weights = rng.uniform(0.1, 2.0, size=4)
    boosted = boost_prejudice_weights(model, mean, weights)
    assert np.all(boosted >= weights - 1e-12)


def test_boost_rejects_mean_length_mismatch():
    model = _two_layer_model()
    with pytest.raises(ValueError, match="ai_mean_parameters"):
        boost_prejudice_weights(model, np.zeros(3), np.ones(4))


def test_boost_rejects_weights_length_mismatch():
    model = _two_layer_model()
    with pytest.raises(ValueError, match="prejudice_weights"):
        boost_prejudice_weights(model, np.zeros(4), np.ones(3))


def test_boost_rejects_negative_weights():
    model = _two_layer_model()
    with pytest.raises(ValueError, match="non-negative"):
        boost_prejudice_weights(model, np.zeros(4), np.array([-1.0, 0.0, 0.0, 0.0]))
