"""
pycsamt.agents._param_schemas
==============================

Reusable ``@validate_params`` constraint atoms and composed per-agent
schemas, built on :mod:`pycsamt.compat.sklearn`'s version-compatible
``Interval``/``StrOptions``/``HasMethods`` primitives.

Status
------
Pilot schemas only — :data:`INV2D_PARAM_CONSTRAINTS` and
:data:`INV3D_PARAM_CONSTRAINTS` are not yet applied to
:class:`~pycsamt.agents.inv2d_agent.Inv2DAgent` /
:class:`~pycsamt.agents.inv3d_agent.Inv3DAgent`. Wiring them in means
decorating each agent's ``__init__`` with::

    from ..compat.sklearn import validate_params
    from ._param_schemas import INV2D_PARAM_CONSTRAINTS

    @validate_params(INV2D_PARAM_CONSTRAINTS, prefer_skip_nested_validation=True)
    def __init__(self, *, ...):

Design notes
------------
* Atoms are defined once and shared across agents (:data:`POSITIVE_INT`,
  :data:`POSITIVE_FLOAT`, ...) so "must be > 0" means the same thing
  everywhere it is used.
* ``sklearn.utils._param_validation.validate_parameter_constraints`` only
  validates parameters that appear as keys in the constraints dict and
  silently skips anything else — schemas here are composed by dict-merging
  smaller fragments (:data:`_BASE_AGENT_PARAMS`,
  :data:`_MAXWELL_TRAINING_PARAMS`, ...) and can be extended incrementally
  without needing to be exhaustive.
* This only validates each parameter in isolation (type/range/membership).
  It cannot express cross-parameter rules such as "``gcn_hidden`` is only
  meaningful when ``physics='mt2d_tri'``" — those stay as small manual
  checks in the agent's own ``__init__``/``execute``, same as today.
"""

from __future__ import annotations

from numbers import Integral, Real

from ._base import _PROVIDERS
from ..compat.sklearn import Interval, StrOptions

__all__ = [
    "POSITIVE_INT",
    "NON_NEGATIVE_INT",
    "POSITIVE_FLOAT",
    "NON_NEGATIVE_FLOAT",
    "OPTIONAL_POSITIVE_INT",
    "OPTIONAL_POSITIVE_FLOAT",
    "OPTIONAL_NON_NEGATIVE_FLOAT",
    "ARRAY_LIKE",
    "ARRAY_LIKE_OR_NONE",
    "OBJECT_OR_NONE",
    "INV2D_PARAM_CONSTRAINTS",
    "INV3D_PARAM_CONSTRAINTS",
]

# ── reusable atoms ──────────────────────────────────────────────────────────

POSITIVE_INT = [Interval(Integral, 1, None, closed="left")]
NON_NEGATIVE_INT = [Interval(Integral, 0, None, closed="left")]
POSITIVE_FLOAT = [Interval(Real, 0, None, closed="neither")]
NON_NEGATIVE_FLOAT = [Interval(Real, 0, None, closed="left")]

OPTIONAL_POSITIVE_INT = [Interval(Integral, 1, None, closed="left"), None]
OPTIONAL_POSITIVE_FLOAT = [Interval(Real, 0, None, closed="neither"), None]
OPTIONAL_NON_NEGATIVE_FLOAT = [Interval(Real, 0, None, closed="left"), None]

ARRAY_LIKE = ["array-like"]
ARRAY_LIKE_OR_NONE = ["array-like", None]
OBJECT_OR_NONE = "no_validation"  # duck-typed, deep-forwarded objects

VERBOSE = ["boolean", Interval(Integral, 0, None, closed="left"), str]

# ── shared fragments ─────────────────────────────────────────────────────────

#: Shared by every :class:`~pycsamt.agents._base.BaseAgent` subclass.
_BASE_AGENT_PARAMS = {
    "api_key": [str, None],
    "model": [str, None],
    "llm_provider": [StrOptions(_PROVIDERS)],
    "verbose": VERBOSE,
}

#: Correlated-field synthetic-training-physics knobs shared between
#: ``Inv2DAgent(physics="mt2d"/"mt2d_tri")`` and ``Inv3DAgent(physics="mt3d")``
#: — same :class:`~pycsamt.compat.sklearn.Interval` objects so "positive"
#: means the same thing in both agents.
_MAXWELL_TRAINING_PARAMS = {
    "log_resistivity_mean": [Real],
    "log_resistivity_std": POSITIVE_FLOAT,
    "mesh_safety_factor": POSITIVE_FLOAT,
    "max_mesh_cells": POSITIVE_INT,
}

#: Shared "how many synthetic profiles / how long to train" knobs.
_TRAINING_LOOP_PARAMS = {
    "n_train_profiles": POSITIVE_INT,
    "epochs": POSITIVE_INT,
    "patience": OPTIONAL_POSITIVE_INT,
    "freqs": ARRAY_LIKE_OR_NONE,
    "depth_max": OPTIONAL_POSITIVE_FLOAT,
}

# ── Inv2DAgent ────────────────────────────────────────────────────────────────

INV2D_PARAM_CONSTRAINTS = {
    **_BASE_AGENT_PARAMS,
    **_MAXWELL_TRAINING_PARAMS,
    **_TRAINING_LOOP_PARAMS,
    "physics": [StrOptions({"mt1d", "mt2d", "mt2d_tri"})],
    "n_depth": POSITIVE_INT,
    "n_freqs": POSITIVE_INT,
    "n_components": POSITIVE_INT,
    "arch": [str],
    "n_stations_per_profile": POSITIVE_INT,
    "station_spacing_m": POSITIVE_FLOAT,
    "correlation_length_x_m": ARRAY_LIKE,
    "correlation_length_z_m": ARRAY_LIKE,
    "lambda_x": NON_NEGATIVE_FLOAT,
    "lambda_z": NON_NEGATIVE_FLOAT,
    "lambda_tv": NON_NEGATIVE_FLOAT,
    "mesh_target_cell_m": POSITIVE_FLOAT,
    "field_grid_cell_m": POSITIVE_FLOAT,
    "topo_x_m": ARRAY_LIKE_OR_NONE,
    "topo_z_m": ARRAY_LIKE_OR_NONE,
    "gcn_hidden": ARRAY_LIKE,
    "gcn_adjacency_radius_m": POSITIVE_FLOAT,
    "mare2dem_adapter": OBJECT_OR_NONE,
}

# ── Inv3DAgent ────────────────────────────────────────────────────────────────

INV3D_PARAM_CONSTRAINTS = {
    **_BASE_AGENT_PARAMS,
    **_MAXWELL_TRAINING_PARAMS,
    **_TRAINING_LOOP_PARAMS,
    "physics": [StrOptions({"mt1d", "mt3d"})],
    "n_layers": POSITIVE_INT,
    "n_freqs": POSITIVE_INT,
    "radius": POSITIVE_FLOAT,
    "hidden": ARRAY_LIKE,
    "dropout": [Interval(Real, 0, 1, closed="left")],
    "n_mc": NON_NEGATIVE_INT,
    "correlation_length_x_m": ARRAY_LIKE,
    "correlation_length_y_m": ARRAY_LIKE,
    "correlation_length_z_m": ARRAY_LIKE,
    "cells_per_skin_depth": OPTIONAL_POSITIVE_FLOAT,
    "geology_grid_nx_ny": POSITIVE_INT,
    "geology_grid_nz": OPTIONAL_POSITIVE_INT,
}
