# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Central parameter contracts for the public Occam1D API.

Schemas in this module validate inexpensive, top-level properties at public
method boundaries. Scientific array contents and cross-parameter invariants
remain the responsibility of the owning object. This separation produces
clear type errors without duplicating numerical validation throughout the
subpackage.

The schemas are compatible with :func:`pycsamt.compat.sklearn.validate_params`
and intentionally contain no parsing or scientific processing logic.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from numbers import Integral, Real
from os import PathLike

from ...compat.sklearn import HasMethods, Interval, StrOptions

__all__ = [
    "BATCH_BUILD_SCHEMA",
    "DATA_SCHEMA",
    "ITERATION_SELECTOR_SCHEMA",
    "LOGGER_SCHEMA",
    "LOG_CONVERGENCE_SCHEMA",
    "LOG_ITERATION_SCHEMA",
    "LOG_READ_SCHEMA",
    "LOG_WRITE_SCHEMA",
    "MODE_OPTIONS",
    "MODEL_BUILD_SCHEMA",
    "MODEL_SCHEMA",
    "PATH_SCHEMA",
    "PLOT_SCHEMA",
    "RUNNER_RUN_SCHEMA",
    "RUNNER_SCHEMA",
    "RESPONSE_READ_SCHEMA",
    "RESPONSE_SELECT_SCHEMA",
    "RESPONSE_WRITE_SCHEMA",
    "STARTUP_SCHEMA",
]

MODE_OPTIONS = frozenset({"xy", "yx", "te", "tm", "determinant"})
"""Canonical lowercase response modes accepted at public boundaries."""

PATH_SCHEMA = {"path": [str, PathLike]}
"""Path contract shared by native readers and writers."""

LOG_READ_SCHEMA = {
    "path": [str, PathLike],
    "strict": ["boolean"],
}
"""Contract for :meth:`Occam1DLog.read`."""

LOG_WRITE_SCHEMA = {
    "path": [str, PathLike],
    "include_missing": ["boolean"],
}
"""Contract for tabular convergence-history export."""

LOG_CONVERGENCE_SCHEMA = {
    "target": [Interval(Real, 0, None, closed="neither")],
}
"""Contract for target-RMS convergence queries."""

LOG_ITERATION_SCHEMA = {
    "iteration": [Interval(Integral, 0, None, closed="left")],
    "all_matches": ["boolean"],
}
"""Contract for selecting a parsed iteration number."""

RESPONSE_READ_SCHEMA = {
    "path": [str, PathLike],
    "strict": ["boolean"],
}
"""Contract for :meth:`Occam1DResponse.read`."""

RESPONSE_WRITE_SCHEMA = {
    "path": [str, PathLike],
    "physical": ["boolean"],
}
"""Contract for normalized response-table export."""

RESPONSE_SELECT_SCHEMA = {
    "type_code": [
        Interval(Integral, 1, None, closed="left"),
        None,
    ],
    "frequency_index": [
        Interval(Integral, 1, None, closed="left"),
        None,
    ],
}
"""Contract for response-row selection."""

DATA_SCHEMA = {
    "frequency": [Iterable],
    "resistivity": [Iterable],
    "phase": [Iterable],
    "resistivity_error": [Iterable],
    "phase_error": [Iterable],
    "mode": [StrOptions(MODE_OPTIONS)],
    "station": [str],
}
"""Top-level contract for one sounding data container."""

MODEL_SCHEMA = {
    "depth": [Iterable],
    "resistivity": [Iterable],
    "penalty": [Iterable, None],
    "preference": [Iterable, None],
    "weight": [Iterable, None],
}
"""Top-level contract for an explicit vertical layer model."""

MODEL_BUILD_SCHEMA = {
    "n_layers": [Interval(Integral, 2, None, closed="left")],
    "first_thickness": [Interval(Real, 0, None, closed="neither")],
    "depth_max": [Interval(Real, 0, None, closed="neither")],
    "resistivity": [Interval(Real, 0, None, closed="neither")],
}
"""Boundary contract for geometric layer-model generation."""

STARTUP_SCHEMA = {
    "model_file": [str, PathLike],
    "data_file": [str, PathLike],
    "parameters": [Iterable],
    "iteration": [Interval(Integral, 0, None, closed="left")],
    "target_misfit": [Interval(Real, 0, None, closed="neither")],
    "max_iterations": [Interval(Integral, 1, None, closed="left")],
    "roughness_type": [Interval(Integral, 1, 2, closed="both")],
    "stepsize_cut_count": [Interval(Integral, 0, None, closed="left")],
    "debug_level": [Interval(Integral, 0, None, closed="left")],
    "lagrange_start": [Interval(Real, 0, None, closed="neither")],
    "extra_fields": [Mapping, None],
}
"""Boundary contract for native startup and iteration objects."""

RUNNER_RUN_SCHEMA = {
    "max_iterations": [
        Interval(Integral, 1, None, closed="left"),
        None,
    ],
    "target_misfit": [
        Interval(Real, 0, None, closed="neither"),
        None,
    ],
    "timeout": [Interval(Real, 0, None, closed="neither"), None],
    "check": ["boolean"],
    "extra_args": [Iterable, None],
}
"""Contract for synchronous external solver execution."""

RUNNER_SCHEMA = {
    "workdir": [str, PathLike],
    "binary_path": [str, PathLike, None],
    "startup_file": [str],
    "binary_name": [str],
}
"""Constructor contract for the Occam1D process runner."""

ITERATION_SELECTOR_SCHEMA = {
    "iteration": [
        Interval(Integral, 0, None, closed="left"),
        StrOptions({"best", "latest"}),
    ]
}
"""Contract for selecting result iterations."""

PLOT_SCHEMA = {
    "dpi": [Interval(Integral, 1, None, closed="left")],
    "figsize": [Iterable, None],
}
"""Shared lightweight plot-construction contract."""

BATCH_BUILD_SCHEMA = {
    "overrides": [Mapping],
}
"""Documentary contract for validated batch configuration overrides."""

LOGGER_SCHEMA = {
    "logger": [
        HasMethods(["debug", "info", "warning", "error", "exception"]),
        None,
    ]
}
"""Logger-like protocol used by package object constructors."""
