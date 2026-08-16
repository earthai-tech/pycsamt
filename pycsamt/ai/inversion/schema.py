# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Parameter-validation schemas for AI inversion interfaces.

These schemas describe public call-boundary constraints for DUHI and
its geometry-aware mapper. Numerical array contents, broadcasting,
geometry consistency, and cross-parameter relationships are validated
by the implementation after the compatibility decorator accepts the
top-level argument types and scalar ranges.
"""

from __future__ import annotations

from numbers import Real

from ...compat.sklearn import Interval

__all__ = [
    "DUHI_INIT_SCHEMA",
    "DUHI_PREPARE_SCHEMA",
    "DIMENSIONALITY_RELIABILITY_SCHEMA",
    "GRID_MAPPING_SCHEMA",
    "OBSERVATION_RELIABILITY_SCHEMA",
    "RELIABILITY_COMBINATION_SCHEMA",
]


OBSERVATION_RELIABILITY_SCHEMA = {
    "errors": ["array-like"],
    "reliability": ["array-like"],
    "reliability_floor": [
        Interval(Real, 0, 1, closed="right"),
    ],
}
"""Constraints for reliability-weighted observation errors."""


DIMENSIONALITY_RELIABILITY_SCHEMA = {
    "beta_deg": ["array-like"],
    "beta_scale_deg": [Interval(Real, 0, None, closed="neither")],
    "minimum": [Interval(Real, 0, 1, closed="both")],
}
"""Constraints for phase-tensor dimensionality reliability."""


RELIABILITY_COMBINATION_SCHEMA = {
    "measurement": ["array-like"],
    "dimensionality": ["array-like"],
    "minimum": [Interval(Real, 0, 1, closed="both")],
}
"""Constraints for combined observation reliability."""


DUHI_INIT_SCHEMA = {
    "lambda_ai": [Interval(Real, 0, None, closed="left")],
    "sigma_ai_floor": [Interval(Real, 0, None, closed="neither")],
    "reliability_floor": [
        Interval(Real, 0, 1, closed="right"),
    ],
    "prejudice_filename": [str],
    "grid_mapper": [callable, None],
    "verbose": ["verbose"],
}
"""Constraints for :class:`~pycsamt.ai.inversion.DUHIInverter2D`."""


DUHI_PREPARE_SCHEMA = {
    "builder": [object],
    "ai_mean": ["array-like"],
    "ai_std": ["array-like"],
    "observation_reliability": ["array-like"],
    "ai_initialize": ["boolean"],
    "ai_x": ["array-like", None],
    "ai_z": ["array-like", None],
}
"""Constraints for ``DUHIInverter2D.prepare``."""


GRID_MAPPING_SCHEMA = {
    "grid": ["array-like"],
    "model": [object],
    "mesh": [object],
    "x_coordinates": ["array-like", None],
    "z_coordinates": ["array-like", None],
}
"""Constraints for geometry-aware AI-to-Occam grid mapping."""
