# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Parameter-validation schemas for the Occam2D public API.

The constants in this module are consumed by
:func:`pycsamt.compat.sklearn.validate_params`. Keeping the schemas
separate from the file containers makes the accepted public argument
types and scalar ranges reusable and straightforward to inspect.

Array contents and relationships between parameters remain the
responsibility of the corresponding Occam2D object. The compatibility
decorator validates only the public call boundary.
"""

from __future__ import annotations

from collections.abc import Iterable
from numbers import Integral
from os import PathLike

from ...compat.sklearn import Interval
from .config import OccamConfig

__all__ = [
    "OCCAM_PREJUDICE_DENSE_SCHEMA",
    "OCCAM_PREJUDICE_INIT_SCHEMA",
    "OCCAM_PREJUDICE_PARAMETER_COUNT_SCHEMA",
    "OCCAM_PREJUDICE_READ_SCHEMA",
    "OCCAM_PREJUDICE_WRITE_SCHEMA",
    "OCCAM_RUNNER_FORWARD_SCHEMA",
]


OCCAM_PREJUDICE_INIT_SCHEMA = {
    "parameter_indices": [Iterable, None],
    "target_values": [Iterable, None],
    "weights": [Iterable, None],
    "config": [OccamConfig, None],
}
"""Constraints for :class:`~pycsamt.models.occam2d.OccamPrejudice`."""


OCCAM_PREJUDICE_DENSE_SCHEMA = {
    "target_values": [Iterable],
    "weights": [Iterable],
    "include_zero_weight": ["boolean"],
    "config": [OccamConfig, None],
}
"""Constraints for ``OccamPrejudice.from_dense``."""


OCCAM_PREJUDICE_PARAMETER_COUNT_SCHEMA = {
    "n_params": [Interval(Integral, 1, None, closed="left")],
}
"""Constraints for model-size checks on prejudice records."""


OCCAM_PREJUDICE_READ_SCHEMA = {
    "path": [str, PathLike],
    "config": [OccamConfig, None],
}
"""Constraints for reading an Occam prejudice file."""


OCCAM_PREJUDICE_WRITE_SCHEMA = {
    "path": [str, PathLike],
}
"""Constraints for writing an Occam prejudice file."""


OCCAM_RUNNER_FORWARD_SCHEMA = {
    "output_root": [str],
    "auto_compile": ["boolean"],
}
"""Constraints for an Occam2D forward-only solver run."""
