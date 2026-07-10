# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Public API interface for pycsamt topography configuration.

This module re-exports the global singleton and its convenience
helpers so users can configure topography from the canonical
``pycsamt.api`` namespace::

    from pycsamt.api.topo import configure_topo, reset_topo, PYCSAMT_TOPO

See :mod:`pycsamt.topo.config` for the full :class:`TopoConfig`
dataclass documentation.
"""

from pycsamt.topo.config import (  # noqa: F401
    PYCSAMT_TOPO,
    TopoConfig,
    configure_topo,
    reset_topo,
)

__all__ = [
    "TopoConfig",
    "PYCSAMT_TOPO",
    "configure_topo",
    "reset_topo",
]
