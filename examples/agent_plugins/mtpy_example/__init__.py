"""Example pyCSAMT-Agent domain plugin wrapping MTpy-style workflows.

Importing this package never imports MTpy; see :mod:`.plugin`.
"""
from .plugin import (
    MTPY_WORKFLOWS,
    MTpyPlanValidator,
    MTpyRegistry,
    make_plugin,
    mtpy_available,
)

__all__ = [
    "MTPY_WORKFLOWS",
    "MTpyPlanValidator",
    "MTpyRegistry",
    "make_plugin",
    "mtpy_available",
]
