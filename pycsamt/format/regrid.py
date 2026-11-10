# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regridding of unstructured meshes onto 2D grids.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.algorithms.regrid`), the reference implementation of the
format. This module re-exports it under its pyCSAMT name, so that
existing imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.algorithms.regrid : Implementation.
"""

from pcsf.algorithms.regrid import (
    mesh_to_grid2d,
)

__all__ = [
    "mesh_to_grid2d",
]
