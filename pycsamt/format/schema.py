# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""In-memory PCSF schema: geometries, stations and topography.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.core.schema`), the reference implementation of the format.
This module re-exports it under its pyCSAMT name, so that existing
imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.core.schema : Implementation.
"""

from pcsf.core.schema import (
    DERIVATION_METHODS,
    GEOMETRY_KINDS,
    PCSF_VERSION,
    RESISTIVITY_ENCODINGS,
    RESISTIVITY_UNIT,
    TOPOGRAPHY_KINDS,
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    StationTable,
    TopographyPerStation,
    TopographyRaster,
    UnstructuredMeshGeometry,
)

__all__ = [
    "DERIVATION_METHODS",
    "DerivedVolume",
    "GEOMETRY_KINDS",
    "Grid2DGeometry",
    "Grid3DGeometry",
    "LineEntry",
    "MultilineGeometry",
    "PCSFModel",
    "PCSF_VERSION",
    "RESISTIVITY_ENCODINGS",
    "RESISTIVITY_UNIT",
    "StationTable",
    "TOPOGRAPHY_KINDS",
    "TopographyPerStation",
    "TopographyRaster",
    "UnstructuredMeshGeometry",
]
