# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Spatially correlated geological priors for synthetic EM datasets.

Generators in this package will create complete 2-D sections and 3-D volumes
before responses are simulated.  This prevents production training labels
from being assembled by tiling unrelated 1-D models.

Planned primitives include correlated random fields, layers, dipping
interfaces, faults, lenses, resistive/conductive bodies, and topography.  All
generators must accept explicit random-number state and expose provenance.
"""

from .benchmark import (
    ID_BENCHMARK_FAMILIES,
    OOD_BENCHMARK_FAMILIES,
    BenchmarkGeology,
    generate_benchmark_geology,
)
from .fields import (
    CorrelatedField,
    DirectionalVariogram,
    GaussianCorrelation,
    GeologyGrid,
    directional_variogram,
    generate_gaussian_field,
)
from .layers import ElectricalLayer, LayeredGeology, generate_layered_geology
from .lenses import EllipsoidalLens, LensGeology, insert_lenses
from .topography import (
    TopographicSurface,
    interpolate_topography,
    topography_from_sites,
)

__all__ = [
    "GeologyGrid",
    "GaussianCorrelation",
    "CorrelatedField",
    "DirectionalVariogram",
    "generate_gaussian_field",
    "directional_variogram",
    "ID_BENCHMARK_FAMILIES",
    "OOD_BENCHMARK_FAMILIES",
    "BenchmarkGeology",
    "generate_benchmark_geology",
    "ElectricalLayer",
    "LayeredGeology",
    "generate_layered_geology",
    "EllipsoidalLens",
    "LensGeology",
    "insert_lenses",
    "TopographicSurface",
    "interpolate_topography",
    "topography_from_sites",
]
