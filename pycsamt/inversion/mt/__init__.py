# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Magnetotelluric inversion convenience API."""

from .inv1d import MT1DInversion
from .inv2d import MT2DInversion
from .inv3d import MT3DInversion, SimPEGMT3DInversion

__all__ = [
    "MT1DInversion",
    "MT2DInversion",
    "MT3DInversion",
    "SimPEGMT3DInversion",
]
