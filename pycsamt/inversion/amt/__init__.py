# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Audio-magnetotelluric inversion convenience API."""

from .inv1d import AMT1DInversion
from .inv2d import AMT2DInversion

__all__ = ["AMT1DInversion", "AMT2DInversion"]
