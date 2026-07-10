# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Time-domain EM inversion convenience API."""

from .inv1d import TDEM1DInversion
from .inv2d import TDEM2DInversion

__all__ = ["TDEM1DInversion", "TDEM2DInversion"]
