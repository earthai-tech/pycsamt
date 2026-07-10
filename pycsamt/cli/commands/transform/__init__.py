# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.transform
==============================

Click command group for all raw-data → EDI conversion sub-commands.

Sub-command modules
-------------------
_base      Root Click group.
spectra    ``pycsamt transform spectra`` — Spectra EDI → Impedance EDI.
avg        ``pycsamt transform avg``     — Zonge AVG   → Impedance EDI.
j          ``pycsamt transform j``       — Jones J     → Impedance EDI.
"""

from . import (
    avg,  # noqa: F401  (registers @transform.command("avg"))
    j,  # noqa: F401  (registers @transform.command("j"))
    spectra,  # noqa: F401  (registers @transform.command("spectra"))
)
from ._base import transform  # noqa: F401

__all__ = ["transform"]
