# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.tdem
=========================

Click command group for all time-domain EM (TEMAVG) sub-commands.

Sub-command modules
-------------------
_base      Root Click group + ``_get_survey()`` shared helper.
info       ``pycsamt tdem info``    — inspect a TEMAVG survey folder.
convert    ``pycsamt tdem convert`` — convert TEMAVG → impedance EDI.
plot       ``pycsamt tdem plot``    — visualise decay / section / map.
"""

from . import (
    convert,  # noqa: F401  (registers @tdem.command("convert"))
    info,  # noqa: F401  (registers @tdem.command("info"))
    plot,  # noqa: F401  (registers @tdem.command("plot"))
)
from ._base import tdem  # noqa: F401

__all__ = ["tdem"]
