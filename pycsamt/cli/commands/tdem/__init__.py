# -*- coding: utf-8 -*-
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

from ._base import tdem  # noqa: F401

from . import info     # noqa: F401  (registers @tdem.command("info"))
from . import convert  # noqa: F401  (registers @tdem.command("convert"))
from . import plot     # noqa: F401  (registers @tdem.command("plot"))

__all__ = ["tdem"]
