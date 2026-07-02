# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.map
========================

Click command group for station-location and map-display workflows.

Importing this package registers every sub-command on the ``map`` group.
External code only needs::

    from pycsamt.cli.commands.map import map

Sub-command modules
-------------------
_base       Root Click group and shared station-coordinate helpers.
stations    ``pycsamt map stations`` - list/export station coordinates.
plot        ``pycsamt map plot`` - save a static station map figure.
"""

from ._base import map  # noqa: F401,A001  (re-exported as package public API)

from . import stations  # noqa: F401
from . import plot      # noqa: F401

__all__ = ["map"]
