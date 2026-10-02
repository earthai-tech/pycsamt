# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.airborne
==============================

Click command group for airborne EM datasets (ZTEM, MobileMT, AFMAG).

Importing this package registers every sub-command on the ``airborne``
group.  External code only needs::

    from pycsamt.cli.commands.airborne import airborne

Sub-command modules
--------------------
_base       Root Click group (``@click.group("airborne")``).
info        ``pycsamt airborne info``     — dataset/site summary report.
diagnose    ``pycsamt airborne diagnose`` — technology-aware diagnostic
            table (ZTEM divergence, MobileMT admittance, AFMAG/AirMt
            tilt).
"""

# Importing each sub-module triggers its @airborne.command decorator,
# which registers the sub-command on the group defined in _base.
from . import (
    diagnose,  # noqa: F401
    info,  # noqa: F401
)
from ._base import (
    airborne,  # noqa: F401  (re-exported as package public API)
)

__all__ = ["airborne"]
