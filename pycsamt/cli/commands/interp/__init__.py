# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.interp
============================

Click command group for geological interpretation of 2-D EM resistivity
models.

Importing this package registers every sub-command on the ``interp``
group.  External code only needs::

    from pycsamt.cli.commands.interp import interp

Sub-command modules
-------------------
_base       Root Click group (``@click.group("interp")``).
rocks       ``pycsamt interp rocks``    — browse / query the built-in
            rock-physics database.
classify    ``pycsamt interp classify`` — load a 2-D EM inversion result
            and build pseudo-stratigraphic logs.
export      ``pycsamt interp export``   — export an interpreted model to
            Oasis Montaj XYZ, LAS 2.0, CSV, or VTK.
"""

from ._base import interp  # noqa: F401  (re-exported as package public API)

# Importing each sub-module triggers its @interp.command / @interp.group
# decorator, which registers the sub-command on the group defined in _base.
from . import rocks     # noqa: F401
from . import classify  # noqa: F401
from . import export    # noqa: F401

__all__ = ["interp"]
