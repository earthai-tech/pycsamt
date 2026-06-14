# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.site
=========================

Click command group for all site inspection, selection, editing, export
and computation sub-commands.

Importing this package registers every sub-command on the ``site`` group.
External code only needs::

    from pycsamt.cli.commands.site import site

Sub-command modules
-------------------
_base       Root Click group + _get_sites() shared resolver.
info        ``pycsamt site info``    — rich survey / single-station report.
select      ``pycsamt site select``  — filter by name, freq, bbox, quality.
edit        ``pycsamt site edit``    — rotate, slice, fill, re-coordinate.
export      ``pycsamt site export``  — write EDIs or zip archive.
compute     ``pycsamt site compute`` — strike, resistivity, phase-slope, tipper.
"""

from ._base import site  # noqa: F401  (re-exported as package public API)

# Importing each sub-module triggers @site.command / @site.group decorators,
# which register the commands on the group defined in _base.
from . import info     # noqa: F401
from . import select   # noqa: F401
from . import edit     # noqa: F401
from . import export   # noqa: F401
from . import compute  # noqa: F401

__all__ = ["site"]
