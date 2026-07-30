# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.build
===========================

Click command group for compiling the external EM solvers (ModEM,
Occam2D, MARE2DEM). Importing this package registers every sub-command
on the ``build`` group. External code only needs::

    from pycsamt.cli.commands.build import build

Sub-command modules
--------------------
_base       Root Click group (``@click.group("build")``) and the shared
            script-locate-and-run dispatcher.
commands    Leaf commands: modem2d, modem3d, occam2d, mare2dem.
"""

# Importing the commands module triggers its @build.add_command() calls,
# which register each sub-command on the group defined in _base.
from . import commands  # noqa: F401
from ._base import build  # noqa: F401  (re-exported as package public API)

__all__ = ["build"]
