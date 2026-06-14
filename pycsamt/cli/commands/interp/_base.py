# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.interp._base
==================================

Root Click group for all geophysical-geological interpretation
sub-commands.

Registered on the root ``main`` group via
:func:`~pycsamt.cli._base._register_commands`.

Sub-commands
------------
rocks       Browse the built-in rock-physics database or classify
            a single resistivity value  (``pycsamt interp rocks``).
classify    Load a 2-D EM inversion result, build pseudo-stratigraphic
            logs and print a lithology summary
            (``pycsamt interp classify``).
export      Export an interpreted model to industry formats
            (Oasis Montaj XYZ, LAS 2.0, CSV, VTK)
            (``pycsamt interp export``).
"""

from __future__ import annotations

import click


@click.group("interp")
@click.pass_context
def interp(ctx: click.Context) -> None:
    """Geological interpretation of 2-D EM resistivity models.

    \b
    Sub-commands:
      rocks     Browse / query the built-in rock-physics database
      classify  Build pseudo-stratigraphic logs from an inversion workdir
      export    Export interpretation results to XYZ, LAS, CSV, or VTK

    \b
    Examples:
      pycsamt interp rocks
      pycsamt interp rocks --rho 250
      pycsamt interp classify occam_run/ --format csv
      pycsamt interp classify occam_run/ --output logs.csv
      pycsamt interp export  occam_run/ --format xyz --output-dir exports/
    """
    ctx.ensure_object(dict)
