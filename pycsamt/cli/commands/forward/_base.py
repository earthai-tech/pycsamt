# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.forward._base
===================================

Root Click group for all forward-modelling sub-commands.

Registered on the root ``main`` group via
:func:`~pycsamt.cli._base._register_commands`.

Sub-commands
------------
model       Inspect geological scenarios (``pycsamt forward model geology``).
run         Run a single 1-D EM forward model (``pycsamt forward run``).
generate    Batch-generate synthetic training datasets
            (``pycsamt forward generate``).
"""

from __future__ import annotations

import click


@click.group("forward")
@click.pass_context
def forward(ctx: click.Context) -> None:
    """Run 1-D EM forward models and generate synthetic datasets.

    \b
    Sub-commands:
      model       Inspect / list geological model scenarios
      run         Run a single 1-D forward model and print the response
      generate    Batch-generate (model, response) training datasets

    \b
    Examples:
      pycsamt forward model geology
      pycsamt forward model geology --name sedimentary
      pycsamt forward run --resistivities 100,10,500 --thicknesses 300,800
      pycsamt forward run --geology geothermal --solver mt --freqs 0.01:1000
      pycsamt forward generate --n-samples 5000 --output mt1d_train.npz
    """
    ctx.ensure_object(dict)
