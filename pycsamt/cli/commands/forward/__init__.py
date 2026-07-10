# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.forward
============================

Click command group for all forward-modelling sub-commands.

Importing this package registers every sub-command on the ``forward``
group.  External code only needs::

    from pycsamt.cli.commands.forward import forward

Sub-command modules
-------------------
_base       Root Click group (``@click.group("forward")``).
model       ``pycsamt forward model geology`` — list / inspect geological
            scenarios from the built-in catalog.
run         ``pycsamt forward run``  — run a single 1-D EM forward model
            and print the response table.
generate    ``pycsamt forward generate`` — batch-generate synthetic
            (model, response) training datasets.
"""

# Importing each sub-module triggers its @forward.command / @forward.group
# decorator, which registers the sub-command on the group defined in _base.
from . import (
    generate,  # noqa: F401
    model,  # noqa: F401
    run,  # noqa: F401
)
from ._base import (
    forward,  # noqa: F401  (re-exported as package public API)
)

__all__ = ["forward"]
