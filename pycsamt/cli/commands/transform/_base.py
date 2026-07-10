# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt transform — root Click group.

Each sub-module registers its commands via ``@transform.command(...)``
and is imported by ``__init__.py`` to trigger the registrations.
"""

from __future__ import annotations

import click


@click.group("transform")
@click.pass_context
def transform(ctx: click.Context) -> None:
    """Convert raw MT data files to impedance EDI format.

    \b
    Supported conversions:
      spectra    Spectra EDI  →  Impedance EDI  (with Tipper)
      avg        Zonge AVG    →  Impedance EDI
      j          Jones J-file →  Impedance EDI

    \b
    Examples:
      pycsamt transform spectra site.edi --output-dir imp/
      pycsamt transform spectra spectra_dir/ --output-dir imp/ --estimate-error
      pycsamt transform avg survey.avg    --output-dir imp/
      pycsamt transform j   station.j    --output-dir imp/
    """
    ctx.ensure_object(dict)
