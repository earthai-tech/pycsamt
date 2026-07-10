# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt site — root Click group and shared helpers.

Every sub-module imports ``site`` from here and registers its commands
via ``@site.command(...)`` or ``@site.group(...)``.  The ``__init__.py``
imports each sub-module to trigger the registrations.

Shared utility
--------------
_get_sites(explicit, survey_path, fresh, verbose) → Sites
    Single call that resolves any combination of positional path,
    --survey flag, or active context to a ``Sites`` object.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

if TYPE_CHECKING:
    from pycsamt.site.base import Sites


def _get_sites(
    explicit: Path | None,
    survey_path: Path | None,
    fresh: bool,
    verbose: int,
) -> Sites:
    """Resolve EDI source to a ``Sites`` object.

    Priority: explicit path > --survey flag > active context.
    Raises :class:`click.UsageError` when nothing is set.
    """
    from pycsamt.cli.survey import (
        resolve_survey,  # noqa: PLC0415
    )
    return resolve_survey(explicit or survey_path, fresh=fresh, verbose=verbose)


def _unwrap_edis(sites: Any) -> list:
    """Return a list of underlying EDIFile objects from a Sites collection.

    Compute functions in ``pycsamt.site.compute`` expect EDI-like objects
    that directly expose a ``.Z`` section.  ``Site`` objects wrap the EDI;
    unwrapping them gives the compute layer what it needs.
    """
    edis = []
    for s in sites:
        edis.append(s.edi if hasattr(s, "edi") else s)
    return edis


@click.group("site")
@click.pass_context
def site(ctx: click.Context) -> None:
    """Inspect, filter, edit, export and analyse MT sites.

    \b
    Typical workflow:
      pycsamt survey set   data/willy/
      pycsamt site info
      pycsamt site select  --stations S01,S05 --output-dir subset/
      pycsamt site edit    --rotate 30 --output-dir rotated/
      pycsamt site compute strike
      pycsamt site export  --output-dir final_edis/
    """
    ctx.ensure_object(dict)
