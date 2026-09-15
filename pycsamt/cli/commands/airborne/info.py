# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt airborne info — report on an airborne EM dataset or one site."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import _get_asites, airborne


@airborne.command("info")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="SOURCE",
)
@click.option(
    "--site",
    "-s",
    "site_name",
    default=None,
    metavar="NAME",
    help="Show detailed report for a single site only.",
)
@click.option(
    "--no-recursive",
    is_flag=True,
    default=False,
    help="When SOURCE is a directory, do not search it recursively.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def info(
    ctx: click.Context,
    source: Path,
    site_name: str | None,
    no_recursive: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Report on an airborne EM dataset (ZTEM, MobileMT, AFMAG).

    SOURCE is a single EMTF-XML file or a directory of them.

    \b
    Examples:
      pycsamt airborne info survey_xml/
      pycsamt airborne info survey_xml/ --site L100_0500
      pycsamt airborne info survey_xml/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    try:
        sites = _get_asites(source, not no_recursive, verbose)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error loading airborne sites: {exc}", err=True)
        sys.exit(1)

    if len(sites) == 0:
        click.echo(f"No airborne sites found under {source}.", err=True)
        sys.exit(1)

    # --- single site ---
    if site_name is not None:
        match = sites.get(site_name)
        if match is None:
            names = ", ".join(s.name for s in sites)
            click.echo(
                f"Site {site_name!r} not found.  Available: {names}",
                err=True,
            )
            sys.exit(1)

        summary = match.summary()
        if output_format == "json":
            click.echo(json.dumps(summary, indent=2, default=str))
        elif output_format == "csv":
            import pandas as pd  # noqa: PLC0415

            click.echo(pd.DataFrame([summary]).to_csv(index=False))
        else:
            click.echo(repr(match))
            for key, value in summary.items():
                click.echo(f"  {key:<12} {value}")
        return

    # --- full collection ---
    summaries = [s.summary() for s in sites]

    if output_format == "json":
        payload = {
            "n_sites": len(sites),
            "technologies": list(sites.technologies),
            "line_ids": list(sites.line_ids),
            "sites": summaries,
        }
        click.echo(json.dumps(payload, indent=2, default=str))
        return

    import pandas as pd  # noqa: PLC0415

    df = pd.DataFrame(summaries)
    if output_format == "csv":
        click.echo(df.to_csv(index=False))
        return

    click.echo(f"Sites:        {len(sites)}")
    click.echo(f"Technologies: {', '.join(sites.technologies) or '-'}")
    click.echo(f"Lines:        {', '.join(sites.line_ids) or '-'}")
    click.echo()
    click.echo(df.to_string(index=False))
