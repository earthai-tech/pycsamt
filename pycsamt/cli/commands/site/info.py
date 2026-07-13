# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt site info — display a rich report for a survey or single station."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    fresh_option,
    no_color_option,
    survey_option,
    verbose_option,
)
from ._base import _get_sites, site


@site.command("info")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--station",
    "-s",
    default=None,
    metavar="NAME",
    help="Show detailed report for a single station only.",
)
@click.option(
    "--top",
    "-n",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Limit the per-station table to the first N stations.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def info(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    station: str | None,
    top: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Display a rich report for a survey or a single station.

    SOURCE is optional when an active survey is set.
    Priority: SOURCE > --survey > active context.

    \b
    Examples:
      # collection report (active survey):
      pycsamt site info

      # explicit directory:
      pycsamt site info survey/

      # single station detail:
      pycsamt site info --station HBH03

      # limit table to first 20 stations:
      pycsamt site info --top 20

      # machine-readable output:
      pycsamt site info --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    sites = _get_sites(edi_source, survey_path, fresh, verbose)

    from pycsamt.site.report import (  # noqa: PLC0415
        SiteReport,
        SitesReport,
    )

    # --- single station ---
    if station is not None:
        match = None
        for s in sites:
            if s.name.upper() == station.upper():
                match = s
                break
        if match is None:
            click.echo(
                f"Station {station!r} not found.  "
                f"Available: {', '.join(s.name for s in sites)}",
                err=True,
            )
            sys.exit(1)

        rep = SiteReport(match)
        if output_format == "json":
            click.echo(json.dumps(rep.to_dict(), indent=2, default=str))
        elif output_format == "csv":
            df = rep.to_dataframe("resphase")
            click.echo(df.to_csv())
        else:
            rep.report()
        return

    # --- full collection ---
    rep = SitesReport(sites)
    if output_format == "json":
        click.echo(json.dumps(rep.to_dict(), indent=2, default=str))
    elif output_format == "csv":
        click.echo(rep.to_dataframe().to_csv(index=False))
    else:
        rep.report(top=top)
