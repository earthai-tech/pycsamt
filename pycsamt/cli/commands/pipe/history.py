# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe history — list previously logged pipeline runs."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import pipe


@pipe.command("history")
@click.option(
    "--file",
    "history_file",
    type=click.Path(path_type=Path),
    default=None,
    metavar="FILE",
    help="Run-history log to read.  Defaults to ~/.pycsamt/pipeline_history.jsonl.",
)
@click.option(
    "--last",
    type=click.IntRange(min=1),
    default=None,
    metavar="N",
    help="Show only the N most recent runs.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def history(
    ctx: click.Context,
    history_file: Path | None,
    last: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """List runs previously logged with  pycsamt pipe run --history.

    Each run is logged as one compact JSON line — pipeline name, overall
    status, timing, site counts, and a per-step summary — appended by
    ``Pipeline.run(..., history=True)`` or ``pycsamt pipe run --history``.
    Logging is opt-in; this command shows nothing until a run has opted in.

    \b
    Examples:
      # List all logged runs
      pycsamt pipe history

      # Just the 5 most recent
      pycsamt pipe history --last 5

      # Machine-readable
      pycsamt pipe history --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.pipeline import load_history  # noqa: PLC0415

    records = load_history(history_file, last=last)

    if output_format == "json":
        click.echo(json.dumps(records, indent=2))
        return

    if output_format == "csv":
        click.echo("timestamp,pipeline_name,ok,n_errors,elapsed_sec,n_sites_in,n_sites_out")
        for r in records:
            click.echo(
                f"{r['timestamp']},{r['pipeline_name']},{r['ok']},"
                f"{r['n_errors']},{r['elapsed_sec']},"
                f"{r['n_sites_in']},{r['n_sites_out']}"
            )
        return

    # text
    if not records:
        click.echo("No pipeline runs logged yet.")
        click.echo("Run  pycsamt pipe run --history  to start logging.")
        return

    click.echo(f"Logged {len(records)} pipeline run(s):")
    for r in records:
        status = "OK" if r["ok"] else f"FAILED ({r['n_errors']} err)"
        click.echo(
            f"  {r['timestamp']}  {r['pipeline_name']:<20} {status:<14} "
            f"{r['elapsed_sec']:>7.2f}s  sites {r['n_sites_in']}→{r['n_sites_out']}"
        )
