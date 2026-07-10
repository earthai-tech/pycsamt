# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt transform avg — convert Zonge AVG files to impedance EDI."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    output_dir_option,
    verbose_option,
)
from ._base import transform


@transform.command("avg")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="FILE_OR_DIR",
)
@output_dir_option
@click.option(
    "--station-name",
    default=None,
    metavar="NAME",
    help="Override the station name (single-file input only).",
)
@click.option(
    "--format", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Validate inputs without writing files.",
)
@verbose_option
@no_color_option
@click.pass_context
def avg(
    ctx: click.Context,
    source: Path,
    output_dir: Path,
    station_name: str | None,
    output_format: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Convert a Zonge AVG file (or directory) to Impedance EDI.

    \b
    Examples:
      pycsamt transform avg survey.avg      --output-dir imp/
      pycsamt transform avg avg_dir/        --output-dir imp/
      pycsamt transform avg survey.avg      --dry-run
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    files = sorted(source.glob("*.avg")) if source.is_dir() else [source]
    if not files:
        click.echo(f"No .avg files found in {source}.", err=True)
        sys.exit(1)

    if dry_run:
        click.echo(f"Dry run — {len(files)} AVG file(s) would be converted:")
        for f in files:
            click.echo(f"  {f.name}")
        return

    if str(output_dir) == ".":
        raise click.UsageError("--output-dir is required")

    from pycsamt.transformers import AVGtoEDI  # noqa: PLC0415

    output_dir.mkdir(parents=True, exist_ok=True)
    edis = []
    failures = []

    for f in files:
        try:
            col = AVGtoEDI().transform(f, name=station_name)
            for ed in col:
                try:
                    ed.write(savepath=str(output_dir))
                except Exception as we:  # noqa: BLE001
                    failures.append({"source": str(f), "error": str(we)})
                    continue
                edis.append({"station": getattr(ed, "station", "?"),
                             "source": f.name})
        except Exception as exc:  # noqa: BLE001
            failures.append({"source": str(f), "error": str(exc)})

    summary = {
        "n_input": len(files),
        "n_ok": len(edis),
        "n_fail": len(failures),
        "converted": edis,
        "failures": failures,
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2))
    else:
        click.echo(
            f"Converted: {len(edis)}/{len(files)}  "
            f"| Errors: {len(failures)}  "
            f"| Output: {output_dir}/"
        )
        if failures:
            for r in failures:
                click.echo(f"  ✗  {r['source']}: {r['error']}", err=True)

    if failures and not edis:
        sys.exit(1)
