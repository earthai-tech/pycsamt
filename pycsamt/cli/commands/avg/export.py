# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt avg export — write a Zonge AVG file to modern, legacy, or netCDF."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    output_dir_option,
    overwrite_option,
    verbose_option,
)
from ._base import avg, _get_avg


@avg.command("export")
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="AVG_FILE",
)
@click.option(
    "--format", "export_format",
    type=click.Choice(["modern", "legacy", "xarray"], case_sensitive=False),
    default="modern",
    show_default=True,
    help=(
        "Output format.  "
        "modern: CSV-based kind-2 AVG (default).  "
        "legacy: fixed-width kind-1 AVG.  "
        "xarray: netCDF4 via xarray (requires xarray + netcdf4)."
    ),
)
@click.option(
    "--stem",
    default=None,
    metavar="NAME",
    help=(
        "Output file base name without extension.  "
        "Default: source stem + '_exported'."
    ),
)
@output_dir_option
@overwrite_option
@click.option(
    "--format-out", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
    help="CLI output format (not the AVG file format).",
)
@verbose_option
@no_color_option
@click.pass_context
def export(
    ctx: click.Context,
    source: Path,
    export_format: str,
    stem: str | None,
    output_dir: Path,
    overwrite: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Write a Zonge AVG file to a different format.

    Use ``--format modern`` (default) to produce a clean CSV-based kind-2
    AVG file.  Use ``--format legacy`` for fixed-width kind-1 output.
    Use ``--format xarray`` to write a netCDF4 dataset (requires *xarray*
    and *netcdf4*).

    \b
    Examples:
      # re-export to modern format in a new directory:
      pycsamt avg export data/avg/K2.AVG --output-dir clean/

      # convert modern to legacy:
      pycsamt avg export data/avg/K2.AVG --format legacy --output-dir legacy/

      # write netCDF4 (xarray required):
      pycsamt avg export data/avg/K2.AVG --format xarray --output-dir nc/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    obj = _get_avg(source, verbose=verbose)

    out_stem = stem or (source.stem + "_exported")
    ext_map  = {"modern": ".avg", "legacy": ".avg", "xarray": ".nc"}
    ext      = ext_map[export_format]
    out_path = output_dir / (out_stem + ext)

    output_dir.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not overwrite:
        raise click.UsageError(
            f"{out_path} already exists.  Use --overwrite to replace it."
        )

    try:
        if export_format == "xarray":
            ds = obj.to_xarray()
            ds.to_netcdf(str(out_path))
        else:
            fmt_arg = "modern" if export_format == "modern" else "legacy"
            obj.write(out_path, fmt=fmt_arg)
    except ImportError as exc:
        click.echo(
            f"Export to {export_format!r} requires an optional dependency:\n"
            f"  {exc}\n\n"
            "Install it with:  pip install xarray netcdf4",
            err=True,
        )
        sys.exit(1)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Export failed: {exc}", err=True)
        sys.exit(1)

    summary = {
        "source":  str(source),
        "output":  str(out_path),
        "format":  export_format,
        "size_kb": round(out_path.stat().st_size / 1024, 1) if out_path.exists() else None,
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2, default=str))
    else:
        sz = summary["size_kb"]
        click.echo(
            f"Exported  {export_format}  →  {out_path}"
            + (f"  ({sz} KB)" if sz else "")
        )
