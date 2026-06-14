# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt site select — filter sites and optionally write the subset to disk."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    fresh_option,
    no_color_option,
    output_dir_option,
    overwrite_option,
    survey_option,
    verbose_option,
)
from ....api.cli.params import FreqRange, StationList

from ._base import site, _get_sites


@site.command("select")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
# --- filter options ---
@click.option(
    "--stations",
    type=StationList(),
    default=None,
    metavar="S01,S02,…",
    help="Keep only these station IDs (comma-separated).",
)
@click.option(
    "--freq", "freq_range",
    type=FreqRange(),
    default=None,
    metavar="FMIN:FMAX",
    help="Keep only frequencies in this Hz range.",
)
@click.option(
    "--bbox",
    default=None,
    metavar="LAT_MIN,LON_MIN,LAT_MAX,LON_MAX",
    help="Keep only stations inside this geographic bounding box.",
)
@click.option(
    "--drop-empty",
    is_flag=True,
    default=False,
    help="Remove stations that carry no data arrays.",
)
@click.option(
    "--keep-finite",
    is_flag=True,
    default=False,
    help="Remove stations whose Z tensor contains no finite values.",
)
@click.option(
    "--phase-err-thresh",
    type=float,
    default=None,
    metavar="DEGREES",
    help="Mask frequencies where phase error exceeds this threshold.",
)
# --- output ---
@output_dir_option
@overwrite_option
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show what would be selected without writing any files.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def select(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    stations: list[str] | None,
    freq_range: tuple[float, float] | None,
    bbox: str | None,
    drop_empty: bool,
    keep_finite: bool,
    phase_err_thresh: float | None,
    output_dir: Path,
    overwrite: bool,
    dry_run: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Filter sites by name, frequency range, bounding box, or data quality.

    Without --output-dir, prints a summary of the selected subset.
    With --output-dir, writes the filtered EDIs to that directory.

    \b
    Examples:
      # filter by station name:
      pycsamt site select --stations S01,S02,S03

      # filter by frequency range:
      pycsamt site select --freq 0.1:1000 --output-dir filtered/

      # geographic box (lat_min,lon_min,lat_max,lon_max):
      pycsamt site select --bbox 27.0,101.0,29.0,103.0

      # combined filters, write output:
      pycsamt site select --stations S01,S05 --freq 1:100 \\
              --drop-empty --output-dir subset/

      # dry run — show what would be selected:
      pycsamt site select --freq 0.1:1000 --dry-run
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    from pycsamt.site import selection as sel  # noqa: PLC0415

    # --- apply filters in order ---
    result = sites_obj

    if stations is not None:
        result = sel.by_names(result, stations)
        if verbose >= 1:
            click.echo(
                f"  → {len(result)} station(s) after name filter", err=True
            )

    if freq_range is not None:
        result = sel.by_freq(result, freq_range[0], freq_range[1])
        if verbose >= 1:
            click.echo(
                f"  → {len(result)} station(s) after freq filter", err=True
            )

    if bbox is not None:
        try:
            parts = [float(x.strip()) for x in bbox.split(",")]
            lat_min, lon_min, lat_max, lon_max = parts
        except (ValueError, TypeError):
            raise click.BadParameter(
                "Expected LAT_MIN,LON_MIN,LAT_MAX,LON_MAX",
                param_hint="--bbox",
            )
        result = sel.by_bbox(result, lat_min, lon_min, lat_max, lon_max)
        if verbose >= 1:
            click.echo(
                f"  → {len(result)} station(s) after bbox filter", err=True
            )

    if drop_empty:
        result = sel.drop_empty(result)
        if verbose >= 1:
            click.echo(
                f"  → {len(result)} station(s) after drop_empty", err=True
            )

    if keep_finite:
        result = sel.keep_finite_z(result)
        if verbose >= 1:
            click.echo(
                f"  → {len(result)} station(s) after keep_finite", err=True
            )

    if phase_err_thresh is not None:
        result = sel.mask_large_phase_err(result, phase_err_thresh)
        if verbose >= 1:
            click.echo(
                f"  → {len(result)} station(s) after phase-err mask", err=True
            )

    n_out = len(result)
    n_in  = len(sites_obj)

    if n_out == 0:
        click.echo(
            "No stations remain after filtering.  "
            "Review your filter criteria.",
            err=True,
        )
        sys.exit(1)

    # --- summary ---
    names = [s.name for s in result]
    summary: dict[str, Any] = {
        "n_input":    n_in,
        "n_selected": n_out,
        "stations":   names,
    }

    if dry_run:
        if output_format == "json":
            click.echo(json.dumps(summary, indent=2))
        else:
            click.echo(f"Dry run — {n_out}/{n_in} station(s) would be selected:")
            for nm in names:
                click.echo(f"  {nm}")
        return

    # --- write to disk if requested ---
    if str(output_dir) != ".":
        output_dir.mkdir(parents=True, exist_ok=True)
        from pycsamt.site.export import write_sites  # noqa: PLC0415
        written = write_sites(result, output_dir, exist_ok=overwrite)
        summary["written"] = [str(p) for p in written]
        if verbose >= 1:
            click.echo(
                f"  Wrote {len(written)} EDI file(s) to {output_dir}/",
                err=True,
            )

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2, default=str))
    elif output_format == "csv":
        import pandas as pd  # noqa: PLC0415
        df = pd.DataFrame({"station": names})
        click.echo(df.to_csv(index=False))
    else:
        click.echo(
            f"Selected {n_out}/{n_in} station(s): "
            + ", ".join(names[:10])
            + ("…" if n_out > 10 else "")
        )
        if "written" in summary:
            click.echo(f"Written to: {output_dir}/")
