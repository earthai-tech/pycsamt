# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt tdem convert — convert a TEMAVG survey folder to impedance EDI."""

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
from ._base import tdem


@tdem.command("convert")
@click.argument(
    "survey_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="SURVEY_DIR",
)
@output_dir_option
@click.option(
    "--method",
    type=click.Choice(["late_time", "fourier"], case_sensitive=False),
    default="late_time",
    show_default=True,
    help="Transform method: late-time approximation (fast) or Fourier cosine transform (rigorous).",
)
@click.option(
    "--component",
    default="Hz",
    show_default=True,
    metavar="COMP",
    help="TEMAVG component column to use for the sounding decay.",
)
@click.option(
    "--pattern",
    default="*.AVG",
    show_default=True,
    metavar="GLOB",
    help="Glob pattern to select processed-average files inside SURVEY_DIR.",
)
@click.option(
    "--stems",
    default=None,
    metavar="STEM[,STEM…]",
    help=(
        "Comma-separated file stems to convert (e.g. TEM100,TEM200).  "
        "When omitted, all discovered AVG files are converted."
    ),
)
@click.option(
    "--coord-file",
    "coordinate_file",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    metavar="FILE",
    help=(
        "Coordinate table to attach to each sounding.  "
        "When omitted, common filenames are discovered automatically."
    ),
)
@click.option(
    "--freq-convention",
    "freq_convention",
    type=click.Choice(["skin_depth", "fourier_time"], case_sensitive=False),
    default="skin_depth",
    show_default=True,
    help="Pseudo-frequency convention used by the late-time transform.",
)
@click.option(
    "--phase-mode",
    "phase_mode",
    type=click.Choice(["homogeneous", "weidelt"], case_sensitive=False),
    default="homogeneous",
    show_default=True,
    help="Phase model used when synthesising impedance from apparent resistivity.",
)
@click.option(
    "--no-geometry-correction",
    "loop_geometry_correction",
    is_flag=True,
    default=True,
    flag_value=False,
    help="Disable the in-loop Biot-Savart geometry correction.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Validate inputs and count soundings without writing EDI files.",
)
@verbose_option
@no_color_option
@click.pass_context
def convert(
    ctx: click.Context,
    survey_dir: Path,
    output_dir: Path,
    method: str,
    component: str,
    pattern: str,
    stems: str | None,
    coordinate_file: Path | None,
    freq_convention: str,
    phase_mode: str,
    loop_geometry_correction: bool,
    output_format: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Convert a TEMAVG survey folder to frequency-domain EDI files.

    SURVEY_DIR must contain Zonge TEMAVG processed-average files.
    The command reads every matching AVG file, builds one
    :class:`TEMSounding` per station, transforms each sounding to the
    frequency domain, and writes synthetic impedance EDI files.

    \b
    Examples:
      # convert all AVG files in a folder:
      pycsamt tdem convert data/TEMAVG/JIANGSU --output-dir imp/

      # only convert two profiles using the Fourier method:
      pycsamt tdem convert data/TEMAVG/JIANGSU \\
          --stems TEM100,TEM200 --method fourier --output-dir imp/

      # dry run — report soundings without writing:
      pycsamt tdem convert data/TEMAVG/JIANGSU --dry-run

      # machine-readable summary:
      pycsamt tdem convert data/TEMAVG/JIANGSU \\
          --output-dir imp/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    stems_list = [s.strip() for s in stems.split(",") if s.strip()] if stems else None

    if dry_run:
        from pycsamt.tdem.workflow import (
            read_temavg_soundings,  # noqa: PLC0415
        )

        soundings = read_temavg_soundings(
            survey_dir,
            pattern=pattern,
            coordinate_file=coordinate_file,
            stems=stems_list,
            component=component,
            verbose=verbose,
        )
        click.echo(f"Dry run — {len(soundings)} sounding(s) would be converted:")
        for snd in soundings:
            click.echo(f"  {snd.station_name}")
        return

    if str(output_dir) == ".":
        raise click.UsageError(
            "--output-dir is required (or use --dry-run to inspect soundings)"
        )

    from pycsamt.tdem.workflow import (
        transform_temavg_survey,  # noqa: PLC0415
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = transform_temavg_survey(
            survey_dir,
            pattern=pattern,
            coordinate_file=coordinate_file,
            stems=stems_list,
            component=component,
            method=method,
            freq_convention=freq_convention,
            phase_mode=phase_mode,
            loop_geometry_correction=loop_geometry_correction,
            return_collection=False,
            savepath=output_dir,
            verbose=verbose,
        )
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Conversion failed: {exc}", err=True)
        sys.exit(1)

    n_soundings = result.n_soundings
    written = result.written_paths
    n_written = len(written)

    summary = {
        "survey_dir": str(survey_dir),
        "output_dir": str(output_dir),
        "method": method,
        "n_soundings": n_soundings,
        "n_written": n_written,
        "written": written,
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2))
    else:
        click.echo(
            f"Soundings : {n_soundings}"
            f"  | Written: {n_written}"
            f"  | Output : {output_dir}/"
        )
        if verbose >= 1:
            for p in written:
                click.echo(f"  {p}")

    if n_written == 0 and n_soundings > 0:
        click.echo(
            "Warning: soundings were transformed but no EDI files were written.",
            err=True,
        )
        sys.exit(1)
