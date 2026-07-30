# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt interp classify — build pseudo-stratigraphic logs from an inversion
workdir and print a lithology summary.

Workflow
--------
1. Auto-detect the inversion solver from workdir signature files.
2. Load the last (or requested) iteration via ``InversionResult``.
3. Adapt to :class:`~pycsamt.interp.ResistivityModel`.
4. Fit :class:`~pycsamt.interp.ModelCalibrator` (database-only, no boreholes).
5. Build :class:`~pycsamt.interp.StratigraphicLog` for every station.
6. Print a layer table or write it to CSV.

Usage
-----
::

    # basic run — auto-detects Occam2D from workdir files
    pycsamt interp classify occam_run/

    # request a specific iteration
    pycsamt interp classify occam_run/ --iteration 10

    # tighter autolayer merging, save to CSV
    pycsamt interp classify occam_run/ --merge-tol 0.1 --output logs.csv

    # JSON for downstream processing
    pycsamt interp classify occam_run/ --format json
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import interp

# Occam2D workdir signature files
_OCCAM_SIGNATURES = {
    "Occam2DMesh",
    "Occam2DModel",
    "OccamData.dat",
    "OccamStartup",
}


def _detect_solver(workdir: Path) -> str:
    """Return 'occam2d' if workdir has Occam2D files, else 'unknown'."""
    files = {f.name for f in workdir.iterdir() if f.is_file()}
    if files & _OCCAM_SIGNATURES:
        return "occam2d"
    return "unknown"


@interp.command("classify")
@click.argument(
    "workdir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="WORKDIR",
)
@click.option(
    "--solver",
    type=click.Choice(["occam2d", "auto"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Inversion solver.  'auto' detects from workdir files.",
)
@click.option(
    "--iteration",
    "-i",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Iteration number to load (default: last available).",
)
@click.option(
    "--ptol",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.10,
    show_default=True,
    metavar="FLOAT",
    help="Borehole-calibration tolerance as a fraction of true resistivity.",
)
@click.option(
    "--merge-tol",
    type=click.FloatRange(min=0.0),
    default=0.2,
    show_default=True,
    metavar="FLOAT",
    help=(
        "Log10-rho difference threshold for merging adjacent cells "
        "into one stratigraphic layer."
    ),
)
@click.option(
    "--max-depth",
    type=float,
    default=None,
    metavar="METRES",
    help="Clip the model at this depth (metres) before classifying.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(writable=True, path_type=Path),
    default=None,
    metavar="FILE",
    help="Save the layer table to this CSV file.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def classify(
    ctx: click.Context,
    workdir: Path,
    solver: str,
    iteration: int | None,
    ptol: float,
    merge_tol: float,
    max_depth: float | None,
    output: Path | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Build pseudo-stratigraphic logs from an inversion workdir.

    Loads the model, classifies each depth cell against the built-in
    rock-physics database, and merges adjacent cells that share the same
    lithology into discrete :class:`~pycsamt.interp.lithology.Layer` objects.

    \b
    Output columns:
      station   depth_top_m   depth_bot_m   thickness_m
      rock      rho_min       rho_max       description

    \b
    Examples:
      pycsamt interp classify occam_run/
      pycsamt interp classify occam_run/ --iteration 10
      pycsamt interp classify occam_run/ --merge-tol 0.1 --output logs.csv
      pycsamt interp classify occam_run/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    # ── detect solver ─────────────────────────────────────────────────────────
    effective_solver = solver.lower()
    if effective_solver == "auto":
        effective_solver = _detect_solver(workdir)
        if effective_solver == "unknown":
            raise click.UsageError(
                f"Could not detect an inversion solver from {workdir}.  "
                "Pass --solver explicitly (e.g. --solver occam2d)."
            )
        if verbose >= 1:
            click.echo(f"Auto-detected solver: {effective_solver}", err=True)

    if effective_solver != "occam2d":
        raise click.UsageError(
            f"Solver '{effective_solver}' is not yet supported by "
            "'pycsamt interp classify'.  Only occam2d is currently available."
        )

    # ── load model ────────────────────────────────────────────────────────────
    try:
        from ....interp import (  # noqa: PLC0415
            ModelCalibrator,
            ResistivityModel,
        )
        from ....models.occam2d.results import (
            InversionResult,  # noqa: PLC0415
        )

        if verbose >= 1:
            click.echo(
                f"Loading Occam2D result from {workdir}/"
                + (
                    f"  iteration={iteration}"
                    if iteration
                    else "  (last iteration)"
                ),
                err=True,
            )

        result = InversionResult(workdir, iteration=iteration)
        model = ResistivityModel.from_occam2d(result)

        if max_depth is not None:
            mask = model.z_centers <= max_depth
            model = ResistivityModel(
                x_centers=model.x_centers,
                z_centers=model.z_centers[mask],
                rho_2d=model.rho_2d[mask, :],
                station_x=model.station_x,
                station_names=model.station_names,
                method=model.method,
                rms=model.rms,
            )

    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error loading inversion result: {exc}", err=True)
        sys.exit(1)

    # ── classify ──────────────────────────────────────────────────────────────
    try:
        cal = ModelCalibrator(ptol=ptol, verbose=verbose >= 2).fit(model)
        logs = cal.stratigraphic_logs(merge_tolerance=merge_tol)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error during classification: {exc}", err=True)
        sys.exit(1)

    if verbose >= 1:
        click.echo(
            f"Classified {len(logs)} station(s)  "
            f"RMS={model.rms:.3f}  solver={effective_solver}",
            err=True,
        )

    # ── emit ──────────────────────────────────────────────────────────────────
    df = _logs_to_dataframe(logs)

    if output is not None:
        df.to_csv(output, index=False)
        click.echo(f"Layer table saved → {output}")

    _emit(df, output_format)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _logs_to_dataframe(logs):
    import pandas as pd  # noqa: PLC0415

    rows = []
    for log in logs:
        for layer in log.layers:
            rows.append(
                {
                    "station": log.station_name,
                    "depth_top_m": round(layer.depth_top, 2),
                    "depth_bot_m": round(layer.depth_bot, 2),
                    "thickness_m": round(layer.depth_bot - layer.depth_top, 2),
                    "rock": layer.rock.name,
                    "rho_min": layer.rock.rho_min,
                    "rho_max": layer.rock.rho_max,
                    "description": layer.rock.description,
                }
            )
    return pd.DataFrame(rows)


def _emit(df, output_format: str) -> None:
    if output_format == "json":
        click.echo(df.to_json(orient="records", indent=2, default_handler=str))
    elif output_format == "csv":
        click.echo(df.to_csv(index=False))
    else:
        click.echo(df.to_string(index=False))
