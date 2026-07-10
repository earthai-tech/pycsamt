# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt interp export — export an interpreted EM model to an industry format.

Supported output formats
------------------------
* **xyz**  — Oasis Montaj XYZ point-data profile
* **las**  — CWLS LAS 2.0 well-log ASCII (one file per station)
* **csv**  — flat layer table (station, depth, rock, resistivity)
* **vtk**  — ASCII rectilinear-grid format for ParaView / QGIS

Workflow
--------
1. Detect / confirm solver from workdir.
2. Load the Occam2D ``InversionResult``.
3. Adapt to ``ResistivityModel``, calibrate, build stratigraphic logs.
4. Write output files to ``--output-dir``.

Usage
-----
::

    pycsamt interp export occam_run/ --format xyz --output-dir exports/
    pycsamt interp export occam_run/ --format csv
    pycsamt interp export occam_run/ --format las --output-dir las_logs/
    pycsamt interp export occam_run/ --format vtk --output-dir vtk_out/
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    overwrite_option,
    verbose_option,
)
from ._base import interp
from .classify import _detect_solver


@interp.command("export")
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
    "--iteration", "-i",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Iteration to load (default: last available).",
)
@click.option(
    "--format", "-f", "export_format",
    type=click.Choice(["xyz", "las", "csv", "vtk"], case_sensitive=False),
    required=True,
    help="Output format.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    default=Path("."),
    show_default=True,
    help="Directory where output files are written.",
)
@click.option(
    "--ptol",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.10,
    show_default=True,
    metavar="FLOAT",
    help="Borehole-calibration tolerance.",
)
@click.option(
    "--merge-tol",
    type=click.FloatRange(min=0.0),
    default=0.2,
    show_default=True,
    metavar="FLOAT",
    help="Log10-rho threshold for merging adjacent cells.",
)
@click.option(
    "--max-depth",
    type=float,
    default=None,
    metavar="METRES",
    help="Clip model at this depth before exporting.",
)
@overwrite_option
@verbose_option
@no_color_option
@click.pass_context
def export(
    ctx: click.Context,
    workdir: Path,
    solver: str,
    iteration: int | None,
    export_format: str,
    output_dir: Path,
    ptol: float,
    merge_tol: float,
    max_depth: float | None,
    overwrite: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Export an interpreted EM model to an industry-standard format.

    Loads the inversion result, builds pseudo-stratigraphic logs, then
    writes one or more output files into ``--output-dir``.

    \b
    Format details:
      xyz   Single profile file — Oasis Montaj / GIS-compatible XYZ
      las   One LAS 2.0 file per station (CWLS well-log format)
      csv   Flat layer table — all stations in a single file
      vtk   ASCII rectilinear VTK grid — for ParaView / QGIS

    \b
    Examples:
      pycsamt interp export occam_run/ --format xyz --output-dir exports/
      pycsamt interp export occam_run/ --format las --output-dir las_logs/
      pycsamt interp export occam_run/ --format csv
      pycsamt interp export occam_run/ --format vtk --output-dir vtk_out/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    # ── detect solver ─────────────────────────────────────────────────────────
    effective_solver = solver.lower()
    if effective_solver == "auto":
        effective_solver = _detect_solver(workdir)
        if effective_solver == "unknown":
            raise click.UsageError(
                f"Could not detect an inversion solver from {workdir}.  "
                "Pass --solver explicitly."
            )

    if effective_solver != "occam2d":
        raise click.UsageError(
            f"Solver '{effective_solver}' is not yet supported.  "
            "Only occam2d is currently available."
        )

    # ── load and classify ─────────────────────────────────────────────────────
    try:
        from ....interp import (  # noqa: PLC0415
            ModelCalibrator,
            ResistivityModel,
        )
        from ....models.occam2d.results import (
            InversionResult,  # noqa: PLC0415
        )

        result = InversionResult(workdir, iteration=iteration)
        model  = ResistivityModel.from_occam2d(result)

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

        cal  = ModelCalibrator(ptol=ptol).fit(model)
        logs = cal.stratigraphic_logs(merge_tolerance=merge_tol)

    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error loading / classifying model: {exc}", err=True)
        sys.exit(1)

    if verbose >= 1:
        click.echo(
            f"Classified {len(logs)} station(s)  "
            f"format={export_format}  output={output_dir}",
            err=True,
        )

    # ── write output ──────────────────────────────────────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)
    fmt = export_format.lower()

    try:
        from ....interp import (
            export as _export,  # noqa: PLC0415
        )

        if fmt == "xyz":
            out = output_dir / "profile.xyz"
            if out.exists() and not overwrite:
                raise click.UsageError(
                    f"{out} already exists.  Pass --overwrite to replace it."
                )
            _export.to_oasis_montaj_xyz(logs, str(out))
            click.echo(f"Written → {out}")

        elif fmt == "las":
            for log in logs:
                name = log.station_name.replace(" ", "_")
                out  = output_dir / f"{name}.las"
                if out.exists() and not overwrite:
                    raise click.UsageError(
                        f"{out} already exists.  Pass --overwrite to replace it."
                    )
                _export.to_las(log, str(out))
            click.echo(
                f"Written {len(logs)} LAS file(s) → {output_dir}/"
            )

        elif fmt == "csv":
            out = output_dir / "layers.csv"
            if out.exists() and not overwrite:
                raise click.UsageError(
                    f"{out} already exists.  Pass --overwrite to replace it."
                )
            _export.to_csv(logs, str(out))
            click.echo(f"Written → {out}")

        elif fmt == "vtk":
            out = output_dir / "model.vtk"
            if out.exists() and not overwrite:
                raise click.UsageError(
                    f"{out} already exists.  Pass --overwrite to replace it."
                )
            _export.to_vtk(model, str(out))
            click.echo(f"Written → {out}")

    except (click.exceptions.UsageError, click.exceptions.BadParameter):
        raise
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Export error: {exc}", err=True)
        sys.exit(1)
