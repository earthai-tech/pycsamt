# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt tdem plot — visualise TEMAVG decay curves, sections, and maps."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import no_color_option, verbose_option

from ._base import tdem, _get_survey

_KINDS = [
    "decay",
    "rho",
    "section",
    "z-section",
    "map",
    "overview",
    "gate-profile",
    "elevation",
    "dashboard",
]


@tdem.command("plot")
@click.argument(
    "survey_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="SURVEY_DIR",
)
@click.option(
    "--kind", "-k",
    type=click.Choice(_KINDS, case_sensitive=False),
    default="decay",
    show_default=True,
    help=(
        "Plot type to produce.  "
        "decay: time-domain decay curves.  "
        "rho: apparent-resistivity section (requires transformation).  "
        "section: TEMAVG gate-value pseudo-section.  "
        "z-section: Z-contour pseudo-section (requires .Z files).  "
        "map: station location map.  "
        "overview: multi-panel survey overview.  "
        "gate-profile: gate values vs station along a profile.  "
        "elevation: station elevation profile.  "
        "dashboard: four-panel TDEM dashboard."
    ),
)
@click.option(
    "--stems",
    default=None,
    metavar="STEM[,STEM…]",
    help="Comma-separated AVG file stems to include (default: all).",
)
@click.option(
    "--pattern",
    default="*.AVG",
    show_default=True,
    metavar="GLOB",
    help="Glob pattern to locate AVG files.",
)
@click.option(
    "--component",
    default="Hz",
    show_default=True,
    metavar="COMP",
    help="TEMAVG component column for sounding extraction.",
)
@click.option(
    "--output-dir", "-o",
    "output_dir",
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    default=None,
    help=(
        "Directory where the figure is saved.  "
        "When omitted the figure is shown interactively."
    ),
)
@click.option(
    "--fmt",
    type=click.Choice(["png", "pdf", "svg", "jpg"], case_sensitive=False),
    default="png",
    show_default=True,
    help="Output file format (only relevant when --output-dir is set).",
)
@click.option(
    "--dpi",
    type=click.IntRange(min=50),
    default=150,
    show_default=True,
    help="Figure resolution in dots per inch.",
)
@verbose_option
@no_color_option
@click.pass_context
def plot(
    ctx: click.Context,
    survey_dir: Path,
    kind: str,
    stems: str | None,
    pattern: str,
    component: str,
    output_dir: Path | None,
    fmt: str,
    dpi: int,
    verbose: int,
    no_color: bool,
) -> None:
    """Visualise TEMAVG decay data, apparent-ρ sections, or station maps.

    SURVEY_DIR must contain Zonge TEMAVG processed-average files.  The
    ``--kind`` option selects which visualisation to produce; use
    ``--output-dir`` to save instead of displaying interactively.

    \b
    Available kinds
    ---------------
    decay        Time-domain decay curves for each station.
    rho          Apparent-resistivity vs frequency (per station).
    section      TEMAVG gate-value pseudo-section.
    z-section    Z-contour pseudo-section (needs companion .Z files).
    map          Station location map with profile lines.
    overview     Multi-panel survey overview (map + elevation).
    gate-profile Gate values vs station position along a profile.
    elevation    Elevation profile along the survey traverse.
    dashboard    Four-panel TDEM dashboard (AVG + Z + decay).

    \b
    Examples:
      # show decay curves interactively:
      pycsamt tdem plot data/TEMAVG/JIANGSU --kind decay

      # save apparent-resistivity section as PNG:
      pycsamt tdem plot data/TEMAVG/JIANGSU --kind rho \\
          --output-dir figs/ --dpi 200

      # save dashboard for specific AVG files:
      pycsamt tdem plot data/TEMAVG/JIANGSU --kind dashboard \\
          --stems TEM100,TEM200 --output-dir figs/

      # station location map:
      pycsamt tdem plot data/TEMAVG/JIANGSU --kind map --output-dir figs/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    import matplotlib  # noqa: PLC0415

    if output_dir is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    stems_list = (
        [s.strip() for s in stems.split(",") if s.strip()]
        if stems else None
    )

    survey = _get_survey(survey_dir, pattern=pattern, verbose=verbose)

    if stems_list is not None:
        missing = [s for s in stems_list if s not in survey.avg_files]
        if missing:
            click.echo(
                f"Warning: stems not found in survey: {', '.join(missing)}",
                err=True,
            )

    from pycsamt.tdem import plot as _tdem_plot  # noqa: PLC0415

    ax_or_axes = None

    try:
        if kind == "decay":
            soundings = survey.to_soundings(
                stems=stems_list, component=component, verbose=verbose
            )
            if not soundings:
                click.echo(
                    "No soundings extracted — check --stems / --component.",
                    err=True,
                )
                sys.exit(1)
            ax_or_axes = _tdem_plot.plot_decay(soundings)

        elif kind == "rho":
            soundings = survey.to_soundings(
                stems=stems_list, component=component, verbose=verbose
            )
            if not soundings:
                click.echo("No soundings extracted.", err=True)
                sys.exit(1)
            ax_or_axes = _tdem_plot.plot_transformed_rho(soundings)

        elif kind == "section":
            avg_data = (
                {s: survey.avg_files[s] for s in stems_list
                 if s in survey.avg_files}
                if stems_list else survey.avg_files
            )
            if not avg_data:
                click.echo("No AVG data available for section plot.", err=True)
                sys.exit(1)
            first_avg = next(iter(avg_data.values()))
            ax_or_axes = _tdem_plot.plot_temavg_section(first_avg)

        elif kind == "z-section":
            z_data = (
                {s: survey.z_files[s] for s in stems_list
                 if s in survey.z_files}
                if stems_list else survey.z_files
            )
            if not z_data:
                click.echo(
                    "No companion .Z files found.  "
                    "Z-section requires Zonge TEMAVG .Z contour files.",
                    err=True,
                )
                sys.exit(1)
            first_z = next(iter(z_data.values()))
            ax_or_axes = _tdem_plot.plot_tem_z_section(first_z)

        elif kind == "map":
            ax_or_axes = _tdem_plot.plot_survey_map(survey)

        elif kind == "overview":
            ax_or_axes = _tdem_plot.plot_survey_overview(survey)

        elif kind == "gate-profile":
            avg_data = (
                {s: survey.avg_files[s] for s in stems_list
                 if s in survey.avg_files}
                if stems_list else survey.avg_files
            )
            if not avg_data:
                click.echo("No AVG data available.", err=True)
                sys.exit(1)
            first_avg = next(iter(avg_data.values()))
            ax_or_axes = _tdem_plot.plot_gate_profile(first_avg)

        elif kind == "elevation":
            soundings = survey.to_soundings(
                stems=stems_list, component=component, verbose=verbose
            )
            ax_or_axes = _tdem_plot.plot_elevation_profile(soundings)

        elif kind == "dashboard":
            avg_data = (
                {s: survey.avg_files[s] for s in stems_list
                 if s in survey.avg_files}
                if stems_list else survey.avg_files
            )
            z_data = (
                {s: survey.z_files[s] for s in stems_list
                 if s in survey.z_files}
                if stems_list else survey.z_files
            )
            soundings = survey.to_soundings(
                stems=stems_list, component=component, verbose=verbose
            )
            first_avg = next(iter(avg_data.values()), None)
            first_z   = next(iter(z_data.values()), None)
            ax_or_axes = _tdem_plot.plot_tem_dashboard(
                first_avg, first_z, soundings
            )

    except Exception as exc:  # noqa: BLE001
        click.echo(f"Plot failed ({kind}): {exc}", err=True)
        sys.exit(1)

    if ax_or_axes is None:
        click.echo("No figure produced.", err=True)
        sys.exit(1)

    # Retrieve the figure from the returned axes object
    if hasattr(ax_or_axes, "figure"):
        fig = ax_or_axes.figure
    elif hasattr(ax_or_axes, "__len__") and hasattr(ax_or_axes[0], "figure"):
        fig = ax_or_axes[0].figure
    else:
        fig = plt.gcf()

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = survey_dir.name or "tdem"
        out_path = output_dir / f"{stem}_{kind}.{fmt}"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        click.echo(f"Saved: {out_path}")
    else:
        plt.show()
