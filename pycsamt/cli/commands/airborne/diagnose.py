# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt airborne diagnose — technology-aware diagnostic table.

Auto-detects the airborne technology from the loaded sites'
:attr:`~pycsamt.airborne.site.AirborneSite.technology` and dispatches
to that technology's canonical ``emtools`` table function:

=============== =======================================================
ztem            :func:`~pycsamt.emtools.ztem.total_divergence_table`
                (along-line total-divergence / Peaker table -- one
                flight line at a time; use ``--line`` on a multi-line
                survey)
mobilemt        :func:`~pycsamt.emtools.mobilemt.admittance_table`
afmag_original  :func:`~pycsamt.emtools.afmag.original_afmag_tilt_table`
afmag_airmt     :func:`~pycsamt.emtools.afmag.airmt_tilt_angles`
=============== =======================================================

Usage
-----
::

    pycsamt airborne diagnose survey_xml/
    pycsamt airborne diagnose ztem_survey/ --line L100
    pycsamt airborne diagnose mixed_survey/ --technology mobilemt
    pycsamt airborne diagnose survey_xml/ --format csv --output table.csv
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
from ._base import _get_asites, airborne

_TECHNOLOGIES = ("ztem", "mobilemt", "afmag_original", "afmag_airmt")


def _select_ztem_line(sites, line_id: str | None):
    """Return the single flight line ``total_divergence_table`` needs.

    Generic EMTF-XML ingestion (:func:`~pycsamt.airborne.site
    .ensure_asites`) never populates
    :attr:`~pycsamt.airborne.site.AirborneSite.line_id` -- that field
    is only set when sites are built explicitly (e.g. via
    ``AirborneSites.from_line``). So a CLI ``--line`` has to work even
    when *every* site's ``line_id`` is ``None``: it is tried first as
    real ``line_id`` metadata, then as a 0-based index into flight
    lines detected geometrically from station coordinates -- the same
    heuristic :func:`~pycsamt.emtools.ztem.plot_ztem_map` uses for its
    per-line ``quantity="divergence"`` grouping.
    """
    import numpy as np

    from ....emtools.ztem import (  # noqa: PLC0415
        _detect_line_groups,
        _station_lonlat,
    )

    items, lon, lat = _station_lonlat(sites)
    groups = (
        _detect_line_groups(lon, lat) if items else np.zeros(0, dtype=int)
    )
    n_groups = (
        int(groups.max()) + 1 if groups.size else (1 if len(sites) else 0)
    )
    name_to_group = {
        getattr(it, "name", None): int(g) for it, g in zip(items, groups)
    }

    if line_id is not None:
        by_id = sites.select(predicate=lambda s: s.line_id == line_id)
        if len(by_id) > 0:
            return by_id
        try:
            wanted = int(line_id)
        except ValueError:
            wanted = None
        if wanted is not None:
            picked = sites.select(
                predicate=lambda s: name_to_group.get(s.name) == wanted
            )
            if len(picked) > 0:
                return picked
        raise click.UsageError(
            f"No sites found for --line {line_id!r} (checked both real "
            f"line_id metadata and the {n_groups} geometrically "
            f"detected flight line(s), indexed 0..{max(n_groups - 1, 0)})."
        )

    if n_groups <= 1:
        return sites
    raise click.UsageError(
        f"Detected {n_groups} flight lines geometrically from station "
        "coordinates -- differentiating the total-divergence table "
        "across a line boundary is not physically meaningful.  Pass "
        f"--line with a 0-based group index (0..{n_groups - 1}) to "
        "select one, e.g. --line 0."
    )


@airborne.command("diagnose")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="SOURCE",
)
@click.option(
    "--technology",
    "-t",
    type=click.Choice(_TECHNOLOGIES, case_sensitive=False),
    default=None,
    help="Force the technology instead of auto-detecting it.",
)
@click.option(
    "--line",
    "-l",
    "line_id",
    default=None,
    metavar="LINE_ID",
    help="ztem only: restrict to one flight line -- a real line_id, "
    "or a 0-based geometrically-detected group index (required when "
    "the dataset has more than one line).",
)
@click.option(
    "--component",
    type=click.Choice(["tzx", "tzy"], case_sensitive=False),
    default="tzx",
    show_default=True,
    help="ztem only: tipper component to differentiate.",
)
@click.option(
    "--spacing-m",
    type=float,
    default=200.0,
    show_default=True,
    metavar="METRES",
    help="ztem only: fall-back inter-station spacing when station "
    "coordinates are unavailable.",
)
@click.option(
    "--no-recursive",
    is_flag=True,
    default=False,
    help="When SOURCE is a directory, do not search it recursively.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(writable=True, path_type=Path),
    default=None,
    metavar="FILE",
    help="Save the table to this CSV file.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def diagnose(
    ctx: click.Context,
    source: Path,
    technology: str | None,
    line_id: str | None,
    component: str,
    spacing_m: float,
    no_recursive: bool,
    output: Path | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Build the canonical diagnostic table for an airborne dataset.

    SOURCE is a single EMTF-XML file or a directory of them.  The
    technology (ZTEM, MobileMT, or AFMAG/AirMt) is auto-detected from
    the data unless ``--technology`` is given.

    \b
    Examples:
      pycsamt airborne diagnose survey_xml/
      pycsamt airborne diagnose ztem_survey/ --line L100
      pycsamt airborne diagnose mixed_survey/ --technology mobilemt
      pycsamt airborne diagnose survey_xml/ --format csv -o table.csv
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

    # --- resolve technology ---
    effective_tech = technology.lower() if technology else None
    if effective_tech is None:
        detected = sites.technologies
        if len(detected) == 1:
            effective_tech = detected[0]
        elif len(detected) == 0:
            click.echo(
                "Could not detect a technology for any site under "
                f"{source}.  Pass --technology explicitly.",
                err=True,
            )
            sys.exit(1)
        else:
            raise click.UsageError(
                f"Multiple technologies present ({', '.join(detected)}).  "
                "Pass --technology to select one."
            )
    if verbose >= 1:
        click.echo(f"Technology: {effective_tech}", err=True)

    # --- build the table ---
    try:
        if effective_tech == "ztem":
            from ....emtools.ztem import (  # noqa: PLC0415
                total_divergence_table,
            )

            line_sites = _select_ztem_line(sites, line_id)
            df = total_divergence_table(
                line_sites, spacing_m, component=component, verbose=verbose
            )

        elif effective_tech == "mobilemt":
            from ....emtools.mobilemt import (  # noqa: PLC0415
                admittance_table,
            )

            df = admittance_table(sites)

        elif effective_tech == "afmag_original":
            from ....emtools.afmag import (  # noqa: PLC0415
                original_afmag_tilt_table,
            )

            df = original_afmag_tilt_table(sites, verbose=verbose)

        else:  # afmag_airmt
            from ....emtools.afmag import (  # noqa: PLC0415
                airmt_tilt_angles,
            )

            df = airmt_tilt_angles(sites, verbose=verbose)

    except click.ClickException:
        raise
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error building diagnostic table: {exc}", err=True)
        sys.exit(1)

    if df.empty:
        click.echo("No diagnostic rows produced from this dataset.", err=True)

    if output is not None:
        df.to_csv(output, index=False)
        click.echo(f"Table saved → {output}")

    if output_format == "json":
        click.echo(df.to_json(orient="records", indent=2, default_handler=str))
    elif output_format == "csv":
        click.echo(df.to_csv(index=False))
    else:
        click.echo(df.to_string(index=False))
