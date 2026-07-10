# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt site edit — rotate, slice, fill or re-coordinate sites."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    fresh_option,
    no_color_option,
    output_dir_option,
    overwrite_option,
    survey_option,
    verbose_option,
)
from ....api.cli.params import FreqRange
from ._base import _get_sites, site


@site.command("edit")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
# --- edit operations ---
@click.option(
    "--rotate",
    "rotate_deg",
    type=float,
    default=None,
    metavar="DEGREES",
    help="Rotate impedance tensor by this azimuthal angle (degrees).",
)
@click.option(
    "--freq", "freq_range",
    type=FreqRange(),
    default=None,
    metavar="FMIN:FMAX",
    help="Slice impedance to this frequency range (Hz).",
)
@click.option(
    "--fill-missing",
    is_flag=True,
    default=False,
    help="Linearly interpolate NaN values in impedance arrays.",
)
@click.option(
    "--recompute-rho-phase",
    is_flag=True,
    default=False,
    help="Recompute apparent resistivity and phase from Z after edits.",
)
@click.option(
    "--set-coords",
    default=None,
    metavar="LAT,LON[,ELEV]",
    help=(
        "Override coordinates for every site with a single LAT,LON,ELEV value.  "
        "Use --coords-table for per-site overrides."
    ),
)
@click.option(
    "--coords-table",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    metavar="FILE",
    help=(
        "CSV/Excel file with columns station,lat,lon[,elev] "
        "to set per-site coordinates."
    ),
)
# --- output ---
@output_dir_option
@overwrite_option
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show what would be done without writing files.",
)
@verbose_option
@no_color_option
@click.pass_context
def edit(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    rotate_deg: float | None,
    freq_range: tuple[float, float] | None,
    fill_missing: bool,
    recompute_rho_phase: bool,
    set_coords: str | None,
    coords_table: Path | None,
    output_dir: Path,
    overwrite: bool,
    dry_run: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Apply edits (rotate, slice, fill, re-coordinate) and write results.

    At least one edit operation must be specified.
    --output-dir is required unless --dry-run is used.

    \b
    Examples:
      # rotate all sites by 30 degrees:
      pycsamt site edit --rotate 30 --output-dir rotated/

      # slice to 0.1–1000 Hz and fill missing values:
      pycsamt site edit --freq 0.1:1000 --fill-missing --output-dir edited/

      # set uniform coordinates:
      pycsamt site edit --set-coords 28.5,102.3,200 --output-dir recoord/

      # apply coords from CSV table:
      pycsamt site edit --coords-table coords.csv --output-dir recoord/

      # dry run — show planned operations:
      pycsamt site edit --rotate 15 --freq 1:100 --dry-run
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    # Require at least one edit operation
    ops = [rotate_deg, freq_range, fill_missing,
           recompute_rho_phase, set_coords, coords_table]
    if not any(op is not None and op is not False for op in ops):
        raise click.UsageError(
            "Specify at least one edit operation: "
            "--rotate, --freq, --fill-missing, "
            "--recompute-rho-phase, --set-coords, or --coords-table"
        )

    # Require output-dir unless dry-run
    if not dry_run and str(output_dir) == ".":
        raise click.UsageError(
            "--output-dir is required when not using --dry-run"
        )

    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    from pycsamt.site import edit as edit_mod  # noqa: PLC0415

    ops_applied: list[str] = []
    result = sites_obj

    if rotate_deg is not None:
        ops_applied.append(f"rotate({rotate_deg}°)")
        if not dry_run:
            result = edit_mod.rotate_all(result, rotate_deg)

    if freq_range is not None:
        ops_applied.append(f"select_freq({freq_range[0]}–{freq_range[1]} Hz)")
        if not dry_run:
            result = edit_mod.select_freq_all(result, freq_range[0], freq_range[1])

    if fill_missing:
        ops_applied.append("fill_missing")
        if not dry_run:
            result = edit_mod.fill_missing(result)

    if recompute_rho_phase:
        ops_applied.append("recompute_res_phase")
        if not dry_run:
            result = edit_mod.recompute_res_phase(result)

    if set_coords is not None:
        try:
            parts = [float(x.strip()) for x in set_coords.split(",")]
            lat, lon = parts[0], parts[1]
            elev = parts[2] if len(parts) > 2 else None
        except (ValueError, TypeError):
            raise click.BadParameter(
                "Expected LAT,LON or LAT,LON,ELEV",
                param_hint="--set-coords",
            )
        ops_applied.append(f"set_coords({lat}, {lon})")
        if not dry_run:
            result = edit_mod.set_coords_all(result, lat, lon, elev)

    if coords_table is not None:
        ops_applied.append(f"set_coords_from_table({coords_table.name})")
        if not dry_run:
            result = edit_mod.set_coords_from_table(result, coords_table)

    if dry_run:
        click.echo(
            f"Dry run — {len(sites_obj)} station(s), "
            f"operations: {', '.join(ops_applied)}"
        )
        return

    # --- write modified sites ---
    output_dir.mkdir(parents=True, exist_ok=True)
    from pycsamt.site.export import (
        write_sites,  # noqa: PLC0415
    )
    try:
        written = write_sites(result, output_dir, exist_ok=overwrite)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error writing sites: {exc}", err=True)
        sys.exit(1)

    click.echo(
        f"Edited {len(sites_obj)} station(s) "
        f"[{', '.join(ops_applied)}] → {output_dir}/ "
        f"({len(written)} file(s) written)"
    )
