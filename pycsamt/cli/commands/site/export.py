# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt site export — write Sites to EDI files or a zip archive."""

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

from ._base import site, _get_sites


@site.command("export")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@output_dir_option
@overwrite_option
@click.option(
    "--zip", "pack_zip",
    is_flag=True,
    default=False,
    help="Pack all EDI outputs into a single zip archive.",
)
@click.option(
    "--zip-name",
    default=None,
    metavar="FILENAME",
    help="Name for the zip archive (default: survey.zip).",
)
@click.option(
    "--template",
    default="{name}.edi",
    show_default=True,
    metavar="TEMPLATE",
    help=(
        "Output filename template.  "
        "Available placeholders: {name}, {index}, {lat}, {lon}."
    ),
)
@click.option(
    "--manifest",
    is_flag=True,
    default=False,
    help="Write a manifest.csv alongside the exported files.",
)
@verbose_option
@no_color_option
@click.pass_context
def export(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    output_dir: Path,
    overwrite: bool,
    pack_zip: bool,
    zip_name: str | None,
    template: str,
    manifest: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Write Sites to EDI files (or a zip archive).

    Without --zip, one EDI file per station is written to --output-dir.
    With --zip, all EDIs are packed into a single archive.

    \b
    Examples:
      # export to directory:
      pycsamt site export --output-dir final_edis/

      # custom filename template:
      pycsamt site export --template '{index:03d}_{name}.edi' \\
              --output-dir final_edis/

      # pack as zip:
      pycsamt site export --zip --zip-name survey_export.zip \\
              --output-dir .

      # include manifest CSV:
      pycsamt site export --manifest --output-dir final_edis/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    if str(output_dir) == ".":
        raise click.UsageError("--output-dir is required for export")

    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    output_dir.mkdir(parents=True, exist_ok=True)

    if pack_zip:
        archive_name = zip_name or "survey.zip"
        archive_path = output_dir / archive_name
        if archive_path.exists() and not overwrite:
            click.echo(
                f"Archive {archive_path} already exists.  "
                "Use --overwrite to replace it.",
                err=True,
            )
            sys.exit(1)

        from pycsamt.site.export import pack_zip as _pack_zip  # noqa: PLC0415
        try:
            out = _pack_zip(sites_obj, archive_path)
        except Exception as exc:  # noqa: BLE001
            click.echo(f"Error creating archive: {exc}", err=True)
            sys.exit(1)
        click.echo(f"Packed {len(sites_obj)} station(s) → {out}")
        return

    # --- write individual EDI files ---
    from pycsamt.site.export import write_sites  # noqa: PLC0415
    try:
        written = write_sites(
            sites_obj,
            output_dir,
            template=template,
            exist_ok=overwrite,
            manifest_csv=str(output_dir / "manifest.csv") if manifest else None,
        )
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error exporting sites: {exc}", err=True)
        sys.exit(1)

    click.echo(
        f"Exported {len(written)} / {len(sites_obj)} station(s) → {output_dir}/"
    )
    if verbose >= 1:
        for p in written:
            click.echo(f"  {p.name}", err=True)
