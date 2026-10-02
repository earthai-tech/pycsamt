# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format build-pcgl
==========================

Build a PCGL (pyCSAMT Common Geology Legend, ``.pcgl.json``) document
from a CSV of named resistivity units -- see
:class:`pycsamt.format.geology.GeologyLegend`.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.options import no_color_option, overwrite_option, verbose_option
from ._base import _rich_table, fmt


@fmt.command("build-pcgl")
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="SOURCE.csv",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    metavar="OUT.pcgl.json",
    help="Output path (default: SOURCE with a .pcgl.json extension).",
)
@click.option("--title", default="", help="Legend title.")
@click.option(
    "--document-id", default=None, help="Explicit document id (default: auto)."
)
@click.option(
    "--created-by",
    default="pycsamt format build-pcgl",
    show_default=True,
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@verbose_option
@no_color_option
@overwrite_option
def build_pcgl(
    source: Path,
    output_path: Path | None,
    title: str,
    document_id: str | None,
    created_by: str,
    output_format: str,
    verbose: int,
    no_color: bool,
    overwrite: bool,
) -> None:
    """Build a PCGL geology legend from SOURCE.csv.

    SOURCE columns: ``name, rho_min, rho_max`` required; ``color,
    description, code, source, pattern_id, pattern_source`` optional --
    the same columns :meth:`~pycsamt.geology.lithology.RockDatabase.from_csv`
    accepts.

    \b
    Examples:
      pycsamt format build-pcgl units.csv
      pycsamt format build-pcgl units.csv -o site_a.pcgl.json --title "Site A"
    """
    from pycsamt.format.geology import legend_from_csv, write_legend

    dst = output_path or source.with_suffix("").with_suffix(".pcgl.json")
    if dst.exists() and not overwrite:
        raise click.ClickException(
            f"{dst} exists — pass --overwrite to replace it."
        )

    try:
        legend = legend_from_csv(
            source,
            document_id=document_id,
            created_by=created_by,
            title=title,
        )
    except Exception as exc:  # noqa: BLE001
        raise click.ClickException(f"Failed to build PCGL: {exc}") from exc

    written = write_legend(legend, dst)
    report = {
        "file": str(written),
        "document_id": legend.document_id,
        "title": legend.title,
        "n_entries": len(legend.entries),
    }

    if output_format == "json":
        click.echo(json.dumps(report, indent=2, default=str))
        return
    _rich_table(
        "format build-pcgl",
        [
            ("wrote", report["file"]),
            ("document_id", report["document_id"]),
            ("title", report["title"] or "-"),
            ("entries", str(report["n_entries"])),
        ],
    )
    click.echo(f"\n✓ {written}")
