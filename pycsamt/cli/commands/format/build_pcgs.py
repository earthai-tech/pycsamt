# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format build-pcgs
==========================

Build a PCGS (structural evidence, ``.pcgs.json``) document from one or
more CSVs of planar/linear measurements and fault traces -- see
:class:`pycsamt.format.structure.StructModel`.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.options import no_color_option, overwrite_option, verbose_option
from ._base import _rich_table, fmt


@fmt.command("build-pcgs")
@click.option(
    "--planar",
    "planar_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="CSV of planar measurements (strike/dip).",
)
@click.option(
    "--linear",
    "linear_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="CSV of linear measurements (trend/plunge).",
)
@click.option(
    "--faults",
    "faults_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="CSV of fault traces.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    required=True,
    metavar="OUT.pcgs.json",
)
@click.option("--title", default="", help="Document title.")
@click.option(
    "--document-id", default=None, help="Explicit document id (default: auto)."
)
@click.option(
    "--created-by",
    default="pycsamt format build-pcgs",
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
def build_pcgs(
    planar_path: Path | None,
    linear_path: Path | None,
    faults_path: Path | None,
    output_path: Path,
    title: str,
    document_id: str | None,
    created_by: str,
    output_format: str,
    verbose: int,
    no_color: bool,
    overwrite: bool,
) -> None:
    """Build a PCGS structural-evidence document.

    At least one of --planar / --linear / --faults is required.

    \b
    Examples:
      pycsamt format build-pcgs --planar strikes.csv -o site.pcgs.json
      pycsamt format build-pcgs --planar p.csv --linear l.csv --faults f.csv \\
          -o site.pcgs.json
    """
    if not (planar_path or linear_path or faults_path):
        raise click.UsageError(
            "Pass at least one of --planar / --linear / --faults."
        )
    if output_path.exists() and not overwrite:
        raise click.ClickException(
            f"{output_path} exists — pass --overwrite to replace it."
        )

    from pycsamt.format.structure import structure_from_csv, write_structure

    try:
        model = structure_from_csv(
            planar_path=planar_path,
            linear_path=linear_path,
            faults_path=faults_path,
            document_id=document_id,
            created_by=created_by,
            title=title,
        )
    except Exception as exc:  # noqa: BLE001
        raise click.ClickException(f"Failed to build PCGS: {exc}") from exc

    written = write_structure(model, output_path)
    report = {
        "file": str(written),
        "document_id": model.document_id,
        "title": model.title,
        "n_planar": len(model.model.planar),
        "n_linear": len(model.model.linear),
        "n_faults": len(model.model.faults),
    }

    if output_format == "json":
        click.echo(json.dumps(report, indent=2, default=str))
        return
    _rich_table(
        "format build-pcgs",
        [
            ("wrote", report["file"]),
            ("document_id", report["document_id"]),
            ("planar", str(report["n_planar"])),
            ("linear", str(report["n_linear"])),
            ("faults", str(report["n_faults"])),
        ],
    )
    click.echo(f"\n✓ {written}")
