# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format edi-to-xml / xml-to-edi
=======================================

Losslessly convert between historical SEG EDI and EMTF-XML transfer
functions, one file or a whole directory at a time -- see
:mod:`pycsamt.emtf.converters.edi` and :mod:`pycsamt.emtf.xml`.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.options import (
    no_color_option,
    output_dir_option,
    overwrite_option,
    verbose_option,
)
from ._base import _rich_table, fmt

_ON_LOSS_CHOICES = ["warn", "raise", "ignore"]


def _iter_sources(source: Path, patterns: tuple[str, ...]) -> list[Path]:
    if source.is_file():
        return [source]
    # On case-insensitive filesystems (Windows, default macOS) a lower-
    # and upper-case pattern both match the same files -- dedupe by
    # resolved path rather than relying on the patterns being disjoint.
    found: dict[Path, None] = {}
    for pattern in patterns:
        for path in source.glob(pattern):
            found.setdefault(path.resolve(), None)
    return sorted(found)


@fmt.command("edi-to-xml")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="SOURCE",
)
@click.option(
    "--prefer-spectra/--no-prefer-spectra",
    default=True,
    show_default=True,
    help="Prefer EDI SPECTRA blocks over impedance/tipper when present.",
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
@output_dir_option
@overwrite_option
def edi_to_xml(
    source: Path,
    prefer_spectra: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
    output_dir: Path,
    overwrite: bool,
) -> None:
    """Convert SOURCE (.edi file or directory of .edi files) to EMTF-XML.

    \b
    Examples:
      pycsamt format edi-to-xml station.edi
      pycsamt format edi-to-xml survey_edi/ -o survey_xml/
    """
    from pycsamt.emtf import edi_to_emtf, write_emtf_xml

    sources = _iter_sources(source, ("*.edi", "*.EDI"))
    if not sources:
        raise click.ClickException(f"No .edi files found under {source}.")

    written: list[dict] = []
    for path in sources:
        try:
            document = edi_to_emtf(path, prefer_spectra=prefer_spectra)
            dst = output_dir / f"{path.stem}.xml"
            if dst.exists() and not overwrite:
                raise click.ClickException(
                    f"{dst} exists — pass --overwrite to replace it."
                )
            write_emtf_xml(document, dst)
        except click.ClickException:
            raise
        except Exception as exc:  # noqa: BLE001
            raise click.ClickException(
                f"{path.name} → XML failed: {exc}"
            ) from exc
        written.append({"source": str(path), "target": str(dst)})
        if verbose:
            click.echo(f"  {path.name} → {dst}")

    if output_format == "json":
        click.echo(json.dumps(written, indent=2, default=str))
        return
    _rich_table(
        "format edi-to-xml",
        [(w["source"], w["target"]) for w in written],
    )
    click.echo(f"\n✓ {len(written)} file(s) written to {output_dir}")


@fmt.command("xml-to-edi")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="SOURCE",
)
@click.option(
    "--on-loss",
    type=click.Choice(_ON_LOSS_CHOICES, case_sensitive=False),
    default="warn",
    show_default=True,
    help="What to do when EMTF-XML content has no lossless EDI equivalent.",
)
@click.option(
    "--strict/--permissive",
    default=True,
    show_default=True,
    help="Reject vs. skip malformed XML elements while reading.",
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
@output_dir_option
@overwrite_option
def xml_to_edi(
    source: Path,
    on_loss: str,
    strict: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
    output_dir: Path,
    overwrite: bool,
) -> None:
    """Convert SOURCE (.xml file or directory of EMTF-XML files) to EDI.

    \b
    Examples:
      pycsamt format xml-to-edi station.xml
      pycsamt format xml-to-edi survey_xml/ -o survey_edi/ --on-loss raise
    """
    from pycsamt.emtf import EMTFXMLReader, write_edi

    sources = _iter_sources(source, ("*.xml", "*.XML"))
    if not sources:
        raise click.ClickException(f"No .xml files found under {source}.")

    reader = EMTFXMLReader(strict=strict)
    written: list[dict] = []
    for path in sources:
        try:
            document = reader.read(path)
            name = document.station or path.stem
            dst = output_dir / f"{name}.edi"
            if dst.exists() and not overwrite:
                raise click.ClickException(
                    f"{dst} exists — pass --overwrite to replace it."
                )
            write_edi(document, dst, on_loss=on_loss.lower())
        except click.ClickException:
            raise
        except Exception as exc:  # noqa: BLE001
            raise click.ClickException(
                f"{path.name} → EDI failed: {exc}"
            ) from exc
        written.append({"source": str(path), "target": str(dst)})
        if verbose:
            click.echo(f"  {path.name} → {dst}")

    if output_format == "json":
        click.echo(json.dumps(written, indent=2, default=str))
        return
    _rich_table(
        "format xml-to-edi",
        [(w["source"], w["target"]) for w in written],
    )
    click.echo(f"\n✓ {len(written)} file(s) written to {output_dir}")
