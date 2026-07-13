# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt site recompute - normalize and rewrite EDI files."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    overwrite_option,
    survey_option,
    verbose_option,
)
from ....api.cli.params import FreqRange
from ._base import site


@site.command("recompute")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    default=None,
)
@survey_option
@click.option(
    "-o",
    "--output-dir",
    "output_root",
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    default=None,
    metavar="DIR",
    help=(
        "Root directory for recomputed EDI files. If omitted, "
        "a recomputed_edis folder is created next to the input line "
        "or inside the input survey directory."
    ),
)
@click.option(
    "--output-name",
    default="recomputed_edis",
    show_default=True,
    metavar="NAME",
    help="Default output folder name when --output-dir is omitted.",
)
@click.option(
    "--flatten",
    is_flag=True,
    default=False,
    help="Write all recomputed EDIs directly under the output root.",
)
@click.option(
    "--template",
    default="{source_stem}.edi",
    show_default=True,
    metavar="TEMPLATE",
    help=(
        "Output filename template. Placeholders: {station}, "
        "{index}, {line}, {source_stem}."
    ),
)
@click.option(
    "--rotate",
    "rotate_angle",
    type=float,
    default=None,
    metavar="DEGREES",
    help="Rotate transfer functions before recomputing.",
)
@click.option(
    "--components",
    default="Z,Tip",
    show_default=True,
    metavar="LIST",
    help="Comma-separated components to rotate: Z, Tip, or Z,Tip.",
)
@click.option(
    "--freq",
    "freq_range",
    type=FreqRange(),
    default=None,
    metavar="FMIN:FMAX",
    help="Keep only this frequency range before recomputing.",
)
@click.option(
    "--fill-missing",
    type=click.Choice(["zero", "nan"], case_sensitive=False),
    default=None,
    metavar="HOW",
    help="Fill missing Z/tipper values before recomputing.",
)
@click.option(
    "--skip-rho-phase",
    is_flag=True,
    default=False,
    help="Do not recompute apparent resistivity and phase from Z.",
)
@click.option(
    "--datatype",
    type=click.Choice(["mt", "emap"], case_sensitive=False),
    default=None,
    help="Force EDI writer datatype instead of auto-detection.",
)
@click.option(
    "--synthesize-spectra",
    is_flag=True,
    default=False,
    help="Ask the EDI writer to synthesize spectra when possible.",
)
@click.option(
    "--no-manifest",
    is_flag=True,
    default=False,
    help="Do not write manifest.csv.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Recompute in memory and report what would be written.",
)
@click.option(
    "--progress",
    is_flag=True,
    default=False,
    help="Show a progress bar while recomputing.",
)
@overwrite_option
@verbose_option
@no_color_option
@click.pass_context
def recompute(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    output_root: Path | None,
    output_name: str,
    flatten: bool,
    template: str,
    rotate_angle: float | None,
    components: str,
    freq_range: tuple[float, float] | None,
    fill_missing: str | None,
    skip_rho_phase: bool,
    datatype: str | None,
    synthesize_spectra: bool,
    no_manifest: bool,
    dry_run: bool,
    progress: bool,
    overwrite: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Recompute foreign or raw EDI files and rewrite pyCSAMT EDIs.

    The command accepts a single EDI file, a line folder, a survey folder
    containing line subfolders, or the active survey. By default output
    is written to a ``recomputed_edis`` directory and line folders are
    preserved.

    \b
    Examples:
      pycsamt site recompute WILLY_DATA/L18PLT --rotate 30
      pycsamt site recompute WILLY_DATA --flatten --template '{line}_{source_stem}.edi'
      pycsamt site recompute --survey WILLY_DATA --output-dir cleaned_edis
    """

    configure_cli(log__level=verbose, log__color=not no_color)
    source = _resolve_source(edi_source, survey_path)
    comps = _parse_components(components)
    fmin = fmax = None
    if freq_range is not None:
        fmin, fmax = freq_range

    from pycsamt.site.recompute import (
        EDIRecomputer,  # noqa: PLC0415
    )

    try:
        result = EDIRecomputer(
            output_root=output_root,
            output_name=output_name,
            preserve_line_dirs=not flatten,
            template=template,
            overwrite=overwrite,
            write=not dry_run,
            manifest_csv=not no_manifest,
            rotate_angle=rotate_angle,
            rotate_components=comps,
            fmin=fmin,
            fmax=fmax,
            fill_missing_values=fill_missing,
            recompute_resphase=not skip_rho_phase,
            datatype=datatype,
            synthesize_spectra=synthesize_spectra,
            progress=progress,
            verbose=verbose,
        ).run(source)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error recomputing EDI files: {exc}", err=True)
        sys.exit(1)

    n_total = len(result.records)
    n_ok = n_total - len(result.failed)

    if dry_run:
        click.echo(
            f"Dry run - recomputed {n_ok}/{n_total} EDI file(s) "
            "in memory; no files written."
        )
    else:
        click.echo(
            f"Recomputed {n_ok}/{n_total} EDI file(s) -> {result.output_root}"
        )
        if not no_manifest and result.output_root is not None:
            click.echo(f"Manifest: {result.output_root / 'manifest.csv'}")

    if verbose >= 1:
        for rec in result.records:
            label = rec.station or (
                rec.source.name if rec.source is not None else "unknown"
            )
            if rec.status == "ok":
                target = rec.output if rec.output is not None else "<memory>"
                click.echo(f"  ok      {label} -> {target}", err=True)
            else:
                click.echo(
                    f"  failed  {label}: {rec.message}",
                    err=True,
                )

    if result.failed:
        sys.exit(1)


def _resolve_source(
    edi_source: Path | None,
    survey_path: Path | None,
) -> Path:
    if edi_source is not None:
        return edi_source
    if survey_path is not None:
        return survey_path

    from pycsamt.cli.survey import (
        survey_summary,  # noqa: PLC0415
    )

    summary = survey_summary()
    if summary and summary.get("survey_path"):
        return Path(str(summary["survey_path"]))
    raise click.UsageError(
        "No EDI source provided. Pass EDI_SOURCE, use --survey, "
        "or set an active survey with 'pycsamt survey set'."
    )


def _parse_components(value: str) -> tuple[str, ...]:
    comps = tuple(
        part.strip() for part in str(value).split(",") if part.strip()
    )
    allowed = {"z", "imp", "impedance", "t", "tip", "tipper"}
    bad = [c for c in comps if c.lower() not in allowed]
    if bad:
        raise click.BadParameter(
            "Unknown component(s): "
            + ", ".join(bad)
            + ". Expected Z, Tip, or Z,Tip.",
            param_hint="--components",
        )
    return comps or ("Z", "Tip")
