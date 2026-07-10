# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt transform spectra — convert Spectra EDI to Impedance EDI."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    output_dir_option,
    overwrite_option,
    verbose_option,
)
from ._base import transform


@transform.command("spectra")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="FILE_OR_DIR",
)
@output_dir_option
@overwrite_option
# --- conversion options ---
@click.option(
    "--station-suffix",
    default="",
    show_default=True,
    metavar="STR",
    help=(
        "Suffix appended to every station name, e.g. '_IMP' "
        "produces HBH03 → HBH03_IMP."
    ),
)
@click.option(
    "--station-name",
    default=None,
    metavar="NAME",
    help="Override the station name (single-file input only).",
)
@click.option(
    "--estimate-error",
    is_flag=True,
    default=False,
    help=(
        "Propagate 1-σ impedance errors into >ZXX.VAR / >ZXY.VAR blocks "
        "using first-order complex-Wishart error propagation."
    ),
)
@click.option(
    "--ridge",
    type=float,
    default=None,
    metavar="FLOAT",
    help=(
        "Tikhonov regularisation for the magnetic block inversion.  "
        "Useful for numerically ill-conditioned spectra."
    ),
)
@click.option(
    "--dof",
    type=float,
    default=None,
    metavar="FLOAT",
    help=(
        "Override effective degrees of freedom (per frequency) "
        "used in error estimation.  Auto-inferred when not provided."
    ),
)
@click.option(
    "--use-remote",
    is_flag=True,
    default=False,
    help=(
        "When both local and remote electric channels are present, "
        "use the remote reference for impedance estimation."
    ),
)
@click.option(
    "--e-labels",
    default="EX,EY",
    show_default=True,
    metavar="EX,EY",
    help="Comma-separated electric channel type labels.",
)
@click.option(
    "--h-labels",
    default="HX,HY",
    show_default=True,
    metavar="HX,HY",
    help="Comma-separated horizontal magnetic channel type labels.",
)
# --- output control ---
@click.option(
    "--format", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
    help="Output format for the conversion summary.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help=(
        "Validate inputs and print what would be converted "
        "without writing any files."
    ),
)
@verbose_option
@no_color_option
@click.pass_context
def spectra(
    ctx: click.Context,
    source: Path,
    output_dir: Path,
    overwrite: bool,
    station_suffix: str,
    station_name: str | None,
    estimate_error: bool,
    ridge: float | None,
    dof: float | None,
    use_remote: bool,
    e_labels: str,
    h_labels: str,
    output_format: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Convert Spectra EDI file(s) to Impedance EDI with optional Tipper.

    SOURCE may be a single .edi file or a directory.  All .edi files
    in a directory are processed in sorted order.

    The conversion recovers the impedance tensor Z (and tipper T when
    a vertical magnetic channel HZ is present) from the cross-spectral
    matrices embedded in the Spectra-format EDI.

    \b
    Examples:
      # single file → output directory
      pycsamt transform spectra HBH03.edi --output-dir imp_edis/

      # append _IMP suffix to station names
      pycsamt transform spectra spectra_dir/ --station-suffix _IMP \\
              --output-dir imp_edis/

      # propagate 1-σ errors
      pycsamt transform spectra spectra_dir/ --estimate-error \\
              --output-dir imp_edis/

      # dry run — validate inputs, show what would be converted
      pycsamt transform spectra spectra_dir/ --dry-run

      # machine-readable summary
      pycsamt transform spectra spectra_dir/ --output-dir imp/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    e_list = [e.strip() for e in e_labels.split(",")]
    h_list = [h.strip() for h in h_labels.split(",")]

    # Collect candidate files
    if source.is_dir():
        files = sorted(source.glob("*.edi"))
    else:
        files = [source]

    if not files:
        click.echo(f"No .edi files found in {source}.", err=True)
        sys.exit(1)

    if dry_run:
        click.echo(
            f"Dry run — {len(files)} file(s) would be converted:\n"
        )
        for f in files:
            name = (station_name if station_name else f.stem) + station_suffix
            click.echo(f"  {f.name}  →  {name}.edi")
        return

    if str(output_dir) == ".":
        raise click.UsageError("--output-dir is required when not using --dry-run")

    from pycsamt.transformers import (
        SpectraToEDI,  # noqa: PLC0415
    )

    transformer = SpectraToEDI(
        e_labels=tuple(e_list),
        h_labels=tuple(h_list),
        ridge=ridge,
        estimate_error=estimate_error,
        dof=dof,
        use_remote=use_remote,
        station_suffix=station_suffix,
        skip_errors=True,
        verbose=verbose,
    )

    if verbose >= 1:
        click.echo(
            f"Converting {len(files)} spectra EDI(s) → {output_dir}/",
            err=True,
        )

    result = transformer.transform_batch(
        source,
        output_dir=output_dir,
        station_name=station_name,
    )

    # --- report ---
    summary = {
        "n_input":   len(files),
        "n_ok":      result.n_ok,
        "n_fail":    result.n_fail,
        "output_dir": str(output_dir),
        "converted": [
            {"station": ed.station, "n_freq": ed.Z.n_freq,
             "has_tipper": ed.has_tipper}
            for ed in result.collection
        ],
        "failures": [
            {"source": r.source, "error": r.error}
            for r in result.failures
        ],
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2))
    else:
        click.echo(
            f"Converted: {result.n_ok}/{len(files)}  "
            f"| Errors: {result.n_fail}  "
            f"| Output: {output_dir}/"
        )
        if verbose >= 1:
            for ed in result.collection:
                tip = "tipper ✓" if ed.has_tipper else "no tipper"
                click.echo(f"  ✓  {ed.station:30}  {ed.Z.n_freq:3d} freq  {tip}")
        if result.failures:
            click.echo("\nFailed conversions:", err=True)
            for r in result.failures:
                click.echo(f"  ✗  {r.source}: {r.error}", err=True)

    if result.n_fail and result.n_ok == 0:
        sys.exit(1)
