# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt edi rotate — rotate Z and Tipper tensors and write output EDIs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    output_dir_option,
    verbose_option,
)
from ._base import _get_collection, _get_edi, edi


@edi.command("rotate")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_OR_DIR",
)
@click.option(
    "--angle",
    "-a",
    "angle_deg",
    type=float,
    required=True,
    help="Rotation angle in degrees (clockwise positive, geophysical convention).",
)
@output_dir_option
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show which files would be rotated without writing anything.",
)
@verbose_option
@no_color_option
@click.pass_context
def rotate(
    ctx: click.Context,
    source: Path,
    angle_deg: float,
    output_dir: Path,
    output_format: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Rotate impedance tensors (and Tipper) by a constant angle.

    Applies  Z' = R·Z·Rᵀ  where R is a 2×2 rotation matrix for the given
    clockwise angle.  Tipper vectors are rotated as  t' = R·t.
    The rotated EDI files are written to ``--output-dir``.

    \b
    Examples:
      pycsamt edi rotate data/3edis/ --angle 30 --output-dir rotated/
      pycsamt edi rotate data/3edis/HBH03.edi --angle -15 --output-dir rotated/
      pycsamt edi rotate data/AMT/ --angle 45 --dry-run
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.seg.ops import (  # noqa: PLC0415
        rotate_impedance,
        rotate_tipper,
    )

    # Collect EDI files
    if source.is_file():
        edis = [_get_edi(source, verbose=verbose)]
        paths = [source]
    else:
        coll = _get_collection(source, verbose=verbose)
        edis = list(coll)
        paths = [
            Path(e.path_str) if e.path_str else source / f"{e.station}.edi"
            for e in edis
        ]

    if not edis:
        click.echo(f"No EDI files found in {source}.", err=True)
        sys.exit(1)

    if dry_run:
        click.echo(f"Dry run — {len(edis)} file(s) would be rotated by {angle_deg}°:")
        for p in paths:
            click.echo(f"  {p.name}")
        return

    if str(output_dir) == ".":
        raise click.UsageError("--output-dir is required when not using --dry-run")

    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    failures: list[dict] = []

    for edi, src_path in zip(edis, paths):
        try:
            z = getattr(edi.Z, "z", None)
            if z is not None and z.size > 0:
                edi.Z.z = rotate_impedance(z, angle_deg)

            if edi.has_tipper:
                t = getattr(edi.Tip, "tipper", None)
                if t is not None and t.size > 0:
                    edi.Tip.tipper = rotate_tipper(t, angle_deg)

            out_name = src_path.name
            out_path = edi.write(
                savepath=str(output_dir),
                new_edifn=out_name,
            )
            written.append(str(out_path))
            if verbose >= 1:
                click.echo(f"  ✓  {out_name}")

        except Exception as exc:  # noqa: BLE001
            failures.append({"source": str(src_path), "error": str(exc)})
            if verbose >= 1:
                click.echo(f"  ✗  {src_path.name}: {exc}", err=True)

    summary = {
        "angle_deg": angle_deg,
        "output_dir": str(output_dir),
        "n_input": len(edis),
        "n_written": len(written),
        "n_failed": len(failures),
        "written": written,
        "failures": failures,
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2))
    else:
        click.echo(
            f"Rotated  {len(written)}/{len(edis)}  by {angle_deg}°"
            f"  | Errors: {len(failures)}"
            f"  | Output: {output_dir}/"
        )
        for f in failures:
            click.echo(f"  ✗  {f['source']}: {f['error']}", err=True)

    if failures and not written:
        sys.exit(1)
