# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt jones select — filter a J-file collection and export a subset."""

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
from ._base import _get_collection, jones


@jones.command("select")
@click.argument(
    "source",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="J_DIR",
)
@output_dir_option
@click.option(
    "--stations",
    "-s",
    "station_ids",
    default=None,
    metavar="ID[,ID…]",
    help="Comma-separated station IDs to keep.",
)
@click.option(
    "--pattern",
    "-p",
    default=None,
    metavar="GLOB",
    help="Glob pattern for station names (e.g. 'S0*', 'KB0-*').",
)
@click.option(
    "--has-z",
    is_flag=True,
    default=False,
    help="Keep only stations that include complex impedance (Z) data.",
)
@click.option(
    "--has-r",
    is_flag=True,
    default=False,
    help="Keep only stations that include resistivity/phase (R) data.",
)
@click.option(
    "--has-t",
    is_flag=True,
    default=False,
    help="Keep only stations that include tipper (T) data.",
)
@click.option(
    "--min-nfreq",
    type=click.IntRange(min=1),
    default=None,
    metavar="N",
    help="Keep only stations with at least N frequency points.",
)
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
    help="List matching stations without writing any files.",
)
@overwrite_option
@verbose_option
@no_color_option
@click.pass_context
def select(
    ctx: click.Context,
    source: Path,
    output_dir: Path,
    station_ids: str | None,
    pattern: str | None,
    has_z: bool,
    has_r: bool,
    has_t: bool,
    min_nfreq: int | None,
    output_format: str,
    dry_run: bool,
    overwrite: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Filter a J-file collection and export the selected subset.

    All filters are applied with logical AND.
    ``--output-dir`` is required unless ``--dry-run`` is used.

    \b
    Filter options:
      --stations   exact station IDs (comma-separated)
      --pattern    glob for station names
      --has-z      must include complex impedance data
      --has-r      must include resistivity/phase data
      --has-t      must include tipper data
      --min-nfreq  minimum number of frequency points

    \b
    Examples:
      pycsamt jones select data/j/ --has-z --output-dir subset/
      pycsamt jones select data/j/ --stations S01,S02 --output-dir subset/
      pycsamt jones select data/j/ --pattern "KB0-*" --min-nfreq 20 --dry-run
      pycsamt jones select data/j/ --has-r --has-t --output-dir subset/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    import fnmatch  # noqa: PLC0415

    coll = _get_collection(source, verbose=verbose)
    if len(coll) == 0:
        click.echo(f"No J-files found in {source}.", err=True)
        sys.exit(1)

    keep = list(coll)

    if station_ids is not None:
        ids = {s.strip() for s in station_ids.split(",") if s.strip()}
        keep = [j for j in keep if (j.station or "") in ids]

    if pattern is not None:
        keep = [j for j in keep if fnmatch.fnmatch(j.station or "", pattern)]

    if has_z:
        keep = [
            j
            for j in keep
            if j.Z is not None and getattr(j.Z, "z", None) is not None
        ]

    if has_r:
        keep = [j for j in keep if j.Res is not None]

    if has_t:
        keep = [j for j in keep if j.Tip is not None]

    if min_nfreq is not None:
        keep = [j for j in keep if int(j.n_freq or 0) >= min_nfreq]

    n_selected = len(keep)
    selected_names = [j.station or "?" for j in keep]

    if dry_run:
        msg = f"Dry run — {n_selected}/{len(coll)} station(s) match:"
        if output_format == "json":
            click.echo(
                json.dumps(
                    {
                        "n_total": len(coll),
                        "n_selected": n_selected,
                        "selected": selected_names,
                    },
                    indent=2,
                )
            )
        else:
            click.echo(msg)
            for name in selected_names:
                click.echo(f"  {name}")
        return

    if str(output_dir) == ".":
        raise click.UsageError(
            "--output-dir is required (use --dry-run to preview)"
        )

    if n_selected == 0:
        click.echo(
            "No stations match the filter — nothing to export.", err=True
        )
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    from pycsamt.jones.collection import (
        JCollection,  # noqa: PLC0415
    )

    sub = JCollection(items=keep, verbose=verbose)
    # Redirect tqdm to stderr so JSON stdout stays clean
    try:
        import sys as _sys

        import tqdm as _tqdm  # noqa: PLC0415

        _orig = _tqdm.tqdm.__init__

        def _patched_init(self, *a, file=None, **kw):
            _orig(self, *a, file=_sys.stderr, **kw)

        _tqdm.tqdm.__init__ = _patched_init
        result = sub.export(output_dir)
        _tqdm.tqdm.__init__ = _orig
    except (ImportError, Exception):
        result = sub.export(output_dir)

    n_ok = len(result.get("successful", []))
    n_fail = len(result.get("failed", []))

    summary = {
        "source": str(source),
        "output_dir": str(output_dir),
        "n_total": len(coll),
        "n_selected": n_selected,
        "n_written": n_ok,
        "n_failed": n_fail,
        "selected": selected_names,
        "written": result.get("successful", []),
        "failed": [
            {"station": s, "error": str(e)}
            for s, e in result.get("failed", [])
        ],
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2, default=str))
    else:
        click.echo(
            f"Selected : {n_selected}/{len(coll)}"
            f"  | Written: {n_ok}"
            f"  | Errors: {n_fail}"
            f"  | Output: {output_dir}/"
        )
        for s, e in result.get("failed", []):
            click.echo(f"  ✗  {s}: {e}", err=True)

    if n_fail and not n_ok:
        sys.exit(1)
