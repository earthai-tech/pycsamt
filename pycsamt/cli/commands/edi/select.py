# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt edi select — filter a collection and export a subset."""

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
from ._base import _get_collection, edi


@edi.command("select")
@click.argument(
    "source",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="EDI_DIR",
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
    help="Glob pattern for station names (e.g. 'S0*', '18-02*').",
)
@click.option(
    "--has-tipper",
    is_flag=True,
    default=False,
    help="Keep only stations that include Tipper data.",
)
@click.option(
    "--min-freq",
    type=float,
    default=None,
    metavar="HZ",
    help="Keep only stations whose minimum frequency ≤ MIN_FREQ.",
)
@click.option(
    "--max-freq",
    type=float,
    default=None,
    metavar="HZ",
    help="Keep only stations whose maximum frequency ≥ MAX_FREQ.",
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
    help="List which stations would be selected without writing any files.",
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
    has_tipper: bool,
    min_freq: float | None,
    max_freq: float | None,
    min_nfreq: int | None,
    output_format: str,
    dry_run: bool,
    overwrite: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Filter an EDI collection and export the selected subset.

    Multiple filters are applied with logical AND.
    ``--output-dir`` is required unless ``--dry-run`` is used.

    \b
    Filter options:
      --stations   exact station IDs (comma-separated)
      --pattern    glob for station names
      --has-tipper only stations with Tipper data
      --min-freq   minimum frequency threshold
      --max-freq   maximum frequency threshold
      --min-nfreq  minimum number of frequency points

    \b
    Examples:
      # export two named stations:
      pycsamt edi select data/AMT/ \\
          --stations HBH03,HBH05 --output-dir subset/

      # all stations whose name starts with "18-02":
      pycsamt edi select data/AMT/WILLY_DATA/L18PLT/ \\
          --pattern "18-02*" --output-dir subset/

      # stations with tipper and at least 40 frequencies:
      pycsamt edi select data/AMT/ \\
          --has-tipper --min-nfreq 40 --output-dir subset/

      # dry run — preview the selection:
      pycsamt edi select data/AMT/ --has-tipper --dry-run
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    import fnmatch  # noqa: PLC0415

    coll = _get_collection(source, verbose=verbose)
    if len(coll) == 0:
        click.echo(f"No EDI files found in {source}.", err=True)
        sys.exit(1)

    # -- build filtering pipeline --
    keep = list(coll)

    if station_ids is not None:
        ids = {s.strip() for s in station_ids.split(",") if s.strip()}
        keep = [e for e in keep if (e.station or "") in ids]

    if pattern is not None:
        keep = [e for e in keep if fnmatch.fnmatch(e.station or "", pattern)]

    if has_tipper:
        keep = [e for e in keep if e.has_tipper]

    if min_freq is not None:

        def _fmin(e) -> float:
            freq = getattr(e.Z, "freq", None)
            return (
                float(freq.min())
                if freq is not None and freq.size
                else float("inf")
            )

        keep = [e for e in keep if _fmin(e) <= min_freq]

    if max_freq is not None:

        def _fmax(e) -> float:
            freq = getattr(e.Z, "freq", None)
            return (
                float(freq.max())
                if freq is not None and freq.size
                else float("-inf")
            )

        keep = [e for e in keep if _fmax(e) >= max_freq]

    if min_nfreq is not None:
        keep = [
            e for e in keep if int(getattr(e.Z, "n_freq", 0) or 0) >= min_nfreq
        ]

    n_selected = len(keep)
    selected_names = [e.station or "?" for e in keep]

    if dry_run:
        msg = (
            f"Dry run — {n_selected}/{len(coll)} station(s) match the filter:"
        )
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

    from pycsamt.seg.collection import (
        EDICollection,  # noqa: PLC0415
    )

    sub = EDICollection(items=keep, verbose=verbose)
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
        click.echo(json.dumps(summary, indent=2))
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
