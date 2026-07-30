# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt avg correct — apply static-shift or capacitive-coupling correction."""

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
from ._base import _get_avg, avg


@avg.command("correct")
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="AVG_FILE",
)
@click.option(
    "--method",
    "-m",
    type=click.Choice(
        ["static-shift", "capacitive", "both"], case_sensitive=False
    ),
    default="static-shift",
    show_default=True,
    help=(
        "Correction to apply.  "
        "static-shift: remove spatial static-shift via spatial filtering.  "
        "capacitive: remedy high-frequency E-field distortion.  "
        "both: apply static-shift then capacitive-coupling."
    ),
)
@click.option(
    "--filter",
    "filter_method",
    type=click.Choice(["tma", "flma", "ama"], case_sensitive=False),
    default="tma",
    show_default=True,
    help=(
        "Spatial filter for static-shift correction.  "
        "tma: Trimmed Moving Average (works on ρa).  "
        "flma: Fixed-Length Moving Average (works on Z).  "
        "ama: Adaptive Moving Average (works on Z)."
    ),
)
@click.option(
    "--ref-freq",
    "ref_freq",
    type=float,
    default=None,
    metavar="HZ",
    help=(
        "Reference frequency (Hz) for static-shift correction.  "
        "When omitted, the highest available frequency is used."
    ),
)
@click.option(
    "--window",
    type=int,
    default=5,
    show_default=True,
    metavar="N",
    help="Window size for the TMA filter (number of stations).",
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
    help="Show what corrections would be applied without writing files.",
)
@verbose_option
@no_color_option
@click.pass_context
def correct(
    ctx: click.Context,
    source: Path,
    method: str,
    filter_method: str,
    ref_freq: float | None,
    window: int,
    output_dir: Path,
    output_format: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Apply physics-based corrections to a Zonge AVG file.

    Static-shift correction estimates a spatially smooth resistivity
    profile at a reference frequency and computes per-station shift
    factors.  Capacitive-coupling correction remedies high-frequency
    E-field distortion due to cable-ground capacitance.

    The corrected file is written to ``--output-dir`` (default: current
    directory) in modern (kind-2) format.

    \b
    Examples:
      # dry run — see what would be corrected:
      pycsamt avg correct data/avg/K2.AVG --dry-run

      # apply TMA static-shift at highest frequency:
      pycsamt avg correct data/avg/K2.AVG --method static-shift \\
          --output-dir corrected/

      # use a specific reference frequency:
      pycsamt avg correct data/avg/K2.AVG --ref-freq 2048 \\
          --filter tma --output-dir corrected/

      # apply both corrections:
      pycsamt avg correct data/avg/K2.AVG --method both \\
          --output-dir corrected/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    obj = _get_avg(source, verbose=verbose)

    from pycsamt.zonge.processing import (
        ASTATIC,  # noqa: PLC0415
    )

    proc = ASTATIC(verbose=bool(verbose)).read(obj)

    # Resolve reference frequency
    df = obj.df
    if df is None:
        click.echo("No data loaded.", err=True)
        sys.exit(1)

    freqs = sorted(df["freq"].unique())
    if ref_freq is None:
        ref_freq = float(freqs[-1])

    if ref_freq not in freqs:
        # pick nearest
        import numpy as np  # noqa: PLC0415

        ref_freq = float(
            freqs[int(np.argmin(np.abs(np.array(freqs) - ref_freq)))]
        )
        click.echo(
            f"Note: adjusted ref-freq to nearest available: {ref_freq} Hz",
            err=True,
        )

    if dry_run:
        click.echo(
            f"Dry run — {source.name}:\n"
            f"  Method       : {method}\n"
            f"  Filter       : {filter_method} (window={window})\n"
            f"  Ref freq     : {ref_freq} Hz\n"
            f"  Stations     : {len(freqs and df['station'].unique())}\n"
            f"  No files written."
        )
        return

    results: list[dict] = []

    # --- Static shift ---
    if method in ("static-shift", "both"):
        kwargs = {"window_size": window} if filter_method == "tma" else {}
        ss_df = proc.correct_static_shift(
            reference_freq=ref_freq,
            filter_method=filter_method,
            update_components=True,
            **kwargs,
        )
        results.append(
            {
                "correction": "static_shift",
                "filter": filter_method,
                "ref_freq": ref_freq,
                "n_stations": len(ss_df),
                "shift_min": float(ss_df["shift_factor"].min()),
                "shift_max": float(ss_df["shift_factor"].max()),
                "shift_mean": float(ss_df["shift_factor"].mean()),
            }
        )

    # --- Capacitive coupling ---
    if method in ("capacitive", "both"):
        try:
            proc.correct_capacitive_coupling(update_components=True)
            results.append({"correction": "capacitive_coupling"})
        except Exception as exc:  # noqa: BLE001
            click.echo(
                f"Warning: capacitive correction failed: {exc}", err=True
            )

    # Write output
    out_name = source.stem + "_corrected.avg"
    out_path = output_dir / out_name
    output_dir.mkdir(parents=True, exist_ok=True)
    proc.avg.write(out_path, fmt="modern")

    summary = {
        "source": str(source),
        "output": str(out_path),
        "corrections": results,
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2, default=str))
    else:
        click.echo(f"Written: {out_path}")
        for r in results:
            if r["correction"] == "static_shift":
                click.echo(
                    f"  Static shift: filter={r['filter']}, "
                    f"ref={r['ref_freq']} Hz, "
                    f"shift={r['shift_min']:.3f}–{r['shift_max']:.3f} "
                    f"(mean {r['shift_mean']:.3f})"
                )
            else:
                click.echo("  Capacitive coupling corrected")
