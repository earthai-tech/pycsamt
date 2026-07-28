# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt site compute
====================

Sub-group for computing derived geophysical quantities from Sites.

Sub-commands
------------
strike          Geoelectric strike angle (Swift / phase-difference methods).
resistivity     Apparent resistivity at a target frequency.
phase-slope     Phase slope (proxy for depth of investigation).
tipper          Tipper magnitude and direction.

Usage
-----
::

    pycsamt site compute strike
    pycsamt site compute strike --method phase_diff --format json

    pycsamt site compute resistivity --freq 1.0
    pycsamt site compute resistivity --freq 1.0 --how interp

    pycsamt site compute phase-slope
    pycsamt site compute tipper
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    fresh_option,
    no_color_option,
    survey_option,
    verbose_option,
)
from ._base import _get_sites, _unwrap_edis, site

# ---------------------------------------------------------------------------
# compute — sub-group
# ---------------------------------------------------------------------------


@site.group("compute")
@click.pass_context
def compute(ctx: click.Context) -> None:
    """Compute derived geophysical quantities from a survey.

    \b
    Sub-commands:
      strike       Geoelectric strike angle
      resistivity  Apparent resistivity at a target frequency
      phase-slope  Phase slope (depth-of-investigation proxy)
      tipper       Tipper magnitude and direction

    \b
    Examples:
      pycsamt site compute strike
      pycsamt site compute resistivity --freq 1.0
      pycsamt site compute phase-slope --format csv
      pycsamt site compute tipper --format json
    """
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# Shared output helper
# ---------------------------------------------------------------------------


def _emit(data: object, output_format: str) -> None:
    """Print *data* (a dict, DataFrame, or APIFrame) in the requested format."""
    import pandas as pd  # noqa: PLC0415

    # site.compute helpers return an APIFrame wrapper when the API view
    # is enabled; unwrap it or the fallback constructor below strips the
    # column names
    if hasattr(data, "df") and isinstance(data.df, pd.DataFrame):
        data = data.df

    if isinstance(data, dict):
        df = pd.DataFrame([data])
    elif isinstance(data, pd.DataFrame):
        df = data
    else:
        df = pd.DataFrame(data)

    if output_format == "json":
        click.echo(df.to_json(orient="records", indent=2, default_handler=str))
    elif output_format == "csv":
        click.echo(df.to_csv(index=False))
    else:
        click.echo(df.to_string(index=False))


# ---------------------------------------------------------------------------
# strike
# ---------------------------------------------------------------------------


@compute.command("strike")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--method",
    type=click.Choice(["swift", "phase_diff"], case_sensitive=False),
    default="swift",
    show_default=True,
    help="Strike estimation method.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def compute_strike(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    method: str,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Compute the geoelectric strike angle for each station.

    Returns one strike angle per station (in degrees from north).

    \b
    Methods:
      swift       Swift (1967) off-diagonal power minimisation (default)
      phase_diff  Phase-difference method

    \b
    Examples:
      pycsamt site compute strike
      pycsamt site compute strike --method phase_diff --format json
      pycsamt site compute strike --format csv > strike.csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    from pycsamt.site.compute import (
        strike_estimate,  # noqa: PLC0415
    )

    try:
        result = strike_estimate(_unwrap_edis(sites_obj), method=method)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _emit(result, output_format)


# ---------------------------------------------------------------------------
# resistivity
# ---------------------------------------------------------------------------


@compute.command("resistivity")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--freq",
    "target_freq",
    type=float,
    required=True,
    metavar="HZ",
    help="Target frequency in Hz.",
)
@click.option(
    "--how",
    type=click.Choice(["nearest", "interp"], case_sensitive=False),
    default="nearest",
    show_default=True,
    help="How to match the target frequency (nearest or interpolate).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def compute_resistivity(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    target_freq: float,
    how: str,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Compute apparent resistivity (ρ_xy, ρ_yx) at a target frequency.

    \b
    Examples:
      pycsamt site compute resistivity --freq 1.0
      pycsamt site compute resistivity --freq 10.0 --how interp
      pycsamt site compute resistivity --freq 0.1 --format csv > rho.csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    from pycsamt.site.compute import (
        res_at_freq,  # noqa: PLC0415
    )

    try:
        result = res_at_freq(_unwrap_edis(sites_obj), target_freq, how=how)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _emit(result, output_format)


# ---------------------------------------------------------------------------
# phase-slope
# ---------------------------------------------------------------------------


@compute.command("phase-slope")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--band",
    type=str,
    default=None,
    metavar="FMIN:FMAX",
    help=(
        "Frequency band in Hz for the slope regression, e.g. 0.1:1000.  "
        "Defaults to the full frequency range of each site."
    ),
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def compute_phase_slope(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    band: str | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Compute the phase slope (depth-of-investigation proxy).

    The phase slope dφ/d(log f) is a depth proxy derived from the
    impedance phase versus log-frequency curve over a frequency band.
    Provide --band FMIN:FMAX to restrict the regression range.

    \b
    Examples:
      pycsamt site compute phase-slope --band 0.1:1000
      pycsamt site compute phase-slope --band 1:100 --format csv > slope.csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    # Resolve band: either from --band flag or full range of first site
    if band is not None:
        try:
            lo, hi = (float(x) for x in band.split(":"))
            freq_band = (lo, hi)
        except (ValueError, TypeError):
            raise click.BadParameter(
                "Expected FMIN:FMAX (e.g. 0.1:1000)", param_hint="--band"
            )
    else:
        # infer full range from the data
        import numpy as np  # noqa: PLC0415

        all_freqs = []
        for s in sites_obj:
            try:
                f = np.asarray(getattr(s, "freq", None) or [], float)
                if f.size:
                    all_freqs.extend([float(f.min()), float(f.max())])
            except Exception:  # noqa: BLE001
                pass
        if not all_freqs:
            raise click.UsageError(
                "Cannot infer frequency range. Pass --band FMIN:FMAX explicitly."
            )
        freq_band = (min(all_freqs), max(all_freqs))

    from pycsamt.site.compute import (
        phase_slope,  # noqa: PLC0415
    )

    try:
        result = phase_slope(_unwrap_edis(sites_obj), freq_band)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _emit(result, output_format)


# ---------------------------------------------------------------------------
# tipper
# ---------------------------------------------------------------------------


@compute.command("tipper")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--freq",
    "target_freq",
    type=float,
    default=None,
    metavar="HZ",
    help="Evaluate tipper at this frequency only (default: all frequencies).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def compute_tipper(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    target_freq: float | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Compute tipper magnitude and direction for each station.

    \b
    Examples:
      pycsamt site compute tipper
      pycsamt site compute tipper --freq 1.0
      pycsamt site compute tipper --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    sites_obj = _get_sites(edi_source, survey_path, fresh, verbose)

    from pycsamt.site.compute import (
        tipper_magnitude,  # noqa: PLC0415
    )

    try:
        result = tipper_magnitude(_unwrap_edis(sites_obj))
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _emit(result, output_format)
