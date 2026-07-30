# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt forward run — execute a single 1-D EM forward model.

Accepts a model either as explicit layer parameters
(``--resistivities`` + ``--thicknesses``) or drawn from a pre-defined
geological scenario (``--geology``).

Usage
-----
::

    # explicit model, MT solver, default 30 log-spaced freqs
    pycsamt forward run \\
        --resistivities 100,10,500 \\
        --thicknesses 300,800

    # geology-based random model
    pycsamt forward run --geology sedimentary --solver mt

    # TEM solver, JSON output
    pycsamt forward run \\
        --resistivities 200,5,1000 \\
        --thicknesses 100,500 \\
        --solver tem \\
        --format json

    # restrict frequency band, CSV for further processing
    pycsamt forward run \\
        --resistivities 50,500,20,800 \\
        --thicknesses 200,600,1500 \\
        --freqs 0.01:1000 \\
        --n-freqs 40 \\
        --format csv > response.csv
"""

from __future__ import annotations

import sys
from typing import Any

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ....api.cli.params import FreqRange
from ._base import forward

# ---------------------------------------------------------------------------
# run command
# ---------------------------------------------------------------------------


@forward.command("run")
@click.option(
    "--resistivities",
    "-r",
    default=None,
    metavar="R1,R2,...",
    help=(
        "Comma-separated layer resistivities in Ω·m, ordered top→bottom.  "
        "Required unless --geology is given."
    ),
)
@click.option(
    "--thicknesses",
    "-t",
    default=None,
    metavar="T1,T2,...",
    help=(
        "Comma-separated layer thicknesses in metres.  "
        "Must have exactly one fewer value than --resistivities.  "
        "Required unless --geology is given."
    ),
)
@click.option(
    "--geology",
    "-g",
    default=None,
    metavar="SCENARIO",
    help=(
        "Draw a random LayeredModel from a pre-defined geological scenario.  "
        "One of: sedimentary, crystalline, geothermal, marine, permafrost.  "
        "Cannot be combined with --resistivities / --thicknesses."
    ),
)
@click.option(
    "--seed",
    type=int,
    default=None,
    metavar="INT",
    help="Random seed for reproducible model generation (--geology only).",
)
@click.option(
    "--solver",
    type=click.Choice(["mt", "tem", "csamt"], case_sensitive=False),
    default="mt",
    show_default=True,
    help="1-D EM solver to use.",
)
@click.option(
    "--freqs",
    "freq_range",
    type=FreqRange(),
    default="0.001:10000",
    show_default=True,
    metavar="FMIN:FMAX",
    help="Frequency range in Hz.",
)
@click.option(
    "--n-freqs",
    type=click.IntRange(min=2),
    default=30,
    show_default=True,
    metavar="INT",
    help="Number of log-spaced frequencies between FMIN and FMAX.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def run(
    ctx: click.Context,
    resistivities: str | None,
    thicknesses: str | None,
    geology: str | None,
    seed: int | None,
    solver: str,
    freq_range: tuple[float, float],
    n_freqs: int,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Run a single 1-D EM forward model and print the response table.

    Define the earth model via ``--resistivities`` + ``--thicknesses``, or
    use ``--geology`` to draw a random model from a geological scenario.
    The ``--geology`` flag uses the n_layers range configured in the chosen
    scenario prior; pass ``--seed`` for reproducibility.

    \b
    Output columns (MT / CSAMT):
      freq_Hz   period_s   rho_a_Ohm_m   phase_deg

    \b
    Output columns (TEM):
      time_s   dBz_dt_T_per_s

    \b
    Examples:
      pycsamt forward run \\
          --resistivities 100,10,500 --thicknesses 300,800

      pycsamt forward run --geology sedimentary --seed 42

      pycsamt forward run \\
          --resistivities 200,5,1000 --thicknesses 100,500 \\
          --solver tem --format json

      pycsamt forward run \\
          --resistivities 50,500,20,800 --thicknesses 200,600,1500 \\
          --freqs 0.01:1000 --n-freqs 40 --format csv > response.csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    import numpy as np  # noqa: PLC0415

    from ....forward import (  # noqa: PLC0415
        GEOLOGY_PRIORS,
        CSAMT1DForward,
        LayeredModel,
        MT1DForward,
        TEM1DForward,
    )

    # ── build model ──────────────────────────────────────────────────────────
    try:
        if geology is not None and (resistivities or thicknesses):
            raise click.UsageError(
                "--geology cannot be combined with "
                "--resistivities / --thicknesses."
            )

        if geology is not None:
            key = geology.lower().strip()
            if key not in GEOLOGY_PRIORS:
                available = ", ".join(sorted(GEOLOGY_PRIORS))
                raise click.BadParameter(
                    f"{geology!r} is not a known scenario.  "
                    f"Available: {available}",
                    param_hint="--geology",
                )
            model = LayeredModel.from_geology(key, seed=seed)

        elif resistivities is not None:
            if thicknesses is None:
                raise click.UsageError(
                    "--thicknesses is required when --resistivities is given."
                )
            try:
                rho = [float(x) for x in resistivities.split(",")]
                thk = [float(x) for x in thicknesses.split(",")]
            except ValueError as exc:
                raise click.UsageError(f"Could not parse layer values: {exc}")

            if len(thk) != len(rho) - 1:
                raise click.UsageError(
                    f"Expected {len(rho) - 1} thickness value(s) "
                    f"for {len(rho)} layers; got {len(thk)}."
                )
            model = LayeredModel(
                resistivity=np.asarray(rho, float),
                thickness=np.asarray(thk, float),
            )

        else:
            raise click.UsageError(
                "Provide --resistivities + --thicknesses, or --geology."
            )

        freqs = np.logspace(
            np.log10(freq_range[0]),
            np.log10(freq_range[1]),
            n_freqs,
        )

        solver_key = solver.lower()
        if solver_key == "mt":
            resp = MT1DForward(freqs).run(model)
        elif solver_key == "tem":
            resp = TEM1DForward(freqs).run(model)
        else:
            resp = CSAMT1DForward(freqs).run(model)

    except (click.exceptions.UsageError, click.exceptions.BadParameter):
        raise
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    if verbose >= 1:
        n = len(model.resistivity)
        click.echo(
            f"Model: {n} layers  solver={solver.upper()}  "
            f"freqs={freq_range[0]:.3g}–{freq_range[1]:.3g} Hz  "
            f"n_freqs={n_freqs}",
            err=True,
        )

    _emit(resp, freqs, solver_key, output_format)


# ---------------------------------------------------------------------------
# Output helper
# ---------------------------------------------------------------------------


def _emit(resp: Any, freqs: Any, solver: str, output_format: str) -> None:
    import numpy as np  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415

    if solver == "tem":
        times = np.atleast_1d(resp.times if resp.times is not None else freqs)
        vals = np.atleast_1d(resp.dBz_dt)
        df = pd.DataFrame(
            {
                "time_s": times,
                "dBz_dt_T_per_s": vals,
            }
        )
    else:
        rho_a = np.atleast_1d(resp.rho_a)
        phase = np.atleast_1d(resp.phase)
        df = pd.DataFrame(
            {
                "freq_Hz": freqs,
                "period_s": np.where(freqs > 0, 1.0 / freqs, float("inf")),
                "rho_a_Ohm_m": rho_a,
                "phase_deg": phase,
            }
        )

    if output_format == "json":
        click.echo(df.to_json(orient="records", indent=2, default_handler=str))
    elif output_format == "csv":
        click.echo(df.to_csv(index=False))
    else:
        click.echo(df.to_string(index=False))
