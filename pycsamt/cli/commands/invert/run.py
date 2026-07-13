# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt invert run — launch an inversion binary."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    verbose_option,
)
from ._base import _resolve_solver, invert


@invert.command("run")
@click.argument(
    "workdir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="WORKDIR",
)
@click.option(
    "--solver",
    type=click.Choice(["occam2d", "modem"], case_sensitive=False),
    default=None,
    help="Inversion code (auto-detected from WORKDIR if omitted).",
)
@click.option(
    "--max-iter",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Override max iterations in the startup/control file.",
)
@click.option(
    "--target-misfit",
    type=float,
    default=None,
    metavar="FLOAT",
    help="Override target RMS misfit.",
)
@click.option(
    "--async",
    "run_async",
    is_flag=True,
    default=False,
    help="Launch non-blocking (returns immediately after starting).",
)
@verbose_option
@no_color_option
@click.pass_context
def run(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    max_iter: int | None,
    target_misfit: float | None,
    run_async: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Launch an inversion binary inside WORKDIR.

    The solver is auto-detected from files present in WORKDIR unless
    --solver is passed explicitly.

    \b
    Examples:
      pycsamt invert run run01/
      pycsamt invert run run01/ --max-iter 100 --target-misfit 1.05
      pycsamt invert run modem_run/ --solver modem --async
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)

    if verbose >= 1:
        click.echo(f"Starting {solver_name.upper()} in {workdir}/", err=True)

    try:
        if solver_name == "occam2d":
            _run_occam2d(workdir, max_iter, target_misfit, run_async, verbose)
        else:
            _run_modem(workdir, max_iter, target_misfit, run_async, verbose)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _run_occam2d(
    workdir: Path,
    max_iter: int | None,
    target_misfit: float | None,
    run_async: bool,
    verbose: int,
) -> None:
    from pycsamt.models.occam2d.runner import (
        OccamRunner,  # noqa: PLC0415
    )

    runner = OccamRunner(workdir=workdir, verbose=verbose)
    if run_async:
        pid = runner.run_async(max_iter=max_iter, target_misfit=target_misfit)
        click.echo(
            f"Occam2D started (PID {pid}).  Logs: {workdir}/occam_stdout.log"
        )
    else:
        code = runner.run(max_iter=max_iter, target_misfit=target_misfit)
        if code != 0:
            click.echo(
                f"Occam2D exited with code {code}.  "
                f"See {workdir}/occam_stderr.log",
                err=True,
            )
            sys.exit(code)
        click.echo("Occam2D finished successfully.")


def _run_modem(
    workdir: Path,
    max_iter: int | None,
    target_misfit: float | None,
    run_async: bool,
    verbose: int,
) -> None:
    from pycsamt.models.modem.runner import (
        ModEmRunner,  # noqa: PLC0415
    )

    runner = ModEmRunner(workdir=workdir, verbose=verbose)
    if run_async:
        pid = runner.run(
            run_async=True, max_iterations=max_iter, target_rms=target_misfit
        )
        click.echo(
            f"ModEM started (PID {pid}).  Logs: {workdir}/modem_stdout.log"
        )
    else:
        code = runner.run(max_iterations=max_iter, target_rms=target_misfit)
        if code != 0:
            click.echo(
                f"ModEM exited with code {code}.  "
                f"See {workdir}/modem_stderr.log",
                err=True,
            )
            sys.exit(code)
        click.echo("ModEM finished successfully.")
