# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt invert results — summarise a completed inversion."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option

from ._base import invert, _load_inversion_result, _resolve_solver, _rich_table


@invert.command("results")
@click.argument(
    "workdir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="WORKDIR",
)
@click.option(
    "--solver",
    type=click.Choice(["occam2d", "modem"], case_sensitive=False),
    default=None,
    help="Inversion code (auto-detected if omitted).",
)
@click.option(
    "--iteration",
    type=int,
    default=None,
    metavar="INT",
    help="Iteration number to load (default: last available).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def results(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Summarise a completed inversion run.

    Loads InversionResult from WORKDIR and reports the best iteration,
    final RMS, model dimensions, and log10-resistivity statistics.

    \b
    Examples:
      pycsamt invert results run01/
      pycsamt invert results run01/ --iteration 50
      pycsamt invert results run01/ --format json
      pycsamt invert results modem_run/ --solver modem
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)

    try:
        info = _results_dict(workdir, solver_name, iteration, verbose)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error loading results: {exc}", err=True)
        sys.exit(1)

    if output_format == "json":
        click.echo(json.dumps(info, indent=2, default=str))
    else:
        _print_results(info, solver_name)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _results_dict(
    workdir: Path, solver: str, iteration: int | None, verbose: int
) -> dict[str, Any]:
    result = _load_inversion_result(workdir, solver, iteration, verbose)

    info: dict[str, Any] = {
        "workdir":      str(workdir),
        "solver":       solver,
        "iteration":    None,
        "rms":          None,
        "model_shape":  None,
        "rho_min":      None,
        "rho_max":      None,
        "rho_mean":     None,
        "n_iter_files": None,
    }

    try:
        info["n_iter_files"] = len(result.iter_files)
    except AttributeError:
        pass

    try:
        best = result.best_iter
        if best is not None:
            info["iteration"] = getattr(best, "iteration", None)
            info["rms"] = (
                getattr(best, "misfit", None) or getattr(best, "rms", None)
            )
    except AttributeError:
        pass

    try:
        import numpy as np  # noqa: PLC0415
        rho = result.rho_2d
        if rho is not None:
            finite = rho[np.isfinite(rho)]
            if finite.size:
                info["model_shape"] = list(rho.shape)
                info["rho_min"]  = round(float(finite.min()), 4)
                info["rho_max"]  = round(float(finite.max()), 4)
                info["rho_mean"] = round(float(finite.mean()), 4)
    except AttributeError:
        pass

    return info


def _print_results(info: dict[str, Any], solver: str) -> None:
    rows: list[tuple[str, str]] = [
        ("Workdir",     info["workdir"]),
        ("Solver",      solver.upper()),
        ("Iter files",  str(info["n_iter_files"] or "—")),
        ("Loaded iter", str(info["iteration"]    or "—")),
        ("Final RMS",   str(info["rms"])          if info["rms"]          is not None else "—"),
        ("Model shape", str(info["model_shape"])  if info["model_shape"]  else "—"),
        ("log10ρ min",  str(info["rho_min"])      if info["rho_min"]      is not None else "—"),
        ("log10ρ max",  str(info["rho_max"])      if info["rho_max"]      is not None else "—"),
        ("log10ρ mean", str(info["rho_mean"])     if info["rho_mean"]     is not None else "—"),
    ]
    _rich_table("Inversion Results", rows, style="green")
