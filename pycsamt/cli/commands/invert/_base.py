# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt invert — root Click group and shared helpers.

Every sub-module in this package imports ``invert`` from here and
decorates its commands with ``@invert.command(...)`` so that
``cli/commands/invert/__init__.py`` only needs to import the sub-modules
to trigger the registrations.

Shared utilities
----------------
_detect_solver      Infer occam2d or modem from workdir file signatures.
_resolve_solver     Return an explicit solver, raising UsageError if unknown.
_rich_table         Print a key/value table via rich (degrades to plain text).
_load_inversion_result  Load InversionResult for the given solver.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

# ---------------------------------------------------------------------------
# Solver fingerprints
# ---------------------------------------------------------------------------

_OCCAM_SIGNATURES = frozenset({
    "OccamDataFile.dat",
    "Occam2DMesh",
    "Occam2DModel",
    "Startup",
    # Backward-compatible aliases from early CLI fixtures and examples.
    "OccamData.dat",
    "Occam2DStartup",
})
_MODEM_SIGNATURES  = frozenset({"ModEMData.dat", "ModEM.inv", "Modular_NLCG.log"})


def _detect_solver(workdir: Path) -> str | None:
    """Return ``'occam2d'``, ``'modem'``, or ``None`` from workdir contents."""
    names = {p.name for p in workdir.iterdir() if p.is_file()}
    if names & _OCCAM_SIGNATURES:
        return "occam2d"
    if names & _MODEM_SIGNATURES:
        return "modem"
    if any(p.suffix in {".iter", ".resp"} for p in workdir.iterdir()):
        return "occam2d"
    return None


def _resolve_solver(solver: str | None, workdir: Path) -> str:
    """Return an explicit solver, auto-detecting from *workdir* when needed."""
    if solver:
        return solver.lower()
    detected = _detect_solver(workdir)
    if detected:
        return detected
    raise click.UsageError(
        f"Cannot detect solver from {workdir}.  "
        "Pass --solver occam2d or --solver modem explicitly."
    )


# ---------------------------------------------------------------------------
# Rich helper
# ---------------------------------------------------------------------------

def _rich_table(
    title: str,
    rows: list[tuple[str, str]],
    style: str = "cyan",
) -> None:
    """Print a two-column key/value table via rich (plain fallback)."""
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table  import Table     # noqa: PLC0415
        console = Console()
        t = Table(title=title, border_style=style, show_header=False)
        t.add_column("Key",   style="bold")
        t.add_column("Value", style="white")
        for k, v in rows:
            t.add_row(k, v)
        console.print(t)
    except ImportError:
        click.echo(f"\n{title}")
        for k, v in rows:
            click.echo(f"  {k:<28} {v}")


# ---------------------------------------------------------------------------
# Shared InversionResult loader
# ---------------------------------------------------------------------------

def _load_inversion_result(
    workdir: Path,
    solver: str,
    iteration: int | None = None,
    verbose: int = 0,
) -> Any:
    """Instantiate and return an InversionResult for *solver*."""
    if solver == "occam2d":
        from pycsamt.models.occam2d.results import InversionResult  # noqa: PLC0415
    else:
        from pycsamt.models.modem.results import InversionResult    # noqa: PLC0415

    kw: dict[str, Any] = {"workdir": workdir, "verbose": verbose}
    if iteration is not None:
        kw["iteration"] = iteration
    return InversionResult(**kw)


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------

@click.group("invert")
@click.pass_context
def invert(ctx: click.Context) -> None:
    """Build, run, and inspect AMT/MT inversions (Occam2D or ModEM).

    \b
    Typical workflow:
      pycsamt invert build   [EDI_DIR]  --solver occam2d --workdir run01/
      pycsamt invert run     run01/     --max-iter 100
      pycsamt invert status  run01/
      pycsamt invert results run01/
      pycsamt invert plot    model run01/
    """
    ctx.ensure_object(dict)
