# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt invert build — build inversion input files from an EDI collection."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    fresh_option,
    no_color_option,
    survey_option,
    verbose_option,
)
from ....api.cli.params import FreqRange
from ...survey import resolve_survey
from ._base import invert


@invert.command("build")
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
    "--solver",
    type=click.Choice(["occam2d", "modem"], case_sensitive=False),
    default="occam2d",
    show_default=True,
    help="Inversion code to build inputs for.",
)
@click.option(
    "--workdir", "-w",
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    default=Path("invert_run"),
    show_default=True,
    help="Working directory for inversion files.",
)
# --- data / frequency ---
@click.option(
    "--modes",
    default="TE,TM",
    show_default=True,
    metavar="TE,TM[,DET]",
    help="Occam2D only: comma-separated EM modes.",
)
@click.option(
    "--freq", "freq_range",
    type=FreqRange(),
    default=None,
    metavar="FMIN:FMAX",
    help="Restrict frequencies to this Hz range (e.g. 0.1:1000).",
)
@click.option(
    "--error-floor-rho",
    type=float,
    default=None,
    metavar="FLOAT",
    help="Relative apparent-resistivity error floor (e.g. 0.05).",
)
@click.option(
    "--error-floor-phase",
    type=float,
    default=None,
    metavar="FLOAT",
    help="Absolute phase error floor in degrees.",
)
# --- mesh / geometry ---
@click.option(
    "--n-layers",
    type=int,
    default=None,
    metavar="INT",
    help="Number of earth layers (Occam2D) or vertical cells (ModEM).",
)
@click.option(
    "--cell-size",
    type=float,
    default=None,
    metavar="METRES",
    help="Target horizontal cell size in metres near stations.",
)
# --- ModEM-specific ---
@click.option(
    "--modem-mode",
    type=click.Choice(["2d", "3d"], case_sensitive=False),
    default="3d",
    show_default=True,
    help="ModEM only: dimensionality of the inversion.",
)
@click.option(
    "--initial-rho",
    type=float,
    default=None,
    metavar="OHM_M",
    help="ModEM only: initial half-space resistivity (Ω·m).",
)
@verbose_option
@no_color_option
@click.pass_context
def build(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    solver: str,
    workdir: Path,
    modes: str,
    freq_range: tuple[float, float] | None,
    error_floor_rho: float | None,
    error_floor_phase: float | None,
    n_layers: int | None,
    cell_size: float | None,
    modem_mode: str,
    initial_rho: float | None,
    verbose: int,
    no_color: bool,
) -> None:
    """Build inversion input files from an EDI collection.

    EDI_DIR is optional when an active survey has been set with
    ``pycsamt survey set``.  Priority: EDI_DIR > --survey > active context.

    \b
    Examples:
      # active survey already set:
      pycsamt invert build --solver occam2d --workdir run01/

      # explicit path:
      pycsamt invert build survey/ --workdir run01/
      pycsamt invert build survey/ --freq 0.1:1000 --n-layers 60 \\
              --error-floor-rho 0.05 --workdir run01/
      pycsamt invert build survey/ --solver modem --modem-mode 3d \\
              --initial-rho 100 --workdir modem_run/

      # force re-parse from disk:
      pycsamt invert build --fresh --workdir run01/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    explicit = edi_source or survey_path
    sites = resolve_survey(explicit, fresh=fresh, verbose=verbose)

    solver = solver.lower()
    workdir.mkdir(parents=True, exist_ok=True)

    if verbose >= 1:
        click.echo(
            f"Building {solver.upper()} inputs "
            f"({len(sites)} stations) → {workdir}/",
            err=True,
        )

    try:
        if solver == "occam2d":
            _build_occam2d(
                sites, workdir,
                modes=[m.strip().upper() for m in modes.split(",")],
                freq_range=freq_range,
                error_floor_rho=error_floor_rho,
                error_floor_phase=error_floor_phase,
                n_layers=n_layers,
                cell_size=cell_size,
                verbose=verbose,
            )
        else:
            _build_modem(
                sites, workdir,
                modem_mode=modem_mode,
                freq_range=freq_range,
                initial_rho=initial_rho,
                n_layers=n_layers,
                cell_size=cell_size,
                verbose=verbose,
            )
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    click.echo(f"Build complete.  Inputs written to: {workdir}")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_occam2d(
    sites: Any,
    workdir: Path,
    modes: list[str],
    freq_range: tuple[float, float] | None,
    error_floor_rho: float | None,
    error_floor_phase: float | None,
    n_layers: int | None,
    cell_size: float | None,
    verbose: int,
) -> None:
    from pycsamt.models.occam2d.builder import (
        InputBuilder,  # noqa: PLC0415
    )
    from pycsamt.models.occam2d.config import (
        OccamConfig,  # noqa: PLC0415
    )

    cfg = OccamConfig(modes=modes)
    if error_floor_rho   is not None: cfg.error_floor_rho       = error_floor_rho
    if error_floor_phase is not None: cfg.error_floor_phase      = error_floor_phase
    if freq_range        is not None: cfg.freq_min, cfg.freq_max = freq_range
    if n_layers          is not None: cfg.n_layers               = n_layers
    if cell_size         is not None: cfg.cell_size_horizontal   = cell_size

    InputBuilder(source=sites, workdir=workdir, config=cfg, verbose=verbose).build()


def _build_modem(
    sites: Any,
    workdir: Path,
    modem_mode: str,
    freq_range: tuple[float, float] | None,
    initial_rho: float | None,
    n_layers: int | None,
    cell_size: float | None,
    verbose: int,
) -> None:
    from pycsamt.models.modem.builder import (
        InputBuilder,  # noqa: PLC0415
    )
    from pycsamt.models.modem.config import (
        ModEmConfig,  # noqa: PLC0415
    )

    cfg = ModEmConfig(mode=modem_mode)
    if freq_range   is not None: cfg.freq_min, cfg.freq_max = freq_range
    if initial_rho  is not None: cfg.initial_rho            = initial_rho
    if cell_size    is not None: cfg.cell_size_h = cfg.cell_size_h_2d = cell_size
    if n_layers     is not None: cfg.nz          = cfg.nz_2d          = n_layers

    InputBuilder(config=cfg).build(sites, workdir=workdir)
