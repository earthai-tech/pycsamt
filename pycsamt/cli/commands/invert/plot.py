# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt invert plot
===================

Visualise inversion results from a working directory.

Sub-commands map to the plot classes in ``pycsamt.models.occam2d.plot``
and ``pycsamt.models.modem.plot``.  Every sub-command:

* auto-detects the solver from the workdir if ``--solver`` is omitted;
* accepts ``--save FILE`` to write the figure to disk;
* accepts ``--show`` to open an interactive window;
* passes ``--dpi``, ``--cmap``, and ``--iteration`` through to the
  underlying plot class.

Sub-commands
------------
model       2-D resistivity section  (PlotModel / PlotModel2D / PlotModel3D)
misfit      Convergence curve        (PlotMisfit)
response    Observed vs predicted    (PlotResponse)
pseudo      Pseudosection            (PlotPseudo)
section     3-D depth slice          (PlotSection — ModEM only)
1d          1-D depth profiles       (PlotSounding1D — Occam2D only)
per-site    Per-site RMS bar chart   (PlotSiteMisfit — Occam2D only)
grid        Response grid all sites  (PlotResponseGrid — Occam2D only)

Usage
-----
::

    pycsamt invert plot model   run01/ --save section.png
    pycsamt invert plot misfit  run01/ --show
    pycsamt invert plot response run01/ --station S01 --save resp.pdf
    pycsamt invert plot pseudo  run01/ --dpi 200 --save pseudo.png
    pycsamt invert plot section modem_run/ --depth 10000 --save slice.png
    pycsamt invert plot 1d      run01/ --save profiles.png
    pycsamt invert plot per-site run01/ --save site_rms.png
    pycsamt invert plot grid    run01/ --save response_grid.png
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import no_color_option, verbose_option

from ._base import invert, _load_inversion_result, _resolve_solver

# ---------------------------------------------------------------------------
# Shared plot options
# ---------------------------------------------------------------------------

_save_option = click.option(
    "--save", "save_path",
    type=click.Path(path_type=Path),
    default=None,
    metavar="FILE",
    help="Save figure to this file path (e.g. fig.png, fig.pdf, fig.svg).",
)

_show_option = click.option(
    "--show",
    is_flag=True,
    default=False,
    help="Display an interactive matplotlib window.",
)

_dpi_option = click.option(
    "--dpi",
    type=int,
    default=100,
    show_default=True,
    help="Figure resolution in dots per inch.",
)

_cmap_option = click.option(
    "--cmap",
    default="jet_r",
    show_default=True,
    metavar="NAME",
    help="Matplotlib colormap name.",
)

_iteration_option = click.option(
    "--iteration",
    type=int,
    default=None,
    metavar="INT",
    help="Iteration number to load (default: last available).",
)

_solver_option = click.option(
    "--solver",
    type=click.Choice(["occam2d", "modem"], case_sensitive=False),
    default=None,
    help="Inversion code (auto-detected if omitted).",
)


# ---------------------------------------------------------------------------
# Helper: save / show figure
# ---------------------------------------------------------------------------

def _handle_figure(fig: Any, save_path: Path | None, show: bool) -> None:
    """Save and/or display *fig*.  Raises if neither is requested."""
    if save_path is None and not show:
        click.echo(
            "Figure created but not saved.  "
            "Pass --save FILE or --show to display it.",
            err=True,
        )
        return
    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=fig.get_dpi(), bbox_inches="tight")
        click.echo(f"Saved: {save_path}")
    if show:
        import matplotlib.pyplot as plt  # noqa: PLC0415
        plt.show()


# ---------------------------------------------------------------------------
# pycsamt invert plot  (sub-group)
# ---------------------------------------------------------------------------

@invert.group("plot")
@click.pass_context
def plot(ctx: click.Context) -> None:
    """Visualise inversion results from a working directory.

    \b
    Sub-commands:
      model     2-D resistivity section
      misfit    Convergence curve
      response  Observed vs predicted curves
      pseudo    Pseudosection
      section   3-D depth slice   (ModEM only)
      1d        1-D depth profiles (Occam2D only)
      per-site  Per-site RMS bar chart (Occam2D only)
      grid      Full response grid     (Occam2D only)

    \b
    Examples:
      pycsamt invert plot model   run01/ --save section.png
      pycsamt invert plot misfit  run01/ --show
      pycsamt invert plot response run01/ --station S01 --save resp.pdf
    """
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# plot model
# ---------------------------------------------------------------------------

@plot.command("model")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@click.option("--rho-min", type=float, default=1.0,   show_default=True,
              help="Color-scale lower bound (Ω·m).")
@click.option("--rho-max", type=float, default=1000.0, show_default=True,
              help="Color-scale upper bound (Ω·m).")
@click.option("--depth-max", type=float, default=None, metavar="METRES",
              help="Maximum display depth in metres.")
@click.option("--no-stations", is_flag=True, default=False,
              help="Hide station markers on the section.")
@_cmap_option
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_model(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    rho_min: float,
    rho_max: float,
    depth_max: float | None,
    no_stations: bool,
    cmap: str,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot the 2-D resistivity model section.

    \b
    Examples:
      pycsamt invert plot model run01/ --save section.png
      pycsamt invert plot model run01/ --rho-min 1 --rho-max 1000 --depth-max 5000
      pycsamt invert plot model modem_run/ --solver modem --show
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        if solver_name == "occam2d":
            from pycsamt.models.occam2d.plot import PlotModel  # noqa: PLC0415
            fig = PlotModel(
                result,
                rho_min=rho_min,
                rho_max=rho_max,
                depth_max=depth_max,
                show_stations=not no_stations,
                cmap=cmap,
                dpi=dpi,
            ).plot()
        else:
            from pycsamt.models.modem.plot import PlotModel2D  # noqa: PLC0415
            fig = PlotModel2D(
                result,
                rho_min=rho_min,
                rho_max=rho_max,
                depth_max=depth_max,
                cmap=cmap,
                dpi=dpi,
            ).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot misfit
# ---------------------------------------------------------------------------

@plot.command("misfit")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@click.option("--no-roughness", is_flag=True, default=False,
              help="Hide the roughness secondary axis.")
@click.option("--lagrange", is_flag=True, default=False,
              help="Add a lower panel for the Lagrange multiplier.")
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_misfit(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    no_roughness: bool,
    lagrange: bool,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot the inversion convergence curve (RMS vs iteration).

    \b
    Examples:
      pycsamt invert plot misfit run01/ --show
      pycsamt invert plot misfit run01/ --save convergence.png --lagrange
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        if solver_name == "occam2d":
            from pycsamt.models.occam2d.plot import PlotMisfit  # noqa: PLC0415
            fig = PlotMisfit(
                result,
                show_roughness=not no_roughness,
                show_lagrange=lagrange,
                dpi=dpi,
            ).plot()
        else:
            from pycsamt.models.modem.plot import PlotMisfit  # noqa: PLC0415
            fig = PlotMisfit(result, dpi=dpi).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot response
# ---------------------------------------------------------------------------

@plot.command("response")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@click.option("--station", default=None, metavar="NAME",
              help="Show only this station (default: all stations).")
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_response(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    station: str | None,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot observed vs predicted response curves.

    \b
    Examples:
      pycsamt invert plot response run01/ --show
      pycsamt invert plot response run01/ --station S01 --save S01_resp.png
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        if solver_name == "occam2d":
            from pycsamt.models.occam2d.plot import PlotResponse  # noqa: PLC0415
            fig = PlotResponse(result, station=station, dpi=dpi).plot()
        else:
            from pycsamt.models.modem.plot import PlotResponse  # noqa: PLC0415
            fig = PlotResponse(result, station=station, dpi=dpi).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot pseudo
# ---------------------------------------------------------------------------

@plot.command("pseudo")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@_cmap_option
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_pseudo(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    cmap: str,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot the apparent-resistivity / phase pseudosection.

    \b
    Examples:
      pycsamt invert plot pseudo run01/ --show
      pycsamt invert plot pseudo run01/ --cmap viridis --save pseudo.png
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        if solver_name == "occam2d":
            from pycsamt.models.occam2d.plot import PlotPseudo  # noqa: PLC0415
            fig = PlotPseudo(result, cmap=cmap, dpi=dpi).plot()
        else:
            from pycsamt.models.modem.plot import PlotPseudo  # noqa: PLC0415
            fig = PlotPseudo(result, cmap=cmap, dpi=dpi).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot section  (ModEM only — 3-D depth slice)
# ---------------------------------------------------------------------------

@plot.command("section")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@click.option("--depth", type=float, default=None, metavar="METRES",
              help="Depth of the horizontal slice in metres.")
@_cmap_option
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_section(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    depth: float | None,
    cmap: str,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot a 3-D depth-slice section (ModEM only).

    \b
    Examples:
      pycsamt invert plot section modem_run/ --depth 5000 --show
      pycsamt invert plot section modem_run/ --depth 10000 --save slice.png
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    if solver_name != "modem":
        raise click.UsageError(
            "'section' is a ModEM-only sub-command.  "
            "For Occam2D depth profiles use 'pycsamt invert plot 1d'."
        )
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        from pycsamt.models.modem.plot import PlotSection  # noqa: PLC0415
        fig = PlotSection(result, depth=depth, cmap=cmap, dpi=dpi).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot 1d  (Occam2D only — sounding depth profiles)
# ---------------------------------------------------------------------------

@plot.command("1d")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@click.option("--stations", default=None, metavar="S01,S02,…",
              help="Comma-separated station names to include (default: all).")
@_cmap_option
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_1d(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    stations: str | None,
    cmap: str,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot 1-D resistivity–depth profiles at each station (Occam2D only).

    \b
    Examples:
      pycsamt invert plot 1d run01/ --show
      pycsamt invert plot 1d run01/ --stations S01,S05,S10 --save profiles.png
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    if solver_name != "occam2d":
        raise click.UsageError(
            "'1d' is an Occam2D-only sub-command.  "
            "For ModEM depth slices use 'pycsamt invert plot section'."
        )
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)
    station_list = (
        [s.strip() for s in stations.split(",")] if stations else None
    )

    try:
        from pycsamt.models.occam2d.plot import PlotSounding1D  # noqa: PLC0415
        fig = PlotSounding1D(
            result, stations=station_list, cmap=cmap, dpi=dpi
        ).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot per-site  (Occam2D — PlotSiteMisfit)
# ---------------------------------------------------------------------------

@plot.command("per-site")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@_cmap_option
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_per_site(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    cmap: str,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot per-site RMS bar chart + normalised-residual pseudosection (Occam2D).

    \b
    Examples:
      pycsamt invert plot per-site run01/ --show
      pycsamt invert plot per-site run01/ --save site_rms.png
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    if solver_name != "occam2d":
        raise click.UsageError(
            "'per-site' is an Occam2D-only sub-command."
        )
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        from pycsamt.models.occam2d.plot import PlotSiteMisfit  # noqa: PLC0415
        fig = PlotSiteMisfit(result, cmap=cmap, dpi=dpi).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)


# ---------------------------------------------------------------------------
# plot grid  (Occam2D — PlotResponseGrid)
# ---------------------------------------------------------------------------

@plot.command("grid")
@click.argument("workdir", type=click.Path(exists=True, file_okay=False,
                path_type=Path), metavar="WORKDIR")
@_solver_option
@_iteration_option
@_dpi_option
@_save_option
@_show_option
@verbose_option
@no_color_option
@click.pass_context
def plot_grid(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    iteration: int | None,
    dpi: int,
    save_path: Path | None,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot a compact response grid for all stations with per-site RMS (Occam2D).

    \b
    Examples:
      pycsamt invert plot grid run01/ --show
      pycsamt invert plot grid run01/ --save response_grid.png
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    if solver_name != "occam2d":
        raise click.UsageError(
            "'grid' is an Occam2D-only sub-command."
        )
    result = _load_inversion_result(workdir, solver_name, iteration, verbose)

    try:
        from pycsamt.models.occam2d.plot import PlotResponseGrid  # noqa: PLC0415
        fig = PlotResponseGrid(result, dpi=dpi).plot()
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    _handle_figure(fig, save_path, show)
