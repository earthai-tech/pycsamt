# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.format._base
=================================

Root Click group and thin Click-facing wrappers for ``pycsamt format`` —
the CLI face of :mod:`pycsamt.format` (PCSF / PCSM inversion-result
format).

The actual conversion logic (detect / resolve-target / build-model /
write-model / report) lives in :mod:`pycsamt.format.convert_engine`,
which is UI-neutral and also used directly by the standalone converter
GUI (:mod:`pycsamt.app.converter`). This module only adds the Click
plumbing: translating the engine's plain ``ValueError``/``FileNotFoundError``
into ``click.ClickException``/``click.UsageError``, and a rich-table
console helper.

Shared helpers
--------------
_rich_table(title, rows)          Two-column key/value table (rich → plain).
_detect(path, solver)             Click-facing wrap of ``convert_engine.detect``.
_resolve_target(...)              Click-facing wrap of ``convert_engine.resolve_target``.
_build_model(source_kind, opts)   Click-facing wrap of ``convert_engine.build_model``.
_write_model(model, dst, fmt)     ``convert_engine.write_model`` passthrough.
_model_report(model, path)        ``convert_engine.model_report`` passthrough.
"""

from __future__ import annotations

from pathlib import Path

import click

from pycsamt.format import convert_engine as _engine

# Re-exported so existing imports (``from ._base import ConvertOptions``)
# keep working unchanged.
ConvertOptions = _engine.ConvertOptions

# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group("format")
@click.pass_context
def fmt(ctx: click.Context) -> None:
    """PCSF / PCSM inversion-result format — convert, inspect, validate.

    PCSF (``.pcsf``, HDF5) is pyCSAMT's backend-neutral container for a
    finished resistivity model; PCSM (``.pcsm``) is its lossless,
    hand-editable ASCII sibling. Every classical solver (Occam2D,
    ModEM 3-D, MARE2DEM) and any AI/DL inversion converts *to* PCSF.

    \b
    Sub-commands:
      convert    Any solver result / AI array bundle / .pcsf / .pcsm
                 → .pcsf or .pcsm (auto-detects the source).
      detect     Report what a file or folder is, without converting.
      info       Full summary of a .pcsf / .pcsm file.
      validate   Structural + round-trip check of a .pcsf / .pcsm file.

    \b
    Examples:
      pycsamt format convert data/occam2D            model.pcsf
      pycsamt format convert data/modem/run01/       -o out/ --to pcsm
      pycsamt format convert unet_prediction.npz     unet.pcsf
      pycsamt format convert model.pcsf              model.pcsm
      pycsamt format detect  data/mare2dem/demo_mt_inversion
      pycsamt format info    model.pcsf
    """
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# rich helper (degrades to plain text)
# ---------------------------------------------------------------------------


def _rich_table(
    title: str,
    rows: list[tuple[str, str]],
    style: str = "cyan",
) -> None:
    """Print a two-column key/value table via rich (plain fallback)."""
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table  # noqa: PLC0415

        console = Console()
        table = Table(title=title, border_style=style, show_header=False)
        table.add_column("Key", style="bold")
        table.add_column("Value", style="white")
        for key, value in rows:
            table.add_row(str(key), str(value))
        console.print(table)
    except ImportError:
        click.echo(f"\n{title}")
        for key, value in rows:
            click.echo(f"  {key:<26} {value}")


# ---------------------------------------------------------------------------
# Click-facing wrappers around pycsamt.format.convert_engine
# ---------------------------------------------------------------------------


def _detect(path: Path, solver: str | None):
    """Call :func:`pycsamt.format.convert_engine.detect`, raising ClickException."""
    try:
        return _engine.detect(path, solver)
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc


def _resolve_target(
    source: Path,
    target: Path | None,
    to_format: str | None,
    output_dir: Path,
) -> tuple[Path, str]:
    try:
        return _engine.resolve_target(source, target, to_format, output_dir)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


def _build_model(source_kind, opts):
    try:
        return _engine.build_model(source_kind, opts)
    except (ValueError, ImportError) as exc:
        raise click.ClickException(str(exc)) from exc


def _write_model(model, dst: Path, canonical_format: str, log10_view: bool):
    return _engine.write_model(model, dst, canonical_format, log10_view)


def _model_report(model, path: Path):
    return _engine.model_report(model, path)
