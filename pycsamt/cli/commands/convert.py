# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt convert
===============

Convert AMT/MT data files to EDI format.

Supported input formats
-----------------------

==========  =======================================================
Extension   Source format
==========  =======================================================
``.j``      Jones J-file (via ``pycsamt.jones``)
``.avg``    Zonge AVG file (via ``pycsamt.zonge``)
``.edi``    EDI pass-through / re-export / header normalisation
==========  =======================================================

Usage
-----
::

    # single file
    pycsamt convert station.j --output-dir edis/

    # whole directory — converts every recognised file
    pycsamt convert survey/ --output-dir edis/

    # force overwrite of existing outputs
    pycsamt convert survey/ --output-dir edis/ --overwrite

    # dry run — list what would be converted
    pycsamt convert survey/ --dry-run -v
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import click

from ...api.cli.config import configure_cli
from ...api.cli.options import (
    format_option,
    no_color_option,
    output_dir_option,
    overwrite_option,
    verbose_option,
)

# ---------------------------------------------------------------------------
# Supported source extensions
# ---------------------------------------------------------------------------

_SUPPORTED_EXTS = {".j", ".avg", ".edi"}


def _collect_inputs(source: Path) -> list[Path]:
    """Return sorted list of convertible files from a file or directory."""
    if source.is_dir():
        return sorted(
            p
            for p in source.iterdir()
            if p.is_file() and p.suffix.lower() in _SUPPORTED_EXTS
        )
    if source.suffix.lower() in _SUPPORTED_EXTS:
        return [source]
    return []


# ---------------------------------------------------------------------------
# Per-format converters
# ---------------------------------------------------------------------------


def _convert_j(src: Path, dst: Path, verbose: int) -> dict[str, Any]:
    """Convert a Jones J-file to EDI."""
    try:
        from pycsamt.jones import JFile  # noqa: PLC0415
        from pycsamt.jones.xa import (
            _meta_from_jfile,  # noqa: PLC0415
        )
    except ImportError as exc:
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "error",
            "message": f"jones package unavailable: {exc}",
        }

    try:
        jf = JFile(src)
        meta = _meta_from_jfile(jf)
        station = meta.get("station_id", src.stem)
        edi_text = _jfile_to_edi_text(jf, meta)
        dst.write_text(edi_text, encoding="utf-8")
        if verbose >= 1:
            click.echo(f"  {src.name}  →  {dst.name}", err=True)
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "ok",
            "station": station,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "error",
            "message": str(exc),
        }


def _jfile_to_edi_text(jf: Any, meta: dict[str, Any]) -> str:
    """Render a minimal valid EDI from a JFile object and its metadata."""
    import numpy as np  # noqa: PLC0415

    station = str(meta.get("station_id", "UNKNOWN"))
    lat = float(meta.get("lat", 0.0))
    lon = float(meta.get("lon", 0.0))
    elev = float(meta.get("elev", 0.0))

    # Extract frequency array from the JFile blocks
    freqs: list[float] = []
    zxx_r = zxx_i = zxy_r = zxy_i = zyx_r = zyx_i = zyy_r = zyy_i = None
    try:
        rblock = jf.rblock
        freqs = list(np.asarray(rblock.freq, dtype=float))
        zxx_r = list(np.asarray(rblock.zxx_r, dtype=float))
        zxx_i = list(np.asarray(rblock.zxx_i, dtype=float))
        zxy_r = list(np.asarray(rblock.zxy_r, dtype=float))
        zxy_i = list(np.asarray(rblock.zxy_i, dtype=float))
        zyx_r = list(np.asarray(rblock.zyx_r, dtype=float))
        zyx_i = list(np.asarray(rblock.zyx_i, dtype=float))
        zyy_r = list(np.asarray(rblock.zyy_r, dtype=float))
        zyy_i = list(np.asarray(rblock.zyy_i, dtype=float))
    except AttributeError:
        pass

    nfreq = len(freqs)

    def _block(tag: str, values: list[float] | None) -> str:
        if values is None:
            return ""
        nums = "  ".join(f"{v: .6E}" for v in values)
        return f">{tag}  // {nfreq}\n{nums}\n"

    lines = [
        ">HEAD",
        f"  DATAID={station!r}",
        f"  LAT={lat:.6f}",
        f"  LONG={lon:.6f}",
        f"  ELEV={elev:.2f}",
        f"  NFREQ={nfreq}",
        ">END_HEAD",
        "",
        ">INFO",
        "  Converted from Jones J-file by pyCSAMT",
        ">END_INFO",
        "",
        _block("FREQ", freqs),
        _block("ZXXR", zxx_r),
        _block("ZXXI", zxx_i),
        _block("ZXYR", zxy_r),
        _block("ZXYI", zxy_i),
        _block("ZYXR", zyx_r),
        _block("ZYXI", zyx_i),
        _block("ZYYR", zyy_r),
        _block("ZYYI", zyy_i),
        ">END",
    ]
    return "\n".join(lines)


def _convert_avg(src: Path, dst: Path, verbose: int) -> dict[str, Any]:
    """Convert a Zonge AVG file to EDI (stub — wired to zonge package)."""
    try:
        from pycsamt.zonge import AVGFile  # noqa: PLC0415

        avg = AVGFile(src)
        avg.to_edi(dst)
        if verbose >= 1:
            click.echo(f"  {src.name}  →  {dst.name}", err=True)
        return {"src": str(src), "dst": str(dst), "status": "ok"}
    except ImportError:
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "error",
            "message": "zonge package unavailable — AVG conversion not supported yet",
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "error",
            "message": str(exc),
        }


def _convert_edi(src: Path, dst: Path, verbose: int) -> dict[str, Any]:
    """Pass-through / re-export an existing EDI file."""
    import shutil  # noqa: PLC0415

    try:
        shutil.copy2(src, dst)
        if verbose >= 1:
            click.echo(f"  {src.name}  →  {dst.name}  (copy)", err=True)
        return {"src": str(src), "dst": str(dst), "status": "ok"}
    except Exception as exc:  # noqa: BLE001
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "error",
            "message": str(exc),
        }


_CONVERTERS = {
    ".j": _convert_j,
    ".avg": _convert_avg,
    ".edi": _convert_edi,
}


# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------


def _fmt_text(results: list[dict[str, Any]]) -> str:
    ok = [r for r in results if r["status"] == "ok"]
    errors = [r for r in results if r["status"] == "error"]
    lines = [
        f"Converted : {len(ok)} file(s)",
        f"Errors    : {len(errors)} file(s)",
    ]
    if errors:
        lines.append("")
        lines.append("Errors:")
        for r in errors:
            lines.append(f"  {r['src']}: {r['message']}")
    return "\n".join(lines)


def _fmt_json(results: list[dict[str, Any]]) -> str:
    return json.dumps(results, indent=2, default=str)


def _fmt_csv(results: list[dict[str, Any]]) -> str:
    if not results:
        return ""
    header = ",".join(results[0].keys())
    rows = [header]
    for r in results:
        rows.append(",".join(str(v) for v in r.values()))
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------


@click.command("convert")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="FILE_OR_DIR",
)
@click.option(
    "--to",
    "target_format",
    type=click.Choice(["edi"], case_sensitive=False),
    default="edi",
    show_default=True,
    help="Target format (currently only EDI is supported).",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="List files that would be converted without writing anything.",
)
@verbose_option
@no_color_option
@format_option
@output_dir_option
@overwrite_option
@click.pass_context
def convert(
    ctx: click.Context,
    source: Path,
    target_format: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
    output_format: str,
    output_dir: Path,
    overwrite: bool,
) -> None:
    """Convert AMT/MT data files to EDI format.

    SOURCE can be a single file (.j, .avg, .edi) or a directory.  When
    SOURCE is a directory every recognised file inside it is converted.

    \b
    Supported input extensions:
      .j    Jones J-file
      .avg  Zonge AVG file
      .edi  EDI pass-through / re-export

    \b
    Examples:
      pycsamt convert station.j --output-dir edis/
      pycsamt convert survey/   --output-dir edis/ --overwrite
      pycsamt convert survey/   --dry-run -v
      pycsamt convert survey/   --output-dir edis/ --format json
    """
    configure_cli(
        log__level=verbose,
        log__color=not no_color,
        output__format=output_format,
        output__dir=output_dir,
        output__overwrite=overwrite,
    )

    inputs = _collect_inputs(source)
    if not inputs:
        click.echo(
            f"No convertible files found in {source}.\n"
            f"Supported extensions: {', '.join(sorted(_SUPPORTED_EXTS))}",
            err=True,
        )
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    if dry_run:
        click.echo(f"Dry run — {len(inputs)} file(s) would be converted:\n")
        for p in inputs:
            dst = output_dir / (p.stem + ".edi")
            click.echo(f"  {p}  →  {dst}")
        return

    if verbose >= 1:
        click.echo(
            f"Converting {len(inputs)} file(s) → {output_dir}/", err=True
        )

    results: list[dict[str, Any]] = []
    for src in inputs:
        dst = output_dir / (src.stem + ".edi")
        if dst.exists() and not overwrite:
            results.append(
                {
                    "src": str(src),
                    "dst": str(dst),
                    "status": "skipped",
                    "message": "output exists (use --overwrite to replace)",
                }
            )
            continue

        converter = _CONVERTERS.get(src.suffix.lower())
        if converter is None:
            results.append(
                {
                    "src": str(src),
                    "dst": str(dst),
                    "status": "error",
                    "message": f"no converter for extension {src.suffix!r}",
                }
            )
            continue

        results.append(converter(src, dst, verbose))

    if output_format == "json":
        click.echo(_fmt_json(results))
    elif output_format == "csv":
        click.echo(_fmt_csv(results))
    else:
        click.echo(_fmt_text(results))

    errors = [r for r in results if r["status"] == "error"]
    if errors:
        sys.exit(1)
