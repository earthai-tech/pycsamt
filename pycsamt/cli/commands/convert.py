# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt convert
===============

Convert electromagnetic transfer-function files.

The preferred Phase-9 interface converts EDI and EMTF XML through the
format-neutral :mod:`pycsamt.emtf` model::

    pycsamt convert station.edi station.xml
    pycsamt convert station.xml station.edi

Historical Jones/AVG/EDI-to-EDI batch workflows remain supported.
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
            "message": (
                "zonge package unavailable — AVG conversion not "
                "supported yet"
            ),
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
# ---------------------------------------------------------------------------
# Format-neutral transfer-function conversion helpers
# ---------------------------------------------------------------------------


def _canonical_tf_format(name: str | None) -> str | None:
    """Return a canonical registered TF format name."""
    if name is None:
        return None
    from ...io import get_tf_format

    return get_tf_format(name).name


def _tf_extension(name: str) -> str:
    """Return the preferred filename extension for a TF format."""
    from ...io import get_tf_format

    spec = get_tf_format(name)
    if not spec.extensions:
        raise click.ClickException(
            f"No output extension is registered for {spec.name!r}."
        )
    return spec.extensions[0]


def _read_as_emtf(src: Path, source_format: str | None):
    """Read any supported TF source into the neutral EMTF model."""
    from ...emtf import EMTF
    from ...io import detect_tf_format

    canonical = (
        _canonical_tf_format(source_format)
        if source_format is not None
        else detect_tf_format(src)
    )
    if canonical == "edi":
        return EMTF.from_edi(src), canonical
    if canonical == "emtf_xml":
        return EMTF.from_xml(src), canonical
    raise click.ClickException(
        f"Format {canonical!r} cannot be converted through EMTF."
    )


def _convert_tf(
    src: Path,
    dst: Path,
    *,
    source_format: str | None,
    target_format: str | None,
    on_loss: str,
    verbose: int,
) -> dict[str, Any]:
    """Convert EDI/EMTF XML through the format-neutral EMTF model."""
    try:
        from ...io import get_tf_format_for_target, write_transfer_function

        doc, src_fmt = _read_as_emtf(src, source_format)
        if target_format is None:
            dst_fmt = get_tf_format_for_target(dst).name
        else:
            dst_fmt = _canonical_tf_format(target_format)
            assert dst_fmt is not None

        kwargs: dict[str, Any] = {}
        if dst_fmt == "edi":
            kwargs["on_loss"] = on_loss
        write_transfer_function(doc, dst, format=dst_fmt, **kwargs)
        if verbose >= 1:
            click.echo(
                f"  {src.name}  →  {dst.name}  "
                f"({src_fmt} → {dst_fmt})",
                err=True,
            )
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "ok",
            "source_format": src_fmt,
            "target_format": dst_fmt,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "src": str(src),
            "dst": str(dst),
            "status": "error",
            "message": str(exc),
        }


def _generic_tf_mode(
    source: Path,
    target: Path | None,
    source_format: str | None,
    target_format: str | None,
) -> bool:
    """Return whether the new format-neutral conversion path is requested.

    ``--to edi`` alone must stay indistinguishable from the historical
    default (EDI was, and remains, the only legacy output format), so a
    directory SOURCE combined only with ``--to edi`` still runs the
    legacy Jones/AVG/EDI batch path. Only an explicit TARGET, an
    explicit ``--from``, or a target format other than EDI opts a
    directory into the single-file transfer-function path.
    """
    if target is not None or source_format is not None:
        return True
    if target_format is not None:
        try:
            canonical = _canonical_tf_format(target_format)
        except Exception:  # noqa: BLE001
            return True
        if canonical != "edi":
            return True
    return source.is_file() and source.suffix.lower() == ".xml"


def _generic_target(
    source: Path,
    target: Path | None,
    output_dir: Path,
    target_format: str | None,
) -> tuple[Path, str]:
    """Resolve output path and target format for one TF conversion."""
    from ...io import get_tf_format_for_target

    if target is not None:
        canonical = (
            _canonical_tf_format(target_format)
            if target_format is not None
            else get_tf_format_for_target(target).name
        )
        assert canonical is not None
        return target, canonical

    if target_format is None:
        # A lone XML source naturally converts to the historical EDI format.
        # A lone EDI source remains on the legacy pass-through path.
        canonical = "edi"
    else:
        canonical = _canonical_tf_format(target_format)
        assert canonical is not None
    return output_dir / f"{source.stem}{_tf_extension(canonical)}", canonical


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------


@click.command("convert")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="SOURCE",
)
@click.argument(
    "target",
    required=False,
    type=click.Path(path_type=Path),
    metavar="[TARGET]",
)
@click.option(
    "--from",
    "source_format",
    default=None,
    metavar="FORMAT",
    help="Explicit source TF format (edi, emtf-xml, xml).",
)
@click.option(
    "--to",
    "target_format",
    default=None,
    metavar="FORMAT",
    help=(
        "Target TF format (edi or emtf-xml). With no TARGET, the output "
        "filename is derived in --output-dir."
    ),
)
@click.option(
    "--on-loss",
    type=click.Choice(["warn", "raise", "ignore"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Policy for EMTF → EDI information that EDI cannot preserve.",
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
    target: Path | None,
    source_format: str | None,
    target_format: str | None,
    on_loss: str,
    dry_run: bool,
    verbose: int,
    no_color: bool,
    output_format: str,
    output_dir: Path,
    overwrite: bool,
) -> None:
    """Convert electromagnetic transfer-function files.

    The preferred modern form is ``pycsamt convert SOURCE TARGET``. EDI and
    EMTF XML are converted through the format-neutral :class:`EMTF` model.
    Existing Jones/AVG-to-EDI directory workflows remain supported.

    \b
    Modern transfer-function examples:
      pycsamt convert station.edi station.xml
      pycsamt convert station.xml station.edi
      pycsamt convert station.edi --to emtf-xml -o converted/
      pycsamt convert input.dat output.xml --from edi --to emtf-xml

    \b
    Legacy input formats (unchanged):
      .j    Jones J-file
      .avg  Zonge AVG file
      .edi  EDI pass-through / re-export

    \b
    Legacy conversion examples:
      pycsamt convert station.j --output-dir edis/
      pycsamt convert survey/ --output-dir edis/ --overwrite
      pycsamt convert survey/ --dry-run -v
    """
    configure_cli(
        log__level=verbose,
        log__color=not no_color,
        output__format=output_format,
        output__dir=output_dir,
        output__overwrite=overwrite,
    )

    # ------------------------------------------------------------------
    # Modern EDI <-> EMTF XML conversion. This path is intentionally
    # single-file so conversion semantics and loss policy stay explicit.
    # ------------------------------------------------------------------
    if _generic_tf_mode(source, target, source_format, target_format):
        if source.is_dir():
            raise click.ClickException(
                "SOURCE must be a file when TARGET, --from, or --to is used."
            )
        try:
            dst, canonical_target = _generic_target(
                source,
                target,
                output_dir,
                target_format,
            )
        except Exception as exc:  # noqa: BLE001
            raise click.ClickException(str(exc)) from exc

        if source.resolve() == dst.resolve():
            raise click.ClickException(
                "SOURCE and TARGET must be different files."
            )
        if dst.exists() and not overwrite:
            result = {
                "src": str(source),
                "dst": str(dst),
                "status": "skipped",
                "message": "output exists (use --overwrite to replace)",
            }
            if output_format == "json":
                click.echo(_fmt_json([result]))
            elif output_format == "csv":
                click.echo(_fmt_csv([result]))
            else:
                click.echo(_fmt_text([result]))
            return

        if dry_run:
            click.echo(
                f"Dry run — {source}  →  {dst} ({canonical_target})"
            )
            return

        dst.parent.mkdir(parents=True, exist_ok=True)
        result = _convert_tf(
            source,
            dst,
            source_format=source_format,
            target_format=canonical_target,
            on_loss=on_loss.lower(),
            verbose=verbose,
        )
        if output_format == "json":
            click.echo(_fmt_json([result]))
        elif output_format == "csv":
            click.echo(_fmt_csv([result]))
        else:
            click.echo(_fmt_text([result]))
        if result["status"] == "error":
            raise click.ClickException(
                result.get("message", "conversion failed")
            )
        return

    # ------------------------------------------------------------------
    # Historical Jones / AVG / EDI-to-EDI batch behavior.
    # ------------------------------------------------------------------
    inputs = _collect_inputs(source)
    if not inputs:
        click.echo(
            f"No convertible files found in {source}.\n"
            f"Supported legacy extensions: "
            f"{', '.join(sorted(_SUPPORTED_EXTS))}",
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
