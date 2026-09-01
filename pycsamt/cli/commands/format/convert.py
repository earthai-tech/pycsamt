# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format convert
======================

Smart, one-shot conversion of *any* inversion result into PCSF or PCSM.

The SOURCE is auto-detected (see ``pycsamt format detect``):

* an Occam2D / ModEM 3-D / MARE2DEM working directory (or a single
  signature file inside one) → the matching
  :mod:`pycsamt.format.adapters` backend converter;
* a ``.npz`` / ``.npy`` array bundle from an AI/DL inversion → the
  solver-agnostic :mod:`pycsamt.format.adapters.generic` path;
* an existing ``.pcsf`` / ``.pcsm`` file → transcoded to the other
  encoding losslessly.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    output_dir_option,
    overwrite_option,
    verbose_option,
)
from ._base import (
    ConvertOptions,
    _build_model,
    _detect,
    _model_report,
    _resolve_target,
    _rich_table,
    _write_model,
    fmt,
)


def _parse_origin(text: str | None) -> tuple[float, ...] | None:
    if not text:
        return None
    try:
        parts = tuple(float(v) for v in text.replace(";", ",").split(","))
    except ValueError as exc:
        raise click.UsageError(
            f"--origin must be comma-separated numbers: {exc}"
        )
    if len(parts) not in (2, 3):
        raise click.UsageError("--origin needs 2 (x,y) or 3 (x,y,z) values.")
    return parts


@fmt.command("convert")
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
    "--to",
    "to_format",
    type=click.Choice(["pcsf", "pcsm", "pcsm.gz"], case_sensitive=False),
    default=None,
    help="Output format (default: pcsf, or inferred from TARGET's suffix).",
)
@click.option(
    "--solver",
    type=click.Choice(["occam2d", "modem", "mare2dem"], case_sensitive=False),
    default=None,
    help="Force the source solver backend instead of fingerprinting it.",
)
@click.option(
    "--iteration",
    type=int,
    default=None,
    help="Occam2D iteration index to convert (default: the final one).",
)
@click.option(
    "--topo",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    metavar="PATH",
    help="Topography source (.bln/.csv/.stn/EDI dir) attached per station.",
)
@click.option(
    "--epsg",
    type=int,
    default=None,
    help="EPSG code for --topo / lon-lat.",
)
@click.option(
    "--utm-zone",
    default=None,
    metavar="ZONE",
    help="UTM zone (e.g. '48N') for --topo / lon-lat projection.",
)
@click.option(
    "--encoding",
    type=click.Choice(["linear", "log10", "ln"], case_sensitive=False),
    default=None,
    help="Encoding of the AI resistivity array (default: auto/linear).",
)
@click.option(
    "--station-z",
    "station_z_convention",
    type=click.Choice(
        ["auto", "elevation", "depth_down"], case_sensitive=False
    ),
    default="auto",
    show_default=True,
    help=(
        "ModEM: how to read the .dat Z column. 'auto' flips a "
        "positive-down depth to an elevation only when the file's "
        "comments say so; 'depth_down' forces it; 'elevation' trusts it."
    ),
)
@click.option(
    "--air-threshold",
    "air_threshold_ohm_m",
    type=float,
    default=1e8,
    show_default=True,
    help=(
        "ModEM: mask cells above this resistivity (ohm.m) as "
        "above-topography air fill. Pass 0 to disable."
    ),
)
@click.option(
    "--poly",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    metavar="PATH",
    help="MARE2DEM .poly PSLG to rebuild the mesh from (default: auto).",
)
@click.option(
    "--origin",
    default=None,
    metavar="X,Y[,Z]",
    help="Real-world grid origin for AI grid2d/grid3d sources.",
)
@click.option(
    "--azimuth",
    "azimuth_deg",
    type=float,
    default=None,
    help="Profile azimuth (grid2d) or volume rotation (grid3d), degrees.",
)
@click.option(
    "--created-by",
    default="pycsamt format convert",
    show_default=True,
    help="Value stored in the PCSF 'created_by' attribute.",
)
@click.option("--description", default="", help="Free-text PCSF description.")
@click.option(
    "--log10-view",
    is_flag=True,
    default=False,
    help="For --to pcsm: write the human-readable block in log10(rho).",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Detect and report the plan without writing anything.",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
    help="Console output format.",
)
@verbose_option
@no_color_option
@output_dir_option
@overwrite_option
@click.pass_context
def convert(
    ctx: click.Context,
    source: Path,
    target: Path | None,
    to_format: str | None,
    solver: str | None,
    iteration: int | None,
    topo: Path | None,
    epsg: int | None,
    utm_zone: str | None,
    encoding: str | None,
    station_z_convention: str,
    air_threshold_ohm_m: float,
    poly: Path | None,
    origin: str | None,
    azimuth_deg: float | None,
    created_by: str,
    description: str,
    log10_view: bool,
    dry_run: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
    output_dir: Path,
    overwrite: bool,
) -> None:
    """Convert SOURCE to PCSF/PCSM.  SOURCE is auto-detected.

    \b
    Examples:
      pycsamt format convert data/occam2D               occam.pcsf
      pycsamt format convert data/modem/run01/ -o out/  --to pcsm
      pycsamt format convert data/mare2dem/demo_mt_inversion  mare.pcsf
      pycsamt format convert unet_prediction.npz  unet.pcsf --encoding log10
      pycsamt format convert occam.pcsf                 occam.pcsm --log10-view
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    to_format = to_format.lower() if to_format else None

    sk = _detect(source, solver)
    if not sk.convertible:
        raise click.ClickException(
            f"Nothing to convert — {source} looks like: {sk.detail}"
        )

    dst, canonical = _resolve_target(source, target, to_format, output_dir)
    same_format_transcode = (
        sk.category in {"pcsf", "pcsm"}
        and canonical.split(".")[0] == sk.category
    )
    if same_format_transcode and (
        target is None or source.resolve() == dst.resolve()
    ):
        raise click.ClickException(
            f"{source.name} is already {sk.category.upper()}; give a "
            "different TARGET or --to the other format."
        )

    plan = {
        "source": str(source),
        "detected": sk.to_dict(),
        "target": str(dst),
        "target_format": canonical,
    }

    if dry_run:
        if output_format == "json":
            click.echo(json.dumps(plan, indent=2, default=str))
        else:
            _rich_table(
                "format convert — dry run",
                [
                    ("source", str(source)),
                    ("category", sk.category),
                    ("backend", sk.backend or "-"),
                    ("geometry", sk.geometry or sk.target_geometry or "-"),
                    ("confidence", sk.confidence),
                    ("detail", sk.detail),
                    ("→ target", str(dst)),
                    ("→ format", canonical),
                ],
            )
        return

    if dst.exists() and not overwrite:
        raise click.ClickException(
            f"{dst} exists — pass --overwrite to replace it."
        )

    # ------------------------------------------------------------------
    # PCSF <-> PCSM direct transcode (no array-level rebuild)
    # ------------------------------------------------------------------
    if sk.category in {"pcsf", "pcsm"}:
        written = _transcode(source, dst, sk.category, canonical, log10_view)
        _emit_result(written, output_format, verbose, transcode=True)
        return

    # ------------------------------------------------------------------
    # Solver / AI → PCSFModel → file
    # ------------------------------------------------------------------
    opts = ConvertOptions(
        iteration=iteration,
        topo=topo,
        epsg=epsg,
        utm_zone=utm_zone,
        encoding=encoding.lower() if encoding else None,
        station_z_convention=station_z_convention.lower(),
        air_threshold_ohm_m=(air_threshold_ohm_m or None),
        poly=poly,
        origin=_parse_origin(origin),
        azimuth_deg=azimuth_deg,
        created_by=created_by,
        description=description,
        verbose=verbose,
    )

    try:
        model = _build_model(sk, opts)
    except click.ClickException:
        raise
    except Exception as exc:  # noqa: BLE001 - surface adapter errors cleanly
        raise click.ClickException(
            f"{sk.backend or sk.category} → PCSF failed: {exc}"
        ) from exc

    written = _write_model(model, dst, canonical, log10_view)
    report = _model_report(model, Path(written))
    _emit_result(report, output_format, verbose)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _transcode(
    source: Path,
    dst: Path,
    src_category: str,
    canonical: str,
    log10_view: bool,
) -> dict:
    dst.parent.mkdir(parents=True, exist_ok=True)
    target_base = canonical.split(".")[0]

    if src_category == "pcsf" and target_base == "pcsm":
        from pycsamt.format.text import pcsf_to_pcsm

        pcsf_to_pcsm(source, dst, log10_view=log10_view)
    elif src_category == "pcsm" and target_base == "pcsf":
        from pycsamt.format.text import pcsm_to_pcsf

        pcsm_to_pcsf(source, dst)
    else:
        # pcsf->pcsf or pcsm->pcsm normalisation via a full model round-trip
        from pycsamt.format.text import read_pcsf_or_pcsm

        model = read_pcsf_or_pcsm(source)
        _write_model(model, dst, target_base, log10_view)

    from pycsamt.format.text import read_pcsf_or_pcsm

    model = read_pcsf_or_pcsm(dst)
    return _model_report(model, dst)


def _emit_result(
    report: dict,
    output_format: str,
    verbose: int,
    *,
    transcode: bool = False,
) -> None:
    if output_format == "json":
        click.echo(json.dumps(report, indent=2, default=str))
        return

    rows = [
        ("wrote", report["file"]),
        ("size", f"{report.get('size_bytes') or 0:,} bytes"),
        ("kind", report.get("kind", "-")),
        ("backend", report.get("source_backend", "-")),
        ("resistivity", str(report.get("resistivity_shape"))),
        ("native encoding", report.get("native_encoding") or "linear"),
        ("stations", str(report.get("n_stations", 0))),
        ("topography", "yes" if report.get("has_topography") else "no"),
    ]
    rho = report.get("rho_ohm_m")
    if rho:
        rng = (
            f"{rho['min']:.4g} … {rho['max']:.4g} "
            f"(median {rho['median']:.4g})"
        )
        rows.append(("rho range (Ω·m)", rng))
    title = (
        "format convert — transcoded"
        if transcode
        else "format convert — done"
    )
    _rich_table(title, rows)
    click.echo(f"\n✓ {report['file']}")
