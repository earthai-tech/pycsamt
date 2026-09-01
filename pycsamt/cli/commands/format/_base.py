# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.format._base
=================================

Root Click group and shared helpers for ``pycsamt format`` — the CLI
face of :mod:`pycsamt.format` (PCSF / PCSM inversion-result format).

Shared helpers
--------------
_rich_table(title, rows)          Two-column key/value table (rich → plain).
_detect(path, solver)             Wrap :func:`pycsamt.format.detect_source`.
_resolve_target(...)              Work out the output path + format.
_build_model(source_kind, opts)   Any source → :class:`PCSFModel`.
_write_model(model, dst, fmt)     PCSFModel → ``.pcsf`` / ``.pcsm`` file.
_model_report(model, path)        JSON-friendly summary of a written file.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import click

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
# Detection wrapper
# ---------------------------------------------------------------------------


def _detect(path: Path, solver: str | None):
    """Call :func:`pycsamt.format.detect_source`, raising ClickException."""
    from pycsamt.format.detect import detect_source

    try:
        return detect_source(path, solver_hint=solver)
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc


# ---------------------------------------------------------------------------
# Conversion options bundle
# ---------------------------------------------------------------------------


@dataclass
class ConvertOptions:
    """Everything ``_build_model`` needs beyond the detected source."""

    iteration: int | None = None
    topo: Path | None = None
    epsg: int | None = None
    utm_zone: str | None = None
    encoding: str | None = None
    poly: Path | None = None
    origin: tuple[float, ...] | None = None
    azimuth_deg: float | None = None
    created_by: str = "pycsamt format convert"
    description: str = ""
    verbose: int = 0
    # ModEM 3-D only.
    station_z_convention: str = "auto"
    air_threshold_ohm_m: float | None = 1e8


# ---------------------------------------------------------------------------
# Target-path resolution
# ---------------------------------------------------------------------------

_FMT_SUFFIX = {"pcsf": ".pcsf", "pcsm": ".pcsm", "pcsm.gz": ".pcsm.gz"}


def _source_stem(path: Path) -> str:
    """Strip a (possibly double) PCSF/array suffix from *path*'s name."""
    name = path.name
    for suffix in (".pcsm.gz", ".pcsf", ".pcsm", ".npz", ".npy"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem if path.is_file() else path.name


def _infer_format_from_name(name: str) -> str | None:
    low = name.lower()
    if low.endswith(".pcsm.gz"):
        return "pcsm.gz"
    if low.endswith(".pcsm"):
        return "pcsm"
    if low.endswith(".pcsf"):
        return "pcsf"
    return None


def _resolve_target(
    source: Path,
    target: Path | None,
    to_format: str | None,
    output_dir: Path,
) -> tuple[Path, str]:
    """Return ``(output_path, canonical_format)``.

    ``canonical_format`` is one of ``"pcsf"``, ``"pcsm"``, ``"pcsm.gz"``.
    """
    if target is not None:
        inferred = _infer_format_from_name(target.name)
        if to_format is not None:
            canonical = to_format
            base_inferred = inferred.split(".")[0] if inferred else None
            if base_inferred is not None and base_inferred != canonical.split(
                "."
            )[0]:
                raise click.UsageError(
                    f"TARGET name {target.name!r} and --to {to_format!r} "
                    "disagree."
                )
        elif inferred is not None:
            canonical = inferred
        else:
            raise click.UsageError(
                f"Cannot tell the output format from {target.name!r}; add "
                "a .pcsf/.pcsm extension or pass --to."
            )
        return target, canonical

    canonical = to_format or "pcsf"
    stem = _source_stem(source)
    return output_dir / f"{stem}{_FMT_SUFFIX[canonical]}", canonical


# ---------------------------------------------------------------------------
# Source → PCSFModel
# ---------------------------------------------------------------------------


def _build_model(source_kind, opts: ConvertOptions):
    """Dispatch on ``source_kind.category`` and return a PCSFModel.

    ``pcsf`` / ``pcsm`` sources return ``None`` — those are handled by a
    direct file-level transcode in ``convert`` (no array round-trip).
    """
    category = source_kind.category
    if category in {"pcsf", "pcsm"}:
        return None
    if category == "solver":
        return _build_from_solver(source_kind, opts)
    if category == "ai_arrays":
        return _build_from_ai_arrays(source_kind, opts)
    raise click.ClickException(
        f"Don't know how to convert a {category!r} source: "
        f"{source_kind.detail}"
    )


def _build_from_solver(sk, opts: ConvertOptions):
    backend = sk.backend
    workdir = sk.path if sk.is_dir else sk.path.parent

    common = {
        "created_by": opts.created_by,
        "description": opts.description or f"{backend} → PCSF via pycsamt CLI",
    }

    if backend == "occam2d":
        from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
        from pycsamt.models.occam2d.results import InversionResult

        result = InversionResult(
            workdir=workdir,
            iteration=opts.iteration,
            verbose=opts.verbose,
        )
        return occam2d_to_pcsf(
            result,
            topo=str(opts.topo) if opts.topo else None,
            epsg=opts.epsg,
            utm_zone=opts.utm_zone,
            origin=opts.origin,
            azimuth_deg=opts.azimuth_deg,
            **common,
        )

    if backend == "modem":
        from pycsamt.format.adapters.modem3d import modem3d_to_pcsf
        from pycsamt.models.modem.results import InversionResult

        result = InversionResult(workdir=workdir, load_data=True)
        return modem3d_to_pcsf(
            result,
            station_z_convention=opts.station_z_convention,
            air_threshold_ohm_m=opts.air_threshold_ohm_m,
            topo=str(opts.topo) if opts.topo else None,
            epsg=opts.epsg,
            utm_zone=opts.utm_zone,
            **common,
        )

    if backend == "mare2dem":
        from pycsamt.format.adapters.mare2dem import mare2dem_to_pcsf
        from pycsamt.models.mare2dem.results import InversionResult

        result = InversionResult(workdir=workdir)
        mesh = _mare2dem_mesh(_resolve_poly(sk, workdir, opts))
        return mare2dem_to_pcsf(result, mesh, **common)

    raise click.ClickException(f"Unsupported solver backend: {backend!r}")


def _resolve_poly(sk, workdir: Path, opts: ConvertOptions) -> Path:
    if opts.poly is not None:
        return opts.poly
    hint = sk.hints.get("poly")
    if hint is not None:
        return Path(hint)
    polys = sorted(workdir.glob("*.poly"))
    if not polys:
        raise click.ClickException(
            f"No .poly PSLG found in {workdir}. MARE2DEM → PCSF needs the "
            "run's polygon mesh file; pass --poly PATH."
        )
    return polys[0]


def _mare2dem_mesh(poly_path: Path):
    """Rebuild the run's TriMesh from its ``.poly`` PSLG (needs *triangle*)."""
    try:
        import triangle  # noqa: PLC0415
    except ImportError as exc:  # pragma: no cover - triangle normally present
        raise click.ClickException(
            "MARE2DEM → PCSF needs the 'triangle' package to rebuild the "
            "mesh from the run's .poly PSLG.  pip install triangle"
        ) from exc

    import numpy as np  # noqa: PLC0415

    from pycsamt.forward.maxwell.contracts_tri import TriMesh  # noqa: PLC0415
    from pycsamt.models.mare2dem.iotools.poly import read_poly  # noqa: PLC0415

    poly = read_poly(poly_path)
    pslg: dict[str, Any] = {
        "vertices": poly.nodes,
        "segments": poly.segments - 1,
    }
    if getattr(poly, "regions", None) is not None and len(poly.regions):
        pslg["regions"] = poly.regions
    triangulated = triangle.triangulate(pslg, "pA")
    attrs = triangulated.get("triangle_attributes")
    if attrs is not None:
        region_ids = np.round(np.asarray(attrs).ravel()).astype(np.int64)
    else:
        region_ids = np.zeros(len(triangulated["triangles"]), dtype=np.int64)
    return TriMesh(
        nodes_m=triangulated["vertices"],
        triangles=triangulated["triangles"],
        region_ids=region_ids,
    )


def _build_from_ai_arrays(sk, opts: ConvertOptions):
    import numpy as np  # noqa: PLC0415

    hints = sk.hints
    encoding = opts.encoding or hints.get("encoding", "linear")

    common = {
        "encoding": encoding,
        "source_backend": "ai",
        "created_by": opts.created_by,
        "description": (
            opts.description or "AI/DL inversion → PCSF via pycsamt CLI"
        ),
        "topo": str(opts.topo) if opts.topo else None,
        "epsg": opts.epsg,
        "utm_zone": opts.utm_zone,
    }

    if sk.path.suffix.lower() == ".npy":
        rho = np.load(sk.path, allow_pickle=False)
        arrays: dict[str, Any] = {}
    else:
        with np.load(sk.path, allow_pickle=False) as npz:
            arrays = {k: npz[k] for k in npz.files}
        rho = arrays[hints["resistivity_key"]]

    def _get(key_hint: str):
        name = hints.get(key_hint)
        return arrays.get(name) if name else None

    uncertainty = _get("uncertainty_key")
    sensitivity = _get("sensitivity_key")

    geom = sk.target_geometry

    if geom == "grid2d":
        from pycsamt.format.adapters.generic import grid2d_to_pcsf

        rho = np.asarray(rho, dtype=float)
        x = _get("x_key")
        z = _get("z_key")
        if x is None or z is None:
            z_n, x_n = rho.shape
            x = np.arange(x_n, dtype=float) if x is None else x
            z = np.arange(z_n, dtype=float) if z is None else z
        return grid2d_to_pcsf(
            rho,
            np.asarray(x, dtype=float),
            np.asarray(z, dtype=float),
            x_nodes=_get("x_nodes_key"),
            z_nodes=_get("z_nodes_key"),
            origin=opts.origin,
            azimuth_deg=opts.azimuth_deg,
            uncertainty=uncertainty,
            sensitivity=sensitivity,
            **common,
        )

    if geom == "grid3d":
        from pycsamt.format.adapters.generic import grid3d_to_pcsf

        rho = np.asarray(rho, dtype=float)
        x = _get("x_key")
        y = _get("y_key")
        z = _get("z_key")
        if x is None or y is None or z is None:
            z_n, y_n, x_n = rho.shape
            x = np.arange(x_n, dtype=float) if x is None else x
            y = np.arange(y_n, dtype=float) if y is None else y
            z = np.arange(z_n, dtype=float) if z is None else z
        return grid3d_to_pcsf(
            rho,
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            x_nodes=_get("x_nodes_key"),
            y_nodes=_get("y_nodes_key"),
            z_nodes=_get("z_nodes_key"),
            origin=opts.origin,
            rotation_deg=opts.azimuth_deg or 0.0,
            uncertainty=uncertainty,
            sensitivity=sensitivity,
            **common,
        )

    if geom == "mesh_unstructured":
        from pycsamt.format.adapters.generic import mesh_to_pcsf

        nodes = arrays[hints["nodes_key"]]
        conn = arrays[hints["connectivity_key"]]
        region_ids = _get("region_ids_key")
        rho = np.asarray(rho, dtype=float)
        n_tri = len(conn)
        n_node = len(nodes)
        kw: dict[str, Any] = {}
        if rho.shape[0] == n_tri and rho.shape[0] != n_node:
            kw["resistivity"] = rho
        elif rho.shape[0] == n_node:
            kw["resistivity_by_node"] = rho
        else:
            kw["resistivity"] = rho
        return mesh_to_pcsf(
            np.asarray(nodes, dtype=float),
            np.asarray(conn),
            region_ids=region_ids,
            uncertainty=uncertainty,
            sensitivity=sensitivity,
            **kw,
            **common,
        )

    raise click.ClickException(f"Unsupported AI geometry: {geom!r}")


# ---------------------------------------------------------------------------
# PCSFModel → file
# ---------------------------------------------------------------------------


def _write_model(model, dst: Path, canonical_format: str, log10_view: bool):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if canonical_format == "pcsf":
        from pycsamt.format.io import write_pcsf

        return write_pcsf(model, dst)
    from pycsamt.format.text import write_pcsm

    return write_pcsm(model, dst, log10_view=log10_view)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _model_report(model, path: Path) -> dict[str, Any]:
    import numpy as np  # noqa: PLC0415

    rho = model.resistivity
    report: dict[str, Any] = {
        "file": str(path),
        "size_bytes": path.stat().st_size if path.exists() else None,
        "kind": model.kind,
        "source_backend": model.source_backend,
        "created_by": model.created_by,
        "created_at": model.created_at,
        "description": model.description,
        "crs": model.crs,
        "resistivity_shape": (
            list(rho.shape) if rho is not None else None
        ),
        "native_encoding": model.resistivity_native_encoding,
        "n_stations": (
            len(model.stations.name) if model.stations is not None else 0
        ),
        "has_topography": model.topography is not None,
        "has_uncertainty": model.uncertainty is not None,
        "has_sensitivity": model.sensitivity is not None,
    }
    if rho is not None and rho.size:
        finite = rho[np.isfinite(rho)]
        if finite.size:
            report["rho_ohm_m"] = {
                "min": float(np.min(finite)),
                "max": float(np.max(finite)),
                "median": float(np.median(finite)),
            }
    try:
        meta = model.metadata_dict()
        if isinstance(meta, dict) and meta.get("model_provenance"):
            report["model_provenance"] = meta["model_provenance"]
    except Exception:  # noqa: BLE001
        pass
    return report
