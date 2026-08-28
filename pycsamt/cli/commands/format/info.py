# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format info
===================

Full, no-nonsense summary of a ``.pcsf`` / ``.pcsm`` / ``.pcsm.gz``
file: geometry, array shapes, resistivity statistics, station table,
topography, provenance, and inversion history.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ._base import _rich_table, fmt


@fmt.command("info")
@click.argument(
    "file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="FILE",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
def info(file: Path, output_format: str) -> None:
    """Summarise a PCSF / PCSM file.

    \b
    Examples:
      pycsamt format info model.pcsf
      pycsamt format info model.pcsm.gz -f json
    """
    import numpy as np

    from pycsamt.format.text import peek_kind, read_pcsf_or_pcsm

    name = file.name.lower()
    if not name.endswith((".pcsf", ".pcsm", ".pcsm.gz")):
        raise click.ClickException(
            f"{file.name} is not a .pcsf / .pcsm / .pcsm.gz file."
        )

    try:
        kind = peek_kind(file)
    except Exception as exc:  # noqa: BLE001
        raise click.ClickException(
            f"Cannot read header of {file}: {exc}"
        ) from exc

    try:
        model = read_pcsf_or_pcsm(file)
    except Exception as exc:  # noqa: BLE001
        raise click.ClickException(f"Cannot read {file}: {exc}") from exc

    geom = model.geometry
    report: dict = {
        "file": str(file),
        "size_bytes": file.stat().st_size,
        "container": (
            "pcsm" if name.endswith((".pcsm", ".pcsm.gz")) else "pcsf"
        ),
        "geometry_kind": kind,
        "source_backend": model.source_backend,
        "created_by": model.created_by,
        "created_at": model.created_at,
        "description": model.description,
        "crs": model.crs,
    }

    # -- geometry detail ------------------------------------------------
    gdetail: dict = {}
    for axis in ("x", "y", "z"):
        arr = getattr(geom, axis, None)
        if arr is not None:
            arr = np.asarray(arr)
            gdetail[axis] = {
                "n": int(arr.size),
                "min": float(arr.min()) if arr.size else None,
                "max": float(arr.max()) if arr.size else None,
            }
    for attr in ("nodes", "connectivity", "region_ids", "lines"):
        val = getattr(geom, attr, None)
        if val is not None:
            gdetail[attr] = (
                list(np.asarray(val).shape)
                if hasattr(val, "shape")
                else len(val)
            )
    for attr in ("origin", "rotation_deg", "azimuth_deg", "n_air", "plane"):
        val = getattr(geom, attr, None)
        if val is not None:
            gdetail[attr] = (
                np.asarray(val).tolist() if hasattr(val, "tolist") else val
            )
    report["geometry"] = gdetail

    # -- resistivity stats -------------------------------------------
    for label, arr in (
        ("resistivity", model.resistivity),
        ("resistivity_native", model.resistivity_native),
        ("resistivity_by_region", model.resistivity_by_region),
        ("resistivity_by_node", model.resistivity_by_node),
        ("uncertainty", model.uncertainty),
        ("sensitivity", model.sensitivity),
    ):
        if arr is None:
            continue
        arr = np.asarray(arr, dtype=float)
        finite = arr[np.isfinite(arr)]
        entry = {
            "shape": list(arr.shape),
            "n_nan": int(arr.size - finite.size),
        }
        if finite.size:
            entry.update(
                min=float(finite.min()),
                max=float(finite.max()),
                median=float(np.median(finite)),
            )
        report[label] = entry
    report["native_encoding"] = model.resistivity_native_encoding

    # -- stations / topography -------------------------------------
    if model.stations is not None:
        st = model.stations
        report["stations"] = {
            "n": len(st.name),
            "has_lonlat": st.lon is not None and st.lat is not None,
            "first": st.name[:3],
            "last": st.name[-3:],
        }
    if model.topography is not None:
        topo = model.topography
        ids = getattr(topo, "station_id", None)
        report["topography"] = {
            "kind": "per_station" if ids is not None else "raster",
            "n": len(ids) if ids is not None else None,
        }

    # -- history / provenance -------------------------------------
    hist = getattr(model, "history", None)
    if hist:
        report["history"] = {
            k: list(np.asarray(v).shape) for k, v in hist.items()
        }
    try:
        meta = model.metadata_dict()
        if isinstance(meta, dict) and meta.get("model_provenance"):
            report["model_provenance"] = meta["model_provenance"]
    except Exception:  # noqa: BLE001
        pass

    if output_format == "json":
        click.echo(json.dumps(report, indent=2, default=str))
        return

    rows = [
        ("file", report["file"]),
        ("size", f"{report['size_bytes']:,} bytes"),
        ("container", report["container"]),
        ("geometry kind", kind),
        ("source backend", report["source_backend"]),
        ("created by", report["created_by"] or "-"),
        ("created at", report["created_at"] or "-"),
        ("CRS", report["crs"] or "-"),
        ("native encoding", report["native_encoding"] or "linear"),
    ]
    if "resistivity" in report:
        r = report["resistivity"]
        rng = (
            f"{r['min']:.4g} … {r['max']:.4g} (median {r['median']:.4g})"
            if "min" in r
            else "all NaN"
        )
        rows.append(("resistivity (Ω·m)", f"{r['shape']}  {rng}"))
        if r["n_nan"]:
            rows.append(("resistivity NaNs", str(r["n_nan"])))
    if "stations" in report:
        rows.append(("stations", str(report["stations"]["n"])))
        rows.append(
            ("lon/lat", "yes" if report["stations"]["has_lonlat"] else "no")
        )
    if "topography" in report:
        rows.append(
            (
                "topography",
                f"{report['topography']['kind']} "
                f"({report['topography'].get('n', '-')})",
            )
        )
    if "history" in report:
        rows.append(("history keys", ", ".join(report["history"]) or "-"))
    if "model_provenance" in report:
        prov = report["model_provenance"]
        rows.append(
            (
                "model provenance",
                f"{prov.get('framework', '?')}/"
                f"{prov.get('architecture', '?')}",
            )
        )
    if report["description"]:
        rows.append(("description", report["description"]))
    _rich_table("format info", rows)

    if gdetail:
        grows = []
        for key, val in gdetail.items():
            if isinstance(val, dict) and "n" in val:
                grows.append(
                    (
                        key,
                        f"n={val['n']}  "
                        f"[{val['min']:.4g}, {val['max']:.4g}]",
                    )
                )
            else:
                grows.append((key, str(val)))
        _rich_table("geometry", grows, style="green")
