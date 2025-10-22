# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable,  Dict
from tempfile import TemporaryDirectory
import csv
import zipfile

import numpy as np

from .utils import (
    iter_edifiles, station_name, ensure_head,
)


__all__ = [
    "write_site",
    "write_sites",
    "pack_zip",
]

class _SafeDict(dict):
    def __missing__(self, key):
        return ""


def _context_for(ed: Any, index: int) -> Dict[str, Any]:
    h = ensure_head(ed)
    nm = station_name(ed)
    lat = getattr(h, "lat", np.nan)
    lon = getattr(h, "long", np.nan)  
    elv = getattr(h, "elev", np.nan)
    ch = getattr(ed, "chainage", np.nan)
    return {
        "station": nm or "",
        "index": int(index),
        "lat": float(lat) if lat is not None else np.nan,
        "lon": float(lon) if lon is not None else np.nan,
        "elev": float(elv) if elv is not None else np.nan,
        "chainage": float(ch) if ch is not None else np.nan,
    }


def _render_name(template: str, ctx: Dict[str, Any]) -> str:
    return str(template).format_map(_SafeDict(ctx))


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _write_via_backend(ed: Any, path: Path) -> None:
    """
    Best-effort write using common backend spellings.
    """
    # 1) most pycsamt EDI backends accept `write(new_edifn=...)`
    try:
        ed.write(new_edifn=str(path))
        return
    except Exception:
        pass

    # 2) sometimes `write(path)` exists
    try:
        ed.write(str(path))
        return
    except Exception:
        pass

    # 3) or `to_file(path)` is provided
    try:
        ed.to_file(str(path))
        return
    except Exception:
        pass

    # 4) last chance: `save(path)`
    try:
        ed.save(str(path))
        return
    except Exception:
        pass

    raise RuntimeError("Cannot write EDI file for site")


def _rows_for_manifest(
    items: Iterable[tuple[int, Any, Path]]
) -> list[dict]:
    rows: list[dict] = []
    for idx, ed, fp in items:
        c = _context_for(ed, idx)
        rows.append(
            {
                "index": c["index"],
                "station": c["station"],
                "lat": c["lat"],
                "lon": c["lon"],
                "elev": c["elev"],
                "chainage": c["chainage"],
                "filename": str(fp.name),
                "path": str(fp),
            }
        )
    return rows


def _write_manifest_csv(
    rows: list[dict], csv_path: Path
) -> None:
    if not rows:
        return
    _ensure_parent(csv_path)
    fields = [
        "index", "station", "lat", "lon",
        "elev", "chainage", "filename", "path",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_site(
    site: Any,
    path: str | Path,
) -> Path:
    """
    Write a single site (EDI) to `path`. Returns `Path`.
    """
    out = Path(path)
    _ensure_parent(out)
    _write_via_backend(site, out)
    return out


def write_sites(
    sites: Any,
    outdir: str | Path,
    *,
    template: str = "{station}.edi",
    exist_ok: bool = False,
    manifest_csv: str | Path | None = None,
) -> list[Path]:
    """
    Write a collection of sites to `outdir` using a template.
    """
    out_root = Path(outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    written: list[tuple[int, Any, Path]] = []
    for i, ed in enumerate(iter_edifiles(sites)):
        ctx = _context_for(ed, i)
        name = _render_name(template, ctx) or "site.edi"
        if not name.lower().endswith(".edi"):
            name = f"{name}.edi"
        dest = out_root / name
        if dest.exists() and not exist_ok:
            raise FileExistsError(dest)
        _write_via_backend(ed, dest)
        written.append((i, ed, dest))

    if manifest_csv:
        _write_manifest_csv(
            _rows_for_manifest(written),
            Path(manifest_csv),
        )

    return [p for _, _, p in written]


def pack_zip(
    sites: Any,
    out_zip: str | Path,
    *,
    template: str = "{station}.edi",
    manifest_csv: str | Path | None = None,
) -> Path:
    """
    Pack sites into a zip file. Filenames use `template`.
    """
    out_zip = Path(out_zip)
    _ensure_parent(out_zip)

    rows: list[dict] = []
    with TemporaryDirectory() as td:
        tmp = Path(td)
        file_map: list[tuple[Path, str]] = []

        for i, ed in enumerate(iter_edifiles(sites)):
            ctx = _context_for(ed, i)
            name = _render_name(template, ctx) or "site.edi"
            if not name.lower().endswith(".edi"):
                name = f"{name}.edi"
            dst = tmp / name
            _ensure_parent(dst)
            _write_via_backend(ed, dst)
            file_map.append((dst, name))
            rows.append(
                {
                    "index": ctx["index"],
                    "station": ctx["station"],
                    "lat": ctx["lat"],
                    "lon": ctx["lon"],
                    "elev": ctx["elev"],
                    "chainage": ctx["chainage"],
                    "filename": name,
                    "path": str(out_zip),
                }
            )

        with zipfile.ZipFile(out_zip, "w",
                             compression=zipfile.ZIP_DEFLATED) as zf:
            for fp, arcname in file_map:
                zf.write(fp, arcname)

    if manifest_csv:
        _write_manifest_csv(rows, Path(manifest_csv))

    return out_zip
