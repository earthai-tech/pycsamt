# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Model sources for the Interpretation Studio (Qt-free).

Interpretation works on a 2-D :class:`~pycsamt.interp.ResistivityModel`
(log10 rho on a profile x depth grid).  :func:`load_model` builds one from
anything the desktop can produce or open:

=====================  ====================================================
``.pcsf`` / ``.pcsm``  ``grid2d``: the full model grid.  ``multiline`` /
                       ``grid3d`` / ``mesh_unstructured``: the per-line
                       station sections (the ones Map View draws), stations
                       placed by their distance along the line.
Occam2D run folder     native ``InversionResult`` (full mesh)
ModEM 2-D run folder   ``ResistivityModel.from_any``
ModEM 3-D / MARE2DEM   converted to PCSF in memory, then sections
Occam1D run folder     stations' final 1-D models stitched along the line
=====================  ====================================================

It raises :class:`ModelSourceError` with a user-facing reason when a
source cannot give a 2-D model.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["ModelInfo", "ModelSourceError", "list_lines", "load_model"]


class ModelSourceError(ValueError):
    """The source cannot provide a 2-D model; the message says why."""


@dataclass
class ModelInfo:
    source: str  # human-readable origin ("ModEM 3-D PCSF, line L1")
    path: str
    kind: str  # "pcsf" | "occam2d" | "modem2d" | "occam1d" | ...
    line: str | None = None
    lines: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def _is_pcsf(p: Path) -> bool:
    name = p.name.lower()
    return name.endswith((".pcsf", ".pcsm", ".pcsm.gz"))


def list_lines(path) -> list[str]:
    """Profile lines a PCSF/PCSM file offers (one for grid2d)."""
    p = Path(path)
    if not _is_pcsf(p):
        return []
    from pycsamt.map import MapView

    view = MapView.from_pcsf(str(p), fetch_elevation=False)
    return [str(k) for k in (view.data.metadata or {}).get("sections", {})]


def _along_line_distance(xy: np.ndarray) -> np.ndarray:
    """Distance of each point along the line's principal direction,
    starting at 0 (stations need not be in order)."""
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 2 or not np.isfinite(xy).all():
        return np.arange(len(xy), dtype=float)
    c = xy - xy.mean(axis=0)
    _u, _s, vt = np.linalg.svd(c, full_matrices=False)
    d = c @ vt[0]
    return d - d.min()


def _station_xy(model) -> dict[str, tuple[float, float]]:
    st = model.stations
    if st is None:
        return {}
    names = [str(n) for n in st.name]
    x = np.asarray(st.x, dtype=float)
    y = np.asarray(st.y, dtype=float) if st.y is not None else \
        np.zeros_like(x)
    return {n: (a, b) for n, a, b in zip(names, x, y)}


def _from_pcsf(path: Path, line: str | None):
    from pycsamt.format.text import read_pcsf_or_pcsm
    from pycsamt.interp import ResistivityModel
    from pycsamt.map import MapView

    pm = read_pcsf_or_pcsm(str(path))
    kind = getattr(pm, "kind", "") or type(pm.geometry).__name__
    backend = getattr(pm, "source_backend", "") or "PCSF"
    if kind == "grid2d":
        g = pm.geometry
        rho = np.asarray(pm.resistivity, dtype=float)
        st = pm.stations
        names = [str(n) for n in st.name] if st is not None else None
        sx = np.asarray(st.x, dtype=float) if st is not None else None
        model = ResistivityModel.from_array(
            np.log10(np.where(rho > 0, rho, np.nan)),
            np.asarray(g.x, float), np.asarray(g.z, float),
            station_x=sx, station_names=names, method=str(backend))
        return model, ModelInfo(f"{backend} 2-D model (PCSF)", str(path),
                                "pcsf", lines=["line"])

    view = MapView.from_pcsf(str(path), fetch_elevation=False)
    sections = (view.data.metadata or {}).get("sections", {})
    if not sections:
        raise ModelSourceError(
            f"{path.name}: no profile sections could be read from this "
            "model (it has no station table).")
    lines = [str(k) for k in sections]
    chosen = line if line in sections else lines[0]
    sec = sections[chosen]
    names = [str(n) for n in sec["stations"]]
    rho = np.asarray(sec["rho"], dtype=float)  # (nz, n_station)
    z = np.asarray(sec["z"], dtype=float)
    xy = _station_xy(pm)
    pts = np.array([xy.get(n, (np.nan, np.nan)) for n in names])
    dist = _along_line_distance(pts) if np.isfinite(pts).all() else \
        np.arange(len(names), dtype=float)
    order = np.argsort(dist)
    model = ResistivityModel.from_array(
        np.log10(np.where(rho > 0, rho, np.nan))[:, order], dist[order], z,
        station_x=dist[order], station_names=[names[i] for i in order],
        method=str(backend))
    label = {"grid3d": "3-D", "multiline": "multi-line",
             "mesh_unstructured": "unstructured-mesh"}.get(kind, kind)
    info = ModelInfo(f"{backend} {label} model (PCSF), line {chosen}",
                     str(path), "pcsf", line=chosen, lines=lines)
    if kind == "grid3d":
        info.notes.append("3-D model: the section follows the stations of "
                          "the line (one column per station).")
    return model, info


def _from_occam1d(path: Path):
    from pycsamt.app.desktop.controllers.inversion_engines import (
        _o1d_load_station,
        engine,
    )
    from pycsamt.interp import ResistivityModel

    run = engine("occam1d").load(path)
    dirs = run.extra["dirs"]
    cols, names, depth = [], [], None
    for name, wd in dirs.items():
        st = _o1d_load_station(wd)
        if st["result"] is None:
            continue
        if depth is None:
            depth = st["depth"]
        if st["depth"].shape != depth.shape:
            continue
        cols.append(np.log10(st["result"].final.resistivity))
        names.append(name)
    if len(cols) < 2:
        raise ModelSourceError("An Occam1D run needs at least two inverted "
                               "stations to make a section.")
    z = depth + np.gradient(depth) / 2.0  # layer mid-depths
    x = np.arange(len(names), dtype=float)
    model = ResistivityModel.from_array(np.column_stack(cols), x, z,
                                        station_x=x, station_names=names,
                                        method="occam1d")
    info = ModelInfo("Occam1D stations stitched along the line", str(path),
                     "occam1d")
    info.notes.append("Stations are evenly spaced (Occam1D run folders "
                      "carry no positions).")
    return model, info


def _from_mare2dem(path: Path, nx: int = 160, nz: int = 80):
    """Resample a MARE2DEM model (triangular regions) onto a regular grid
    over its inverted (free) regions; receivers become the stations."""
    import matplotlib.tri as mtri

    from pycsamt.interp import ResistivityModel
    from pycsamt.models.mare2dem import InversionResult
    from pycsamt.models.mare2dem.plot import PlotModel

    res = InversionResult(path)
    if res.model is None:
        raise ModelSourceError(f"No .resistivity model in {path.name}.")
    pm = PlotModel(res)
    nodes, tris, logrho = pm._load_mesh()
    if nodes is None:
        raise ModelSourceError(
            f"{path.name}: the model mesh (.poly) could not be rebuilt "
            "(needs the 'triangle' package).")
    free = getattr(pm, "_free_tris", None)
    use = free if free is not None and free.any() else         np.isfinite(logrho)
    pts = np.unique(np.asarray(tris)[use].ravel())
    y0, y1 = nodes[pts, 0].min(), nodes[pts, 0].max()
    z0, z1 = max(nodes[pts, 1].min(), 0.0), nodes[pts, 1].max()
    names, sx = None, None
    em = res.data
    mt = getattr(em, "mt", None) if em is not None else None
    if mt is not None and len(getattr(mt, "receivers", [])):
        sx = np.asarray(mt.receivers, float)[:, 1]
        names = list(mt.receiver_name) or None
    if sx is not None and sx.size > 1:
        # the stations' span (+20 %), not the whole inverted domain
        span = max(np.ptp(sx), 1.0)
        y0 = max(y0, sx.min() - 0.2 * span)
        y1 = min(y1, sx.max() + 0.2 * span)
    xs = np.linspace(y0, y1, nx)
    # geometric depth spacing: shallow detail, still reaching the base
    top = max(z0, (z1 - z0) / 2000.0, 1.0)
    zs = np.concatenate([[top], np.geomspace(top * 1.5, z1, nz - 1)])
    tri = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)
    finder = tri.get_trifinder()
    X, Z = np.meshgrid(xs, zs)
    idx = finder(X, Z)
    grid = np.where(idx >= 0, np.asarray(logrho)[np.clip(idx, 0, None)],
                    np.nan)
    model = ResistivityModel.from_array(grid, xs, zs, station_x=sx,
                                        station_names=names,
                                        method="mare2dem")
    info = ModelInfo("MARE2DEM run folder (resampled onto a regular grid)",
                     str(path), "mare2dem")
    info.notes.append(f"Triangular regions sampled on {nx} x {nz} cells "
                      "over the inverted area.")
    return model, info


def load_model(path, *, line: str | None = None):
    """``(ResistivityModel, ModelInfo)`` for *path* (file or run folder)."""
    p = Path(path)
    if not p.exists():
        raise ModelSourceError(f"{p} does not exist.")
    if p.is_file():
        if _is_pcsf(p):
            return _from_pcsf(p, line)
        raise ModelSourceError(
            f"{p.name}: open a .pcsf/.pcsm file or an inversion run folder.")
    from pycsamt.app.desktop.controllers.inversion_engines import (
        detect_engine,
    )

    key = detect_engine(p)
    if key is None:
        raise ModelSourceError(
            f"No Occam1D, Occam2D, ModEM or MARE2DEM run found in {p}.")
    if key == "occam2d":
        from pycsamt.interp import ResistivityModel
        from pycsamt.models.occam2d import InversionResult

        model = ResistivityModel.from_occam2d(InversionResult(str(p)))
        return model, ModelInfo("Occam2D run folder", str(p), "occam2d")
    if key == "modem2d":
        from pycsamt.interp import ResistivityModel
        from pycsamt.models.modem import InversionResult

        model = ResistivityModel.from_any(InversionResult(p))
        return model, ModelInfo("ModEM 2-D run folder", str(p), "modem2d")
    if key == "occam1d":
        return _from_occam1d(p)
    if key == "mare2dem":
        return _from_mare2dem(p)
    # ModEM 3-D: through PCSF (a 3-D grid is not a 2-D section), exactly
    # as "Export PCSF" would write it.
    from pycsamt.format import convert_engine as ce

    try:
        pm = ce.build_model(ce.detect(p, None), ce.ConvertOptions())
    except Exception as exc:
        raise ModelSourceError(f"Could not convert {p.name}: {exc}") from exc
    tmp = Path(tempfile.mkdtemp(prefix="pycsamt_interp_")) / f"{p.name}.pcsf"
    ce.write_model(pm, tmp, "pcsf", False)
    model, info = _from_pcsf(tmp, line)
    info.source = info.source.replace("(PCSF)", f"({p.name})")
    info.path = str(p)
    info.kind = key
    return model, info
