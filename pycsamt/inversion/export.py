# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Export helpers for :mod:`pycsamt.inversion` results."""

from __future__ import annotations

import csv
import json
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Union

import numpy as np

from .results import InversionResult

PathLike = Union[str, Path]

__all__ = ["to_archive", "to_csv", "to_geojson", "to_geotiff", "to_npz", "to_vtk"]


def to_csv(result: InversionResult, path: PathLike, *, log_rho: bool = True) -> Path:
    """Export a recovered resistivity model as long-form CSV.

    The CSV writer stores one row per model cell with profile position, depth,
    resistivity value, and station label. It is the simplest exchange format for
    spreadsheets, quick inspection, and downstream scripts that do not need mesh
    topology.

    Parameters
    ----------
    result : InversionResult
        Inversion result convertible through
        :meth:`pycsamt.inversion.results.InversionResult.to_resistivity_model`.
    path : path-like
        Output CSV path. Parent directories are created automatically.
    log_rho : bool, default True
        If ``True``, write ``log10(rho / ohm m)`` values. If ``False``, write
        linear resistivity in ohm metres.

    Returns
    -------
    pathlib.Path
        Path to the written CSV file.

    Examples
    --------
    >>> from pycsamt.inversion.export import to_csv
    >>> to_csv(result, "profile.csv")  # doctest: +SKIP
    >>> to_csv(result, "profile_ohm_m.csv", log_rho=False)  # doctest: +SKIP

    References
    ----------
    .. [1] Shafranovich, Y. (2005). Common Format and MIME Type for
       Comma-Separated Values (CSV) Files. RFC 4180.
    """
    model = result.to_resistivity_model()
    rho = np.asarray(model.rho_2d, dtype=float)
    values = rho if log_rho else 10.0 ** rho
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        unit = "log10_ohm_m" if log_rho else "ohm_m"
        writer.writerow(["x_m", "z_m", f"rho_{unit}", "station"])
        for ix, x in enumerate(model.x_centers):
            station = ""
            if ix < len(model.station_names):
                station = model.station_names[ix]
            for iz, z in enumerate(model.z_centers):
                writer.writerow([float(x), float(z), float(values[iz, ix]), station])
    return out


def to_npz(result: InversionResult, path: PathLike) -> Path:
    """Export common result arrays to a compressed NumPy archive.

    The NPZ writer preserves numerical arrays used by the inversion API:
    resistivity grid, coordinates, station metadata, RMS, and optional
    uncertainty/history arrays. It is the preferred lightweight format for
    Python workflows because arrays are stored without text parsing.

    Parameters
    ----------
    result : InversionResult
        Inversion result convertible to a resistivity model. Uncertainty and
        convergence-history arrays are exported when present.
    path : path-like
        Output ``.npz`` path. Parent directories are created automatically.

    Returns
    -------
    pathlib.Path
        Path to the written compressed NumPy archive.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.inversion.export import to_npz
    >>> path = to_npz(result, "profile.npz")  # doctest: +SKIP
    >>> arrays = np.load(path)  # doctest: +SKIP
    >>> arrays["rho_2d"].shape  # doctest: +SKIP
    (10, 20)

    References
    ----------
    .. [1] NumPy Developers. ``numpy.savez_compressed`` documentation.
    """
    model = result.to_resistivity_model()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "rho_2d": np.asarray(model.rho_2d, dtype=float),
        "x_centers": np.asarray(model.x_centers, dtype=float),
        "z_centers": np.asarray(model.z_centers, dtype=float),
        "station_x": np.asarray(model.station_x, dtype=float),
        "station_names": np.asarray(model.station_names, dtype=str),
        "method": result.method,
        "dimension": result.dimension,
        "backend": result.backend,
        "rms": result.rms,
    }
    if result.uncertainty is not None:
        for key in (
            "model_std",
            "covariance_diag",
            "sensitivity",
            "confidence",
            "station_confidence",
            "depth_confidence",
        ):
            value = getattr(result.uncertainty, key, None)
            if value is not None:
                payload[f"uncertainty_{key}"] = np.asarray(value, dtype=float)
    if result.history is not None:
        for key, value in result.history.arrays().items():
            payload[f"history_{key}"] = np.asarray(value, dtype=float)
    np.savez_compressed(out, **payload)
    return out


def to_geojson(result: InversionResult, path: PathLike, *, log_rho: bool = True) -> Path:
    """Export a 2-D inversion section as GeoJSON cell polygons.

    Coordinates are profile-distance/depth pairs in metres. Depth is positive
    downward, matching :class:`pycsamt.interp.ResistivityModel`.

    Parameters
    ----------
    result : InversionResult
        Inversion result convertible to a 2-D resistivity model.
    path : path-like
        Output GeoJSON path. Parent directories are created automatically.
    log_rho : bool, default True
        If ``True``, each feature contains ``rho_log10_ohm_m``. If ``False``,
        each feature contains ``rho_ohm_m``.

    Returns
    -------
    pathlib.Path
        Path to the written GeoJSON file.

    Notes
    -----
    Each model cell is written as one polygon feature. Cell properties include
    ``ix``, ``iz``, cell-center coordinates, resistivity, and matching
    uncertainty maps such as ``uncertainty_confidence`` when available.

    Examples
    --------
    >>> from pycsamt.inversion.export import to_geojson
    >>> to_geojson(result, "profile.geojson")  # doctest: +SKIP

    References
    ----------
    .. [1] Butler, H. et al. (2016). The GeoJSON Format. RFC 7946.
    """
    model = result.to_resistivity_model()
    rho = np.asarray(model.rho_2d, dtype=float)
    values = rho if log_rho else 10.0 ** rho
    x_edges = _cell_edges(model.x_centers)
    z_edges = _cell_edges(model.z_centers)
    value_key = _rho_key(log_rho)
    uncertainty = _uncertainty_arrays(result, rho.shape)
    features = []
    for iz in range(rho.shape[0]):
        z0, z1 = float(z_edges[iz]), float(z_edges[iz + 1])
        for ix in range(rho.shape[1]):
            x0, x1 = float(x_edges[ix]), float(x_edges[ix + 1])
            properties = {
                "ix": ix,
                "iz": iz,
                "x_center_m": float(model.x_centers[ix]),
                "z_center_m": float(model.z_centers[iz]),
                value_key: float(values[iz, ix]),
            }
            properties.update({
                key: float(array[iz, ix])
                for key, array in uncertainty.items()
            })
            features.append({
                "type": "Feature",
                "properties": properties,
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [x0, z0],
                        [x1, z0],
                        [x1, z1],
                        [x0, z1],
                        [x0, z0],
                    ]],
                },
            })

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    collection = {
        "type": "FeatureCollection",
        "features": features,
        "properties": _metadata(result),
    }
    with out.open("w", encoding="utf-8") as fh:
        json.dump(collection, fh, indent=2, sort_keys=True)
    return out


def to_vtk(result: InversionResult, path: PathLike, *, log_rho: bool = True) -> Path:
    """Export a 2-D inversion section as legacy ASCII VTK.

    The file is a ``RECTILINEAR_GRID`` with one cell in the cross-line
    direction, suitable for ParaView and other lightweight model viewers.

    Parameters
    ----------
    result : InversionResult
        Inversion result convertible to a 2-D resistivity model.
    path : path-like
        Output ``.vtk`` path. Parent directories are created automatically.
    log_rho : bool, default True
        If ``True``, write ``rho_log10_ohm_m`` as the primary cell scalar. If
        ``False``, write ``rho_ohm_m``.

    Returns
    -------
    pathlib.Path
        Path to the written legacy VTK file.

    Notes
    -----
    The profile coordinate is written on the VTK X axis, a dummy cross-line
    axis is written on Y, and positive-down depth is written on Z. Matching
    2-D uncertainty arrays are included as additional cell scalars.

    Examples
    --------
    >>> from pycsamt.inversion.export import to_vtk
    >>> to_vtk(result, "profile.vtk")  # doctest: +SKIP

    References
    ----------
    .. [1] Schroeder, W., Martin, K. and Lorensen, B. (2006). *The
       Visualization Toolkit*, 4th edition.
    """
    model = result.to_resistivity_model()
    rho = np.asarray(model.rho_2d, dtype=float)
    values = rho if log_rho else 10.0 ** rho
    x_edges = _cell_edges(model.x_centers)
    z_edges = _cell_edges(model.z_centers)
    uncertainty = _uncertainty_arrays(result, rho.shape)

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write("pycsamt inversion result\n")
        fh.write("ASCII\n")
        fh.write("DATASET RECTILINEAR_GRID\n")
        fh.write(f"DIMENSIONS {len(x_edges)} 2 {len(z_edges)}\n")
        _write_vtk_coords(fh, "X_COORDINATES", x_edges)
        _write_vtk_coords(fh, "Y_COORDINATES", np.array([0.0, 1.0]))
        _write_vtk_coords(fh, "Z_COORDINATES", z_edges)
        fh.write(f"CELL_DATA {rho.size}\n")
        _write_vtk_scalar(fh, _rho_key(log_rho), values)
        for key, array in uncertainty.items():
            _write_vtk_scalar(fh, key, array)
    return out


def to_geotiff(
    result: InversionResult,
    path: PathLike,
    *,
    log_rho: bool = True,
    crs: Any = None,
) -> Path:
    """Export a 2-D inversion section as a single-band GeoTIFF.

    This writer requires :mod:`rasterio`. The raster axes are profile distance
    and depth in metres; pass *crs* only when those coordinates are already tied
    to a projected coordinate reference system.

    Parameters
    ----------
    result : InversionResult
        Inversion result convertible to a 2-D resistivity model.
    path : path-like
        Output GeoTIFF path. Parent directories are created automatically.
    log_rho : bool, default True
        If ``True``, write log10 resistivity. If ``False``, write linear
        resistivity in ohm metres.
    crs : object, optional
        Coordinate reference system passed to :mod:`rasterio`. Use this only
        when profile/depth coordinates are already in a projected CRS.

    Returns
    -------
    pathlib.Path
        Path to the written GeoTIFF.

    Raises
    ------
    ImportError
        If :mod:`rasterio` is not installed.

    Examples
    --------
    >>> from pycsamt.inversion.export import to_geotiff
    >>> to_geotiff(result, "profile.tif")  # doctest: +SKIP

    References
    ----------
    .. [1] Ritter, N. and Ruth, M. (1997). GeoTIFF Format Specification.
    .. [2] Gillies, S. et al. Rasterio documentation.
    """
    try:
        import rasterio
        from rasterio.transform import Affine
    except ImportError as exc:
        raise ImportError("GeoTIFF export requires rasterio.") from exc

    model = result.to_resistivity_model()
    rho = np.asarray(model.rho_2d, dtype=float)
    values = rho if log_rho else 10.0 ** rho
    x_edges = _cell_edges(model.x_centers)
    z_edges = _cell_edges(model.z_centers)
    dx = _mean_spacing(x_edges)
    dz = _mean_spacing(z_edges)
    transform = Affine(dx, 0.0, float(x_edges[0]), 0.0, dz, float(z_edges[0]))

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        out,
        "w",
        driver="GTiff",
        height=values.shape[0],
        width=values.shape[1],
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(values.astype("float32"), 1)
        dst.update_tags(
            backend=result.backend,
            method=result.method,
            dimension=result.dimension,
            unit=_rho_key(log_rho),
        )
    return out


def to_archive(
    result: InversionResult,
    path: PathLike,
    *,
    include_native: bool = True,
    log_rho: bool = True,
) -> Path:
    """Write a portable ZIP snapshot of an inversion result.

    The archive always contains common products (metadata, NPZ, CSV). Existing
    backend-native files referenced by ``result.files`` are included under
    ``native_files/`` when *include_native* is true.

    Parameters
    ----------
    result : InversionResult
        Inversion result to snapshot.
    path : path-like
        Output ``.zip`` path. Parent directories are created automatically.
    include_native : bool, default True
        Include existing backend-native files listed in ``result.files`` under
        ``native_files/``.
    log_rho : bool, default True
        Passed to the CSV export stored inside the archive.

    Returns
    -------
    pathlib.Path
        Path to the written ZIP archive.

    Notes
    -----
    The archive contains:

    * ``metadata.json`` with backend, RMS, warnings, file references, and
      history/uncertainty metadata.
    * ``result.npz`` with numerical arrays.
    * ``model.csv`` with the long-form resistivity model.
    * ``native_files/`` entries when available and requested.

    Examples
    --------
    >>> from pycsamt.inversion.export import to_archive
    >>> to_archive(result, "run_snapshot.zip")  # doctest: +SKIP

    References
    ----------
    .. [1] PKWARE Inc. ZIP File Format Specification.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pycsamt-inv-export-") as tmp_name:
        tmp = Path(tmp_name)
        metadata_path = tmp / "metadata.json"
        with metadata_path.open("w", encoding="utf-8") as fh:
            json.dump(_metadata(result), fh, indent=2, sort_keys=True)
        npz_path = to_npz(result, tmp / "result.npz")
        csv_path = to_csv(result, tmp / "model.csv", log_rho=log_rho)

        with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.write(metadata_path, "metadata.json")
            zf.write(npz_path, "result.npz")
            zf.write(csv_path, "model.csv")
            if include_native:
                for key, value in sorted(result.files.items()):
                    source = Path(value)
                    if source.exists() and source.is_file():
                        zf.write(source, f"native_files/{key}_{source.name}")
    return out


def _cell_edges(centers: Any) -> np.ndarray:
    centers = np.asarray(centers, dtype=float).reshape(-1)
    if centers.size == 0:
        raise ValueError("cannot export model with empty cell coordinates.")
    if centers.size == 1:
        half_width = max(abs(float(centers[0])) * 0.05, 0.5)
        return np.array([centers[0] - half_width, centers[0] + half_width], dtype=float)
    mids = 0.5 * (centers[:-1] + centers[1:])
    first = centers[0] - (mids[0] - centers[0])
    last = centers[-1] + (centers[-1] - mids[-1])
    return np.r_[first, mids, last].astype(float)


def _mean_spacing(edges: np.ndarray) -> float:
    spacing = np.diff(np.asarray(edges, dtype=float))
    if spacing.size == 0 or not np.all(np.isfinite(spacing)):
        return 1.0
    return float(np.nanmean(np.abs(spacing))) or 1.0


def _rho_key(log_rho: bool) -> str:
    return "rho_log10_ohm_m" if log_rho else "rho_ohm_m"


def _uncertainty_arrays(result: InversionResult, shape: tuple[int, int]) -> dict[str, np.ndarray]:
    if result.uncertainty is None:
        return {}
    out: dict[str, np.ndarray] = {}
    for attr in ("model_std", "covariance_diag", "sensitivity", "confidence"):
        value = getattr(result.uncertainty, attr, None)
        if value is None:
            continue
        array = np.asarray(value, dtype=float)
        if array.shape == shape:
            out[f"uncertainty_{attr}"] = array
    return out


def _write_vtk_coords(fh: Any, name: str, values: np.ndarray) -> None:
    values = np.asarray(values, dtype=float).reshape(-1)
    fh.write(f"{name} {values.size} float\n")
    fh.write(" ".join(f"{float(value):.12g}" for value in values))
    fh.write("\n")


def _write_vtk_scalar(fh: Any, name: str, values: np.ndarray) -> None:
    array = np.asarray(values, dtype=float)
    fh.write(f"SCALARS {name} float 1\n")
    fh.write("LOOKUP_TABLE default\n")
    for value in array.reshape(-1):
        fh.write(f"{float(value):.12g}\n")


def _metadata(result: InversionResult) -> dict[str, Any]:
    payload = {
        "method": result.method,
        "dimension": result.dimension,
        "backend": result.backend,
        "status": result.status,
        "rms": result.rms,
        "objective": result.objective,
        "n_iter": result.n_iter,
        "workdir": result.workdir,
        "files": dict(result.files),
        "warnings": list(result.warnings),
        "metadata": dict(result.metadata),
    }
    if result.uncertainty is not None:
        payload["uncertainty"] = {"metadata": dict(result.uncertainty.metadata)}
    if result.history is not None:
        payload["history"] = {
            "n_records": result.history.n_records,
            "metadata": dict(result.history.metadata),
        }
    return _json_ready(payload)


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    try:
        json.dumps(value)
    except TypeError:
        return str(value)
    return value
