# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""HDF5 reader/writer for the pyCSAMT Common Subsurface Format (PCSF).

This module is intentionally serialization-only, mirroring
:mod:`pycsamt.io.formats`'s philosophy: it knows how to turn a
:class:`~pycsamt.format.schema.PCSFModel` into a ``.pcsf`` file and
back, but carries no backend-specific (Occam2D/ModEM/MARE2DEM/DUHI)
conversion logic — that lives in :mod:`pycsamt.format.adapters`
(Phases 2-4 of ``PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md``).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from os import PathLike
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from ._version import check_pcsf_version
from .schema import (
    PCSF_VERSION,
    RESISTIVITY_UNIT,
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    StationTable,
    TopographyPerStation,
    TopographyRaster,
    UnstructuredMeshGeometry,
)

__all__ = ["write_pcsf", "read_pcsf"]

# Below this element count, gzip's chunking overhead outweighs the
# space it saves — small arrays (a handful of stations, a short
# history series) are stored uncompressed.
_COMPRESS_MIN_SIZE = 64


# ---------------------------------------------------------------------
# Low-level HDF5 helpers
# ---------------------------------------------------------------------


def _write_arr(
    group: h5py.Group, name: str, arr: np.ndarray | None, *, dtype: Any = None
) -> None:
    if arr is None:
        return
    arr = np.asarray(arr, dtype=dtype)
    kwargs: dict[str, Any] = {}
    if arr.size >= _COMPRESS_MIN_SIZE:
        kwargs = {"compression": "gzip", "compression_opts": 4}
    group.create_dataset(name, data=arr, **kwargs)


def _read_arr(group: h5py.Group | None, name: str) -> np.ndarray | None:
    if group is None or name not in group:
        return None
    return np.asarray(group[name][()])


def _write_str_list(group: h5py.Group, name: str, values: list[str]) -> None:
    str_dtype = h5py.string_dtype(encoding="utf-8")
    group.create_dataset(
        name, data=np.asarray(list(values), dtype=object), dtype=str_dtype
    )


def _read_str_list(group: h5py.Group | None, name: str) -> list[str] | None:
    if group is None or name not in group:
        return None
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in group[name][()]
    ]


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    raise TypeError(f"object of type {type(value)!r} is not JSON serializable")


def _write_json(container: h5py.Group | h5py.File, name: str, obj: Any) -> None:
    text = json.dumps(obj, default=_json_default)
    container.create_dataset(
        name, data=text, dtype=h5py.string_dtype(encoding="utf-8")
    )


def _read_json(container: h5py.Group | h5py.File, name: str) -> Any:
    if name not in container:
        return {}
    raw = container[name][()]
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    return json.loads(raw)


def _attr_str(attrs: Any, name: str, default: str = "") -> str:
    value = attrs.get(name, default)
    return default if value is None else str(value)


def _attr_float(attrs: Any, name: str) -> float | None:
    value = attrs.get(name)
    return None if value is None else float(value)


# ---------------------------------------------------------------------
# Geometry (de)serialization
# ---------------------------------------------------------------------


def _write_grid2d(group: h5py.Group, geo: Grid2DGeometry) -> None:
    group.attrs["kind"] = geo.kind
    _write_arr(group, "x", geo.x)
    _write_arr(group, "z", geo.z)
    _write_arr(group, "x_nodes", geo.x_nodes)
    _write_arr(group, "z_nodes", geo.z_nodes)
    _write_arr(group, "origin", geo.origin)
    if geo.azimuth_deg is not None:
        group.attrs["azimuth_deg"] = float(geo.azimuth_deg)


def _read_grid2d(group: h5py.Group) -> Grid2DGeometry:
    return Grid2DGeometry(
        x=_read_arr(group, "x"),
        z=_read_arr(group, "z"),
        x_nodes=_read_arr(group, "x_nodes"),
        z_nodes=_read_arr(group, "z_nodes"),
        origin=_read_arr(group, "origin"),
        azimuth_deg=_attr_float(group.attrs, "azimuth_deg"),
    )


def _write_grid3d(group: h5py.Group, geo: Grid3DGeometry) -> None:
    group.attrs["kind"] = geo.kind
    _write_arr(group, "x", geo.x)
    _write_arr(group, "y", geo.y)
    _write_arr(group, "z", geo.z)
    _write_arr(group, "x_nodes", geo.x_nodes)
    _write_arr(group, "y_nodes", geo.y_nodes)
    _write_arr(group, "z_nodes", geo.z_nodes)
    _write_arr(group, "origin", geo.origin)
    group.attrs["rotation_deg"] = float(geo.rotation_deg)
    group.attrs["n_air"] = int(geo.n_air)


def _read_grid3d(group: h5py.Group) -> Grid3DGeometry:
    return Grid3DGeometry(
        x=_read_arr(group, "x"),
        y=_read_arr(group, "y"),
        z=_read_arr(group, "z"),
        x_nodes=_read_arr(group, "x_nodes"),
        y_nodes=_read_arr(group, "y_nodes"),
        z_nodes=_read_arr(group, "z_nodes"),
        origin=_read_arr(group, "origin"),
        rotation_deg=float(group.attrs.get("rotation_deg", 0.0)),
        n_air=int(group.attrs.get("n_air", 0)),
    )


def _write_mesh(group: h5py.Group, geo: UnstructuredMeshGeometry) -> None:
    group.attrs["kind"] = geo.kind
    _write_arr(group, "nodes", geo.nodes)
    _write_arr(group, "connectivity", geo.connectivity, dtype=np.int64)
    _write_arr(group, "region_ids", geo.region_ids, dtype=np.int32)
    group.attrs["plane"] = geo.plane


def _read_mesh(group: h5py.Group) -> UnstructuredMeshGeometry:
    return UnstructuredMeshGeometry(
        nodes=_read_arr(group, "nodes"),
        connectivity=_read_arr(group, "connectivity"),
        region_ids=_read_arr(group, "region_ids"),
        plane=_attr_str(group.attrs, "plane", "xz"),
    )


def _write_multiline(group: h5py.Group, geo: MultilineGeometry) -> None:
    group.attrs["kind"] = geo.kind
    line_ids = [line.line_id for line in geo.lines]
    _write_str_list(group, "line_order", line_ids)

    lines_group = group.create_group("lines")
    for line in geo.lines:
        line_group = lines_group.create_group(line.line_id)
        _write_grid2d(line_group.create_group("geometry"), line.geometry)
        _write_arr(line_group, "resistivity", line.resistivity)
        line_group.attrs["offset_y"] = float(line.offset_y)
        line_group.attrs["offset_kind"] = line.offset_kind
        if line.azimuth_deg is not None:
            line_group.attrs["azimuth_deg"] = float(line.azimuth_deg)

    if geo.derived_volume is not None:
        dv = geo.derived_volume
        dv_group = group.create_group("derived_volume")
        _write_grid3d(dv_group.create_group("grid"), dv.grid)
        _write_arr(dv_group, "resistivity", dv.resistivity)
        dv_group.attrs["derivation_method"] = dv.derivation_method
        dv_group.attrs["synthesized"] = bool(dv.synthesized)
        _write_str_list(dv_group, "derived_from", dv.derived_from)


def _read_multiline(group: h5py.Group) -> MultilineGeometry:
    line_ids = _read_str_list(group, "line_order") or []
    lines_group = group["lines"]
    lines = []
    for line_id in line_ids:
        line_group = lines_group[line_id]
        lines.append(
            LineEntry(
                line_id=line_id,
                geometry=_read_grid2d(line_group["geometry"]),
                resistivity=_read_arr(line_group, "resistivity"),
                offset_y=float(line_group.attrs.get("offset_y", 0.0)),
                offset_kind=_attr_str(
                    line_group.attrs, "offset_kind", "synthetic"
                ),
                azimuth_deg=_attr_float(line_group.attrs, "azimuth_deg"),
            )
        )

    derived_volume = None
    if "derived_volume" in group:
        dv_group = group["derived_volume"]
        derived_volume = DerivedVolume(
            grid=_read_grid3d(dv_group["grid"]),
            resistivity=_read_arr(dv_group, "resistivity"),
            derivation_method=_attr_str(
                dv_group.attrs, "derivation_method", "linear_interp"
            ),
            derived_from=_read_str_list(dv_group, "derived_from") or [],
            synthesized=bool(dv_group.attrs.get("synthesized", True)),
        )

    return MultilineGeometry(lines=lines, derived_volume=derived_volume)


_GEOMETRY_WRITERS = {
    "grid2d": _write_grid2d,
    "grid3d": _write_grid3d,
    "mesh_unstructured": _write_mesh,
    "multiline": _write_multiline,
}
_GEOMETRY_READERS = {
    "grid2d": _read_grid2d,
    "grid3d": _read_grid3d,
    "mesh_unstructured": _read_mesh,
    "multiline": _read_multiline,
}


# ---------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------


def write_pcsf(model: PCSFModel, path: str | PathLike) -> Path:
    """Write a :class:`PCSFModel` to a ``.pcsf`` (HDF5) file.

    Parameters
    ----------
    model : PCSFModel
        The model to serialize. Validated before anything is written.
    path : path-like
        Destination file. Parent directories are created if missing.

    Returns
    -------
    pathlib.Path
        The path written to.

    Raises
    ------
    ValueError
        If *model* fails :meth:`PCSFModel.validate`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format import Grid2DGeometry, PCSFModel, write_pcsf, read_pcsf
    >>> geometry = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
    >>> model = PCSFModel(
    ...     geometry=geometry,
    ...     resistivity=np.array([[100.0, 120.0], [50.0, 60.0]]),
    ...     source_backend="occam2d",
    ... )
    >>> path = write_pcsf(model, "example.pcsf")  # doctest: +SKIP
    >>> round_tripped = read_pcsf(path)  # doctest: +SKIP
    """
    if not isinstance(model, PCSFModel):
        raise TypeError(f"model must be a PCSFModel, got {type(model)!r}")
    model.validate()

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(path, "w") as fh:
        fh.attrs["pcsf_version"] = PCSF_VERSION
        fh.attrs["source_backend"] = model.source_backend
        fh.attrs["created_by"] = model.created_by
        fh.attrs["created_at"] = (
            model.created_at or datetime.now(timezone.utc).isoformat()
        )
        fh.attrs["resistivity_unit"] = RESISTIVITY_UNIT
        if model.crs:
            fh.attrs["crs"] = model.crs
        if model.description:
            fh.attrs["description"] = model.description

        writer = _GEOMETRY_WRITERS[model.geometry.kind]
        writer(fh.create_group("geometry"), model.geometry)

        model_group = fh.create_group("model")
        _write_arr(model_group, "resistivity", model.resistivity)
        _write_arr(model_group, "resistivity_native", model.resistivity_native)
        if model.resistivity_native_encoding is not None:
            model_group.attrs["encoding"] = model.resistivity_native_encoding
        _write_arr(model_group, "uncertainty", model.uncertainty)
        _write_arr(model_group, "sensitivity", model.sensitivity)
        _write_arr(
            model_group, "resistivity_by_region", model.resistivity_by_region
        )
        _write_arr(
            model_group, "resistivity_by_node", model.resistivity_by_node
        )

        if model.stations is not None:
            st_group = fh.create_group("stations")
            _write_str_list(st_group, "name", model.stations.name)
            _write_arr(st_group, "x", model.stations.x)
            _write_arr(st_group, "y", model.stations.y)
            _write_arr(st_group, "z", model.stations.z)
            if model.stations.line_id is not None:
                _write_str_list(st_group, "line_id", model.stations.line_id)
            if model.stations.lon is not None:
                _write_arr(st_group, "lon", model.stations.lon)
                _write_arr(st_group, "lat", model.stations.lat)

        if model.topography is not None:
            topo_group = fh.create_group("topography")
            topo_group.attrs["kind"] = model.topography.kind
            if model.topography.kind == "raster":
                r_group = topo_group.create_group("raster")
                _write_arr(r_group, "x", model.topography.x)
                _write_arr(r_group, "y", model.topography.y)
                _write_arr(r_group, "elevation", model.topography.elevation)
            else:
                ps_group = topo_group.create_group("per_station")
                _write_str_list(
                    ps_group, "station_id", model.topography.station_id
                )
                _write_arr(ps_group, "elevation", model.topography.elevation)

        if model.survey:
            _write_json(fh, "survey_json", model.survey)

        if model.history:
            hist_group = fh.create_group("history")
            for key, arr in model.history.items():
                _write_arr(hist_group, key, arr)

        if model.metadata:
            _write_json(fh, "metadata_json", model.metadata)

        if model.boreholes is not None:
            from .borehole.pcsf import association_to_dict

            boreholes = fh.create_group("boreholes")
            boreholes.attrs["kind"] = "pcbh"
            boreholes.attrs["version"] = "0.1.0"
            _write_json(
                boreholes,
                "association_json",
                association_to_dict(model.boreholes),
            )

    return path


def read_pcsf(path: str | PathLike) -> PCSFModel:
    """Read a :class:`PCSFModel` back from a ``.pcsf`` (HDF5) file.

    Parameters
    ----------
    path : path-like
        Source file.

    Returns
    -------
    PCSFModel
        Fully reconstructed and re-validated model.

    Raises
    ------
    ValueError
        If ``pcsf_version`` is missing, malformed, or names an
        unrecognised MAJOR version (see ``pycsamt/format/SPEC.md``
        section 5); if ``geometry/kind`` is missing or not a
        recognised value; or if the reconstructed model fails
        :meth:`PCSFModel.validate`.

    Warns
    -----
    UserWarning
        If the file's ``pcsf_version`` MINOR component is newer than
        this reader's — fields added since then are silently ignored
        rather than causing a hard failure.
    """
    path = Path(path)

    with h5py.File(path, "r") as fh:
        if "pcsf_version" not in fh.attrs:
            raise ValueError(
                f"{path}: missing required 'pcsf_version' root attribute"
            )
        check_pcsf_version(str(fh.attrs["pcsf_version"]), PCSF_VERSION)

        if "geometry" not in fh:
            raise ValueError(f"{path}: missing required 'geometry' group")
        geometry_group = fh["geometry"]
        kind = _attr_str(geometry_group.attrs, "kind")
        reader = _GEOMETRY_READERS.get(kind)
        if reader is None:
            raise ValueError(
                f"{path}: unknown geometry kind {kind!r}; expected one "
                f"of {tuple(_GEOMETRY_READERS)}"
            )
        geometry = reader(geometry_group)

        model_group = fh.get("model")
        encoding = None
        if model_group is not None and "encoding" in model_group.attrs:
            encoding = str(model_group.attrs["encoding"])

        stations = None
        if "stations" in fh:
            st_group = fh["stations"]
            stations = StationTable(
                name=_read_str_list(st_group, "name") or [],
                x=_read_arr(st_group, "x"),
                y=_read_arr(st_group, "y"),
                z=_read_arr(st_group, "z"),
                line_id=_read_str_list(st_group, "line_id"),
                lon=_read_arr(st_group, "lon"),
                lat=_read_arr(st_group, "lat"),
            )

        topography = None
        if "topography" in fh:
            topo_group = fh["topography"]
            topo_kind = _attr_str(topo_group.attrs, "kind", "per_station")
            if topo_kind == "raster":
                r_group = topo_group["raster"]
                topography = TopographyRaster(
                    x=_read_arr(r_group, "x"),
                    y=_read_arr(r_group, "y"),
                    elevation=_read_arr(r_group, "elevation"),
                )
            elif topo_kind == "per_station":
                ps_group = topo_group["per_station"]
                topography = TopographyPerStation(
                    station_id=_read_str_list(ps_group, "station_id") or [],
                    elevation=_read_arr(ps_group, "elevation"),
                )

        history = {}
        if "history" in fh:
            hist_group = fh["history"]
            history = {key: _read_arr(hist_group, key) for key in hist_group}

        boreholes = None
        if "boreholes" in fh:
            from .borehole.pcsf import association_from_dict

            boreholes = association_from_dict(
                _read_json(fh["boreholes"], "association_json")
            )

        model = PCSFModel(
            geometry=geometry,
            resistivity=_read_arr(model_group, "resistivity"),
            resistivity_native=_read_arr(model_group, "resistivity_native"),
            resistivity_native_encoding=encoding,
            uncertainty=_read_arr(model_group, "uncertainty"),
            sensitivity=_read_arr(model_group, "sensitivity"),
            resistivity_by_region=_read_arr(
                model_group, "resistivity_by_region"
            ),
            resistivity_by_node=_read_arr(
                model_group, "resistivity_by_node"
            ),
            stations=stations,
            topography=topography,
            survey=_read_json(fh, "survey_json"),
            history=history,
            source_backend=_attr_str(fh.attrs, "source_backend", "generic"),
            created_by=_attr_str(fh.attrs, "created_by"),
            created_at=_attr_str(fh.attrs, "created_at"),
            crs=(str(fh.attrs["crs"]) if "crs" in fh.attrs else None),
            description=_attr_str(fh.attrs, "description"),
            metadata=_read_json(fh, "metadata_json"),
            boreholes=boreholes,
        )

    model.validate()
    return model
