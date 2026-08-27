# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Generic array-based -> PCSF adapter, for any AI/DL inversion result.

Every other module in :mod:`pycsamt.format.adapters` converts a *specific*
solver's own in-memory result object (Occam2D's ``InversionResult``, ModEM's
``ModEmModel3D``, MARE2DEM's ``TriMesh``). This module needs none of that:
it builds a :class:`~pycsamt.format.schema.PCSFModel` straight from plain
``numpy`` arrays, so a UNet, a GCN, a ResNet, or any third-party AI/DL
inversion tool -- with zero dependency on pycsamt's own
:mod:`pycsamt.ai` subpackage -- can publish its result as a citable,
reproducible ``.pcsf``/``.pcsm`` file others can reload and check.

Three functions, one per non-multiline geometry kind:

- :func:`grid2d_to_pcsf` -- a regular 2-D section (e.g. a UNet slice).
- :func:`grid3d_to_pcsf` -- a regular 3-D volume (e.g. a 3-D CNN/ResNet).
- :func:`mesh_to_pcsf` -- a triangular mesh, per-node *or* per-cell (e.g. a
  GCN predicting one value per graph vertex).

An AI model producing several 2-D lines with real cross-strike offsets
already has a home: :func:`pycsamt.format.multiline.build_multiline_pcsf`,
whose ``source_backend`` is already free text -- no changes needed there.

All three functions share the same trailing keyword surface as the
solver-specific adapters (station identity/position, elevation/topography
via the same smart ``topo=`` resolver, an ``origin``/rotation "offset" for
placing the result in real-world space, uncertainty/sensitivity, iteration
history, and free-form survey/model metadata), so an AI-produced file gets
the same spatial richness a real solver's file gets -- see
:mod:`pycsamt.format.provenance` for the accompanying ML-provenance
convention (architecture/framework/checkpoint/hyperparameters).
"""

from __future__ import annotations

import warnings
from typing import Any, Mapping, Sequence

import numpy as np

from ..provenance import ModelProvenance
from ..schema import (
    RESISTIVITY_ENCODINGS,
    Grid2DGeometry,
    Grid3DGeometry,
    PCSFModel,
    StationTable,
    TopographyPerStation,
    TopographyRaster,
    UnstructuredMeshGeometry,
)
from ..topo_source import resolve_topo

__all__ = ["grid2d_to_pcsf", "grid3d_to_pcsf", "mesh_to_pcsf"]


def _survey_to_dict(survey: Any | Mapping[str, Any] | None) -> dict[str, Any]:
    if survey is None:
        return {}
    to_dict = getattr(survey, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())
    return dict(survey)


def _provenance_to_metadata(
    provenance: ModelProvenance | Mapping[str, Any] | None,
    metadata: Mapping[str, Any] | None,
) -> dict[str, Any]:
    merged = dict(metadata or {})
    if provenance is not None:
        merged["model_provenance"] = (
            provenance.to_dict()
            if isinstance(provenance, ModelProvenance)
            else dict(provenance)
        )
    return merged


def _apply_encoding(
    values: Any, encoding: str
) -> tuple[np.ndarray, np.ndarray | None, str | None]:
    """Returns ``(canonical_linear, native, native_encoding)``.

    *native*/*native_encoding* are ``None`` when *encoding* is already
    ``"linear"`` -- matching every other adapter's convention of only
    populating ``resistivity_native`` when a real conversion happened.
    """
    arr = np.asarray(values, dtype=float)
    if encoding == "linear":
        return arr, None, None
    if encoding == "log10":
        return 10.0**arr, arr, "log10"
    if encoding == "ln":
        return np.exp(arr), arr, "ln"
    raise ValueError(
        f"encoding must be one of {RESISTIVITY_ENCODINGS}, got {encoding!r}"
    )


def _build_stations(
    names: Sequence[str],
    x: Sequence[float] | None,
    y: Sequence[float] | None,
    z: Sequence[float] | None,
    *,
    station_elevations: Mapping[str, float] | None,
    station_lonlat: Mapping[str, tuple[float, float]] | None,
    topo: Any,
    epsg: int | None,
    utm_zone: Any | None,
    latlon: bool,
    on_mismatch: str,
    caller: str,
) -> tuple[StationTable, TopographyPerStation | None]:
    names = list(names)
    n = len(names)
    x_arr = np.zeros(n) if x is None else np.asarray(x, dtype=float)
    y_arr = np.zeros(n) if y is None else np.asarray(y, dtype=float)
    if z is not None:
        z_arr = np.asarray(z, dtype=float)
    elif station_elevations:
        z_arr = np.array(
            [float(station_elevations.get(name, np.nan)) for name in names]
        )
    else:
        z_arr = np.full(n, np.nan)

    lon = lat = None
    if station_lonlat:
        lonlat = np.array(
            [station_lonlat.get(name, (np.nan, np.nan)) for name in names],
            dtype=float,
        )
        if not np.all(np.isnan(lonlat)):
            lon, lat = lonlat[:, 0], lonlat[:, 1]

    elevation_overrides = dict(station_elevations or {})
    if topo is not None:
        if station_lonlat or station_elevations:
            warnings.warn(
                f"{caller}: 'topo' takes precedence over "
                "'station_lonlat'/'station_elevations' for every station "
                "it resolves.",
                UserWarning,
                stacklevel=3,
            )
        attr = resolve_topo(
            topo, names, epsg=epsg, utm_zone=utm_zone, latlon=latlon,
            on_mismatch=on_mismatch,
        )
        lon = np.array(
            [attr.lon.get(n, lon[i] if lon is not None else np.nan) for i, n in enumerate(names)]
        )
        lat = np.array(
            [attr.lat.get(n, lat[i] if lat is not None else np.nan) for i, n in enumerate(names)]
        )
        if np.all(np.isnan(lon)):
            lon = lat = None
        z_arr = np.array(
            [attr.elevation.get(n, z_arr[i]) for i, n in enumerate(names)]
        )
        elevation_overrides.update(attr.elevation)

    stations = StationTable(
        name=names, x=x_arr, y=y_arr, z=z_arr, lon=lon, lat=lat
    )
    topography = None
    if elevation_overrides:
        matched = [
            (name, elevation_overrides[name])
            for name in names
            if name in elevation_overrides
        ]
        if matched:
            topography = TopographyPerStation(
                station_id=[name for name, _ in matched],
                elevation=np.array([elev for _, elev in matched], dtype=float),
            )
    return stations, topography


def _resolve_spatial(
    *,
    stations: StationTable | None,
    topography: TopographyPerStation | TopographyRaster | None,
    station_names: Sequence[str] | None,
    station_x: Sequence[float] | None,
    station_y: Sequence[float] | None,
    station_z: Sequence[float] | None,
    station_elevations: Mapping[str, float] | None,
    station_lonlat: Mapping[str, tuple[float, float]] | None,
    topo: Any,
    epsg: int | None,
    utm_zone: Any | None,
    latlon: bool,
    on_mismatch: str,
    caller: str,
) -> tuple[StationTable | None, TopographyPerStation | TopographyRaster | None]:
    if stations is not None:
        return stations, topography
    if station_names is None or len(station_names) == 0:
        return None, topography
    built_stations, built_topography = _build_stations(
        station_names, station_x, station_y, station_z,
        station_elevations=station_elevations,
        station_lonlat=station_lonlat,
        topo=topo, epsg=epsg, utm_zone=utm_zone, latlon=latlon,
        on_mismatch=on_mismatch, caller=caller,
    )
    return built_stations, (topography or built_topography)


def grid2d_to_pcsf(
    resistivity: Any,
    x: Any,
    z: Any,
    *,
    encoding: str = "linear",
    x_nodes: Any = None,
    z_nodes: Any = None,
    origin: Any = None,
    azimuth_deg: float | None = None,
    uncertainty: Any = None,
    sensitivity: Any = None,
    stations: StationTable | None = None,
    topography: TopographyPerStation | TopographyRaster | None = None,
    station_names: Sequence[str] | None = None,
    station_x: Sequence[float] | None = None,
    station_elevations: Mapping[str, float] | None = None,
    station_lonlat: Mapping[str, tuple[float, float]] | None = None,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
    history: Mapping[str, Any] | None = None,
    provenance: ModelProvenance | Mapping[str, Any] | None = None,
    survey: Any | Mapping[str, Any] | None = None,
    source_backend: str = "ai",
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
    metadata: Mapping[str, Any] | None = None,
) -> PCSFModel:
    r"""Build a ``grid2d`` :class:`PCSFModel` from a bare AI/DL prediction.

    The natural fit for a 2-D model (e.g. a UNet trained on inversion
    sections): give it its predicted array plus the coordinates it was
    predicted on, and it becomes the exact same on-disk artifact
    :func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf` produces for a
    real Occam2D run -- readable by the same :func:`pycsamt.format.read_pcsf`,
    ``app/mapview``, and the web 3-D view.

    Parameters
    ----------
    resistivity : array-like, shape (n_z, n_x)
        The model's predicted resistivity, in whatever *encoding* it was
        trained/predicted in.
    x, z : array-like
        Cell-centre coordinates (metres) matching
        :class:`~pycsamt.format.schema.Grid2DGeometry`.
    encoding : {"linear", "log10", "ln"}, default "linear"
        Encoding of *resistivity*. Many DL models predict log-resistivity
        for training stability; when not ``"linear"``, the given array is
        kept verbatim as ``resistivity_native`` and the canonical linear
        ``resistivity`` is derived automatically.
    x_nodes, z_nodes, origin, azimuth_deg
        Forwarded to :class:`~pycsamt.format.schema.Grid2DGeometry`.
        *origin*/*azimuth_deg* are this geometry's own "offset" -- the
        real-world placement of an otherwise locally-referenced section.
    uncertainty, sensitivity : array-like, optional
        Same shape as *resistivity* -- a direct fit for a predictive
        standard deviation from an MC-dropout/Bayesian model.
    stations : StationTable, optional
        A pre-built table, used as-is when given (bypasses every
        ``station_*``/``topo`` parameter below).
    topography : TopographyPerStation or TopographyRaster, optional
        A pre-built topography, used as-is when given -- e.g. a DEM raster
        from :func:`pycsamt.format.topography.topography_from_grid`, which
        no ``station_elevations``/``topo`` combination below can produce.
    station_names, station_x : sequence, optional
        Station identity and along-profile chainage (same frame as *x*).
    station_elevations, station_lonlat, topo, epsg, utm_zone, latlon, on_mismatch
        Identical to :func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf`'s
        own parameters of the same name -- see that function's docstring
        for the full description of the smart ``topo=`` resolver.
    history : mapping of str to array-like, optional
        Per-iteration/per-epoch series (e.g. training/validation loss).
    provenance : ModelProvenance or mapping, optional
        Folded into ``metadata["model_provenance"]`` -- see
        :mod:`pycsamt.format.provenance`.
    survey, source_backend, created_by, crs, description, metadata
        Passed through to :class:`PCSFModel`; *source_backend* defaults
        to ``"ai"`` but any free-text string is accepted (e.g.
        ``"unet"``, ``"gcn"``, ``"resnet"``, a third party's own name).

    Returns
    -------
    PCSFModel

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format.adapters.generic import grid2d_to_pcsf
    >>> model = grid2d_to_pcsf(
    ...     resistivity=np.array([[2.0, 2.1], [1.7, 1.8]]),
    ...     x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]),
    ...     encoding="log10", source_backend="unet",
    ... )
    >>> model.kind, model.source_backend
    ('grid2d', 'unet')
    """
    geometry = Grid2DGeometry(
        x=x, z=z, x_nodes=x_nodes, z_nodes=z_nodes,
        origin=None if origin is None else np.asarray(origin, dtype=float),
        azimuth_deg=azimuth_deg,
    )
    resistivity_lin, resistivity_native, native_encoding = _apply_encoding(
        resistivity, encoding
    )
    stations, topography = _resolve_spatial(
        stations=stations, topography=topography,
        station_names=station_names, station_x=station_x, station_y=None,
        station_z=None, station_elevations=station_elevations,
        station_lonlat=station_lonlat, topo=topo, epsg=epsg,
        utm_zone=utm_zone, latlon=latlon, on_mismatch=on_mismatch,
        caller="grid2d_to_pcsf",
    )

    return PCSFModel(
        geometry=geometry,
        resistivity=resistivity_lin,
        resistivity_native=resistivity_native,
        resistivity_native_encoding=native_encoding,
        uncertainty=uncertainty,
        sensitivity=sensitivity,
        stations=stations,
        topography=topography,
        survey=_survey_to_dict(survey),
        history=dict(history or {}),
        source_backend=source_backend,
        created_by=created_by,
        crs=crs,
        description=description,
        metadata=_provenance_to_metadata(provenance, metadata),
    )


def grid3d_to_pcsf(
    resistivity: Any,
    x: Any,
    y: Any,
    z: Any,
    *,
    encoding: str = "linear",
    x_nodes: Any = None,
    y_nodes: Any = None,
    z_nodes: Any = None,
    origin: Any = None,
    rotation_deg: float = 0.0,
    n_air: int = 0,
    uncertainty: Any = None,
    sensitivity: Any = None,
    stations: StationTable | None = None,
    topography: TopographyPerStation | TopographyRaster | None = None,
    station_names: Sequence[str] | None = None,
    station_x: Sequence[float] | None = None,
    station_y: Sequence[float] | None = None,
    station_z: Sequence[float] | None = None,
    station_elevations: Mapping[str, float] | None = None,
    station_lonlat: Mapping[str, tuple[float, float]] | None = None,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
    history: Mapping[str, Any] | None = None,
    provenance: ModelProvenance | Mapping[str, Any] | None = None,
    survey: Any | Mapping[str, Any] | None = None,
    source_backend: str = "ai",
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
    metadata: Mapping[str, Any] | None = None,
) -> PCSFModel:
    r"""Build a ``grid3d`` :class:`PCSFModel` from a bare AI/DL prediction.

    The natural fit for a volumetric model (a 3-D CNN/ResNet-style
    architecture predicting a full tensor volume): the resulting file is
    the same artifact :func:`pycsamt.format.adapters.modem3d.modem3d_to_pcsf`
    produces for a real ModEM 3-D run.

    Parameters
    ----------
    resistivity : array-like, shape (n_z, n_y, n_x)
        The model's predicted resistivity (see
        :attr:`~pycsamt.format.schema.Grid3DGeometry.resistivity_shape` for
        this axis order -- ModEM's own native convention).
    x, y, z : array-like
        Cell-centre coordinates (metres).
    encoding : {"linear", "log10", "ln"}, default "linear"
        See :func:`grid2d_to_pcsf`'s identical parameter.
    x_nodes, y_nodes, z_nodes, origin, rotation_deg, n_air
        Forwarded to :class:`~pycsamt.format.schema.Grid3DGeometry`.
        *origin*/*rotation_deg* are this geometry's "offset" -- real-world
        placement and bearing of an otherwise locally-referenced volume.
    uncertainty, sensitivity, stations, topography, station_names,
    station_x, station_y, station_z, station_elevations, station_lonlat,
    topo, epsg, utm_zone, latlon, on_mismatch, history, provenance, survey,
    source_backend, created_by, crs, description, metadata
        Same meaning as :func:`grid2d_to_pcsf`'s identical parameters
        (``station_y``/``station_z`` additionally accepted here, matching
        a native 3-D station table).

    Returns
    -------
    PCSFModel

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format.adapters.generic import grid3d_to_pcsf
    >>> model = grid3d_to_pcsf(
    ...     resistivity=np.full((2, 2, 2), 100.0),
    ...     x=np.array([0.0, 100.0]), y=np.array([0.0, 100.0]),
    ...     z=np.array([10.0, 50.0]), source_backend="resnet",
    ... )
    >>> model.kind
    'grid3d'
    """
    geometry = Grid3DGeometry(
        x=x, y=y, z=z, x_nodes=x_nodes, y_nodes=y_nodes, z_nodes=z_nodes,
        origin=None if origin is None else np.asarray(origin, dtype=float),
        rotation_deg=float(rotation_deg), n_air=int(n_air),
    )
    resistivity_lin, resistivity_native, native_encoding = _apply_encoding(
        resistivity, encoding
    )
    stations, topography = _resolve_spatial(
        stations=stations, topography=topography,
        station_names=station_names, station_x=station_x,
        station_y=station_y, station_z=station_z,
        station_elevations=station_elevations,
        station_lonlat=station_lonlat, topo=topo, epsg=epsg,
        utm_zone=utm_zone, latlon=latlon, on_mismatch=on_mismatch,
        caller="grid3d_to_pcsf",
    )

    return PCSFModel(
        geometry=geometry,
        resistivity=resistivity_lin,
        resistivity_native=resistivity_native,
        resistivity_native_encoding=native_encoding,
        uncertainty=uncertainty,
        sensitivity=sensitivity,
        stations=stations,
        topography=topography,
        survey=_survey_to_dict(survey),
        history=dict(history or {}),
        source_backend=source_backend,
        created_by=created_by,
        crs=crs,
        description=description,
        metadata=_provenance_to_metadata(provenance, metadata),
    )


def _node_values_to_triangles(
    connectivity: np.ndarray, node_values: np.ndarray
) -> np.ndarray:
    """Per-triangle value as the arithmetic mean of its 3 vertex values.

    A deliberate, documented, reproducible convention -- not a hidden
    default -- for turning a graph-model's native per-node output into
    PCSF's canonical per-cell ``resistivity``.
    """
    return node_values[connectivity].mean(axis=1)


def _expand_region_resistivity(
    resistivity_by_region: np.ndarray, region_ids: np.ndarray
) -> np.ndarray:
    """0-based region-id expansion.

    Unlike :func:`pycsamt.format.adapters.mare2dem.mare2dem_to_pcsf`'s own
    equivalent (which follows MARE2DEM's native 1-based region numbering),
    ``mesh_to_pcsf`` has no external file format to match, so its
    ``region_ids`` are plain 0-based indices into *resistivity_by_region*.
    """
    ids = np.asarray(region_ids)
    n_regions = np.asarray(resistivity_by_region).shape[0]
    if ids.size and (ids.min() < 0 or ids.max() >= n_regions):
        raise ValueError(
            f"mesh region_ids span [{ids.min()}, {ids.max()}], outside "
            f"the resistivity_by_region table's [0, {n_regions - 1}] range"
        )
    return np.asarray(resistivity_by_region)[ids.astype(np.int64)]


def mesh_to_pcsf(
    nodes: Any,
    connectivity: Any,
    *,
    region_ids: Any = None,
    resistivity: Any = None,
    resistivity_by_node: Any = None,
    resistivity_by_region: Any = None,
    encoding: str = "linear",
    plane: str = "xz",
    uncertainty: Any = None,
    sensitivity: Any = None,
    stations: StationTable | None = None,
    topography: TopographyPerStation | TopographyRaster | None = None,
    station_names: Sequence[str] | None = None,
    station_x: Sequence[float] | None = None,
    station_y: Sequence[float] | None = None,
    station_z: Sequence[float] | None = None,
    station_elevations: Mapping[str, float] | None = None,
    station_lonlat: Mapping[str, tuple[float, float]] | None = None,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
    history: Mapping[str, Any] | None = None,
    provenance: ModelProvenance | Mapping[str, Any] | None = None,
    survey: Any | Mapping[str, Any] | None = None,
    source_backend: str = "ai",
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
    metadata: Mapping[str, Any] | None = None,
) -> PCSFModel:
    r"""Build a ``mesh_unstructured`` :class:`PCSFModel` from a bare
    AI/DL prediction on a triangular mesh.

    The natural fit for a graph-based model (a GCN predicting one value
    per mesh vertex/graph node), but also works for any model that already
    predicts per-triangle or per-region values on a mesh it did not itself
    generate (e.g. a caller-supplied mesh from
    :func:`pycsamt.models.mare2dem.tri_mesh.tri_mesh_from_poly`).

    Parameters
    ----------
    nodes : array-like, shape (n, 2) or (n, 3)
        Mesh node coordinates, metres.
    connectivity : array-like, shape (m, 3), int
        Triangle node indices.
    region_ids : array-like, shape (m,), int, optional
        Per-triangle region id. Defaults to all-zeros (a single, generic
        region -- "no region structure") when omitted.
    resistivity, resistivity_by_node, resistivity_by_region : array-like, optional
        At least one is required. *resistivity* is per-triangle/per-cell
        (shape ``(m,)``); *resistivity_by_node* is per-vertex (shape
        ``(n,)``, the natural GCN output); *resistivity_by_region* is a
        compact per-region table (shape ``(n_regions,)``). Precedence when
        more than one is given: *resistivity* > *resistivity_by_node* >
        *resistivity_by_region* -- the first present is treated as the
        canonical source (*encoding* applies to it), and the canonical
        per-triangle ``resistivity`` is derived from it: taken as-is, or
        averaged from its 3 vertex values per triangle (documented
        arithmetic mean, see :func:`_node_values_to_triangles`), or
        expanded via 0-based *region_ids* (see
        :func:`_expand_region_resistivity`). Any of the other two tables
        also given alongside the canonical source is stored as-is,
        assumed already linear ohm.m.
    encoding : {"linear", "log10", "ln"}, default "linear"
        Encoding of whichever of the three resistivity arrays above is
        the canonical source (see precedence above).
    plane : {"xz", "xy", "3d"}, default "xz"
        Forwarded to :class:`~pycsamt.format.schema.UnstructuredMeshGeometry`.
    uncertainty, sensitivity : array-like, optional
        Same shape as the canonical per-triangle ``resistivity``.
    stations, topography, station_names, station_x, station_y, station_z,
    station_elevations, station_lonlat, topo, epsg, utm_zone, latlon,
    on_mismatch, history, provenance, survey, source_backend, created_by,
    crs, description, metadata
        Same meaning as :func:`grid3d_to_pcsf`'s identical parameters. A
        mesh's own node coordinates already carry real position -- these
        are for naming/geo-referencing discrete *receivers* distinct from
        the mesh nodes (mirrors
        :func:`pycsamt.format.adapters.mare2dem.mare2dem_to_pcsf`'s own
        caller-supplied ``stations``).

    Returns
    -------
    PCSFModel

    Raises
    ------
    ValueError
        If none of *resistivity*/*resistivity_by_node*/*resistivity_by_region*
        is given, or *region_ids* falls outside the resistivity-by-region
        table's range.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format.adapters.generic import mesh_to_pcsf
    >>> nodes = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
    >>> connectivity = np.array([[0, 1, 2], [1, 3, 2]])
    >>> node_values = np.array([1.0, 2.0, 3.0, 4.0])  # log10, one per node
    >>> model = mesh_to_pcsf(
    ...     nodes, connectivity, resistivity_by_node=node_values,
    ...     encoding="log10", source_backend="gcn",
    ... )
    >>> model.kind
    'mesh_unstructured'
    """
    nodes_arr = np.asarray(nodes, dtype=float)
    connectivity_arr = np.asarray(connectivity, dtype=np.int64)
    n_triangles = connectivity_arr.shape[0]
    region_ids_arr = (
        np.zeros(n_triangles, dtype=np.int32)
        if region_ids is None
        else np.asarray(region_ids, dtype=np.int32)
    )

    geometry = UnstructuredMeshGeometry(
        nodes=nodes_arr, connectivity=connectivity_arr,
        region_ids=region_ids_arr, plane=plane,
    )

    sources = (
        ("resistivity", resistivity),
        ("resistivity_by_node", resistivity_by_node),
        ("resistivity_by_region", resistivity_by_region),
    )
    given = [name for name, value in sources if value is not None]
    if not given:
        raise ValueError(
            "mesh_to_pcsf needs at least one of 'resistivity', "
            "'resistivity_by_node', or 'resistivity_by_region'"
        )
    primary = given[0]
    primary_values = dict(sources)[primary]
    primary_lin, resistivity_native, native_encoding = _apply_encoding(
        primary_values, encoding
    )

    if primary == "resistivity":
        resistivity_lin = primary_lin
    elif primary == "resistivity_by_node":
        resistivity_lin = _node_values_to_triangles(connectivity_arr, primary_lin)
    else:
        resistivity_lin = _expand_region_resistivity(primary_lin, region_ids_arr)

    resistivity_by_node_out = (
        primary_lin if primary == "resistivity_by_node"
        else (np.asarray(resistivity_by_node, dtype=float) if resistivity_by_node is not None else None)
    )
    resistivity_by_region_out = (
        primary_lin if primary == "resistivity_by_region"
        else (np.asarray(resistivity_by_region, dtype=float) if resistivity_by_region is not None else None)
    )

    stations, topography = _resolve_spatial(
        stations=stations, topography=topography,
        station_names=station_names, station_x=station_x,
        station_y=station_y, station_z=station_z,
        station_elevations=station_elevations,
        station_lonlat=station_lonlat, topo=topo, epsg=epsg,
        utm_zone=utm_zone, latlon=latlon, on_mismatch=on_mismatch,
        caller="mesh_to_pcsf",
    )

    return PCSFModel(
        geometry=geometry,
        resistivity=resistivity_lin,
        resistivity_native=resistivity_native,
        resistivity_native_encoding=native_encoding,
        resistivity_by_node=resistivity_by_node_out,
        resistivity_by_region=resistivity_by_region_out,
        uncertainty=uncertainty,
        sensitivity=sensitivity,
        stations=stations,
        topography=topography,
        survey=_survey_to_dict(survey),
        history=dict(history or {}),
        source_backend=source_backend,
        created_by=created_by,
        crs=crs,
        description=description,
        metadata=_provenance_to_metadata(provenance, metadata),
    )
