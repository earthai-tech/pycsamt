# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM -> PCSF adapter (Phase 4 of the PCSF format plan).

Converts a completed :class:`pycsamt.models.mare2dem.results.InversionResult`
into a backend-neutral :class:`~pycsamt.format.schema.PCSFModel` with a
``mesh_unstructured`` geometry — MARE2DEM's real triangular FEM mesh is
preserved as-is, never forced onto a rectilinear grid (unlike
``pycsamt.interp._base.ResistivityModel.from_any``, which explicitly
raises ``NotImplementedError`` for MARE2DEM today).

Unlike Occam2D/ModEM, ``InversionResult`` never loads a mesh at all —
only the per-region resistivity table (see
``pycsamt/models/mare2dem/results.py: InversionResult._scan``, which
scans for ``.log``/``.resistivity``/``.emdata`` files but no
``.poly``/``.node``/``.ele`` mesh). Building the actual
:class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh` is therefore a
separate, already-existing concern
(:func:`pycsamt.models.mare2dem.tri_mesh.build_survey_mesh` /
:func:`~pycsamt.models.mare2dem.tri_mesh.tri_mesh_from_poly`, both of
which need a real Triangle run) — this adapter's job is only to fold a
mesh the caller already has together with the result's resistivity
table, matching the geometry/resistivity split every other adapter in
this package follows.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Mapping

import numpy as np

from ..schema import PCSFModel, StationTable, UnstructuredMeshGeometry

if TYPE_CHECKING:
    from pycsamt.forward.maxwell.contracts_tri import TriMesh
    from pycsamt.models.mare2dem.results import InversionResult

__all__ = ["mare2dem_to_pcsf"]


def _survey_to_dict(survey: Any | Mapping[str, Any] | None) -> dict[str, Any]:
    if survey is None:
        return {}
    to_dict = getattr(survey, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())
    return dict(survey)


def _expand_region_resistivity(
    resistivity_by_region: np.ndarray, region_ids: np.ndarray
) -> np.ndarray:
    """Map each mesh triangle's 1-based region id to its resistivity.

    Parameters
    ----------
    resistivity_by_region : ndarray (n_regions,)
        Row ``i`` is the resistivity for region ``i + 1`` (MARE2DEM's
        ``.resistivity`` file is a 1-based row-per-region table).
    region_ids : ndarray (n_triangles,)
        Each mesh triangle's region id, as written by Triangle's own
        per-element region attribute.

    Raises
    ------
    ValueError
        If a region id falls outside ``[1, n_regions]``.
    """
    ids = np.asarray(region_ids)
    n_regions = resistivity_by_region.shape[0]
    if ids.size and (ids.min() < 1 or ids.max() > n_regions):
        raise ValueError(
            f"mesh region_ids span [{ids.min()}, {ids.max()}], outside "
            f"the resistivity table's [1, {n_regions}] range"
        )
    return resistivity_by_region[ids.astype(np.int64) - 1]


def mare2dem_to_pcsf(
    result: InversionResult,
    mesh: TriMesh,
    *,
    stations: StationTable | None = None,
    survey: Any | Mapping[str, Any] | None = None,
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
) -> PCSFModel:
    r"""Convert a MARE2DEM :class:`InversionResult` to a :class:`PCSFModel`.

    Parameters
    ----------
    result : InversionResult
        A loaded MARE2DEM working directory with a readable
        ``.resistivity`` file (``result.model``).
    mesh : TriMesh
        The real triangular mesh paired with *result* — the caller
        provides it (e.g. from
        :func:`pycsamt.models.mare2dem.tri_mesh.tri_mesh_from_poly`, or
        an in-process constrained triangulation of the run's own
        ``.poly`` PSLG) since ``InversionResult`` does not load one
        itself. ``mesh.region_ids`` must be populated and must match
        the 1-based region numbering of ``result.model.resistivity``.
    stations : StationTable, optional
        MARE2DEM's ``EMData``/``EMDataFile`` receiver geometry has no
        single reliable per-point name across its MT/CSEM/DC variants,
        so station identity is not auto-derived here — pass a
        pre-built table when the caller already has one (e.g. via
        :func:`pycsamt.models.mare2dem.geom.area_of_interest.survey_points`
        plus the caller's own naming). Set its ``lon``/``lat`` (e.g.
        from the same EDI headers the receiver positions were derived
        from) to make the file self-sufficiently geo-referenced too —
        the same field :func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf`'s
        ``station_lonlat`` and :func:`pycsamt.format.adapters.modem3d.modem3d_to_pcsf`'s
        own ``GG_Lat``/``GG_Lon`` passthrough populate.
    survey : SurveyMeta or mapping, optional
        Survey-level metadata, stored the same way as in the other
        adapters in this package.
    created_by, crs, description : str, optional
        Passed straight through to :class:`PCSFModel`.

    Returns
    -------
    PCSFModel
        ``geometry.kind == "mesh_unstructured"``. Canonical
        ``resistivity`` is the per-triangle array expanded from
        ``result.model.resistivity`` via ``mesh.region_ids``; the
        compact per-region table itself is kept in
        :attr:`PCSFModel.resistivity_by_region`. MARE2DEM's
        ``.resistivity`` file is already linear ohm.m (confirmed
        against a real compiled binary — see
        :class:`~pycsamt.models.mare2dem.iotools.resistivity.ResistivityFile`'s
        own docstring), so no ``resistivity_native``/encoding
        conversion applies here.

    Raises
    ------
    ValueError
        If *result* has no resistivity model, *mesh* has no
        ``region_ids``, a region id falls outside the resistivity
        table's range, or the source is anisotropic (only isotropic
        MARE2DEM models are supported by this adapter today).

    Examples
    --------
    >>> from pycsamt.models.mare2dem.results import InversionResult
    >>> from pycsamt.models.mare2dem.tri_mesh import tri_mesh_from_poly
    >>> from pycsamt.format.adapters.mare2dem import mare2dem_to_pcsf
    >>> from pycsamt.format import write_pcsf
    >>> result = InversionResult("mare2dem_run")  # doctest: +SKIP
    >>> mesh = tri_mesh_from_poly("mare2dem_run/mesh.1.poly")  # doctest: +SKIP
    >>> model = mare2dem_to_pcsf(result, mesh)  # doctest: +SKIP
    >>> write_pcsf(model, "mare2dem_run.pcsf")  # doctest: +SKIP
    """
    rf = result.model
    if rf is None or rf.resistivity is None or rf.resistivity.size == 0:
        raise ValueError(
            "InversionResult has no resistivity model — ensure the "
            "workdir contains a readable .resistivity file."
        )
    if mesh.region_ids is None:
        raise ValueError(
            "mesh.region_ids is required to fold result.model's "
            "per-region resistivity onto mesh triangles."
        )
    if rf.anisotropy.lower().strip() != "isotropic":
        raise ValueError(
            f"mare2dem_to_pcsf only supports isotropic models today, "
            f"got anisotropy={rf.anisotropy!r}"
        )

    resistivity_by_region = np.asarray(rf.resistivity[:, 0], dtype=float)
    resistivity = _expand_region_resistivity(
        resistivity_by_region, mesh.region_ids
    )

    geometry = UnstructuredMeshGeometry(
        nodes=mesh.nodes_m,
        connectivity=mesh.triangles,
        region_ids=mesh.region_ids,
        plane="xz",
    )

    history: dict[str, np.ndarray] = {}
    if result.log is not None and result.log.iterations:
        recs = result.log.iterations
        history = {
            "iteration": np.array([r.iteration for r in recs], dtype=float),
            "rms": np.array([r.rms for r in recs], dtype=float),
            "roughness": np.array([r.roughness for r in recs], dtype=float),
            "lambda": np.array([r.lambda_ for r in recs], dtype=float),
        }

    metadata: dict[str, Any] = {
        "workdir": str(result.workdir),
        "final_rms": result.final_rms,
        "n_iterations": result.n_iterations,
        "converged": result.converged,
        "n_regions": int(rf.num_regions),
        "anisotropy": rf.anisotropy,
    }

    return PCSFModel(
        geometry=geometry,
        resistivity=resistivity,
        resistivity_by_region=resistivity_by_region,
        stations=stations,
        survey=_survey_to_dict(survey),
        history=history,
        source_backend="mare2dem",
        created_by=created_by,
        crs=crs,
        description=description,
        metadata=metadata,
    )
