# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regridded ``grid2d`` convenience export for ``mesh_unstructured`` models.

Design principle 6 (``SPEC.md`` section 2) deliberately keeps a
``mesh_unstructured`` file's own triangular mesh un-regridded by
default -- MARE2DEM's mesh stays a mesh. A reader wanting a
rectilinear approximation instead (a quick-look plot, a pipeline stage
that only understands ``grid2d``) previously had to write that
interpolation itself; this module closes that gap with one function,
:func:`mesh_to_grid2d`.

It is not a second interpolation scheme: it reuses the exact
:class:`matplotlib.tri.Triangulation` point-location machinery
:func:`pycsamt.map.inversion.load_pcsf_lines` already uses to slice a
per-station curtain through a ``mesh_unstructured`` model, evaluated on
a regular grid instead of at station positions. A query point outside
the mesh returns ``nan``, never an extrapolated or fabricated value --
the same convention that curtain slicer already follows.

The result is explicitly marked synthesized (in ``description`` and
``metadata``), matching design principle 5's rule that a derived view
must never be mistaken for a native inversion output: no
``resistivity_native``, no ``resistivity_by_region`` survive the
regrid, since neither concept applies to an interpolated value.
"""

from __future__ import annotations

import numpy as np

from .schema import Grid2DGeometry, PCSFModel

__all__ = ["mesh_to_grid2d"]


def mesh_to_grid2d(
    model: PCSFModel,
    *,
    nx: int = 200,
    nz: int = 150,
) -> PCSFModel:
    """Regrid a ``mesh_unstructured`` PCSF model onto a rectilinear
    ``grid2d`` approximation.

    Parameters
    ----------
    model : PCSFModel
        A model with ``geometry.kind == "mesh_unstructured"`` and
        per-cell resistivity (shape ``(n_triangles,)`` -- every
        ``pycsamt.format.adapters`` writer already produces this; a
        region-collapsed ``(n_regions,)`` array must be expanded onto
        each triangle's region id first).
    nx, nz : int, default 200, 150
        Grid resolution. Chosen independently of the source mesh's own
        resolution -- there is no natural rectilinear resolution to
        inherit from an unstructured mesh -- the same kind of
        pragmatic, documented default the rest of this format already
        makes (cf. PCSM's row-width cap, the per-station curtain's
        ``n_z`` default) rather than leaving it unbounded.

    Returns
    -------
    PCSFModel
        A new ``grid2d`` model. Grid points outside the source mesh's
        triangulation are ``nan``. ``source_backend``, ``created_by``,
        and ``crs`` are carried over from *model*; ``description`` and
        ``metadata`` record that this is a synthesized regrid, not a
        native inversion output.

    Raises
    ------
    ValueError
        If *model* is not ``mesh_unstructured``, or its resistivity is
        not already per-cell.
    NotImplementedError
        If the mesh's ``plane`` is not ``"xz"`` -- the only plane any
        current adapter produces; a general 3-D mesh has no single
        rectilinear plane to regrid onto.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format.schema import PCSFModel, UnstructuredMeshGeometry
    >>> from pycsamt.format.regrid import mesh_to_grid2d
    >>> geometry = UnstructuredMeshGeometry(
    ...     nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
    ...     connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
    ...     region_ids=np.array([0, 1]),
    ... )
    >>> model = PCSFModel(geometry=geometry, resistivity=np.array([10.0, 20.0]))
    >>> regridded = mesh_to_grid2d(model, nx=5, nz=5)
    >>> regridded.kind
    'grid2d'
    >>> regridded.resistivity.shape
    (5, 5)
    """
    if model.kind != "mesh_unstructured":
        raise ValueError(
            "mesh_to_grid2d expects a mesh_unstructured model, got "
            f"geometry kind {model.kind!r}"
        )
    geo = model.geometry
    if geo.plane != "xz":
        raise NotImplementedError(
            "mesh_to_grid2d only supports plane='xz' meshes (the only "
            f"kind any current adapter produces), got plane={geo.plane!r}"
        )
    connectivity = np.asarray(geo.connectivity, dtype=np.int64)
    if model.resistivity.shape != (connectivity.shape[0],):
        raise ValueError(
            "mesh_to_grid2d needs per-cell mesh_unstructured resistivity "
            f"(shape ({connectivity.shape[0]},)), got shape "
            f"{model.resistivity.shape} -- expand resistivity_by_region "
            "onto each triangle's region id first (every "
            "pycsamt.format.adapters writer already does this)."
        )

    from matplotlib.tri import Triangulation

    nodes = np.asarray(geo.nodes, dtype=float)
    resistivity = np.asarray(model.resistivity, dtype=float)

    triangulation = Triangulation(nodes[:, 0], nodes[:, 1], triangles=connectivity)
    trifinder = triangulation.get_trifinder()

    x = np.linspace(float(nodes[:, 0].min()), float(nodes[:, 0].max()), nx)
    z = np.linspace(float(nodes[:, 1].min()), float(nodes[:, 1].max()), nz)
    x_grid, z_grid = np.meshgrid(x, z)  # both (nz, nx)
    tri_idx = trifinder(x_grid.ravel(), z_grid.ravel())
    grid_rho = np.full(tri_idx.shape, np.nan)
    inside = tri_idx >= 0
    grid_rho[inside] = resistivity[tri_idx[inside]]
    grid_rho = grid_rho.reshape(nz, nx)

    note = "regridded from mesh_unstructured via triangulation point-location (synthesized, not a native inversion output)"
    description = f"{model.description} -- {note}" if model.description else note

    return PCSFModel(
        geometry=Grid2DGeometry(x=x, z=z),
        resistivity=grid_rho,
        source_backend=model.source_backend,
        created_by=model.created_by,
        crs=model.crs,
        description=description,
        survey=dict(model.survey),
        metadata={
            **model.metadata,
            "synthesized": True,
            "regridded_from": "mesh_unstructured",
            "regrid_nx": nx,
            "regrid_nz": nz,
        },
    )
