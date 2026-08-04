# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Build a topography-aware triangular FEM mesh for a survey.

Completes the piece ``builder.py``'s own module docstring flags as
missing: *"The full FEM mesh generation (Triangle ``.poly`` file) is
not yet implemented."* This module composes the already-ported
geometry primitives in :mod:`pycsamt.models.mare2dem.geom` with
:func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle` to produce
a real Triangle-refined mesh, then converts it to the solver-neutral
:class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh` contract so
downstream code (a future ``Mare2DEMAdapter``, AI training, rendering)
never has to touch MARE2DEM's native Triangle file formats directly.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np

from .geom.area_of_interest import estimate_area_of_interest, survey_points
from .geom.simplify import dp_simplify
from .geom.topo import parse_topo
from .iotools.poly import PolyFile, read_triangulation, write_poly
from .triangle_exec import run_triangle

if TYPE_CHECKING:
    from pycsamt.forward.maxwell.contracts_tri import TriMesh

    from .config import Mare2DEMConfig
    from .iotools.emdata import EMDataFile

PathLike = Union[str, Path]

__all__ = ["build_pslg", "build_survey_mesh", "tri_mesh_from_poly"]


def tri_mesh_from_poly(refined_poly_path: PathLike) -> "TriMesh":
    """Convert a Triangle-refined ``.poly`` mesh to a :class:`TriMesh`.

    Parameters
    ----------
    refined_poly_path : path-like
        The ``.1.poly`` file produced by
        :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle` (or
        any ``.poly``/``.node`` path with a companion ``.node``/``.ele``
        pair of the same stem).

    Returns
    -------
    pycsamt.forward.maxwell.contracts_tri.TriMesh
        Validated triangular mesh; ``region_ids`` are Triangle's own
        per-element region attribute (all zero when the mesh was built
        without region constraints).

    Raises
    ------
    FileNotFoundError
        If the companion ``.node``/``.ele`` files are missing.

    Examples
    --------
    Normally called on :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle`'s
    own return value::

        >>> from pycsamt.models.mare2dem.triangle_exec import run_triangle
        >>> refined = run_triangle("mesh.poly")  # doctest: +SKIP
        >>> mesh = tri_mesh_from_poly(refined)  # doctest: +SKIP
    """
    from pycsamt.forward.maxwell.contracts_tri import TriMesh as _TriMesh

    node_path = Path(refined_poly_path).with_suffix(".node")
    nodes, triangles, region_attrs = read_triangulation(node_path)
    return _TriMesh(nodes, triangles, region_ids=region_attrs)


def build_pslg(
    boundary,
    *,
    interior_points=None,
    region_seeds=None,
) -> PolyFile:
    """Assemble a closed-polygon PSLG from boundary and interior points.

    Parameters
    ----------
    boundary : array-like, shape (n, 2)
        Ordered ``(y, z)`` vertices of a single closed outer boundary
        polygon (do not repeat the first point at the end). Consecutive
        vertices become PSLG boundary segments, including the closing
        edge from the last vertex back to the first.
    interior_points : array-like, shape (m, 2), optional
        Extra points Triangle must place a mesh node at exactly (e.g.
        receiver locations), added with no connecting segments.
    region_seeds : array-like, shape (k, 2), optional
        One interior point per resistivity region, written as Triangle
        region-constraint rows (attribute ``i + 1``, no area cap). Omit
        for a single uniform region (``region_ids`` will be all zero).

    Returns
    -------
    PolyFile
        PSLG ready for :func:`~pycsamt.models.mare2dem.iotools.poly.write_poly`.

    Examples
    --------
    >>> box = [[0, 0], [1000, 0], [1000, 500], [0, 500]]
    >>> pf = build_pslg(box)
    >>> pf.n_nodes, pf.n_segments
    (4, 4)
    """
    boundary = np.asarray(boundary, dtype=float)
    if boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 3:
        raise ValueError("boundary must have shape (n, 2) with n >= 3.")
    extra = (
        np.zeros((0, 2))
        if interior_points is None
        else np.asarray(interior_points, dtype=float).reshape(-1, 2)
    )
    nodes = np.vstack([boundary, extra]) if len(extra) else boundary
    n_boundary = len(boundary)

    pf = PolyFile()
    pf.nodes = nodes
    pf.segments = np.column_stack(
        [
            np.arange(1, n_boundary + 1),
            np.r_[np.arange(2, n_boundary + 1), 1],
        ]
    )
    if region_seeds is not None:
        seeds = np.asarray(region_seeds, dtype=float).reshape(-1, 2)
        pf.regions = np.column_stack(
            [seeds, np.arange(1, len(seeds) + 1), np.zeros(len(seeds))]
        )
    return pf


def build_survey_mesh(
    em: "EMDataFile",
    *,
    topo: float | np.ndarray = 0.0,
    depth_below_topo_m: float | None = None,
    simplify_tolerance_m: float | None = None,
    min_angle: float = 30.0,
    max_area: float | None = None,
    workdir: PathLike = ".",
    stem: str = "mesh",
    executable: PathLike | None = None,
    config: "Mare2DEMConfig | None" = None,
) -> "TriMesh":
    """Build a topography-aware triangular mesh for a survey.

    Composes :func:`~pycsamt.models.mare2dem.geom.area_of_interest.estimate_area_of_interest`,
    :func:`~pycsamt.models.mare2dem.geom.topo.parse_topo`,
    :func:`~pycsamt.models.mare2dem.geom.simplify.dp_simplify`,
    :func:`build_pslg`, and
    :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle`, then
    converts the refined mesh via :func:`tri_mesh_from_poly`.

    Parameters
    ----------
    em : EMDataFile
        Survey data supplying receiver/transmitter positions, used both
        for the area-of-interest estimate and as required PSLG interior
        points so Triangle places a node exactly at every receiver.
    topo : float or array-like, shape (n, 2), default=0.0
        Topography, forwarded to
        :func:`~pycsamt.models.mare2dem.geom.topo.parse_topo`. A scalar
        gives a flat surface at that depth.
    depth_below_topo_m : float, optional
        Bottom-boundary depth. Defaults to the estimated area-of-interest
        ``zlim[1]``.
    simplify_tolerance_m : float, optional
        Douglas-Peucker tolerance applied to the sampled topography
        polyline before meshing. Left unsimplified when ``None``.
    min_angle, max_area : float, optional
        Forwarded to :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle`.
    workdir : path-like, default="."
        Directory the ``.poly``/``.node``/``.ele`` files are written to.
    stem : str, default="mesh"
        Base filename (without extension) for the written mesh files.
    executable, config :
        Forwarded to :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle`.

    Returns
    -------
    pycsamt.forward.maxwell.contracts_tri.TriMesh
        Triangular mesh spanning the survey's area of interest, with a
        node at every receiver position.

    Raises
    ------
    ValueError
        If ``em`` has no receiver/transmitter geometry to mesh.
    pycsamt.models.mare2dem.triangle_exec.TriangleRunError
        If no Triangle executable is found or refinement fails.

    Examples
    --------
    Requires a real Triangle executable, so this is illustrative only::

        >>> mesh = build_survey_mesh(em, topo=0.0)  # doctest: +SKIP
    """
    receivers = survey_points(em)
    if len(receivers) == 0:
        raise ValueError("em has no receiver/transmitter geometry to mesh.")
    ylim, zlim = estimate_area_of_interest(em)

    n_profile = max(len(receivers) * 4, 40)
    y_profile = np.linspace(ylim[0], ylim[1], n_profile)
    z_top, _, _ = parse_topo(topo, y_profile)
    topo_line = np.column_stack([y_profile, z_top])
    if simplify_tolerance_m is not None:
        topo_line = dp_simplify(topo_line, simplify_tolerance_m)

    bottom = (
        float(zlim[1])
        if depth_below_topo_m is None
        else float(depth_below_topo_m)
    )
    boundary = np.vstack(
        [topo_line, [[ylim[1], bottom], [ylim[0], bottom]]]
    )

    pf = build_pslg(boundary, interior_points=receivers)
    workdir_path = Path(workdir)
    workdir_path.mkdir(parents=True, exist_ok=True)
    poly_path = write_poly(pf, workdir_path / f"{stem}.poly")

    refined = run_triangle(
        poly_path,
        executable=executable,
        min_angle=min_angle,
        max_area=max_area,
        config=config,
    )
    return tri_mesh_from_poly(refined)
