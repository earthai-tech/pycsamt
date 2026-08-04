# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reader and writer for Triangle ``.poly`` PSLG mesh files.

Port of:
  * ``m2d_readPoly.m``
  * ``m2d_writePoly.m``

The ``.poly`` format is J. R. Shewchuk's Planar Straight-Line Graph
(PSLG) format used by the Triangle mesh generator that is bundled
inside MARE2DEM.  A MARE2DEM ``.poly`` file describes the boundary
geometry, interior features, holes, and area constraints for the 2-D
finite-element mesh.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "PolyFile",
    "read_poly",
    "read_triangulation",
    "write_poly",
    "write_triangulation",
]


def _read_int(tokens: list[str], idx: int, default: int = 0) -> int:
    try:
        return int(tokens[idx])
    except (IndexError, ValueError):
        return default


def _read_float(tokens: list[str], idx: int, default: float = 0.0) -> float:
    try:
        return float(tokens[idx])
    except (IndexError, ValueError):
        return default


class PolyFile:
    """Contents of one Triangle ``.poly`` PSLG file.

    Attributes
    ----------
    nodes : numpy.ndarray, shape (n_nodes, 2)
        Node (x, y) coordinates.  A MARE2DEM ``.poly`` file uses a
        right-handed coordinate system where the second coordinate is
        depth (z positive down).
    node_attributes : numpy.ndarray, shape (n_nodes, n_attr) or shape (0,)
        Optional per-node attribute values.
    node_boundary_markers : numpy.ndarray, shape (n_nodes,) or shape (0,)
        Optional per-node boundary marker integers.
    segments : numpy.ndarray, shape (n_segs, 2) with int dtype
        Segment endpoint node-index pairs (1-based).
    segment_markers : numpy.ndarray, shape (n_segs,) or shape (0,)
        Optional per-segment boundary marker integers.
    holes : numpy.ndarray, shape (n_holes, 2) or shape (0, 2)
        Hole point coordinates.
    regions : numpy.ndarray, shape (n_reg, 4) or shape (0, 4)
        Region attribute and area-constraint points.
        Columns: x, y, attribute, max_area.
    """

    def __init__(self) -> None:
        self.nodes: np.ndarray = np.empty((0, 2), dtype=float)
        self.node_attributes: np.ndarray = np.empty((0,), dtype=float)
        self.node_boundary_markers: np.ndarray = np.empty((0,), dtype=int)
        self.segments: np.ndarray = np.empty((0, 2), dtype=int)
        self.segment_markers: np.ndarray = np.empty((0,), dtype=int)
        self.holes: np.ndarray = np.empty((0, 2), dtype=float)
        self.regions: np.ndarray = np.empty((0, 4), dtype=float)

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_segments(self) -> int:
        return len(self.segments)

    @property
    def n_holes(self) -> int:
        return len(self.holes)

    @property
    def n_regions(self) -> int:
        return len(self.regions)

    def __repr__(self) -> str:
        return (
            f"PolyFile(n_nodes={self.n_nodes}, "
            f"n_segments={self.n_segments}, "
            f"n_holes={self.n_holes}, "
            f"n_regions={self.n_regions})"
        )


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------


def _next_data_line(lines: list[str], pos: int) -> tuple[list[str], int]:
    """Return the next non-comment, non-empty token list and updated pos."""
    while pos < len(lines):
        line = lines[pos].strip()
        pos += 1
        if line and not line.startswith("#"):
            return line.split(), pos
    return [], pos


def read_poly(path: str | Path) -> PolyFile:
    """Read a Triangle ``.poly`` PSLG file.

    Port of ``m2d_readPoly.m``.

    Parameters
    ----------
    path : path-like
        File to read.  If the node count is zero (triangulation already
        done), the function automatically looks for the companion
        ``.node`` and ``.ele`` files in the same directory.

    Returns
    -------
    PolyFile
        Parsed PSLG contents.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.poly import read_poly
    >>> poly = read_poly("mare2dem.poly")
    >>> poly.n_nodes
    4812
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    pf = PolyFile()
    lines = path.read_text(errors="replace").splitlines()
    pos = 0

    # --- Node header ---
    toks, pos = _next_data_line(lines, pos)
    n_nodes = _read_int(toks, 0)
    n_attr = _read_int(toks, 2)
    n_bmark = _read_int(toks, 3)

    if n_nodes == 0:
        # nodes are in a companion .node file; triangulation done
        stem = path.with_suffix("")
        node_path = stem.with_suffix(".node")
        ele_path = stem.with_suffix(".ele")

        if node_path.exists():
            nlines = node_path.read_text(errors="replace").splitlines()
            ntoks, npos = _next_data_line(nlines, 0)
            n_v = _read_int(ntoks, 0)
            n_attr2 = _read_int(ntoks, 2)
            n_bm2 = _read_int(ntoks, 3)
            ncols = 3 + n_attr2 + n_bm2
            node_rows: list[list[float]] = []
            while len(node_rows) < n_v:
                t, npos = _next_data_line(nlines, npos)
                if t:
                    node_rows.append([float(v) for v in t[:ncols]])
            arr = np.array(node_rows, dtype=float)
            pf.nodes = arr[:, 1:3]
            if n_attr2 > 0:
                pf.node_attributes = arr[:, 3 : 3 + n_attr2]
            if n_bm2 > 0:
                pf.node_boundary_markers = arr[:, 3 + n_attr2].astype(int)

        if ele_path.exists():
            elines = ele_path.read_text(errors="replace").splitlines()
            etoks, epos = _next_data_line(elines, 0)
            n_tri = _read_int(etoks, 0)
            n_dpt = _read_int(etoks, 1)
            n_ea = _read_int(etoks, 2)
            ele_rows: list[list[int]] = []
            while len(ele_rows) < n_tri:
                t, epos = _next_data_line(elines, epos)
                if t:
                    ele_rows.append([int(v) for v in t[: 1 + n_dpt + n_ea]])
            if ele_rows:
                arr = np.array(ele_rows, dtype=int)
                # store as segments placeholder (not standard, but useful)
                pf.segments = arr[:, 1 : 1 + n_dpt]
    else:
        # normal .poly file — read nodes
        ncols = 3 + n_attr + n_bmark
        node_rows = []
        while len(node_rows) < n_nodes:
            t, pos = _next_data_line(lines, pos)
            if t:
                node_rows.append([float(v) for v in t[:ncols]])
        arr = np.array(node_rows, dtype=float)
        pf.nodes = arr[:, 1:3]
        if n_attr > 0:
            pf.node_attributes = arr[:, 3 : 3 + n_attr]
        if n_bmark > 0:
            pf.node_boundary_markers = arr[:, 3 + n_attr].astype(int)

    # --- Segments ---
    toks, pos = _next_data_line(lines, pos)
    n_segs = _read_int(toks, 0)
    n_seg_bmark = _read_int(toks, 1)
    seg_ncols = 3 + n_seg_bmark
    seg_rows: list[list[int]] = []
    while len(seg_rows) < n_segs:
        t, pos = _next_data_line(lines, pos)
        if t:
            seg_rows.append([int(float(v)) for v in t[:seg_ncols]])
    if seg_rows:
        arr_s = np.array(seg_rows, dtype=int)
        pf.segments = arr_s[:, 1:3]
        if n_seg_bmark > 0:
            pf.segment_markers = arr_s[:, 3]

    # --- Holes ---
    toks, pos = _next_data_line(lines, pos)
    n_holes = _read_int(toks, 0)
    if n_holes > 0:
        hole_rows: list[list[float]] = []
        while len(hole_rows) < n_holes:
            t, pos = _next_data_line(lines, pos)
            if t:
                hole_rows.append([float(v) for v in t[1:3]])
        pf.holes = np.array(hole_rows, dtype=float)

    # --- Region constraints ---
    toks, pos = _next_data_line(lines, pos)
    n_regions = _read_int(toks, 0) if toks else 0
    if n_regions > 0:
        reg_rows: list[list[float]] = []
        while len(reg_rows) < n_regions:
            t, pos = _next_data_line(lines, pos)
            if t:
                reg_rows.append([float(v) for v in t[1:5]])
        pf.regions = np.array(reg_rows, dtype=float)

    return pf


def read_triangulation(
    node_path: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a Triangle ``.node``/``.ele`` pair into plain arrays.

    Shared by :func:`~pycsamt.models.mare2dem.tri_mesh.tri_mesh_from_poly`
    and :class:`~pycsamt.models.mare2dem.plot.PlotModel` so there is one
    parser for Triangle's element-file format, not two.

    Parameters
    ----------
    node_path : path-like
        The ``.node`` file. Its companion ``.ele`` file (same stem, same
        directory) is read alongside it.

    Returns
    -------
    nodes : numpy.ndarray, shape (n_nodes, 2)
        Node ``(x, y)`` coordinates.
    triangles : numpy.ndarray of int, shape (n_triangles, 3)
        0-based node-index connectivity, normalized to 0-based regardless
        of whether the source files used Triangle's ``-z`` (0-based) or
        default (1-based) numbering -- detected from the ``.node`` file's
        own first index.
    region_attrs : numpy.ndarray of int, shape (n_triangles,)
        Triangle's per-element region attribute (its 4th ``.ele`` column,
        present for ``-A``/region-constrained runs), or all zeros when the
        ``.ele`` file carries no region attribute column.

    Raises
    ------
    FileNotFoundError
        If ``node_path`` or its companion ``.ele`` file is missing.

    Examples
    --------
    Normally called on files written by
    :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle`, not
    constructed by hand.
    """
    node_path = Path(node_path)
    ele_path = node_path.with_suffix(".ele")
    if not node_path.exists():
        raise FileNotFoundError(f"node file not found: {node_path}")
    if not ele_path.exists():
        raise FileNotFoundError(f"companion .ele file not found: {ele_path}")

    nlines = node_path.read_text(errors="replace").splitlines()
    ntoks, npos = _next_data_line(nlines, 0)
    n_nodes = _read_int(ntoks, 0)
    node_rows: list[list[float]] = []
    node_indices: list[int] = []
    while len(node_rows) < n_nodes:
        t, npos = _next_data_line(nlines, npos)
        if t:
            node_indices.append(_read_int(t, 0))
            node_rows.append([float(v) for v in t[1:3]])
    nodes = np.array(node_rows, dtype=float).reshape(-1, 2)
    base = min(node_indices) if node_indices else 0

    elines = ele_path.read_text(errors="replace").splitlines()
    etoks, epos = _next_data_line(elines, 0)
    n_tri = _read_int(etoks, 0)
    n_dpt = _read_int(etoks, 1, 3)
    n_ea = _read_int(etoks, 2)
    tri_rows: list[list[int]] = []
    region_rows: list[int] = []
    while len(tri_rows) < n_tri:
        t, epos = _next_data_line(elines, epos)
        if t:
            tri_rows.append(
                [int(float(v)) - base for v in t[1 : 1 + n_dpt]]
            )
            region_rows.append(
                int(float(t[1 + n_dpt]))
                if n_ea > 0 and len(t) > 1 + n_dpt
                else 0
            )
    triangles = np.array(tri_rows, dtype=int).reshape(-1, n_dpt if n_tri else 3)
    region_attrs = np.array(region_rows, dtype=int)
    return nodes, triangles, region_attrs


def _triangle_neighbors(triangles: np.ndarray) -> np.ndarray:
    """Return ``(n_triangles, 3)`` 0-based neighbor indices, ``-1`` at boundaries.

    Column ``k`` is the triangle sharing the edge *opposite* vertex ``k``
    (Triangle's own ``.neigh`` convention -- neighbor 0 is opposite corner
    0, etc.), matching ``triangle.c``'s ``writeneighbors``.

    Raises
    ------
    ValueError
        If any edge is shared by more than two triangles (non-manifold).
    """
    n = len(triangles)
    neighbors = np.full((n, 3), -1, dtype=int)
    edge_map: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for ti, tri in enumerate(triangles):
        v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
        for local_idx, (a, b) in enumerate(((v1, v2), (v2, v0), (v0, v1))):
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append((ti, local_idx))
    for key, entries in edge_map.items():
        if len(entries) == 2:
            (t1, l1), (t2, l2) = entries
            neighbors[t1, l1] = t2
            neighbors[t2, l2] = t1
        elif len(entries) > 2:
            raise ValueError(
                f"non-manifold mesh: edge {key} is shared by "
                f"{len(entries)} triangles (expected at most 2)."
            )
    return neighbors


def write_triangulation(
    nodes: np.ndarray,
    triangles: np.ndarray,
    region_attrs: np.ndarray | None,
    node_path: str | Path,
) -> Path:
    """Write a ``.node``/``.ele``/``.neigh`` triple plus a companion stub ``.poly``.

    Inverse of :func:`read_triangulation`: hands an already-triangulated
    mesh (e.g. a :class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh`)
    to MARE2DEM directly, without re-invoking Triangle. Node/element
    indices are written 1-based, matching the convention MARE2DEM itself
    reads (Triangle's default, distinct from the ``-z`` 0-based output
    :func:`~pycsamt.models.mare2dem.triangle_exec.run_triangle` requests
    for its own refined mesh).

    The ``.neigh`` file (triangle-triangle adjacency, normally produced by
    Triangle's own ``-n`` switch) is computed here directly, since we are
    bypassing Triangle entirely. Confirmed against a real compiled
    MARE2DEM binary that this file is required, not optional -- MARE2DEM's
    ``readPoly`` reads it unconditionally and stops with "Did you remember
    to include the 'n' flag when calling Triangle?" when it is missing.

    Parameters
    ----------
    nodes : array-like, shape (n_nodes, 2)
        Node ``(x, y)`` coordinates.
    triangles : array-like of int, shape (n_triangles, 3)
        0-based node-index connectivity.
    region_attrs : array-like of int, shape (n_triangles,), optional
        Per-triangle region attribute. Defaults to all ``1`` when omitted.
    node_path : path-like
        Destination ``.node`` file. The companion ``.ele``/``.neigh``/
        ``.poly`` files are written alongside it with the same stem.

    Returns
    -------
    pathlib.Path
        The written ``.node`` path.

    Raises
    ------
    ValueError
        If the connectivity is non-manifold (an edge shared by more than
        two triangles), since no valid ``.neigh`` file can describe that.

    Examples
    --------
    >>> nodes = [[0, 0], [1, 0], [0, 1]]
    >>> triangles = [[0, 1, 2]]
    >>> path = write_triangulation(nodes, triangles, None, "mesh.node")
    >>> read_triangulation(path)[1].tolist()
    [[0, 1, 2]]
    """
    nodes = np.asarray(nodes, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    attrs = (
        np.ones(len(triangles), dtype=int)
        if region_attrs is None
        else np.asarray(region_attrs, dtype=int)
    )
    neighbors = _triangle_neighbors(triangles)
    node_path = Path(node_path)
    ele_path = node_path.with_suffix(".ele")
    neigh_path = node_path.with_suffix(".neigh")
    poly_path = node_path.with_suffix(".poly")
    node_path.parent.mkdir(parents=True, exist_ok=True)

    with node_path.open("w") as fh:
        fh.write(f"{len(nodes)} 2 0 0\n")
        for i, (x, y) in enumerate(nodes):
            fh.write(f"{i + 1} {x:.16g} {y:.16g}\n")

    with ele_path.open("w") as fh:
        fh.write(f"{len(triangles)} 3 1\n")
        for i, (tri, attr) in enumerate(zip(triangles, attrs)):
            fh.write(
                f"{i + 1} {tri[0] + 1} {tri[1] + 1} {tri[2] + 1} {int(attr)}\n"
            )

    with neigh_path.open("w") as fh:
        fh.write(f"{len(triangles)} 3\n")
        for i, row in enumerate(neighbors):
            vals = [(-1 if v < 0 else v + 1) for v in row]
            fh.write(f"{i + 1} {vals[0]} {vals[1]} {vals[2]}\n")

    # Zero-node header: nodes/elements live in the companion files, the
    # same convention Triangle's own "-p" refined output uses and
    # read_poly()/read_triangulation() already understand.
    with poly_path.open("w") as fh:
        fh.write("0 2 0 0\n0 0\n0\n0\n")

    return node_path


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------


def write_poly(
    pf: PolyFile,
    path: str | Path,
) -> Path:
    """Write a :class:`PolyFile` to *path*.

    Port of ``m2d_writePoly.m``.

    Parameters
    ----------
    pf : PolyFile
        PSLG data to write.
    path : path-like
        Destination ``.poly`` file.

    Returns
    -------
    pathlib.Path
        Path of the written file.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.poly import write_poly
    >>> write_poly(pf, "mare2dem.poly")
    PosixPath('mare2dem.poly')
    """
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)

    nodes = pf.nodes
    n_nodes = len(nodes)
    has_attr = pf.node_attributes is not None and len(pf.node_attributes) > 0
    pf.node_attributes.shape[
        1
    ] if has_attr and pf.node_attributes.ndim == 2 else 0
    n_bmark = (
        1
        if pf.node_boundary_markers is not None
        and len(pf.node_boundary_markers)
        else 0
    )

    with dest.open("w") as fh:
        # --- Nodes ---
        fh.write(f"{n_nodes} 2 0 {n_bmark}\n")
        for i, nd in enumerate(nodes):
            row = f"{i + 1} {nd[0]:.16g} {nd[1]:.16g}"
            if n_bmark:
                bm = (
                    int(pf.node_boundary_markers[i])
                    if i < len(pf.node_boundary_markers)
                    else 0
                )
                row += f" {bm}"
            fh.write(row + "\n")

        # --- Segments ---
        segs = pf.segments
        n_segs = len(segs)
        has_smk = (
            pf.segment_markers is not None
            and len(pf.segment_markers) == n_segs
        )
        fh.write(f"{n_segs} {1 if has_smk else 0}\n")
        for i, seg in enumerate(segs):
            row = f"{i + 1} {int(seg[0])} {int(seg[1])}"
            if has_smk:
                row += f" {int(pf.segment_markers[i])}"
            fh.write(row + "\n")

        # --- Holes ---
        holes = pf.holes
        n_holes = len(holes)
        fh.write(f"{n_holes}\n")
        for i, h in enumerate(holes):
            fh.write(f"{i + 1} {h[0]:.16g} {h[1]:.16g}\n")

        # --- Regions ---
        regions = pf.regions
        n_regions = len(regions)
        fh.write(f"{n_regions}\n")
        for i, reg in enumerate(regions):
            fh.write(
                f"{i + 1} {reg[0]:.16g} {reg[1]:.16g} {reg[2]:g} {reg[3]:g}\n"
            )

    return dest
