"""Global mesh-display controls for pyCSAMT figures.

This module centralizes how a computational or inversion mesh is drawn,
for both mesh families pyCSAMT supports:

* a **structured (rectilinear)** tensor grid --
  :class:`~pycsamt.forward.maxwell.contracts.MaxwellMesh` or
  :class:`~pycsamt.ai.geology.fields.GeologyGrid` -- drawn by
  :func:`draw_mesh` from 1-D edge arrays;
* an **unstructured triangular** mesh --
  :class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh` -- drawn by
  :func:`draw_tri_mesh` from its node/triangle connectivity.

Both share the same three-preset display model: a plain color fill (the
historical default, identical to a bare ``pcolormesh``/``tripcolor``
call), a color fill with cell boundaries drawn on top for QC review, or a
bare wireframe diagram with no fill at all -- e.g. to inspect a solver
mesh's discretization before any inversion has run -- via one shared
:class:`MeshStyle`/:class:`PyCSAMTMesh` configuration object.

A :class:`TriMesh` built by
:func:`~pycsamt.models.mare2dem.tri_mesh.build_survey_mesh` already has
its node z-coordinates following real topography (the mesh's own top
boundary *is* the topo line), so unlike the rectilinear
:func:`~pycsamt.topo.section.build_topo_section` path -- which warps a
flat regular grid onto a terrain surface via interpolation --
:func:`draw_tri_mesh` needs no separate drape step: plotting the mesh's
own node coordinates on a depth-positive-down (inverted) y-axis already
shows a topography-draped section.

Examples
--------
Draw a plain filled section (current default behaviour, no edges)::

    from pycsamt.api.mesh import draw_mesh

    fill, edges = draw_mesh(ax, x_edges, z_edges, values)
    assert edges is None

Overlay mesh cell boundaries on top of the color fill for review::

    fill, edges = draw_mesh(ax, x_edges, z_edges, values, preset="review")

Show only the mesh wireframe, no color fill — useful to QC a solver mesh
before running an inversion::

    draw_mesh(ax, x_edges, z_edges, preset="diagram")

Draw a triangular-mesh section the same way, sharing the same presets::

    from pycsamt.api.mesh import draw_tri_mesh

    fill, edges = draw_tri_mesh(ax, tri_mesh, per_triangle_values)
    fill, edges = draw_tri_mesh(ax, tri_mesh, values, preset="review")

Customize one preset globally (applies to both mesh families)::

    from pycsamt.api import configure_mesh

    configure_mesh(review__edge__alpha=0.4, review__edge__linewidth=0.25)
"""

from __future__ import annotations

import copy
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.tri import Triangulation

    from pycsamt.ai.geology.fields import GeologyGrid
    from pycsamt.forward.maxwell.contracts import MaxwellMesh
    from pycsamt.forward.maxwell.contracts_tri import TriMesh

_PRESETS = {"filled", "review", "diagram"}


def cell_edges_from_centres(centres: np.ndarray) -> np.ndarray:
    """Return an ``(n+1,)`` cell-boundary array from an ``(n,)`` centre array.

    Assumes ``centres`` is monotonically increasing; boundaries fall at the
    midpoint between neighbouring centres, with the first/last boundary
    extrapolated by the first/last spacing. A single-centre input returns
    two boundaries half a unit either side.
    """
    centres = np.asarray(centres, dtype=float)
    if len(centres) < 2:
        half = 0.5
        return np.array([centres[0] - half, centres[0] + half])
    d = np.diff(centres)
    mids = centres[:-1] + d / 2.0
    lo = centres[0] - d[0] / 2.0
    hi = centres[-1] + d[-1] / 2.0
    return np.concatenate([[lo], mids, [hi]])


@dataclass
class MeshFillStyle:
    """Color-fill layer of a mesh display."""

    show: bool = True
    alpha: float = 1.0
    shading: str = "flat"


@dataclass
class MeshEdgeStyle:
    """Cell-boundary line layer of a mesh display."""

    show: bool = False
    color: str = "#2a2a2a"
    linewidth: float = 0.4
    linestyle: str = "-"
    alpha: float = 0.6
    zorder: int = 3


@dataclass
class MeshStyle:
    """Complete visual-control bundle for one mesh-display preset."""

    fill: MeshFillStyle = field(default_factory=MeshFillStyle)
    edge: MeshEdgeStyle = field(default_factory=MeshEdgeStyle)

    def copy(self, **kw: Any) -> MeshStyle:
        """Return a modified deep copy of this mesh style."""
        new = copy.deepcopy(self)
        for key, value in kw.items():
            setattr(new, key, value)
        return new


def draw_mesh(
    ax: Axes,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    values: np.ndarray | None = None,
    *,
    style: MeshStyle | None = None,
    preset: str = "filled",
    cmap: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Any | None = None,
    **pcolormesh_kw: Any,
) -> tuple[Any | None, Any | None]:
    """Draw a structured mesh: an optional color fill, optional cell edges.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    x_edges, y_edges : ndarray
        1-D, strictly monotonic cell-boundary arrays (``n+1`` entries for
        ``n`` cells along that axis). Use :func:`cell_edges_from_centres`,
        :func:`edges_from_maxwell_mesh`, or :func:`edges_from_geology_grid`
        to build these from cell centres or a mesh contract object.
    values : ndarray, shape ``(len(y_edges) - 1, len(x_edges) - 1)``, optional
        Cell values for the color fill. Required when the resolved style's
        fill is enabled; omit entirely for a pure wireframe diagram
        (``preset="diagram"``).
    style : MeshStyle, optional
        Explicit style. Overrides ``preset`` when supplied.
    preset : str, default "filled"
        One of ``"filled"`` (color only, no edges — matches a bare
        ``pcolormesh`` call), ``"review"`` (color fill with edges on top),
        or ``"diagram"`` (edges only, no fill). Looked up on
        :data:`PYCSAMT_MESH` when ``style`` is not given.
    cmap, vmin, vmax, norm : optional
        Forwarded to the fill layer's ``pcolormesh`` call.
    **pcolormesh_kw
        Additional keywords forwarded to the fill layer's ``pcolormesh``
        call only (never the edge layer, which fixes its own face/edge
        keywords).

    Returns
    -------
    (fill_mappable, edge_collection) : tuple
        Either element is ``None`` when that layer was not drawn. The fill
        mappable can still be handed to a colorbar helper, e.g.
        :func:`pycsamt.api.plot.add_colorbar`.
    """
    resolved = style or PYCSAMT_MESH.style_for(preset)

    fill_im = None
    if resolved.fill.show:
        if values is None:
            msg = (
                "draw_mesh: 'values' is required when the fill layer is "
                "enabled (style.fill.show=True / preset != 'diagram')."
            )
            raise ValueError(msg)
        fill_im = ax.pcolormesh(
            x_edges,
            y_edges,
            values,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            norm=norm,
            shading=resolved.fill.shading,
            alpha=resolved.fill.alpha,
            edgecolor="none",
            **pcolormesh_kw,
        )

    edge_coll = None
    if resolved.edge.show:
        # matplotlib.pcolormesh ties a single `alpha` to both face and edge
        # on one QuadMesh, so independent fill/edge alpha needs a second,
        # face-transparent layer rather than edge kwargs on the fill call.
        # x_edges/y_edges may be 1-D (rectilinear grid) or 2-D (a draped/
        # topographic mesh, e.g. from pycsamt.topo.drape_section) -- the
        # dummy array's cell count comes from the trailing two axes either
        # way, since a 2-D edge array's own shape already is (ny+1, nx+1).
        x_arr = np.asarray(x_edges)
        y_arr = np.asarray(y_edges)
        if x_arr.ndim == 1 and y_arr.ndim == 1:
            dummy_shape = (y_arr.shape[0] - 1, x_arr.shape[0] - 1)
        else:
            dummy_shape = (y_arr.shape[-2] - 1, y_arr.shape[-1] - 1)
        dummy = np.zeros(dummy_shape)
        edge_coll = ax.pcolormesh(
            x_edges,
            y_edges,
            dummy,
            facecolor="none",
            edgecolor=resolved.edge.color,
            linewidth=resolved.edge.linewidth,
            linestyle=resolved.edge.linestyle,
            alpha=resolved.edge.alpha,
            zorder=resolved.edge.zorder,
            shading=resolved.fill.shading,
        )
        edge_coll.set_array(None)

    return fill_im, edge_coll


def edges_from_maxwell_mesh(
    mesh: MaxwellMesh,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(x_edges, z_edges)`` for a 2-D :class:`MaxwellMesh`.

    ``MaxwellMesh`` already stores cell boundaries (``x_edges_m``,
    ``z_edges_m``), so this is a thin, allocation-free unpack. Raises
    :class:`NotImplementedError` for a 3-D mesh (``y_edges_m is not
    None``) since :func:`draw_mesh` renders one 2-D slice at a time —
    select a slice's ``x``/``z`` (or ``x``/``y``) edges before calling.
    """
    if getattr(mesh, "y_edges_m", None) is not None:
        msg = (
            "edges_from_maxwell_mesh: draw_mesh renders one 2-D (x, z) or "
            "(x, y) slice at a time; a 3-D MaxwellMesh has no single slice "
            "to return. Pass mesh.x_edges_m paired with mesh.z_edges_m or "
            "mesh.y_edges_m directly for the slice you want."
        )
        raise NotImplementedError(msg)
    return np.asarray(mesh.x_edges_m), np.asarray(mesh.z_edges_m)


def edges_from_geology_grid(
    grid: GeologyGrid,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(x_edges, z_edges)`` for a 2-D :class:`GeologyGrid`.

    ``GeologyGrid`` stores regularly-spaced cell *centres* (``x_m``,
    ``z_m``), not edges, so this reconstructs boundaries with
    :func:`cell_edges_from_centres`. Raises :class:`NotImplementedError`
    for a 3-D grid (``y_m is not None``), matching
    :func:`edges_from_maxwell_mesh`.
    """
    if getattr(grid, "y_m", None) is not None:
        msg = (
            "edges_from_geology_grid: draw_mesh renders one 2-D (x, z) or "
            "(x, y) slice at a time; a 3-D GeologyGrid has no single slice "
            "to return. Build edges from grid.x_m paired with grid.z_m or "
            "grid.y_m directly for the slice you want."
        )
        raise NotImplementedError(msg)
    return (
        cell_edges_from_centres(grid.x_m),
        cell_edges_from_centres(grid.z_m),
    )


def triangulation_from_tri_mesh(mesh: TriMesh) -> Triangulation:
    """Return a :class:`matplotlib.tri.Triangulation` for a 2-D :class:`TriMesh`.

    Parameters
    ----------
    mesh : TriMesh
        Solver-neutral unstructured triangular mesh.

    Returns
    -------
    matplotlib.tri.Triangulation
        Built directly from ``mesh.nodes_m`` and ``mesh.triangles``, ready
        for :func:`draw_tri_mesh` or a direct ``ax.tripcolor``/``ax.triplot``
        call.

    Examples
    --------
    >>> from pycsamt.forward.maxwell import TriMesh
    >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
    >>> tri = triangulation_from_tri_mesh(mesh)
    >>> tri.triangles.tolist()
    [[0, 1, 2]]
    """
    from matplotlib.tri import Triangulation

    return Triangulation(
        mesh.nodes_m[:, 0], mesh.nodes_m[:, 1], mesh.triangles
    )


def draw_tri_mesh(
    ax: Axes,
    tri_mesh: TriMesh,
    values: np.ndarray | None = None,
    *,
    style: MeshStyle | None = None,
    preset: str = "filled",
    cmap: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Any | None = None,
    **tripcolor_kw: Any,
) -> tuple[Any | None, Any | None]:
    """Draw an unstructured triangular mesh: color fill and/or cell edges.

    The triangular-mesh counterpart of :func:`draw_mesh`, sharing the same
    :class:`MeshStyle`/:data:`PYCSAMT_MESH` preset system so
    ``"filled"``/``"review"``/``"diagram"`` behave identically across both
    mesh families.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    tri_mesh : TriMesh
        Solver-neutral unstructured triangular mesh (see
        :func:`triangulation_from_tri_mesh`).
    values : ndarray, shape ``(tri_mesh.n_triangles,)``, optional
        Per-triangle cell values for the color fill (e.g.
        ``log10(resistivity_ohm_m)``). Required when the resolved
        style's fill is enabled; omit for a pure wireframe diagram
        (``preset="diagram"``).
    style : MeshStyle, optional
        Explicit style. Overrides ``preset`` when supplied.
    preset : str, default "filled"
        One of ``"filled"``, ``"review"``, or ``"diagram"`` -- see
        :func:`draw_mesh`.
    cmap, vmin, vmax, norm : optional
        Forwarded to the fill layer's ``ax.tripcolor`` call.
    **tripcolor_kw
        Additional keywords forwarded to the fill layer's
        ``ax.tripcolor`` call only.

    Returns
    -------
    (fill_mappable, edge_lines) : tuple
        ``fill_mappable`` is the ``ax.tripcolor`` return value (or
        ``None``); ``edge_lines`` is the 2-tuple of ``Line2D`` objects
        ``ax.triplot`` returns (or ``None``). The fill mappable can
        still be handed to a colorbar helper, e.g.
        :func:`pycsamt.api.plot.add_colorbar`.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> import matplotlib.pyplot as plt
    >>> from pycsamt.forward.maxwell import TriMesh
    >>> mesh = TriMesh(
    ...     [[0, 0], [1000, 0], [1000, 500], [0, 500]],
    ...     [[0, 1, 2], [0, 2, 3]],
    ... )
    >>> fig, ax = plt.subplots()
    >>> fill, edges = draw_tri_mesh(ax, mesh, [1.0, 2.0])
    >>> fill is not None, edges is None
    (True, True)
    >>> ax.invert_yaxis()  # depth positive down
    """
    resolved = style or PYCSAMT_MESH.style_for(preset)
    triangulation = triangulation_from_tri_mesh(tri_mesh)

    fill_im = None
    if resolved.fill.show:
        if values is None:
            msg = (
                "draw_tri_mesh: 'values' is required when the fill layer "
                "is enabled (style.fill.show=True / preset != 'diagram')."
            )
            raise ValueError(msg)
        fill_im = ax.tripcolor(
            triangulation,
            facecolors=np.asarray(values),
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            norm=norm,
            alpha=resolved.fill.alpha,
            edgecolor="none",
            **tripcolor_kw,
        )

    edge_lines = None
    if resolved.edge.show:
        edge_lines = ax.triplot(
            triangulation,
            color=resolved.edge.color,
            linewidth=resolved.edge.linewidth,
            linestyle=resolved.edge.linestyle,
            alpha=resolved.edge.alpha,
            zorder=resolved.edge.zorder,
        )

    return fill_im, edge_lines


class PyCSAMTMesh:
    """Package-wide mesh-display control container."""

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.filled = MeshStyle(
            fill=MeshFillStyle(show=True),
            edge=MeshEdgeStyle(show=False),
        )
        self.review = MeshStyle(
            fill=MeshFillStyle(show=True, alpha=0.92),
            edge=MeshEdgeStyle(
                show=True,
                alpha=0.55,
                linewidth=0.35,
                color="#1a1a1a",
            ),
        )
        self.diagram = MeshStyle(
            fill=MeshFillStyle(show=False),
            edge=MeshEdgeStyle(
                show=True,
                alpha=0.9,
                linewidth=0.6,
                color="#2a2a2a",
            ),
        )

    def style_for(self, preset: str = "filled") -> MeshStyle:
        """Return the mesh style associated with ``preset``."""
        preset = str(preset).lower()
        if preset not in _PRESETS:
            msg = f"mesh preset must be one of {sorted(_PRESETS)}."
            raise ValueError(msg)
        return getattr(self, preset)

    def configure(self, **kw: Any) -> None:
        """Configure mesh styles using ``preset__layer__attribute`` paths."""
        for path, value in kw.items():
            parts = path.split("__")
            obj = self
            for part in parts[:-1]:
                obj = getattr(obj, part)
            setattr(obj, parts[-1], value)

    @contextmanager
    def context(
        self,
        preset: str | None = None,
        **kw: Any,
    ) -> Generator[PyCSAMTMesh, None, None]:
        """Temporarily override mesh styles, then restore them."""
        snapshot = self._snapshot()
        try:
            if preset:
                self.use_preset(preset)
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self._restore(snapshot)

    def use_preset(self, preset: str) -> None:
        """Copy one preset into the filled slot."""
        self.filled = copy.deepcopy(self.style_for(preset))

    def reset(self) -> None:
        """Reset mesh controls to package defaults."""
        self._init_defaults()

    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = ["PyCSAMTMesh"]
        for name in sorted(_PRESETS):
            style = self.style_for(name)
            lines.append(
                "  "
                f"{name}: fill={style.fill.show!r}, "
                f"edge={style.edge.show!r}, "
                f"edge_alpha={style.edge.alpha}",
            )
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

    def _snapshot(self) -> dict[str, MeshStyle]:
        return {name: copy.deepcopy(self.style_for(name)) for name in _PRESETS}

    def _restore(self, snapshot: dict[str, MeshStyle]) -> None:
        for name, value in snapshot.items():
            setattr(self, name, value)


PYCSAMT_MESH = PyCSAMTMesh()


def configure_mesh(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_MESH`."""
    PYCSAMT_MESH.configure(**kw)


def reset_mesh() -> None:
    """Reset :data:`PYCSAMT_MESH` to package defaults."""
    PYCSAMT_MESH.reset()


__all__ = [
    "PYCSAMT_MESH",
    "MeshEdgeStyle",
    "MeshFillStyle",
    "MeshStyle",
    "PyCSAMTMesh",
    "cell_edges_from_centres",
    "configure_mesh",
    "draw_mesh",
    "draw_tri_mesh",
    "edges_from_geology_grid",
    "edges_from_maxwell_mesh",
    "reset_mesh",
    "triangulation_from_tri_mesh",
]
