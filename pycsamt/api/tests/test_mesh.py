"""Tests for global mesh-display controls."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.ai.geology.fields import GeologyGrid
from pycsamt.api import (
    PYCSAMT_MESH,
    MeshEdgeStyle,
    MeshFillStyle,
    MeshStyle,
    cell_edges_from_centres,
    configure_mesh,
    draw_mesh,
    edges_from_geology_grid,
    edges_from_maxwell_mesh,
    reset_mesh,
)
from pycsamt.forward.maxwell.contracts import MaxwellMesh


def test_mesh_presets_are_available():
    """Mesh presets should expose consistent fill/edge defaults."""
    reset_mesh()

    filled = PYCSAMT_MESH.style_for("filled")
    review = PYCSAMT_MESH.style_for("review")
    diagram = PYCSAMT_MESH.style_for("diagram")

    assert filled.fill.show is True
    assert filled.edge.show is False  # back-compat: current pcolormesh look

    assert review.fill.show is True
    assert review.edge.show is True

    assert diagram.fill.show is False
    assert diagram.edge.show is True


def test_mesh_configure_and_context_restore():
    """Nested mesh configuration should restore after context."""
    reset_mesh()
    original = PYCSAMT_MESH.review.edge.alpha

    configure_mesh(review__edge__alpha=0.9)
    assert PYCSAMT_MESH.review.edge.alpha == 0.9

    with PYCSAMT_MESH.context(review__edge__alpha=0.2):
        assert PYCSAMT_MESH.review.edge.alpha == 0.2

    assert PYCSAMT_MESH.review.edge.alpha == 0.9
    reset_mesh()
    assert PYCSAMT_MESH.review.edge.alpha == original


def test_draw_mesh_filled_preset_has_no_edge_collection():
    """The 'filled' preset should draw only the color fill."""
    fig, ax = plt.subplots()
    x_edges = np.array([0.0, 1.0, 2.0])
    y_edges = np.array([0.0, 1.0])
    values = np.array([[1.0, 2.0]])

    fill, edge = draw_mesh(ax, x_edges, y_edges, values, preset="filled")

    assert fill is not None
    assert edge is None
    assert len(ax.collections) == 1
    plt.close(fig)


def test_draw_mesh_review_preset_returns_both_layers():
    """The 'review' preset should draw fill and edges with independent alpha."""
    fig, ax = plt.subplots()
    x_edges = np.array([0.0, 1.0, 2.0])
    y_edges = np.array([0.0, 1.0])
    values = np.array([[1.0, 2.0]])
    style = MeshStyle(
        fill=MeshFillStyle(show=True, alpha=1.0),
        edge=MeshEdgeStyle(show=True, alpha=0.3),
    )

    fill, edge = draw_mesh(ax, x_edges, y_edges, values, style=style)

    assert fill is not None
    assert edge is not None
    assert len(ax.collections) == 2
    assert fill.get_alpha() != edge.get_alpha()
    plt.close(fig)


def test_draw_mesh_accepts_2d_draped_edges():
    """A 2-D (e.g. topography-draped) edge array should render both layers."""
    fig, ax = plt.subplots()
    x_edges = np.array([0.0, 1.0, 2.0, 3.0])
    x_2d = np.broadcast_to(x_edges[None, :], (3, 4)).astype(float)
    y_2d = np.random.default_rng(0).random((3, 4)) + np.arange(3)[:, None]
    values = np.random.default_rng(1).random((2, 3))

    fill, edge = draw_mesh(ax, x_2d, y_2d, values, preset="review")

    assert fill is not None
    assert edge is not None
    plt.close(fig)


def test_draw_mesh_diagram_preset_has_no_fill():
    """The 'diagram' preset should draw edges only, values optional."""
    fig, ax = plt.subplots()
    x_edges = np.array([0.0, 1.0, 2.0])
    y_edges = np.array([0.0, 1.0])

    fill, edge = draw_mesh(ax, x_edges, y_edges, preset="diagram")

    assert fill is None
    assert edge is not None
    assert len(ax.collections) == 1
    plt.close(fig)


def test_draw_mesh_toggle_off_via_style():
    """A style with both layers off should draw nothing."""
    fig, ax = plt.subplots()
    x_edges = np.array([0.0, 1.0, 2.0])
    y_edges = np.array([0.0, 1.0])
    style = MeshStyle(
        fill=MeshFillStyle(show=False),
        edge=MeshEdgeStyle(show=False),
    )

    fill, edge = draw_mesh(ax, x_edges, y_edges, style=style)

    assert fill is None
    assert edge is None
    assert len(ax.collections) == 0
    plt.close(fig)


def test_draw_mesh_requires_values_for_fill():
    """Requesting a fill layer without values should raise, not silently skip."""
    fig, ax = plt.subplots()
    x_edges = np.array([0.0, 1.0, 2.0])
    y_edges = np.array([0.0, 1.0])

    with pytest.raises(ValueError, match="values"):
        draw_mesh(ax, x_edges, y_edges, values=None, preset="filled")
    plt.close(fig)


def test_edges_from_maxwell_mesh_2d():
    """The MaxwellMesh adapter should return its own edges unchanged."""
    mesh = MaxwellMesh([0, 100, 250], [0, 50, 150])

    x_edges, z_edges = edges_from_maxwell_mesh(mesh)

    assert np.array_equal(x_edges, mesh.x_edges_m)
    assert np.array_equal(z_edges, mesh.z_edges_m)


def test_edges_from_maxwell_mesh_3d_raises():
    """A 3-D MaxwellMesh has no single 2-D slice to return."""
    mesh = MaxwellMesh([0, 100, 250], [0, 50, 150], y_edges_m=[0, 10, 20])

    with pytest.raises(NotImplementedError):
        edges_from_maxwell_mesh(mesh)


def test_edges_from_geology_grid_matches_cell_edges_from_centres():
    """The GeologyGrid adapter should match direct centre-to-edge conversion."""
    grid = GeologyGrid(x_m=[50, 150, 250], z_m=[25, 75])

    x_edges, z_edges = edges_from_geology_grid(grid)

    assert np.array_equal(x_edges, cell_edges_from_centres(grid.x_m))
    assert np.array_equal(z_edges, cell_edges_from_centres(grid.z_m))
    assert np.allclose(x_edges, [0.0, 100.0, 200.0, 300.0])
    assert np.allclose(z_edges, [0.0, 50.0, 100.0])


def test_edges_from_geology_grid_3d_raises():
    """A 3-D GeologyGrid has no single 2-D slice to return."""
    grid = GeologyGrid(x_m=[50, 150], z_m=[25, 75], y_m=[10, 20])

    with pytest.raises(NotImplementedError):
        edges_from_geology_grid(grid)


def test_mesh_reset_restores_defaults():
    """reset_mesh() should undo any configure()/context() drift."""
    configure_mesh(diagram__edge__linewidth=5.0)
    assert PYCSAMT_MESH.diagram.edge.linewidth == 5.0

    reset_mesh()

    assert PYCSAMT_MESH.diagram.edge.linewidth != 5.0
