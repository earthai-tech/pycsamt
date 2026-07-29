"""Unit tests for :mod:`pycsamt.forward.em2d`'s finite-difference assembly.

Focused on ``_assemble_tm``'s interface-resistivity averaging: it was
changed from a plain arithmetic mean to a thickness-weighted harmonic
mean, matching the parallel-conductance combination
``_surface_impedance_tm`` already used for the analogous surface
extraction. See the module docstring of
:mod:`pycsamt.ai.training.dataset2d` for the investigation this fix
came out of, and why the earlier "TM mode does not converge" finding
was a measurement artifact (independently varying x/z resolution and
regenerating a new random realization at every resolution) rather
than a defect in the solver itself.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.em2d import _assemble_tm, _node_index
from pycsamt.forward.grid2d import Grid2D


def _two_layer_grid(rho_top: float, rho_bottom: float) -> Grid2D:
    """A 2x2-cell grid, uniform in x, with a resistivity contrast in z
    only -- isolates the x-direction coefficient (which combines the
    two z-neighbour cells) from the z-direction one.
    """
    dx = np.array([100.0, 150.0])
    dz = np.array([50.0, 80.0])
    resistivity = np.array([[rho_top, rho_top], [rho_bottom, rho_bottom]])
    return Grid2D(dx, dz, resistivity, np.array([50.0]), n_pad=0)


def _two_column_grid(rho_left: float, rho_right: float) -> Grid2D:
    """A 2x2-cell grid, uniform in z, with a resistivity contrast in x
    only -- isolates the z-direction coefficient (which combines the
    two x-neighbour cells).
    """
    dx = np.array([100.0, 150.0])
    dz = np.array([50.0, 80.0])
    resistivity = np.array([[rho_left, rho_right], [rho_left, rho_right]])
    return Grid2D(dx, dz, resistivity, np.array([50.0]), n_pad=0)


def test_x_direction_coefficient_uses_thickness_weighted_harmonic_mean():
    rho_top, rho_bottom = 10.0, 1000.0
    grid = _two_layer_grid(rho_top, rho_bottom)
    omega = 2.0 * np.pi * 1.0
    A = _assemble_tm(grid, omega)

    # The single interior node (i=1, j=1) sits exactly on the z-interface
    # between the top row (rho_top) and bottom row (rho_bottom); its
    # x-direction neighbours are (1, 2) [right] and (1, 0) [left].
    nx1 = grid.nx + 1
    k = _node_index(1, 1, nx1)
    k_right = _node_index(1, 2, nx1)
    dz_u, dz_d = grid.dz[0], grid.dz[1]
    dx_r = grid.dx[1]
    dx_avg = 0.5 * (grid.dx[0] + grid.dx[1])

    expected_rho_xr = (dz_u + dz_d) / (dz_u / rho_top + dz_d / rho_bottom)
    expected_c_xr = expected_rho_xr / (dx_r * dx_avg)
    naive_arithmetic_rho_xr = 0.5 * (rho_top + rho_bottom)

    assert A[k, k_right] == pytest.approx(expected_c_xr)
    # The two averaging conventions must actually differ for this
    # contrast, or the test would not discriminate between them.
    assert expected_rho_xr != pytest.approx(naive_arithmetic_rho_xr)
    assert A[k, k_right] != pytest.approx(
        naive_arithmetic_rho_xr / (dx_r * dx_avg)
    )


def test_z_direction_coefficient_uses_thickness_weighted_harmonic_mean():
    rho_left, rho_right = 5.0, 2000.0
    grid = _two_column_grid(rho_left, rho_right)
    omega = 2.0 * np.pi * 1.0
    A = _assemble_tm(grid, omega)

    nx1 = grid.nx + 1
    k = _node_index(1, 1, nx1)
    k_down = _node_index(2, 1, nx1)
    dx_l, dx_r = grid.dx[0], grid.dx[1]
    dz_d = grid.dz[1]
    dz_avg = 0.5 * (grid.dz[0] + grid.dz[1])

    expected_rho_zd = (dx_l + dx_r) / (dx_l / rho_left + dx_r / rho_right)
    expected_c_zd = expected_rho_zd / (dz_d * dz_avg)
    naive_arithmetic_rho_zd = 0.5 * (rho_left + rho_right)

    assert A[k, k_down] == pytest.approx(expected_c_zd)
    assert expected_rho_zd != pytest.approx(naive_arithmetic_rho_zd)


def test_assembly_unaffected_by_averaging_choice_for_uniform_medium():
    """A change in averaging scheme must be a no-op when the two
    values being averaged are equal -- the reason the already-passing
    half-space/layered-earth benchmarks (laterally and, respectively,
    vertically uniform) are untouched by this fix.
    """
    grid = _two_layer_grid(75.0, 75.0)
    omega = 2.0 * np.pi * 2.0
    A = _assemble_tm(grid, omega)
    nx1 = grid.nx + 1
    k = _node_index(1, 1, nx1)
    k_right = _node_index(1, 2, nx1)
    dx_avg = 0.5 * (grid.dx[0] + grid.dx[1])
    expected = 75.0 / (grid.dx[1] * dx_avg)
    assert A[k, k_right] == pytest.approx(expected)


def _resample(rho_fine, fine_x_edges, fine_z_edges, nx, nz, x_max, z_max):
    x_edges = np.linspace(0, x_max, nx + 1)
    z_edges = np.linspace(0, z_max, nz + 1)
    x_c = 0.5 * (x_edges[:-1] + x_edges[1:])
    z_c = 0.5 * (z_edges[:-1] + z_edges[1:])
    xi = np.clip(
        np.searchsorted(fine_x_edges, x_c, side="right") - 1, 0, len(fine_x_edges) - 2
    )
    zi = np.clip(
        np.searchsorted(fine_z_edges, z_c, side="right") - 1, 0, len(fine_z_edges) - 2
    )
    return rho_fine[np.ix_(zi, xi)]


def test_tm_mode_converges_under_joint_mesh_refinement():
    """Resample one fixed, laterally-varying resistivity field onto
    increasingly fine, *jointly* refined meshes (x and z together) and
    solve TM mode at each resolution. A true refinement study of the
    same physical model must show shrinking increments, not the
    oscillation a prior investigation saw when it varied resolution in
    only one direction at a time, or regenerated an independent random
    field at every resolution instead of resampling one fixed field
    (both of which look like non-convergence without being it).
    """
    from pycsamt.ai.geology import GeologyGrid, GaussianCorrelation, generate_gaussian_field
    from pycsamt.forward.em2d import MT2DForward

    x_max, z_max = 12000.0, 4000.0
    nx_fine, nz_fine = 160, 96
    fine_grid = GeologyGrid.regular_2d(
        nx=nx_fine, nz=nz_fine, dx_m=x_max / nx_fine, dz_m=z_max / nz_fine
    )
    field = generate_gaussian_field(
        fine_grid,
        GaussianCorrelation(length_x_m=1200.0, length_z_m=500.0),
        seed=13,
    )
    rho_fine = 10.0 ** (2.0 + 0.4 * field.values)
    fine_x_edges = np.linspace(0, x_max, nx_fine + 1)
    fine_z_edges = np.linspace(0, z_max, nz_fine + 1)

    def rho_a_tm(nx, nz):
        resistivity = _resample(
            rho_fine, fine_x_edges, fine_z_edges, nx, nz, x_max, z_max
        )
        dx = np.full(nx, x_max / nx)
        dz = np.full(nz, z_max / nz)
        grid = Grid2D(dx, dz, resistivity, np.array([x_max / 2.0]), n_pad=0)
        response = MT2DForward(np.array([0.5]), grid, verbose=False).run()
        return float(response.rho_a_tm[0, 0])

    coarse = rho_a_tm(20, 12)
    medium = rho_a_tm(40, 24)
    fine = rho_a_tm(80, 48)

    coarse_to_medium = abs(medium - coarse)
    medium_to_fine = abs(fine - medium)
    assert medium_to_fine < coarse_to_medium
