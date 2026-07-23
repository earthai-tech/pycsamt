# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Behaviour-contract tests for :mod:`pycsamt.forward.grid3d`.

Covers ``Grid3D`` construction/validation, coordinate/property math,
the XZ/YZ 2-D slice extraction used by the quasi-3D solver, the
station-to-cell lookup helpers, and the ``halfspace``/``block_anomaly``/
``random_layered`` constructors.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from pycsamt.forward.grid2d import Grid2D  # noqa: E402
from pycsamt.forward.grid3d import Grid3D, _ensure_rng  # noqa: E402


# ─────────────────────────────────────────────────────────────────────────────
# _ensure_rng
# ─────────────────────────────────────────────────────────────────────────────


def test_ensure_rng_passes_through_existing_generator():
    rng = np.random.default_rng(123)
    assert _ensure_rng(rng) is rng


def test_ensure_rng_builds_generator_from_seed():
    rng = _ensure_rng(7)
    assert isinstance(rng, np.random.Generator)


# ─────────────────────────────────────────────────────────────────────────────
# A small hand-built grid with distinguishable cell values, used for the
# slicing / property / station-lookup tests below.
# ─────────────────────────────────────────────────────────────────────────────


def _make_grid(n_pad=0):
    dx = np.array([10.0, 20.0, 30.0])
    dy = np.array([5.0, 5.0, 5.0, 5.0])
    dz = np.array([100.0, 200.0])
    nz, ny, nx = 2, 4, 3
    rho = np.empty((nz, ny, nx))
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                rho[iz, iy, ix] = 100 * iz + 10 * iy + ix + 1
    stations = np.array(
        [[15.0, 7.0], [25.0, 12.0], [5.0, 2.0], [50.0, 18.0]]
    )
    return Grid3D(
        dx=dx, dy=dy, dz=dz, resistivity=rho, stations_xy=stations,
        n_pad=n_pad,
    )


# ─────────────────────────────────────────────────────────────────────────────
# __post_init__ / properties
# ─────────────────────────────────────────────────────────────────────────────


def test_post_init_shapes_dtypes_and_coordinates():
    g = _make_grid()
    assert g.dx.dtype == np.float64
    assert g.nx == 3 and g.ny == 4 and g.nz == 2
    assert g.n_stations == 4
    assert np.allclose(g.x_nodes, [0, 10, 30, 60])
    assert np.allclose(g.y_nodes, [0, 5, 10, 15, 20])
    assert np.allclose(g.z_nodes, [0, 100, 300])
    assert np.allclose(g.x_centers, [5, 20, 45])
    assert np.allclose(g.y_centers, [2.5, 7.5, 12.5, 17.5])
    assert np.allclose(g.z_centers, [50, 200])
    assert g.x_extent == pytest.approx(60.0)
    assert g.y_extent == pytest.approx(20.0)
    assert g.z_extent == pytest.approx(300.0)


def test_post_init_reshapes_1d_stations_xy():
    dx = np.array([10.0, 10.0])
    dy = np.array([10.0, 10.0])
    dz = np.array([10.0])
    rho = np.full((1, 2, 2), 100.0)
    g = Grid3D(dx=dx, dy=dy, dz=dz, resistivity=rho, stations_xy=[5.0, 5.0])
    assert g.stations_xy.shape == (1, 2)
    assert np.allclose(g.stations_xy, [[5.0, 5.0]])


def test_post_init_rejects_bad_resistivity_shape():
    dx = np.array([10.0, 10.0])
    dy = np.array([10.0, 10.0])
    dz = np.array([10.0])
    bad_rho = np.full((1, 2, 3), 100.0)  # nx should be 2, not 3
    with pytest.raises(ValueError, match="does not match"):
        Grid3D(
            dx=dx, dy=dy, dz=dz, resistivity=bad_rho,
            stations_xy=[[5.0, 5.0]],
        )


def test_post_init_rejects_nonpositive_resistivity():
    dx = np.array([10.0, 10.0])
    dy = np.array([10.0, 10.0])
    dz = np.array([10.0])
    rho = np.full((1, 2, 2), 100.0)
    rho[0, 0, 0] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        Grid3D(
            dx=dx, dy=dy, dz=dz, resistivity=rho,
            stations_xy=[[5.0, 5.0]],
        )


def test_post_init_rejects_station_outside_x_extent():
    dx = np.array([10.0, 10.0])
    dy = np.array([10.0, 10.0])
    dz = np.array([10.0])
    rho = np.full((1, 2, 2), 100.0)
    with pytest.raises(ValueError, match="outside the grid extent"):
        Grid3D(
            dx=dx, dy=dy, dz=dz, resistivity=rho,
            stations_xy=[[999.0, 5.0]],
        )


def test_post_init_rejects_station_outside_y_extent():
    dx = np.array([10.0, 10.0])
    dy = np.array([10.0, 10.0])
    dz = np.array([10.0])
    rho = np.full((1, 2, 2), 100.0)
    with pytest.raises(ValueError, match="outside the grid extent"):
        Grid3D(
            dx=dx, dy=dy, dz=dz, resistivity=rho,
            stations_xy=[[5.0, 999.0]],
        )


# ─────────────────────────────────────────────────────────────────────────────
# station cell lookup
# ─────────────────────────────────────────────────────────────────────────────


def test_station_cell_lookup_matches_expected_cells():
    g = _make_grid()
    # x_nodes = [0, 10, 30, 60]; stations x = [15, 25, 5, 50]
    assert np.array_equal(g._station_x_cells(), [1, 1, 0, 2])
    # y_nodes = [0, 5, 10, 15, 20]; stations y = [7, 12, 2, 18]
    assert np.array_equal(g._station_y_cells(), [1, 2, 0, 3])


# ─────────────────────────────────────────────────────────────────────────────
# xz_slice / yz_slice
# ─────────────────────────────────────────────────────────────────────────────


def test_xz_slice_matches_3d_array_and_station_indices():
    g = _make_grid()
    y_cells = g._station_y_cells()
    for yi in range(g.ny):
        g2d, idx = g.xz_slice(yi)
        assert isinstance(g2d, Grid2D)
        assert g2d.resistivity.shape == (g.nz, g.nx)
        assert np.array_equal(g2d.resistivity, g.resistivity[:, yi, :])
        assert np.array_equal(g2d.dx, g.dx)
        assert np.array_equal(g2d.dz, g.dz)
        assert g2d.n_pad == g.n_pad

        expected_idx = np.where(y_cells == yi)[0]
        assert np.array_equal(idx, expected_idx)
        if expected_idx.size:
            assert np.array_equal(
                g2d.x_stations, g.stations_xy[expected_idx, 0]
            )


def test_yz_slice_matches_3d_array_and_station_indices():
    g = _make_grid()
    x_cells = g._station_x_cells()
    for xi in range(g.nx):
        g2d, idx = g.yz_slice(xi)
        assert isinstance(g2d, Grid2D)
        assert g2d.resistivity.shape == (g.nz, g.ny)
        assert np.array_equal(g2d.resistivity, g.resistivity[:, :, xi])
        assert np.array_equal(g2d.dx, g.dy)
        assert np.array_equal(g2d.dz, g.dz)
        assert g2d.n_pad == g.n_pad

        expected_idx = np.where(x_cells == xi)[0]
        assert np.array_equal(idx, expected_idx)
        if expected_idx.size:
            assert np.array_equal(
                g2d.x_stations, g.stations_xy[expected_idx, 1]
            )


def test_xz_slice_falls_back_to_grid_midpoint_when_row_empty():
    dx = np.array([10.0, 20.0, 30.0])
    dy = np.array([5.0, 5.0, 5.0, 5.0])
    dz = np.array([100.0, 200.0])
    rho = np.full((2, 4, 3), 50.0)
    stations = np.array([[15.0, 7.0]])  # only y-cell 1 is populated
    g = Grid3D(dx=dx, dy=dy, dz=dz, resistivity=rho, stations_xy=stations)

    g2d, idx = g.xz_slice(3)  # empty row
    assert idx.size == 0
    assert g2d.x_stations.shape == (1,)
    assert g2d.x_stations[0] == pytest.approx(
        g.x_nodes[g.n_pad + g.nx // 2]
    )


def test_yz_slice_falls_back_to_grid_midpoint_when_column_empty():
    dx = np.array([10.0, 20.0, 30.0])
    dy = np.array([5.0, 5.0, 5.0, 5.0])
    dz = np.array([100.0, 200.0])
    rho = np.full((2, 4, 3), 50.0)
    stations = np.array([[15.0, 7.0]])  # only x-cell 1 is populated
    g = Grid3D(dx=dx, dy=dy, dz=dz, resistivity=rho, stations_xy=stations)

    g2d, idx = g.yz_slice(2)  # empty column
    assert idx.size == 0
    assert g2d.x_stations.shape == (1,)
    assert g2d.x_stations[0] == pytest.approx(
        g.y_nodes[g.n_pad + g.ny // 2]
    )


# ─────────────────────────────────────────────────────────────────────────────
# column_profile_3d / conductivity
# ─────────────────────────────────────────────────────────────────────────────


def test_column_profile_3d_matches_3d_array():
    g = _make_grid()
    rho, thick = g.column_profile_3d(xi=2, yi=1)
    assert np.array_equal(rho, g.resistivity[:, 1, 2])
    assert np.array_equal(thick, g.dz[:-1])


def test_conductivity_is_inverse_of_resistivity():
    g = _make_grid()
    assert np.allclose(g.conductivity, 1.0 / g.resistivity)


# ─────────────────────────────────────────────────────────────────────────────
# core (_cx/_cy/_cz) slices
# ─────────────────────────────────────────────────────────────────────────────


def test_core_slices_exclude_padding():
    g = Grid3D.halfspace(
        nx=6, ny=8, nz=5, n_pad=3,
        x_max=600.0, y_max=800.0, z_max=500.0,
        nx_stations=2, ny_stations=2,
    )
    assert g._cx == slice(3, g.nx - 3)
    assert g._cy == slice(3, g.ny - 3)
    assert g._cz == slice(None, g.nz - 3)

    assert len(np.arange(g.nx)[g._cx]) == 6
    assert len(np.arange(g.ny)[g._cy]) == 8
    assert len(np.arange(g.nz)[g._cz]) == 5


def test_core_slices_are_full_slices_when_no_padding():
    g = _make_grid(n_pad=0)
    assert g._cx == slice(None)
    assert g._cy == slice(None)
    assert g._cz == slice(None)


# ─────────────────────────────────────────────────────────────────────────────
# plot()
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_returns_figure_and_three_axes():
    g = Grid3D.halfspace(
        nx=6, ny=6, nz=5, n_pad=2,
        x_max=600.0, y_max=600.0, z_max=500.0,
        nx_stations=3, ny_stations=3,
    )
    fig, axs = g.plot()
    try:
        assert axs.shape == (3,)
        assert fig is not None
    finally:
        plt.close(fig)


def test_plot_kwargs_branches_linear_scale_no_clip_no_stations():
    g = Grid3D.block_anomaly(
        nx=6, ny=6, nz=5, n_pad=2,
        x_max=600.0, y_max=600.0, z_max=500.0,
        bounds=(100.0, 400.0, 100.0, 400.0, 50.0, 300.0),
        nx_stations=2, ny_stations=2,
    )
    fig, axs = g.plot(
        log_scale=False, clip_core=False, show_stations=False,
        vmin=1.0, vmax=1000.0,
    )
    try:
        assert axs.shape == (3,)
    finally:
        plt.close(fig)


def test_plot_with_no_padding_and_no_stations():
    g = Grid3D.halfspace(
        nx=4, ny=4, nz=4, n_pad=0,
        x_max=400.0, y_max=400.0, z_max=400.0,
        nx_stations=0, ny_stations=0,
    )
    assert g.n_stations == 0
    fig, axs = g.plot()
    try:
        assert axs.shape == (3,)
    finally:
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# __repr__
# ─────────────────────────────────────────────────────────────────────────────


def test_repr_contains_key_info():
    g = _make_grid()
    r = repr(g)
    assert "Grid3D(nx=3" in r
    assert "n_stations=4" in r
    assert "n_pad=0" in r


def test_repr_includes_name_when_set():
    g = Grid3D.halfspace(
        name="myhalfspace",
        nx=4, ny=4, nz=4, n_pad=0,
        x_max=400.0, y_max=400.0, z_max=400.0,
        nx_stations=2, ny_stations=2,
    )
    assert "name='myhalfspace'" in repr(g)


# ─────────────────────────────────────────────────────────────────────────────
# halfspace()
# ─────────────────────────────────────────────────────────────────────────────


def test_halfspace_uniform_resistivity_and_sizing():
    g = Grid3D.halfspace(
        rho=250.0, nx=10, ny=12, nz=8,
        x_max=1000.0, y_max=1200.0, z_max=800.0,
        n_pad=4, nx_stations=3, ny_stations=3,
    )
    assert np.all(g.resistivity == 250.0)
    # padding is added on both sides in x/y, only at the bottom in z
    assert g.nx == 10 + 2 * 4
    assert g.ny == 12 + 2 * 4
    assert g.nz == 8 + 4
    assert g.n_stations == 9


def test_halfspace_station_grid_positions():
    g = Grid3D.halfspace(
        nx=10, ny=10, nz=5,
        x_max=1000.0, y_max=1000.0, z_max=500.0,
        n_pad=4, pad_factor=1.3, nx_stations=3, ny_stations=3,
    )
    x_pad_off = g.x_nodes[g.n_pad]  # cumulative left-padding width
    y_pad_off = g.y_nodes[g.n_pad]
    xs = np.linspace(0.0, 1000.0, 3) + x_pad_off
    ys = np.linspace(0.0, 1000.0, 3) + y_pad_off
    gx, gy = np.meshgrid(xs, ys)
    expected = np.column_stack([gx.ravel(), gy.ravel()])
    assert np.allclose(g.stations_xy, expected)


# ─────────────────────────────────────────────────────────────────────────────
# block_anomaly()
# ─────────────────────────────────────────────────────────────────────────────


def test_block_anomaly_inside_and_outside_values():
    # nx=10 core cells over x_max=1000 -> dx=100/cell; likewise y and z.
    g = Grid3D.block_anomaly(
        bg_rho=300.0, anomaly_rho=3.0,
        bounds=(200.0, 600.0, 200.0, 600.0, 100.0, 400.0),
        nx=10, ny=10, nz=8,
        x_max=1000.0, y_max=1000.0, z_max=800.0,
        n_pad=3, nx_stations=2, ny_stations=2,
    )
    # A cell clearly inside the block: core column/row 4 (x,y in
    # [400, 500]), core depth row 2 (z in [200, 300]); global indices
    # shift x/y by n_pad (no shift needed in z, padding is at the bottom).
    assert g.resistivity[2, g.n_pad + 4, g.n_pad + 4] == pytest.approx(3.0)
    # A cell clearly outside in x (core column 0, x in [0, 100]):
    assert g.resistivity[2, g.n_pad + 4, g.n_pad + 0] == pytest.approx(
        300.0
    )
    # A cell clearly outside in z (core depth row 6, z in [600, 700]):
    assert g.resistivity[6, g.n_pad + 4, g.n_pad + 4] == pytest.approx(
        300.0
    )
    # Padding cells are always background.
    assert np.all(g.resistivity[:, : g.n_pad, :] == 300.0)
    assert np.all(g.resistivity[:, :, : g.n_pad] == 300.0)


def test_block_anomaly_custom_name():
    g = Grid3D.block_anomaly(
        name="custom",
        nx=4, ny=4, nz=4, n_pad=0,
        x_max=400.0, y_max=400.0, z_max=400.0,
        nx_stations=2, ny_stations=2,
    )
    assert g.name == "custom"


def test_block_anomaly_default_name_mentions_both_resistivities():
    g = Grid3D.block_anomaly(
        bg_rho=500.0, anomaly_rho=5.0,
        nx=4, ny=4, nz=4, n_pad=0,
        x_max=400.0, y_max=400.0, z_max=400.0,
        nx_stations=2, ny_stations=2,
    )
    assert "500" in g.name
    assert "5" in g.name


# ─────────────────────────────────────────────────────────────────────────────
# random_layered()
# ─────────────────────────────────────────────────────────────────────────────


def test_random_layered_reproducible_with_same_seed():
    kwargs = dict(
        nx=6, ny=6, nz=6, n_pad=2,
        x_max=600.0, y_max=600.0, z_max=600.0,
        nx_stations=2, ny_stations=2,
    )
    g1 = Grid3D.random_layered(seed=42, **kwargs)
    g2 = Grid3D.random_layered(seed=42, **kwargs)
    assert np.array_equal(g1.resistivity, g2.resistivity)


def test_random_layered_differs_across_seeds():
    kwargs = dict(
        nx=6, ny=6, nz=6, n_pad=2,
        x_max=600.0, y_max=600.0, z_max=600.0,
        nx_stations=2, ny_stations=2,
    )
    g1 = Grid3D.random_layered(seed=1, **kwargs)
    g2 = Grid3D.random_layered(seed=2, **kwargs)
    assert not np.array_equal(g1.resistivity, g2.resistivity)


def test_random_layered_resistivity_within_bounds():
    g = Grid3D.random_layered(
        seed=7, nx=8, ny=8, nz=8, n_pad=2,
        rho_min=1.0, rho_max=10_000.0,
        x_max=800.0, y_max=800.0, z_max=800.0,
        nx_stations=2, ny_stations=2,
    )
    assert np.all(g.resistivity >= 1.0)
    assert np.all(g.resistivity <= 10_000.0)


def test_random_layered_without_lateral_variation_is_uniform_per_layer():
    g = Grid3D.random_layered(
        seed=3, n_layers=3, nx=6, ny=6, nz=9, n_pad=0,
        lateral_variation=False,
        x_max=600.0, y_max=600.0, z_max=900.0,
        nx_stations=2, ny_stations=2,
    )
    for iz in range(g.nz):
        layer = g.resistivity[iz]
        assert np.all(layer == layer[0, 0])
