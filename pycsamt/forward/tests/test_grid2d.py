# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Behaviour-contract tests for :mod:`pycsamt.forward.grid2d`.

Covers the ``make_padding`` helper and the :class:`Grid2D` dataclass:
validation in ``__post_init__``, coordinate/shape properties, padding
slices, station indexing, conductivity/harmonic-mean helpers, the
classmethod constructors (``halfspace``, ``from_1d_layers``,
``with_anomaly``, ``random``), ``plot`` and ``__repr__``.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.forward.grid2d import Grid2D, make_padding
from pycsamt.forward.synthetic import LayeredModel


# ─────────────────────────────────────────────────────────────────────────
# make_padding
# ─────────────────────────────────────────────────────────────────────────


def test_make_padding_cell_count_and_monotonic_growth():
    pad = make_padding(core_spacing=100.0, n_pad=5, factor=1.3)

    assert pad.shape == (5,)
    # strictly increasing away from the core edge
    assert np.all(np.diff(pad) > 0.0)
    # first cell is core_spacing * factor, subsequent ones grow by `factor`
    assert pad[0] == pytest.approx(100.0 * 1.3)
    assert np.allclose(pad[1:] / pad[:-1], 1.3)


def test_make_padding_default_factor():
    pad = make_padding(50.0, 3)
    assert np.allclose(pad, 50.0 * 1.3 ** np.arange(1, 4))


# ─────────────────────────────────────────────────────────────────────────
# __post_init__ validation
# ─────────────────────────────────────────────────────────────────────────


def test_post_init_rejects_bad_resistivity_shape():
    with pytest.raises(ValueError, match="does not match"):
        Grid2D(
            dx=np.full(4, 10.0),
            dz=np.full(3, 10.0),
            resistivity=np.full((3, 5), 100.0),  # wrong nx
            x_stations=np.array([0.0]),
        )


def test_post_init_rejects_nonpositive_resistivity():
    rho = np.full((3, 4), 100.0)
    rho[1, 1] = 0.0
    with pytest.raises(ValueError, match="strictly positive"):
        Grid2D(
            dx=np.full(4, 10.0),
            dz=np.full(3, 10.0),
            resistivity=rho,
            x_stations=np.array([0.0]),
        )


def test_post_init_rejects_nonpositive_spacing():
    dx = np.full(4, 10.0)
    dx[0] = -1.0
    with pytest.raises(ValueError, match="strictly positive"):
        Grid2D(
            dx=dx,
            dz=np.full(3, 10.0),
            resistivity=np.full((3, 4), 100.0),
            x_stations=np.array([0.0]),
        )


def test_post_init_rejects_stations_outside_grid():
    with pytest.raises(ValueError, match="lie outside the grid"):
        Grid2D(
            dx=np.full(4, 10.0),
            dz=np.full(3, 10.0),
            resistivity=np.full((3, 4), 100.0),
            x_stations=np.array([-5.0, 20.0]),
        )


def test_post_init_accepts_valid_grid():
    g = Grid2D(
        dx=np.full(4, 10.0),
        dz=np.full(3, 10.0),
        resistivity=np.full((3, 4), 100.0),
        x_stations=np.array([0.0, 20.0, 40.0]),
    )
    assert g.nx == 4
    assert g.nz == 3


# ─────────────────────────────────────────────────────────────────────────
# Shape / coordinate properties on a small hand-built grid
# ─────────────────────────────────────────────────────────────────────────


def _small_grid(n_pad=0):
    dx = np.array([10.0, 20.0, 30.0])
    dz = np.array([5.0, 15.0])
    rho = np.array(
        [
            [1.0, 2.0, 4.0],
            [8.0, 16.0, 32.0],
        ]
    )
    return Grid2D(
        dx=dx,
        dz=dz,
        resistivity=rho,
        x_stations=np.array([0.0, 30.0, 60.0]),
        n_pad=n_pad,
    )


def test_basic_shape_properties():
    g = _small_grid()
    assert g.nx == 3
    assert g.nz == 2
    assert g.n_stations == 3
    assert g.n_nodes == (3 + 1) * (2 + 1)


def test_x_and_z_nodes():
    g = _small_grid()
    assert np.allclose(g.x_nodes, [0.0, 10.0, 30.0, 60.0])
    assert np.allclose(g.z_nodes, [0.0, 5.0, 20.0])
    assert len(g.x_nodes) == g.nx + 1
    assert len(g.z_nodes) == g.nz + 1


def test_x_and_z_centers_are_midpoints():
    g = _small_grid()
    assert np.allclose(g.x_centers, [5.0, 20.0, 45.0])
    assert np.allclose(g.z_centers, [2.5, 12.5])


def test_x_and_z_extent():
    g = _small_grid()
    assert g.x_extent == pytest.approx(60.0)
    assert g.z_extent == pytest.approx(20.0)


def test_core_slices_no_padding():
    g = _small_grid(n_pad=0)
    assert g.core_x_slice == slice(None)
    assert g.core_z_slice == slice(None)
    assert np.array_equal(g.core_resistivity, g.resistivity)


def test_core_slices_with_padding():
    g = _small_grid(n_pad=1)
    assert g.core_x_slice == slice(1, 2)  # nx=3, p=1 -> [1:2]
    assert g.core_z_slice == slice(None, 1)  # nz=2, p=1 -> [:1]
    core = g.core_resistivity
    assert core.shape == (1, 1)
    assert core[0, 0] == g.resistivity[0, 1]


def test_column_profile_returns_rho_and_thickness():
    g = _small_grid()
    rho, thick = g.column_profile(1)
    assert np.allclose(rho, [2.0, 16.0])
    assert np.allclose(thick, g.dz[:-1])
    assert thick.shape == (g.nz - 1,)


def test_conductivity_is_inverse_of_resistivity():
    g = _small_grid()
    assert np.allclose(g.conductivity, 1.0 / g.resistivity)


def test_harmonic_mean_x_formula():
    g = _small_grid()
    s = 1.0 / g.resistivity
    expected = 2.0 * s[:, :-1] * s[:, 1:] / (s[:, :-1] + s[:, 1:])
    hmx = g.harmonic_mean_x()
    assert hmx.shape == (g.nz, g.nx - 1)
    assert np.allclose(hmx, expected)
    # hand-check one entry: row 0, cols 0/1 -> sigma = 1, 0.5
    assert hmx[0, 0] == pytest.approx(2.0 * 1.0 * 0.5 / (1.0 + 0.5))


def test_harmonic_mean_z_formula():
    g = _small_grid()
    s = 1.0 / g.resistivity
    expected = 2.0 * s[:-1, :] * s[1:, :] / (s[:-1, :] + s[1:, :])
    hmz = g.harmonic_mean_z()
    assert hmz.shape == (g.nz - 1, g.nx)
    assert np.allclose(hmz, expected)
    # hand-check column 0: sigma = 1.0 (row0) and 0.125 (row1)
    assert hmz[0, 0] == pytest.approx(2.0 * 1.0 * 0.125 / (1.0 + 0.125))


def test_station_cell_indices():
    g = _small_grid()
    # x_nodes = [0, 10, 30, 60]; stations at 0, 30, 60
    idx = g.station_cell_indices()
    assert idx.tolist() == [0, 2, 2]  # station at 60 clipped to last cell


def test_station_cell_indices_midpoint():
    g = Grid2D(
        dx=np.array([10.0, 10.0]),
        dz=np.array([5.0]),
        resistivity=np.array([[1.0, 2.0]]),
        x_stations=np.array([5.0, 15.0]),
    )
    idx = g.station_cell_indices()
    assert idx.tolist() == [0, 1]


def test_station_node_indices():
    g = _small_grid()
    # x_nodes = [0, 10, 30, 60], stations at [0, 30, 60]
    idx = g.station_node_indices()
    assert idx.tolist() == [0, 2, 3]


def test_repr_contains_key_info():
    g = _small_grid()
    r = repr(g)
    assert isinstance(r, str) and r
    assert "Grid2D" in r
    assert f"nx={g.nx}" in r
    assert f"nz={g.nz}" in r
    assert f"n_stations={g.n_stations}" in r


def test_repr_includes_name_when_set():
    g = _small_grid()
    g.name = "my-model"
    assert "name='my-model'" in repr(g)


def test_repr_omits_name_when_blank():
    g = _small_grid()
    assert g.name == ""
    assert "name=" not in repr(g)


# ─────────────────────────────────────────────────────────────────────────
# Grid2D.halfspace
# ─────────────────────────────────────────────────────────────────────────


def test_halfspace_uniform_resistivity_and_shape():
    g = Grid2D.halfspace(
        rho=100.0,
        nx=10,
        nz=8,
        x_max=1000.0,
        z_max=500.0,
        n_pad=4,
        n_stations=5,
    )
    assert np.all(g.resistivity == 100.0)
    # nx/nz include left+right / bottom padding
    assert g.nx == 10 + 2 * 4
    assert g.nz == 8 + 4
    assert g.n_pad == 4
    assert g.n_stations == 5


def test_halfspace_core_region_matches_requested_extent():
    g = Grid2D.halfspace(
        rho=50.0, nx=10, nz=8, x_max=1000.0, z_max=500.0, n_pad=4
    )
    core_dx = g.dx[g.core_x_slice]
    core_dz = g.dz[g.core_z_slice]
    assert core_dx.sum() == pytest.approx(1000.0)
    assert core_dz.sum() == pytest.approx(500.0)


def test_halfspace_stations_within_core_region():
    g = Grid2D.halfspace(
        rho=100.0,
        nx=10,
        nz=8,
        x_max=1000.0,
        z_max=500.0,
        n_pad=4,
        n_stations=5,
        station_x_max=800.0,
    )
    # station offset should sit inside the padded x-range
    assert g.x_stations.min() >= g.x_nodes[0]
    assert g.x_stations.max() <= g.x_nodes[-1]
    # stations span exactly station_x_max relative to core-left edge
    x_pad_offset = g.x_nodes[g.n_pad]
    assert g.x_stations.min() == pytest.approx(x_pad_offset)
    assert g.x_stations.max() == pytest.approx(x_pad_offset + 800.0)


def test_halfspace_no_padding():
    g = Grid2D.halfspace(rho=10.0, nx=5, nz=4, n_pad=0, n_stations=3)
    assert g.n_pad == 0
    assert g.nx == 5
    assert g.nz == 4
    assert g.core_x_slice == slice(None)


# ─────────────────────────────────────────────────────────────────────────
# Grid2D.with_anomaly
# ─────────────────────────────────────────────────────────────────────────


def test_with_anomaly_inside_and_outside_values():
    g = Grid2D.with_anomaly(
        bg_rho=100.0,
        anomaly_rho=5.0,
        anomaly_bounds=(1000.0, 3000.0, 200.0, 800.0),
        nx=40,
        nz=25,
        x_max=6000.0,
        z_max=4000.0,
        n_pad=8,
        n_stations=12,
    )
    x_pad_offset = g.x_nodes[g.n_pad]
    xc = g.x_centers - x_pad_offset
    zc = g.z_centers

    inside_cols = (xc >= 1000.0) & (xc <= 3000.0)
    inside_rows = (zc >= 200.0) & (zc <= 800.0)
    anomaly_block = g.resistivity[np.ix_(inside_rows, inside_cols)]
    assert np.all(anomaly_block == 5.0)

    # far outside the anomaly (near core edges), background is untouched
    assert g.resistivity[0, g.n_pad] == pytest.approx(100.0)
    assert g.resistivity[-1, g.n_pad] == pytest.approx(100.0)


def test_with_anomaly_touching_grid_boundary():
    # anomaly starting at x=0 (left edge of core) and z=0 (surface)
    g = Grid2D.with_anomaly(
        bg_rho=200.0,
        anomaly_rho=1.0,
        anomaly_bounds=(0.0, 1000.0, 0.0, 500.0),
        nx=20,
        nz=15,
        x_max=6000.0,
        z_max=4000.0,
        n_pad=5,
        n_stations=6,
    )
    x_pad_offset = g.x_nodes[g.n_pad]
    xc = g.x_centers - x_pad_offset
    zc = g.z_centers
    inside_cols = (xc >= 0.0) & (xc <= 1000.0)
    inside_rows = (zc >= 0.0) & (zc <= 500.0)
    assert inside_cols.any() and inside_rows.any()
    block = g.resistivity[np.ix_(inside_rows, inside_cols)]
    assert np.all(block == 1.0)
    # padding columns to the left of the core must remain background
    assert np.all(g.resistivity[:, : g.n_pad] == 200.0)


def test_with_anomaly_default_name_mentions_values():
    g = Grid2D.with_anomaly(bg_rho=50.0, anomaly_rho=2.0, nx=10, nz=8, n_pad=3)
    assert "50.0" in g.name
    assert "2.0" in g.name


# ─────────────────────────────────────────────────────────────────────────
# Grid2D.from_1d_layers
# ─────────────────────────────────────────────────────────────────────────


def test_from_1d_layers_matches_supplied_model_at_every_column():
    model = LayeredModel(resistivity=[100.0, 10.0, 500.0], thickness=[300.0, 800.0])
    g = Grid2D.from_1d_layers(model, nx=15, x_max=6000.0, n_pad=6, n_stations=5)

    # laterally invariant: every column has the same resistivity profile
    core_cols = g.resistivity[:, g.core_x_slice]
    first_col = core_cols[:, 0]
    for j in range(core_cols.shape[1]):
        assert np.allclose(core_cols[:, j], first_col)

    # the profile's core rows equal the model's layer resistivities plus
    # the synthetic halfspace discretisation cell
    assert first_col[0] == pytest.approx(100.0)
    assert first_col[1] == pytest.approx(10.0)
    assert first_col[2] == pytest.approx(500.0)


def test_from_1d_layers_default_name_and_layer_count():
    model = LayeredModel(resistivity=[100.0, 10.0], thickness=[400.0])
    g = Grid2D.from_1d_layers(model, nx=10, n_pad=4, n_stations=4)
    assert "2 layers" in g.name


def test_from_1d_layers_respects_dz_min():
    model = LayeredModel(resistivity=[100.0, 10.0], thickness=[1.0])
    g = Grid2D.from_1d_layers(
        model, nx=10, n_pad=4, n_stations=4, dz_min=50.0
    )
    # the (tiny) 1 m layer thickness must be floored to dz_min
    assert g.dz[0] == pytest.approx(50.0)


# ─────────────────────────────────────────────────────────────────────────
# Grid2D.random
# ─────────────────────────────────────────────────────────────────────────


def test_random_reproducible_with_same_seed():
    g1 = Grid2D.random(nx=12, nz=10, n_pad=3, n_stations=4, seed=42)
    g2 = Grid2D.random(nx=12, nz=10, n_pad=3, n_stations=4, seed=42)
    assert np.array_equal(g1.resistivity, g2.resistivity)


def test_random_differs_with_different_seed():
    g1 = Grid2D.random(nx=12, nz=10, n_pad=3, n_stations=4, seed=1)
    g2 = Grid2D.random(nx=12, nz=10, n_pad=3, n_stations=4, seed=2)
    assert not np.array_equal(g1.resistivity, g2.resistivity)


def test_random_resistivity_within_bounds():
    g = Grid2D.random(
        nx=12,
        nz=10,
        n_pad=3,
        n_stations=4,
        rho_min=1.0,
        rho_max=10_000.0,
        seed=7,
    )
    assert np.all(g.resistivity >= 1.0)
    assert np.all(g.resistivity <= 10_000.0)


def test_random_without_lateral_variation_is_column_invariant():
    g = Grid2D.random(
        nx=12,
        nz=10,
        n_pad=3,
        n_stations=4,
        lateral_variation=False,
        seed=3,
    )
    for i in range(g.nz):
        row = g.resistivity[i, :]
        assert np.allclose(row, row[0])


def test_random_accepts_existing_generator():
    rng = np.random.default_rng(99)
    g = Grid2D.random(nx=8, nz=6, n_pad=2, n_stations=3, seed=rng)
    assert g.resistivity.shape == (g.nz, g.nx)


# ─────────────────────────────────────────────────────────────────────────
# plot()
# ─────────────────────────────────────────────────────────────────────────


def test_plot_returns_axes_default_kwargs():
    g = Grid2D.halfspace(rho=100.0, nx=10, nz=8, n_pad=4, n_stations=5)
    ax = g.plot()
    assert ax is not None
    assert ax.get_xlabel() == "Distance  (m)"
    assert ax.get_ylabel() == "Depth  (m)"
    plt.close(ax.get_figure())


def test_plot_with_existing_axes():
    fig, ax = plt.subplots()
    out = g_for_plot().plot(ax=ax)
    assert out is ax
    plt.close(fig)


def g_for_plot():
    return Grid2D.halfspace(rho=100.0, nx=8, nz=6, n_pad=3, n_stations=4)


def test_plot_linear_scale_no_stations_no_clip():
    g = g_for_plot()
    ax = g.plot(log_scale=False, show_stations=False, clip_core=False)
    assert ax.get_title() == "halfspace"  # g_for_plot() uses the default name
    plt.close(ax.get_figure())


def test_plot_with_vmin_vmax_and_custom_cmap():
    g = g_for_plot()
    ax = g.plot(log_scale=True, vmin=1.0, vmax=3.0, cmap="viridis")
    plt.close(ax.get_figure())


def test_plot_uses_name_as_title():
    g = Grid2D.halfspace(rho=100.0, nx=8, nz=6, n_pad=3, n_stations=4, name="test-model")
    ax = g.plot()
    assert ax.get_title() == "test-model"
    plt.close(ax.get_figure())
