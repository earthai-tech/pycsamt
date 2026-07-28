# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for :mod:`pycsamt.forward.em3d`.

Covers the pure Yee-grid indexing helpers, the ``ForwardResponse3D``
output container, the quasi-3D orthogonal-profile solver, and the full
finite-difference 3-D (FD3D) edge-based solver.  Grids are kept tiny
(a handful of core cells + minimal padding) so the FD3D sparse solves
stay well under a second; physics checks use generous tolerances since
these coarse grids are for code-path coverage, not quantitative
accuracy.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from pycsamt.forward.em3d import (
    ForwardResponse3D,
    MT3DForward,
    _apply_dirichlet_3d,
    _bc_vector_3d,
    _boundary_edge_mask,
    _build_curl_matrix,
    _edge_conductivity,
    _ex_idx,
    _ey_idx,
    _ez_idx,
    _fd3d_edge_counts,
    _fd3d_face_counts,
)
from pycsamt.forward.grid3d import Grid3D

# ─────────────────────────────────────────────────────────────────────────────
# Pure helpers: edge / face counts
# ─────────────────────────────────────────────────────────────────────────────


def test_fd3d_edge_counts_matches_hand_computation():
    nx, ny, nz = 2, 3, 4
    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    assert (n_ex, n_ey, n_ez) == (
        nx * (ny + 1) * (nz + 1),
        (nx + 1) * ny * (nz + 1),
        (nx + 1) * (ny + 1) * nz,
    )
    assert (n_ex, n_ey, n_ez) == (40, 45, 48)


def test_fd3d_face_counts_matches_hand_computation():
    nx, ny, nz = 2, 3, 4
    n_hx, n_hy, n_hz = _fd3d_face_counts(nx, ny, nz)
    assert (n_hx, n_hy, n_hz) == (
        (nx + 1) * ny * nz,
        nx * (ny + 1) * nz,
        nx * ny * (nz + 1),
    )
    assert (n_hx, n_hy, n_hz) == (36, 32, 30)


# ─────────────────────────────────────────────────────────────────────────────
# Edge index helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_edge_index_helpers_are_unique_and_cover_full_range():
    nx, ny, nz = 2, 3, 4
    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    n_edges = n_ex + n_ey + n_ez

    ex_ids = {
        int(_ex_idx(i, j, k, ny, nz))
        for i in range(nx)
        for j in range(ny + 1)
        for k in range(nz + 1)
    }
    ey_ids = {
        int(_ey_idx(i, j, k, nx, ny, nz))
        for i in range(nx + 1)
        for j in range(ny)
        for k in range(nz + 1)
    }
    ez_ids = {
        int(_ez_idx(i, j, k, nx, ny, nz))
        for i in range(nx + 1)
        for j in range(ny + 1)
        for k in range(nz)
    }

    assert len(ex_ids) == n_ex
    assert len(ey_ids) == n_ey
    assert len(ez_ids) == n_ez

    # Each family occupies a disjoint contiguous block.
    assert ex_ids.isdisjoint(ey_ids)
    assert ex_ids.isdisjoint(ez_ids)
    assert ey_ids.isdisjoint(ez_ids)
    assert ex_ids | ey_ids | ez_ids == set(range(n_edges))
    assert min(ex_ids) == 0 and max(ex_ids) == n_ex - 1
    assert min(ey_ids) == n_ex and max(ey_ids) == n_ex + n_ey - 1
    assert min(ez_ids) == n_ex + n_ey and max(ez_ids) == n_edges - 1


# ─────────────────────────────────────────────────────────────────────────────
# Boundary edge mask
# ─────────────────────────────────────────────────────────────────────────────


def test_boundary_edge_mask_interior_count_matches_hand_computation():
    nx, ny, nz = 2, 3, 4
    mask = _boundary_edge_mask(nx, ny, nz)
    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    assert mask.shape == (n_ex + n_ey + n_ez,)
    assert mask.dtype == bool

    # Interior edges are those strictly inside all 6 bounding faces.
    interior_ex = nx * (ny - 1) * (nz - 1)
    interior_ey = (nx - 1) * ny * (nz - 1)
    interior_ez = (nx - 1) * (ny - 1) * nz
    n_interior = interior_ex + interior_ey + interior_ez

    assert (~mask).sum() == n_interior
    assert mask.sum() == mask.size - n_interior

    # Spot-check: an Ex edge at j=0 (front face) must be boundary.
    assert mask[_ex_idx(0, 0, 1, ny, nz)]
    # An Ex edge strictly interior (j=1, k=1..2) must not be boundary
    # for this grid (ny=3 -> interior j in {1}, nz=4 -> interior k in {1,2}).
    assert not mask[_ex_idx(0, 1, 1, ny, nz)]


# ─────────────────────────────────────────────────────────────────────────────
# ForwardResponse3D — direct construction
# ─────────────────────────────────────────────────────────────────────────────


def _make_response3d(nf=3, ns=2, method="quasi3d") -> ForwardResponse3D:
    freqs = np.logspace(-1, 1, nf)
    stations_xy = np.column_stack([np.arange(ns) * 100.0, np.zeros(ns)])
    rng = np.random.default_rng(0)

    def _rand_complex():
        return rng.normal(size=(nf, ns)) + 1j * rng.normal(size=(nf, ns))

    zxy, zyx, zxx, zyy = (_rand_complex() for _ in range(4))
    omega = 2.0 * np.pi * freqs[:, None]
    rho_a_xy = np.abs(zxy) ** 2 / (omega * 4e-7 * np.pi)
    phase_xy = np.angle(zxy, deg=True)
    rho_a_yx = np.abs(zyx) ** 2 / (omega * 4e-7 * np.pi)
    phase_yx = np.angle(zyx, deg=True)
    rho_a_xx = np.abs(zxx) ** 2 / (omega * 4e-7 * np.pi)
    phase_xx = np.angle(zxx, deg=True)
    rho_a_yy = np.abs(zyy) ** 2 / (omega * 4e-7 * np.pi)
    phase_yy = np.angle(zyy, deg=True)

    return ForwardResponse3D(
        freqs=freqs,
        stations_xy=stations_xy,
        zxy=zxy,
        zyx=zyx,
        zxx=zxx,
        zyy=zyy,
        rho_a_xy=rho_a_xy,
        phase_xy=phase_xy,
        rho_a_yx=rho_a_yx,
        phase_yx=phase_yx,
        rho_a_xx=rho_a_xx,
        phase_xx=phase_xx,
        rho_a_yy=rho_a_yy,
        phase_yy=phase_yy,
        method=method,
    )


def test_forward_response3d_basic_properties():
    resp = _make_response3d(nf=5, ns=3)
    assert resp.n_freqs == 5
    assert resp.n_stations == 3
    assert np.allclose(resp.periods, 1.0 / resp.freqs)
    assert "ForwardResponse3D" in repr(resp)
    assert "n_freqs=5" in repr(resp)
    assert "n_stations=3" in repr(resp)
    assert resp.method == "quasi3d"


def test_forward_response3d_station_response_matches_full_arrays():
    resp = _make_response3d(nf=4, ns=3)
    sub = resp.station_response(1)
    assert np.allclose(sub["freqs"], resp.freqs)
    assert np.allclose(sub["periods"], resp.periods)
    assert np.allclose(sub["zxy"], resp.zxy[:, 1])
    assert np.allclose(sub["zyx"], resp.zyx[:, 1])
    assert np.allclose(sub["zxx"], resp.zxx[:, 1])
    assert np.allclose(sub["zyy"], resp.zyy[:, 1])
    assert np.allclose(sub["rho_a_xy"], resp.rho_a_xy[:, 1])
    assert np.allclose(sub["phase_yx"], resp.phase_yx[:, 1])


def test_forward_response3d_to_feature_array_shapes_and_values():
    resp = _make_response3d(nf=6, ns=4)

    default_feat = resp.to_feature_array()  # "xy_yx", log_rho, phase
    assert default_feat.shape == (4, 2 * 2 * 6)

    xy_only_no_phase = resp.to_feature_array(components="xy", include_phase=False)
    assert xy_only_no_phase.shape == (4, 6)
    assert np.allclose(xy_only_no_phase, np.log10(np.maximum(resp.rho_a_xy, 1e-12)).T)

    all_feat = resp.to_feature_array(components="all")
    assert all_feat.shape == (4, 4 * 2 * 6)

    linear_feat = resp.to_feature_array(
        components="xy", log_rho=False, include_phase=False
    )
    assert np.allclose(linear_feat, resp.rho_a_xy.T)


def test_forward_response3d_to_survey_dataset():
    resp = _make_response3d(nf=4, ns=3)
    ds = resp.to_survey_dataset()
    assert ds.n_surveys == 1
    assert ds.n_stations == 3
    assert ds.X.shape == (1, 3, 2 * 2 * 4)
    assert ds.y.shape == (1, 3, 1)
    assert np.allclose(ds.y, 0.0)
    assert np.allclose(ds.coords, resp.stations_xy)
    assert np.allclose(ds.freqs, resp.freqs)

    y_models = np.ones((3, 2), dtype=np.float32)
    ds2 = resp.to_survey_dataset(y_models=y_models)
    assert ds2.y.shape == (1, 3, 2)
    assert np.allclose(ds2.y, 1.0)


def test_forward_response3d_to_survey_dataset_labels_actual_solver_method():
    """``solver`` on the dataset must reflect ``self.method``, not a hardcoded
    'quasi3d' -- regression test for a bug where an fd3d-solved response was
    mislabeled as quasi3d in the exported SurveyDataset3D."""
    resp = _make_response3d(nf=3, ns=2, method="fd3d")
    ds = resp.to_survey_dataset()
    assert ds.solver == "fd3d"


# ─────────────────────────────────────────────────────────────────────────────
# MT3DForward — construction / validation
# ─────────────────────────────────────────────────────────────────────────────


def test_mt3dforward_invalid_method_raises():
    g = Grid3D.halfspace(
        rho=100.0,
        nx=3,
        ny=3,
        nz=3,
        x_max=300.0,
        y_max=300.0,
        z_max=300.0,
        n_pad=1,
        nx_stations=1,
        ny_stations=1,
    )
    with pytest.raises(ValueError):
        MT3DForward([1.0], g, method="bogus")


# ─────────────────────────────────────────────────────────────────────────────
# MT3DForward — quasi3d path
# ─────────────────────────────────────────────────────────────────────────────


def test_quasi3d_halfspace_recovers_background_resistivity():
    rho_bg = 80.0
    g = Grid3D.halfspace(
        rho=rho_bg,
        nx=8,
        ny=8,
        nz=7,
        x_max=3_000.0,
        y_max=3_000.0,
        z_max=2_500.0,
        n_pad=6,
        pad_factor=1.3,
        nx_stations=2,
        ny_stations=2,
    )
    freqs = np.logspace(-1, 1, 4)
    resp = MT3DForward(freqs, g, method="quasi3d", verbose=False).run()

    assert isinstance(resp, ForwardResponse3D)
    assert resp.method == "quasi3d"
    assert resp.rho_a_xy.shape == (4, g.n_stations)
    assert resp.zxy.shape == (4, g.n_stations)

    # Off-diagonal components are identically zero in the 2-D approximation.
    assert np.all(resp.zxx == 0.0)
    assert np.all(resp.zyy == 0.0)

    # Halfspace: both principal apparent resistivities should be within a
    # generous factor of the true background value at every frequency.
    assert np.all(resp.rho_a_xy > 0.3 * rho_bg)
    assert np.all(resp.rho_a_xy < 3.0 * rho_bg)
    assert np.all(resp.rho_a_yx > 0.3 * rho_bg)
    assert np.all(resp.rho_a_yx < 3.0 * rho_bg)

    # Phase should stay in the physical quadrant, near 45 deg for a halfspace.
    assert np.all(resp.phase_xy > 20.0)
    assert np.all(resp.phase_xy < 70.0)


def test_quasi3d_block_anomaly_produces_lateral_variation():
    g = Grid3D.block_anomaly(
        bg_rho=100.0,
        anomaly_rho=5.0,
        bounds=(1_500.0, 2_500.0, 1_500.0, 2_500.0, 300.0, 900.0),
        nx=6,
        ny=6,
        nz=5,
        x_max=4_000.0,
        y_max=4_000.0,
        z_max=2_000.0,
        n_pad=4,
        pad_factor=1.3,
        nx_stations=3,
        ny_stations=3,
    )
    freqs = np.logspace(-1, 1, 4)
    resp = MT3DForward(freqs, g, method="quasi3d", verbose=False).run()

    # A buried conductor beneath a subset of stations must create lateral
    # variation in apparent resistivity across the station array.
    assert np.all(resp.rho_a_xy.std(axis=1) > 1.0)


def test_quasi3d_verbose_prints_progress(capsys):
    # Smoke test for the verbose progress-printing branches (one XZ + one
    # YZ profile is enough to exercise every print statement).
    g = Grid3D.halfspace(
        rho=100.0,
        nx=2,
        ny=2,
        nz=2,
        x_max=400.0,
        y_max=400.0,
        z_max=300.0,
        n_pad=1,
        nx_stations=1,
        ny_stations=1,
    )
    resp = MT3DForward([1.0], g, method="quasi3d", verbose=True).run()
    out = capsys.readouterr().out
    assert "XZ profile" in out
    assert "YZ profile" in out
    assert "done" in out
    assert isinstance(resp, ForwardResponse3D)


# ─────────────────────────────────────────────────────────────────────────────
# MT3DForward — full FD3D path
# ─────────────────────────────────────────────────────────────────────────────


def test_fd3d_verbose_prints_progress(capsys):
    # Smoke test for the verbose branches in the FD3D loop (unknown count
    # + per-frequency progress + done banner).
    g = Grid3D.halfspace(
        rho=100.0,
        nx=2,
        ny=2,
        nz=2,
        x_max=400.0,
        y_max=400.0,
        z_max=300.0,
        n_pad=1,
        nx_stations=1,
        ny_stations=1,
    )
    resp = MT3DForward([1.0], g, method="fd3d", verbose=True).run()
    out = capsys.readouterr().out
    assert "unknowns" in out
    assert "freq" in out
    assert "done" in out
    assert isinstance(resp, ForwardResponse3D)


def test_fd3d_run_tiny_halfspace_grid_is_finite_and_sane():
    rho_bg = 100.0
    g = Grid3D.halfspace(
        rho=rho_bg,
        nx=3,
        ny=3,
        nz=3,
        x_max=1_000.0,
        y_max=1_000.0,
        z_max=800.0,
        n_pad=1,
        pad_factor=1.5,
        nx_stations=2,
        ny_stations=2,
    )
    freqs = np.array([1.0, 0.1])
    resp = MT3DForward(freqs, g, method="fd3d", verbose=False).run()

    assert isinstance(resp, ForwardResponse3D)
    assert resp.method == "fd3d"
    ns = g.n_stations
    assert resp.zxy.shape == (2, ns)
    assert resp.zxx.shape == (2, ns)

    for arr in (resp.zxy, resp.zyx, resp.zxx, resp.zyy):
        assert np.all(np.isfinite(arr))
    for arr in (resp.rho_a_xy, resp.rho_a_yx, resp.rho_a_xx, resp.rho_a_yy):
        assert np.all(np.isfinite(arr))
        assert np.all(arr >= 0.0)

    # FD3D on a very coarse grid is not quantitatively accurate, but it
    # should recover the right order of magnitude for a halfspace.
    assert np.all(resp.rho_a_xy > 0.1 * rho_bg)
    assert np.all(resp.rho_a_xy < 10.0 * rho_bg)
    assert np.all(resp.rho_a_yx > 0.1 * rho_bg)
    assert np.all(resp.rho_a_yx < 10.0 * rho_bg)

    # Off-diagonal terms exist (true 3-D coupling) but should be much
    # smaller than the principal terms for a laterally-uniform halfspace.
    assert np.abs(resp.zxx).max() < np.abs(resp.zxy).max()
    assert np.abs(resp.zyy).max() < np.abs(resp.zyx).max()


# ─────────────────────────────────────────────────────────────────────────────
# Curl matrix
# ─────────────────────────────────────────────────────────────────────────────


def _tiny_uniform_grid(n=2, cell=100.0, n_pad=0):
    return Grid3D.halfspace(
        rho=100.0,
        nx=n,
        ny=n,
        nz=n,
        x_max=n * cell,
        y_max=n * cell,
        z_max=n * cell,
        n_pad=n_pad,
        nx_stations=1,
        ny_stations=1,
    )


def test_build_curl_matrix_shape_and_hx_row0_entries():
    g = _tiny_uniform_grid(n=2, cell=100.0)
    C = _build_curl_matrix(g)

    n_ex, n_ey, n_ez = _fd3d_edge_counts(g.nx, g.ny, g.nz)
    n_hx, n_hy, n_hz = _fd3d_face_counts(g.nx, g.ny, g.nz)
    assert sparse.issparse(C)
    assert C.shape == (n_hx + n_hy + n_hz, n_ex + n_ey + n_ez)

    # Row 0 is the Hx face at (i=0, j=0, k=0): curl_x = dEz/dy - dEy/dz.
    row0 = np.asarray(C.getrow(0).todense()).ravel()
    nz_cols = np.nonzero(row0)[0]
    expected_cols = {
        int(_ez_idx(0, 1, 0, g.nx, g.ny, g.nz)): +1.0 / g.dy[0],
        int(_ez_idx(0, 0, 0, g.nx, g.ny, g.nz)): -1.0 / g.dy[0],
        int(_ey_idx(0, 0, 1, g.nx, g.ny, g.nz)): -1.0 / g.dz[0],
        int(_ey_idx(0, 0, 0, g.nx, g.ny, g.nz)): +1.0 / g.dz[0],
    }
    assert set(nz_cols.tolist()) == set(expected_cols.keys())
    for col in nz_cols:
        assert row0[col] == pytest.approx(expected_cols[int(col)])


# ─────────────────────────────────────────────────────────────────────────────
# Edge conductivity
# ─────────────────────────────────────────────────────────────────────────────


def test_edge_conductivity_uniform_on_halfspace():
    rho = 50.0
    g = Grid3D.halfspace(
        rho=rho,
        nx=3,
        ny=3,
        nz=3,
        x_max=300.0,
        y_max=300.0,
        z_max=300.0,
        n_pad=0,
        nx_stations=1,
        ny_stations=1,
    )
    sig = _edge_conductivity(g)
    n_ex, n_ey, n_ez = _fd3d_edge_counts(g.nx, g.ny, g.nz)
    assert sig.shape == (n_ex + n_ey + n_ez,)
    assert np.allclose(sig, 1.0 / rho)


# ─────────────────────────────────────────────────────────────────────────────
# BC vector
# ─────────────────────────────────────────────────────────────────────────────


def test_bc_vector_3d_x_polarisation_shape_and_pattern():
    g = _tiny_uniform_grid(n=3, cell=100.0)
    omega = 2.0 * np.pi * 1.0
    b = _bc_vector_3d(g, omega, "x")

    n_ex, n_ey, n_ez = _fd3d_edge_counts(g.nx, g.ny, g.nz)
    assert b.shape == (n_ex + n_ey + n_ez,)
    assert np.iscomplexobj(b)

    # Only Ex edges are non-zero for x-polarisation; Ey, Ez fully zero.
    assert np.count_nonzero(b[n_ex:]) == 0

    # Top surface (k=0) is normalised to 1.
    top_vals = b[_ex_idx(np.arange(g.nx), 0, 0, g.ny, g.nz)]
    assert np.allclose(top_vals, 1.0)


def test_bc_vector_3d_y_polarisation_shape_and_pattern():
    g = _tiny_uniform_grid(n=3, cell=100.0)
    omega = 2.0 * np.pi * 1.0
    b = _bc_vector_3d(g, omega, "y")
    n_ex, n_ey, n_ez = _fd3d_edge_counts(g.nx, g.ny, g.nz)

    # Only Ey edges are non-zero for y-polarisation.
    assert np.count_nonzero(b[:n_ex]) == 0
    assert np.count_nonzero(b[n_ex + n_ey :]) == 0

    top_vals = b[_ey_idx(0, np.arange(g.ny), 0, g.nx, g.ny, g.nz)]
    assert np.allclose(top_vals, 1.0)


def test_bc_vector_3d_e_decay_is_monotonically_decaying_with_depth():
    g = _tiny_uniform_grid(n=3, cell=100.0)
    omega = 2.0 * np.pi * 1.0
    b = _bc_vector_3d(g, omega, "x")

    # Side boundary j=0, i=0, over k=0..nz: values are E_decay(z_nodes[k])
    # except k=0 which is fixed to 1 (same starting point anyway).
    mags = np.abs([b[_ex_idx(0, 0, k, g.ny, g.nz)] for k in range(g.nz + 1)])
    assert np.all(np.diff(mags) < 0.0)  # strictly decaying with depth
    assert mags[0] == pytest.approx(1.0)


# ─────────────────────────────────────────────────────────────────────────────
# Dirichlet enforcement
# ─────────────────────────────────────────────────────────────────────────────


def test_apply_dirichlet_3d_overwrites_boundary_rows_and_values():
    g = _tiny_uniform_grid(n=3, cell=100.0)
    omega = 2.0 * np.pi * 1.0

    C = _build_curl_matrix(g)
    sigma = _edge_conductivity(g)
    from pycsamt.forward.em3d import MU0

    M_sig = sparse.diags(sigma, format="csr")
    A = C.T @ C + 1j * omega * MU0 * M_sig
    bc_mask = _boundary_edge_mask(g.nx, g.ny, g.nz)
    bc_idx = np.where(bc_mask)[0]

    b = _bc_vector_3d(g, omega, "x")
    A_mod, b_mod = _apply_dirichlet_3d(A, b, bc_mask, bc_idx)

    assert A_mod.shape == A.shape
    assert b_mod.shape == b.shape

    # Boundary rows become an identity row, and RHS matches the BC value.
    for idx in bc_idx[:5]:
        row = np.asarray(A_mod.getrow(idx).todense()).ravel()
        assert row[idx] == pytest.approx(1.0)
        assert np.count_nonzero(row) == 1
        assert b_mod[idx] == b[idx]

    # Interior rows are untouched in structure (not identity) and the
    # RHS has been updated by subtracting the known boundary contribution.
    int_idx = np.where(~bc_mask)[0]
    assert len(int_idx) > 0
    b_bc = np.zeros_like(b)
    b_bc[bc_idx] = b[bc_idx]
    expected_interior = (b - A @ b_bc)[int_idx]
    assert np.allclose(b_mod[int_idx], expected_interior)
