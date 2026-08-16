# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.models.occam2d.synthetic."""

from __future__ import annotations

import numpy as np

from pycsamt.forward.em2d import MT2DForward
from pycsamt.forward.grid2d import Grid2D
from pycsamt.models.occam2d.data import OccamData
from pycsamt.models.occam2d.synthetic import SyntheticSite, sites_from_response


def _response():
    grid = Grid2D.halfspace(rho=100.0, nx=20, x_max=2000.0, n_stations=5)
    resp = MT2DForward(np.array([200.0, 20.0]), grid, verbose=False).run()
    return grid, resp


def test_sites_from_response_shapes_and_names():
    grid, resp = _response()
    names = [f"S{i}" for i in range(grid.n_stations)]
    sites = sites_from_response(resp, grid.x_stations, names)
    assert len(sites) == grid.n_stations
    assert [s.name for s in sites] == names
    for s in sites:
        assert s.rho.shape == (2, 2, 2)
        assert s.phase.shape == (2, 2, 2)
        assert s.rho_err.shape == (2, 2, 2)
        assert s.phase_err.shape == (2, 2, 2)


def test_sites_from_response_only_off_diagonal_is_populated():
    grid, resp = _response()
    names = [f"S{i}" for i in range(grid.n_stations)]
    sites = sites_from_response(resp, grid.x_stations, names)
    s = sites[0]
    assert np.all(s.rho[:, 0, 0] == 0.0)
    assert np.all(s.rho[:, 1, 1] == 0.0)
    assert np.all(s.rho[:, 0, 1] > 0.0)
    assert np.all(s.rho[:, 1, 0] > 0.0)


def test_sites_from_response_error_fractions():
    grid, resp = _response()
    names = [f"S{i}" for i in range(grid.n_stations)]
    sites = sites_from_response(resp, grid.x_stations, names, rho_err_frac=0.1, phase_err_deg=3.0)
    s = sites[0]
    np.testing.assert_allclose(s.rho_err[:, 0, 1], s.rho[:, 0, 1] * 0.1)
    np.testing.assert_allclose(s.phase_err[:, 0, 1], 3.0)


def test_synthetic_site_coords_recover_chainage():
    site = SyntheticSite("S0", 555.0, [1.0], np.zeros((1, 2, 2)), np.zeros((1, 2, 2)),
                          np.zeros((1, 2, 2)), np.zeros((1, 2, 2)))
    lat, lon, elev = site.coords
    assert lat == 0.0
    assert lon * 111_195.0 == 555.0
    assert elev == 0.0


def test_sites_from_response_accepted_by_occam_data_normalise_source():
    # OccamData._normalise_source falls back to list(source) for anything
    # without ._items/.edic -- confirm the real contract accepts these
    # duck-typed sites without raising.
    grid, resp = _response()
    names = [f"S{i}" for i in range(grid.n_stations)]
    sites = sites_from_response(resp, grid.x_stations, names)
    data = OccamData.from_edi(sites, title="test")
    assert data.n_sites == grid.n_stations
    assert data.n_data > 0
