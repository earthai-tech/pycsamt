from types import SimpleNamespace

import numpy as np

from pycsamt.models.occam2d.results import _build_rho_2d


def test_build_rho_2d_leaves_air_rows_unassigned():
    mesh = SimpleNamespace(n_zcells=3, n_xcells=2, n_airlayers=1)
    model = SimpleNamespace(
        layers=[
            {"n_merge": 1, "params": np.array([2])},
            {"n_merge": 1, "params": np.array([2])},
        ]
    )

    grid = _build_rho_2d(np.array([1.5, 2.5]), model, mesh)

    assert np.isnan(grid[0]).all()
    np.testing.assert_allclose(grid[1], [1.5, 1.5])
    np.testing.assert_allclose(grid[2], [2.5, 2.5])
