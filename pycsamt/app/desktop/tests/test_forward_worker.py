# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ForwardWorker (pycsamt.app.desktop.workers.forward_worker).

Strategy
--------
Unlike the AI inversion worker, the analytic/numeric forward solvers
(``MT1DForward`` / ``MT2DForward`` / ``MT3DForward``) run in well under a
second even for small real grids, so these tests call them for real
(small nx/nz/n_freq/n_stations) rather than mocking — this exercises the
worker's own grid-building, padding, and anomaly-injection arithmetic
against the real solver contracts.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.workers.forward_worker import ForwardWorker


def _run(params):
    w = ForwardWorker(params)
    results = []
    errors = []
    progresses = []
    w.finished.connect(results.append)
    w.error.connect(errors.append)
    w.progress.connect(progresses.append)
    w.run()
    return results, errors, progresses


# ── 1D ────────────────────────────────────────────────────────────────────


class Test1D:
    def test_run_1d_no_noise(self):
        params = {
            "dim": "1D",
            "resistivity": [100.0, 1000.0],
            "thickness": [300.0],
            "n_freq": 6,
            "f_min": 1e-1,
            "f_max": 1e2,
            "noise": "none",
        }
        results, errors, progresses = _run(params)
        assert errors == []
        assert len(results) == 1
        assert progresses[-1] == 100
        assert results[0].rho_a.shape == (6,)

    def test_run_1d_gaussian_noise(self):
        params = {
            "dim": "1D",
            "resistivity": [100.0, 1000.0],
            "thickness": [300.0],
            "n_freq": 5,
            "f_min": 1e-1,
            "f_max": 1e1,
            "noise": "gaussian",
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1

    def test_run_1d_multiplicative_noise(self):
        params = {
            "dim": "1D",
            "resistivity": [100.0, 1000.0],
            "thickness": [300.0],
            "n_freq": 5,
            "f_min": 1e-1,
            "f_max": 1e1,
            "noise": "multiplicative",
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1

    def test_run_1d_missing_resistivity_reports_error(self):
        results, errors, _ = _run({"dim": "1D", "thickness": [300.0]})
        assert results == []
        assert len(errors) == 1

    def test_run_1d_defaults_used_when_params_sparse(self):
        # thickness/resistivity are the only truly required keys
        results, errors, _ = _run(
            {"dim": "1D", "resistivity": [50.0, 500.0], "thickness": [200.0]}
        )
        assert errors == []
        assert len(results) == 1


# ── 2D ────────────────────────────────────────────────────────────────────


class Test2D:
    def test_run_2d_no_anomaly(self):
        params = {
            "dim": "2D",
            "nx": 4,
            "nz": 4,
            "dx": 500.0,
            "dz_min": 50.0,
            "dz_max": 1000.0,
            "n_pad": 2,
            "bg_rho": 100.0,
            "n_freq": 3,
            "f_min": 1e-1,
            "f_max": 1e1,
            "n_stations": 3,
            "anomaly": False,
        }
        results, errors, progresses = _run(params)
        assert errors == []
        assert len(results) == 1
        assert progresses[-1] == 100

    def test_run_2d_with_anomaly(self):
        params = {
            "dim": "2D",
            "nx": 6,
            "nz": 4,
            "dx": 500.0,
            "dz_min": 50.0,
            "dz_max": 1000.0,
            "n_pad": 2,
            "bg_rho": 100.0,
            "n_freq": 3,
            "f_min": 1e-1,
            "f_max": 1e1,
            "n_stations": 3,
            "anomaly": True,
            "anom_x": 1500.0,
            "anom_z": 200.0,
            "anom_w": 1000.0,
            "anom_h": 500.0,
            "anom_rho": 10.0,
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1

    def test_run_2d_anomaly_outside_grid_no_injection(self):
        """Anomaly window entirely outside the grid -> ix/iz empty, the
        `if ix.size and iz.size:` guard skips injection without raising."""
        params = {
            "dim": "2D",
            "nx": 4,
            "nz": 4,
            "dx": 500.0,
            "dz_min": 50.0,
            "dz_max": 1000.0,
            "n_pad": 2,
            "n_freq": 3,
            "f_min": 1e-1,
            "f_max": 1e1,
            "n_stations": 3,
            "anomaly": True,
            "anom_x": -1e6,
            "anom_z": 1e6,
            "anom_w": 10.0,
            "anom_h": 10.0,
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1

    def test_run_2d_exception_reports_error(self, monkeypatch):
        import pycsamt.forward as fwd_pkg

        def _boom(*a, **k):
            raise RuntimeError("grid build failed")

        monkeypatch.setattr(fwd_pkg, "Grid2D", _boom)
        results, errors, _ = _run({"dim": "2D"})
        assert results == []
        assert errors == ["grid build failed"]


# ── 3D ────────────────────────────────────────────────────────────────────


class Test3D:
    def test_run_3d_no_anomaly(self):
        params = {
            "dim": "3D",
            "nx": 4,
            "ny": 4,
            "nz": 3,
            "dx": 1000.0,
            "dy": 1000.0,
            "n_pad": 2,
            "bg_rho": 100.0,
            "n_freq": 2,
            "f_min": 1e-1,
            "f_max": 1e0,
            "n_sx": 2,
            "n_sy": 2,
            "anomaly": False,
        }
        results, errors, progresses = _run(params)
        assert errors == []
        assert len(results) == 1
        assert progresses[-1] == 100

    def test_run_3d_with_anomaly(self):
        params = {
            "dim": "3D",
            "nx": 4,
            "ny": 4,
            "nz": 3,
            "dx": 1000.0,
            "dy": 1000.0,
            "n_pad": 2,
            "bg_rho": 100.0,
            "n_freq": 2,
            "f_min": 1e-1,
            "f_max": 1e0,
            "n_sx": 2,
            "n_sy": 2,
            "anomaly": True,
            "anom_x": 2000.0,
            "anom_y": 2000.0,
            "anom_z": 300.0,
            "anom_rho": 10.0,
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1

    def test_run_3d_anomaly_outside_grid_no_injection(self):
        params = {
            "dim": "3D",
            "nx": 4,
            "ny": 4,
            "nz": 3,
            "dx": 1000.0,
            "dy": 1000.0,
            "n_pad": 2,
            "n_freq": 2,
            "f_min": 1e-1,
            "f_max": 1e0,
            "n_sx": 2,
            "n_sy": 2,
            "anomaly": True,
            "anom_x": -1e6,
            "anom_y": -1e6,
            "anom_z": 1e6,
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1

    def test_run_3d_exception_reports_error(self, monkeypatch):
        import pycsamt.forward as fwd_pkg

        def _boom(*a, **k):
            raise RuntimeError("3D grid failed")

        monkeypatch.setattr(fwd_pkg.Grid3D, "halfspace", staticmethod(_boom))
        results, errors, _ = _run({"dim": "3D"})
        assert results == []
        assert errors == ["3D grid failed"]


class TestDimDispatch:
    def test_unknown_dim_defaults_to_3d_path(self):
        """dim not in {'1D','2D'} falls through to the else (3D) branch."""
        params = {
            "dim": "weird",
            "nx": 4,
            "ny": 4,
            "nz": 3,
            "n_pad": 2,
            "n_freq": 2,
            "n_sx": 2,
            "n_sy": 2,
        }
        results, errors, _ = _run(params)
        assert errors == []
        assert len(results) == 1


# ── Anomaly placement (survey coordinates) ─────────────────────────────────


def _anomaly_centre(nodes, n_pad, idx):
    import numpy as np

    centres = 0.5 * (nodes[:-1] + nodes[1:])
    return float(np.mean(centres[idx])) - nodes[n_pad]


class TestAnomalyPlacement:
    """Regression: the anomaly centre is given in survey coordinates
    (0..nx*dx), but was applied from the outer edge of the padding. In 3-D
    the default anomaly was then never inserted at all; in 2-D it landed
    ~5.9 km left of where it was asked."""

    def test_3d_default_anomaly_is_inserted_where_asked(self):
        import numpy as np

        w = ForwardWorker({
            "dim": "3D", "nx": 12, "ny": 12, "nz": 6, "dx": 1000.0,
            "dy": 1000.0, "n_pad": 4, "n_freq": 1, "n_sx": 2, "n_sy": 2,
            "anomaly": True, "anom_x": 6000.0, "anom_y": 6000.0,
            "anom_z": 400.0, "anom_rho": 5.0,
        })
        w.progress.connect(lambda *_: None)
        g = w._run_3d().grid
        iz, iy, ix = np.where(g.resistivity == 5.0)
        assert iz.size > 0
        assert _anomaly_centre(g.x_nodes, g.n_pad, ix) == pytest.approx(6000.0)
        assert _anomaly_centre(g.y_nodes, g.n_pad, iy) == pytest.approx(6000.0)

    def test_2d_anomaly_is_centred_where_asked(self):
        import numpy as np

        w = ForwardWorker({
            "dim": "2D", "nx": 30, "nz": 10, "dx": 500.0, "n_pad": 5,
            "n_freq": 1, "n_stations": 3, "anomaly": True,
            "anom_x": 7500.0, "anom_z": 500.0, "anom_w": 2000.0,
            "anom_rho": 5.0,
        })
        w.progress.connect(lambda *_: None)
        g = w._run_2d().grid
        iz, ix = np.where(g.resistivity == 5.0)
        xn = np.concatenate([[0.0], np.cumsum(g.dx)])
        assert iz.size > 0
        assert _anomaly_centre(xn, g.n_pad, ix) == pytest.approx(7500.0)

    def test_3d_anomaly_far_outside_grid_is_not_clamped_in(self):
        import numpy as np

        w = ForwardWorker({
            "dim": "3D", "nx": 4, "ny": 4, "nz": 3, "n_pad": 2, "n_freq": 1,
            "n_sx": 2, "n_sy": 2, "anomaly": True, "anom_x": -1e6,
            "anom_y": -1e6, "anom_z": 300.0, "anom_rho": 5.0,
        })
        w.progress.connect(lambda *_: None)
        g = w._run_3d().grid
        assert not np.any(g.resistivity == 5.0)
