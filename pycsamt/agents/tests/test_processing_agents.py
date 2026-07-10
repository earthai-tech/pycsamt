# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Integration tests for the four processing agents exposed through the desktop
AgentWorker (no LLM, no API key required).

  * QC Quicklook       → plot_qc_quicklook
  * Dimensionality     → classify_dimensionality  +  plot_dimensionality_psection
  * Static Shift (fast)→ correct_static_shift     +  plot_rho_phase_bode

Each test loads real EDI data from:
  • TIPPER  (2 stations)  — data/AMT/TIPPER/
  • WILLY   (28 stations) — data/AMT/WILLY_DATA/L18PLT/

Tests skip gracefully when the data directories are absent (CI / minimal env).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")          # headless backend — must come before pyplot

import matplotlib.pyplot as plt
import pytest

# ── project root ──────────────────────────────────────────────────────────────
_HERE         = Path(__file__).resolve().parent       # agents/tests/
_PROJECT_ROOT = _HERE.parents[2]                      # repo root (pycsamt/)
_TIPPER_DIR   = _PROJECT_ROOT / "data" / "AMT" / "TIPPER"
_WILLY_DIR    = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def tipper_dir() -> Path:
    if not _TIPPER_DIR.exists() or not list(_TIPPER_DIR.glob("*.edi")):
        pytest.skip(f"TIPPER EDI dataset not found: {_TIPPER_DIR}")
    return _TIPPER_DIR


@pytest.fixture(scope="session")
def willy_dir() -> Path:
    if not _WILLY_DIR.exists() or not list(_WILLY_DIR.glob("*.edi")):
        pytest.skip(f"WILLY EDI dataset not found: {_WILLY_DIR}")
    return _WILLY_DIR


@pytest.fixture(scope="session")
def tipper_sites(tipper_dir):
    from pycsamt.emtools import ensure_sites
    return ensure_sites(tipper_dir)


@pytest.fixture(scope="session")
def willy_sites(willy_dir):
    from pycsamt.emtools import ensure_sites
    return ensure_sites(willy_dir)


@pytest.fixture(autouse=True)
def _close_figures():
    """Close all open figures after each test to prevent resource leaks."""
    yield
    plt.close("all")


# ══════════════════════════════════════════════════════════════════════════════
# 1. QC Quicklook
# ══════════════════════════════════════════════════════════════════════════════

class TestQCQuicklook:
    """plot_qc_quicklook — must return a matplotlib Figure."""

    def test_returns_figure_tipper(self, tipper_sites):
        import pycsamt.emtools as et
        fig = et.plot_qc_quicklook(tipper_sites, verbose=0)
        assert fig is not None
        assert hasattr(fig, "savefig"), "Expected a matplotlib Figure"

    def test_returns_figure_willy(self, willy_sites):
        import pycsamt.emtools as et
        fig = et.plot_qc_quicklook(willy_sites, verbose=0)
        assert fig is not None
        assert hasattr(fig, "savefig")

    def test_figure_has_axes(self, tipper_sites):
        import pycsamt.emtools as et
        fig = et.plot_qc_quicklook(tipper_sites, verbose=0)
        assert len(fig.axes) >= 2, "Expected at least 2 axes (coverage + SNR)"

    def test_accepts_path_string(self, tipper_dir):
        import pycsamt.emtools as et
        fig = et.plot_qc_quicklook(str(tipper_dir), verbose=0)
        assert hasattr(fig, "savefig")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Dimensionality
# ══════════════════════════════════════════════════════════════════════════════

class TestDimensionality:
    """classify_dimensionality → DataFrame; plot_dimensionality_psection → Axes."""

    def test_classify_returns_dataframe_like(self, tipper_sites):
        import pycsamt.emtools as et
        result = et.classify_dimensionality(tipper_sites, verbose=0)
        assert result is not None
        # result is either a DataFrame or a wrapped frame — both have len()
        assert len(result) >= 0

    def test_classify_has_dim_column(self, willy_sites):
        import pandas as pd

        import pycsamt.emtools as et
        result = et.classify_dimensionality(willy_sites, verbose=0)
        df = getattr(result, "df", result)
        if isinstance(df, pd.DataFrame):
            assert "dim" in df.columns

    def test_classify_skew_threshold(self, willy_sites):
        import pycsamt.emtools as et
        r1 = et.classify_dimensionality(willy_sites, skew_th=2.0, verbose=0)
        r2 = et.classify_dimensionality(willy_sites, skew_th=8.0, verbose=0)
        assert r1 is not None and r2 is not None

    def test_psection_returns_axes(self, willy_sites):
        import pycsamt.emtools as et
        ax = et.plot_dimensionality_psection(willy_sites, verbose=0)
        assert ax is not None
        assert hasattr(ax, "get_figure")

    def test_psection_tipper(self, tipper_sites):
        import pycsamt.emtools as et
        ax = et.plot_dimensionality_psection(tipper_sites, verbose=0)
        assert hasattr(ax, "get_figure")

    def test_psection_has_parent_figure(self, willy_sites):
        import pycsamt.emtools as et
        ax = et.plot_dimensionality_psection(willy_sites, verbose=0)
        fig = ax.get_figure()
        assert fig is not None

    def test_psection_accepts_custom_ax(self, willy_sites):
        import pycsamt.emtools as et
        _, given_ax = plt.subplots()
        returned = et.plot_dimensionality_psection(
            willy_sites, ax=given_ax, verbose=0
        )
        assert returned is given_ax


# ══════════════════════════════════════════════════════════════════════════════
# 3. Static Shift (fast)
# ══════════════════════════════════════════════════════════════════════════════

class TestStaticShift:
    """correct_static_shift → Sites; plot_rho_phase_bode → Figure."""

    def test_correction_returns_sites(self, willy_sites):
        import pycsamt.emtools as et
        from pycsamt.site.base import Sites
        result = et.correct_static_shift(willy_sites, verbose=0)
        assert result is not None
        assert isinstance(result, Sites), (
            f"Expected a Sites object, got {type(result)}"
        )

    def test_correction_not_inplace(self, willy_sites):
        import pycsamt.emtools as et
        corrected = et.correct_static_shift(willy_sites, inplace=False, verbose=0)
        assert corrected is not willy_sites, "inplace=False must return a new object"

    def test_correction_tipper(self, tipper_sites):
        import pycsamt.emtools as et
        corrected = et.correct_static_shift(tipper_sites, verbose=0)
        assert corrected is not None

    def test_bode_plot_original(self, willy_sites):
        import pycsamt.emtools as et
        fig = et.plot_rho_phase_bode(willy_sites, verbose=0)
        assert fig is not None
        assert hasattr(fig, "savefig")

    def test_bode_plot_corrected(self, willy_sites):
        import pycsamt.emtools as et
        corrected = et.correct_static_shift(willy_sites, verbose=0)
        fig = et.plot_rho_phase_bode(corrected, verbose=0)
        assert fig is not None
        assert hasattr(fig, "savefig")

    def test_bode_plot_tipper(self, tipper_sites):
        import pycsamt.emtools as et
        fig = et.plot_rho_phase_bode(tipper_sites, verbose=0)
        assert hasattr(fig, "savefig")

    def test_bode_figure_has_axes(self, willy_sites):
        import pycsamt.emtools as et
        fig = et.plot_rho_phase_bode(willy_sites, verbose=0)
        assert len(fig.axes) >= 2, "Bode plot expects rho + phase axes"

    def test_correction_different_components(self, willy_sites):
        import pycsamt.emtools as et
        for comp in ("det", "xy", "yx"):
            result = et.correct_static_shift(willy_sites, comp=comp, verbose=0)
            assert result is not None, f"comp={comp!r} returned None"

    def test_pipeline_correct_then_plot(self, willy_sites):
        """The full pipeline that Static Shift (fast) processing agent runs."""
        import pycsamt.emtools as et
        corrected = et.correct_static_shift(willy_sites, verbose=0)
        fig = et.plot_rho_phase_bode(corrected, verbose=0)
        assert hasattr(fig, "savefig")


# ══════════════════════════════════════════════════════════════════════════════
# 5. interpolate_grid guard (unit test for the fix)
# ══════════════════════════════════════════════════════════════════════════════

class TestInterpolateGridGuard:
    """Verify the all-NaN guard in interpolate_grid doesn't crash."""

    def test_all_nan_2d_does_not_raise(self):
        import numpy as np

        from pycsamt.utils.arrayops import interpolate_grid

        arr = np.full((5, 8), np.nan)
        result = interpolate_grid(arr, fill_value="auto", view=False)
        assert result is not None
        assert result.shape == arr.shape

    def test_sparse_array_does_not_raise(self):
        import numpy as np

        from pycsamt.utils.arrayops import interpolate_grid

        arr = np.full((6, 10), np.nan)
        arr[2, 3] = 1.0   # only one finite value
        result = interpolate_grid(arr, fill_value="auto", view=False)
        assert result is not None
        assert result.shape == arr.shape

    def test_well_populated_still_interpolates(self):
        import numpy as np

        from pycsamt.utils.arrayops import interpolate_grid

        rng = np.random.default_rng(0)
        arr = rng.uniform(1, 10, (8, 12)).astype(float)
        arr[1, 2] = np.nan
        result = interpolate_grid(arr, fill_value="auto", view=False)
        assert result is not None
        assert np.isfinite(result[1, 2]), "Interpolated NaN should be filled"

    def test_1d_all_nan_does_not_raise(self):
        import numpy as np

        from pycsamt.utils.arrayops import interpolate_grid

        arr = np.full(8, np.nan)
        result = interpolate_grid(arr, fill_value="auto", view=False)
        assert result is not None
        assert result.shape == arr.shape


# ══════════════════════════════════════════════════════════════════════════════
# 6. AgentWorker integration (requires PySide6)
# ══════════════════════════════════════════════════════════════════════════════

_pyside = pytest.importorskip("PySide6", reason="PySide6 required for AgentWorker tests")


@pytest.fixture(scope="session")
def qapp():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    existing = QApplication.instance()
    if existing is not None:
        yield existing
        return
    app = QApplication(["pytest", "-platform", "offscreen"])
    yield app


def _run_worker(agent_name, sites, params=None):
    """Run AgentWorker synchronously; return (result, errors, log_lines)."""
    from pycsamt.app.desktop.workers.agent_worker import (
        AgentWorker,
    )
    results, errors, logs = [], [], []
    w = AgentWorker(agent_name=agent_name, sites=sites, params=params or {})
    w.result.connect(results.append)
    w.error.connect(errors.append)
    w.log_line.connect(logs.append)
    w.run()
    return results, errors, logs


class TestAgentWorkerProcessing:

    def test_qc_quicklook_worker_no_error(self, qapp, willy_sites):
        results, errors, _ = _run_worker("QC Quicklook", willy_sites)
        assert not errors, f"Worker emitted errors: {errors}"
        assert results, "Worker must emit a result"

    def test_qc_quicklook_worker_returns_figure(self, qapp, willy_sites):
        results, errors, _ = _run_worker("QC Quicklook", willy_sites)
        assert not errors
        res = results[0]
        assert hasattr(res, "savefig") or hasattr(res, "get_figure"), (
            f"Expected Figure or Axes, got {type(res)}"
        )

    def test_dimensionality_worker_no_error(self, qapp, willy_sites):
        results, errors, _ = _run_worker("Dimensionality", willy_sites)
        assert not errors, f"Worker emitted errors: {errors}"
        assert results

    def test_dimensionality_worker_emits_axes(self, qapp, willy_sites):
        results, errors, _ = _run_worker("Dimensionality", willy_sites)
        assert not errors
        res = results[0]
        assert hasattr(res, "get_figure") or hasattr(res, "savefig"), (
            f"Expected Figure or Axes from Dimensionality worker, got {type(res)}"
        )

    def test_static_shift_worker_no_error(self, qapp, willy_sites):
        results, errors, _ = _run_worker("Static Shift (fast)", willy_sites)
        assert not errors, f"Worker emitted errors: {errors}"
        assert results

    def test_static_shift_worker_emits_figure(self, qapp, willy_sites):
        """Regression: Static Shift (fast) must emit a renderable figure."""
        results, errors, _ = _run_worker("Static Shift (fast)", willy_sites)
        assert not errors
        res = results[0]
        assert hasattr(res, "savefig") or hasattr(res, "get_figure"), (
            f"Static Shift result must be Figure/Axes, got {type(res)}"
        )

    def test_worker_logs_function_name(self, qapp, tipper_sites):
        _, _, logs = _run_worker("QC Quicklook", tipper_sites)
        assert any("plot_qc_quicklook" in line for line in logs), (
            "Worker must log the emtools function name"
        )
