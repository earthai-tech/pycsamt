# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MARE2DEM log parser, InversionResult scanner, diff, and noise.

Uses the example files in ``data/mare2dem/``.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.models.mare2dem import (
    InversionResult,
    Mare2DEMLog,
    diff_resistivity,
    read_emdata,
)
from pycsamt.models.mare2dem.noise import (
    NoiseConfig,
    add_synthetic_noise,
)

# ==========================================================================
# Mare2DEMLog — OccamLog.2012.0 format
# ==========================================================================


class TestMare2DEMLog:
    """Tests on the real MARE2DEM log file from demo_mt_inversion."""

    def test_read_log_basic(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        assert log.n_iterations > 0

    def test_log_iteration_count(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        # The demo inversion ran 6 iterations
        assert log.n_iterations == 6

    def test_log_first_rms(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        # First iteration misfit should be >> 1 (data not yet fit)
        assert log.iterations[0].rms > 2.0

    def test_log_rms_decreasing(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        rms = log.rms_history()
        # RMS should monotonically decrease (or at least trend down)
        assert rms[0] > rms[-1]

    def test_log_final_rms_near_target(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        # Demo targets RMS=1.0; final should be close
        assert log.final_rms < 1.5

    def test_log_roughness_populated(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        for rec in log.iterations:
            assert rec.roughness >= 0.0

    def test_log_lambda_populated(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        for rec in log.iterations:
            assert rec.lambda_ > 0.0

    def test_log_repr(self, inversion_dir):
        log = Mare2DEMLog(inversion_dir / "demo.logfile")
        r = repr(log)
        assert "n_iterations" in r
        assert "final_rms" in r

    def test_log_nonexistent_file(self, tmp_path):
        log = Mare2DEMLog(tmp_path / "nonexistent.logfile")
        assert log.n_iterations == 0
        assert log.final_rms is None


# ==========================================================================
# InversionResult — scan directory
# ==========================================================================


class TestInversionResult:
    """Tests on InversionResult scanning the demo_mt_inversion directory."""

    def test_scan_inversion_dir(self, inversion_dir):
        result = InversionResult(inversion_dir)
        assert result.log is not None
        assert result.model is not None
        assert result.data is not None

    def test_inversion_result_log_populated(self, inversion_dir):
        result = InversionResult(inversion_dir)
        assert result.n_iterations == 6
        assert result.final_rms is not None
        assert result.final_rms < 1.5

    def test_inversion_result_model_loaded(self, inversion_dir):
        result = InversionResult(inversion_dir)
        # Should load demo.6.resistivity (most recent)
        assert result.model is not None

    def test_inversion_result_response_loaded(self, inversion_dir):
        result = InversionResult(inversion_dir)
        # demo.6.resp should be detected
        assert result.response is not None

    def test_inversion_result_summary(self, inversion_dir):
        result = InversionResult(inversion_dir)
        summary = result.summary()
        assert "converged" in summary.lower() or "rms" in summary.lower()

    def test_inversion_result_nonexistent_dir(self, tmp_path):
        result = InversionResult(tmp_path / "empty_run")
        assert result.log is None
        assert result.model is None

    def test_inversion_result_empty_dir(self, tmp_path):
        result = InversionResult(tmp_path)
        assert result.log is None


# ==========================================================================
# diff_resistivity
# ==========================================================================


class TestDiffResistivity:
    """Tests on the resistivity model diff utility."""

    def test_diff_same_model(self, inversion_dir, tmp_path):
        out = tmp_path / "diff_zero.resistivity"
        dm = diff_resistivity(
            inversion_dir / "demo.0.resistivity",
            inversion_dir / "demo.0.resistivity",
            out,
        )
        assert dm.num_regions > 0
        # log10(x) - log10(x) = 0
        np.testing.assert_allclose(dm.resistivity, 0.0, atol=1e-10)

    def test_diff_initial_vs_final(self, inversion_dir, tmp_path):
        out = tmp_path / "diff_0_vs_6.resistivity"
        dm = diff_resistivity(
            inversion_dir / "demo.0.resistivity",
            inversion_dir / "demo.6.resistivity",
            out,
        )
        assert dm.num_regions > 0
        # Models differ, so diff should be non-zero for at least some regions
        assert not np.allclose(dm.resistivity, 0.0)

    def test_diff_output_file_exists(self, inversion_dir, tmp_path):
        out = tmp_path / "diff_test.resistivity"
        diff_resistivity(
            inversion_dir / "demo.0.resistivity",
            inversion_dir / "demo.6.resistivity",
            out,
        )
        assert out.exists()

    def test_diff_mismatched_regions_raises(self, hill_dir, inversion_dir, tmp_path):
        with pytest.raises(ValueError, match="regions"):
            diff_resistivity(
                hill_dir / "hill.0.resistivity",
                inversion_dir / "demo.0.resistivity",
                tmp_path / "fail.resistivity",
            )

    def test_diff_custom_fn(self, inversion_dir, tmp_path):
        out = tmp_path / "diff_pct.resistivity"
        dm = diff_resistivity(
            inversion_dir / "demo.0.resistivity",
            inversion_dir / "demo.6.resistivity",
            out,
            diff_fn=lambda A, B: np.abs((A - B) / A * 100),
        )
        assert (dm.resistivity >= 0).all()


# ==========================================================================
# add_synthetic_noise
# ==========================================================================


class TestAddSyntheticNoise:
    """Tests on add_synthetic_noise using the hill response file."""

    def test_noise_produces_6col_data(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        nc = NoiseConfig(mt_rel_noise=0.05)
        noisy = add_synthetic_noise(resp, nc, seed=0)
        assert noisy.data.shape[1] == 6
        assert not noisy.is_response

    def test_noise_data_count_le_response(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        nc = NoiseConfig(mt_rel_noise=0.05, mt_abs_noise_tipper=0.01)
        noisy = add_synthetic_noise(resp, nc, seed=0)
        # NaN rows dropped (unhandled codes like 135/136) → count ≤ input
        assert noisy.n_data <= resp.n_data
        assert noisy.n_data > 0

    def test_noise_mt_section_preserved(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        nc = NoiseConfig(mt_rel_noise=0.05)
        noisy = add_synthetic_noise(resp, nc, seed=0)
        assert noisy.mt is not None
        assert len(noisy.mt.frequencies) == len(resp.mt.frequencies)

    def test_noise_std_errors_positive(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        nc = NoiseConfig(mt_rel_noise=0.05)
        noisy = add_synthetic_noise(resp, nc, seed=0)
        assert (noisy.data[:, 5] > 0).all()

    def test_noise_reproducible_with_seed(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        nc = NoiseConfig(mt_rel_noise=0.05)
        n1 = add_synthetic_noise(resp, nc, seed=42)
        n2 = add_synthetic_noise(resp, nc, seed=42)
        np.testing.assert_array_equal(n1.data, n2.data)

    def test_noise_different_each_call_no_seed(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        nc = NoiseConfig(mt_rel_noise=0.05)
        n1 = add_synthetic_noise(resp, nc)
        n2 = add_synthetic_noise(resp, nc)
        assert not np.allclose(n1.data[:, 4], n2.data[:, 4])

    def test_noise_requires_8col_data(self):
        import numpy as np

        from pycsamt.models.mare2dem import EMDataFile

        em = EMDataFile()
        em.data = np.zeros((5, 6))  # 6-column data, not a response file
        nc = NoiseConfig()
        with pytest.raises(ValueError, match="8-column"):
            add_synthetic_noise(em, nc)

    def test_noise_csem_response(self, csem_dir):
        resp = read_emdata(csem_dir / "demo.0.resp")
        nc = NoiseConfig(csem_rel_noise_e=0.05, csem_rel_noise_b=0.05)
        noisy = add_synthetic_noise(resp, nc, seed=1)
        assert noisy.n_data > 0
        assert noisy.data.shape[1] == 6
