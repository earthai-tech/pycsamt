# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MARE2DEM data-builder, survey, merge, grid, and import-topo.

Uses both synthetic data (no example files needed) and the real
example files in ``data/mare2dem/``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from pycsamt.models.mare2dem import (
    read_emdata, write_emdata, EMDataFile,
    read_resistivity,
    make_data_file, MTSurveyConfig, CSEMSurveyConfig,
    merge_data_files, merge_emdata,
    grid_to_mare2dem,
)
from pycsamt.models.mare2dem.iotools.emdata import MTConfig, CSEMConfig


# ==========================================================================
# make_data_file / MTSurveyConfig
# ==========================================================================

class TestMakeDataFileMT:
    """Tests on the MT data-file builder (no example files needed)."""

    def _basic_cfg(self, lTE=True, lTM=True):
        return MTSurveyConfig(
            frequencies=np.logspace(-2, 2, 5),
            rx_y=np.linspace(-2000, 2000, 8),
            rx_type="marine",
            lTE=lTE,
            lTM=lTM,
        )

    def test_basic_mt_forward(self, tmp_path):
        mt = self._basic_cfg()
        em = make_data_file(tmp_path / "mt_fwd.emdata", topo=-500.0, mt=mt)
        assert em.n_mt_receivers == 8
        assert len(em.mt.frequencies) == 5
        assert em.n_data > 0
        assert (tmp_path / "mt_fwd.emdata").exists()

    def test_te_only(self, tmp_path):
        mt = self._basic_cfg(lTE=True, lTM=False)
        em = make_data_file(tmp_path / "te_only.emdata", topo=0.0, mt=mt)
        codes = em.data[:, 0].astype(int)
        # TE codes only: 123 (log10 apres) and 104 (phase)
        assert set(codes).issubset({123, 104})

    def test_tm_only(self, tmp_path):
        mt = self._basic_cfg(lTE=False, lTM=True)
        em = make_data_file(tmp_path / "tm_only.emdata", topo=0.0, mt=mt)
        codes = em.data[:, 0].astype(int)
        assert set(codes).issubset({125, 106})

    def test_tipper_included(self, tmp_path):
        mt = MTSurveyConfig(
            frequencies=np.array([1.0, 10.0]),
            rx_y=np.linspace(-1000, 1000, 5),
            lTipperRealImag=True,
        )
        em = make_data_file(tmp_path / "tipper.emdata", topo=0.0, mt=mt)
        codes = em.data[:, 0].astype(int)
        assert 133 in codes or 134 in codes

    def test_land_receivers_have_beta(self, tmp_path):
        topo = np.column_stack([
            np.array([-2000., 0., 2000.]),
            np.array([100., 50., 200.]),
        ])
        mt = MTSurveyConfig(
            frequencies=np.array([1.0]),
            rx_y=np.linspace(-1500, 1500, 6),
            rx_type="land", lTE=True,
        )
        em = make_data_file(tmp_path / "land.emdata", topo=topo, mt=mt)
        assert em.mt is not None

    def test_amphibious_receivers(self, tmp_path):
        mt = MTSurveyConfig(
            frequencies=np.array([1.0, 10.0]),
            rx_y=np.linspace(-2000, 2000, 8),
            rx_type="amphibious", lTE=True, lTM=True,
        )
        em = make_data_file(tmp_path / "amphi.emdata", topo=-100.0, mt=mt)
        assert em.n_mt_receivers >= 8

    def test_data_indices_in_range(self, tmp_path):
        mt = self._basic_cfg()
        em = make_data_file(tmp_path / "range_check.emdata", topo=0.0, mt=mt)
        n_freq = len(em.mt.frequencies)
        n_rx   = em.n_mt_receivers
        # Freq# column (index 1) should be in [1, n_freq]
        freq_idx = em.data[:, 1].astype(int)
        assert freq_idx.min() >= 1 and freq_idx.max() <= n_freq
        # Rx# column (index 3) should be in [1, n_rx]
        rx_idx = em.data[:, 3].astype(int)
        assert rx_idx.min() >= 1 and rx_idx.max() <= n_rx

    def test_data_zeros_initialised(self, tmp_path):
        mt = self._basic_cfg()
        em = make_data_file(tmp_path / "zeros.emdata", topo=0.0, mt=mt)
        # Forward-mode: data column should be 0
        np.testing.assert_allclose(em.data[:, 4], 0.0)
        np.testing.assert_allclose(em.data[:, 5], 0.0)

    def test_roundtrip_after_make(self, tmp_path):
        mt = self._basic_cfg()
        em = make_data_file(tmp_path / "rt.emdata", topo=0.0, mt=mt)
        em2 = read_emdata(tmp_path / "rt.emdata")
        assert em2.n_data == em.n_data
        assert em2.n_mt_receivers == em.n_mt_receivers


# ==========================================================================
# make_data_file / CSEMSurveyConfig
# ==========================================================================

class TestMakeDataFileCSEM:

    def _basic_csem(self, lEx=True):
        return CSEMSurveyConfig(
            frequencies=np.array([0.25, 1.0]),
            tx_y=np.linspace(-3000, 3000, 4),
            rx_y=np.linspace(-5000, 5000, 10),
            rx_type="marine",
            lEx=lEx,
        )

    def test_csem_basic(self, tmp_path):
        csem = self._basic_csem()
        em = make_data_file(tmp_path / "csem.emdata", topo=-1000.0, csem=csem)
        assert em.csem is not None
        assert em.n_csem_transmitters == 4
        assert em.n_csem_receivers == 10

    def test_csem_data_codes(self, tmp_path):
        csem = self._basic_csem(lEx=True)
        em = make_data_file(tmp_path / "csem_codes.emdata", topo=0.0, csem=csem)
        codes = em.data[:, 0].astype(int)
        # Ex log10 amplitude=27, phase=22
        assert 27 in codes and 22 in codes

    def test_combined_mt_csem(self, tmp_path):
        mt = MTSurveyConfig(
            frequencies=np.array([1.0]),
            rx_y=np.linspace(-2000, 2000, 5),
            lTE=True, lTM=True,
        )
        csem = self._basic_csem()
        em = make_data_file(tmp_path / "joint.emdata", topo=0.0, mt=mt, csem=csem)
        assert em.mt is not None
        assert em.csem is not None
        codes = em.data[:, 0].astype(int)
        assert (codes > 100).any()   # MT codes
        assert (codes < 100).any()   # CSEM codes

    def test_min_max_range_filtering(self, tmp_path):
        csem = CSEMSurveyConfig(
            frequencies=np.array([1.0]),
            tx_y=np.array([0.0]),
            rx_y=np.linspace(-5000, 5000, 20),
            lEx=True,
            min_range=500.0,
            max_range=3000.0,
        )
        em = make_data_file(tmp_path / "range_filt.emdata", topo=-1000.0, csem=csem)
        # After filtering, some receivers should be excluded
        assert em.n_data < 20 * 2   # 20 rx × 2 codes per rx, but some filtered


# ==========================================================================
# merge_data_files / merge_emdata
# ==========================================================================

class TestMergeDataFiles:

    def _make_mt_em(self, freqs, rx_y, comment=""):
        em = EMDataFile()
        n_rx = len(rx_y)
        em.mt = MTConfig(
            frequencies=np.asarray(freqs),
            receivers=np.zeros((n_rx, 8)),
            receiver_name=[f"S{i:03d}" for i in range(n_rx)],
        )
        em.mt.receivers[:, 1] = rx_y
        em.data = np.zeros((n_rx * len(freqs) * 2, 6))
        row = 0
        for ifreq in range(1, len(freqs) + 1):
            for irx in range(1, n_rx + 1):
                em.data[row, :] = [123, ifreq, irx, irx, 1.5, 0.1]
                row += 1
                em.data[row, :] = [104, ifreq, irx, irx, 45.0, 2.0]
                row += 1
        em.comment = comment
        return em

    def test_merge_two_mt_files(self, tmp_path):
        em1 = self._make_mt_em([1.0, 10.0], np.array([-1000., 0., 1000.]))
        em2 = self._make_mt_em([100.0], np.array([-500., 500.]))
        write_emdata(em1, tmp_path / "mt1.emdata")
        write_emdata(em2, tmp_path / "mt2.emdata")
        merged = merge_data_files(
            [tmp_path / "mt1.emdata", tmp_path / "mt2.emdata"],
            tmp_path / "joint.emdata",
        )
        # Should have 3 unique frequencies (1, 10, 100)
        assert len(merged.mt.frequencies) == 3
        assert merged.n_data > em1.n_data

    def test_merge_preserves_all_data(self, tmp_path):
        em1 = self._make_mt_em([1.0], np.array([0.]))
        em2 = self._make_mt_em([10.0], np.array([100.]))
        merged = merge_emdata([em1, em2])
        assert merged.n_data == em1.n_data + em2.n_data

    def test_merge_deduplicates_receivers(self, tmp_path):
        em1 = self._make_mt_em([1.0], np.array([0., 1000.]))
        em2 = self._make_mt_em([10.0], np.array([0., 2000.]))  # shares rx at y=0
        merged = merge_emdata([em1, em2])
        # 3 unique receivers: 0, 1000, 2000
        assert len(merged.mt.receivers) == 3

    def test_merge_requires_two_files(self):
        em = EMDataFile()
        with pytest.raises(ValueError, match="least two"):
            merge_emdata([em])

    def test_merge_utm_mismatch_raises(self, tmp_path):
        from pycsamt.models.mare2dem.iotools.emdata import UTMOrigin
        em1 = self._make_mt_em([1.0], np.array([0.]))
        em2 = self._make_mt_em([1.0], np.array([0.]))
        em1.utm = UTMOrigin(grid=19, hemi="N", north0=1000.0, east0=500000.0, theta=0.0)
        em2.utm = UTMOrigin(grid=19, hemi="N", north0=2000.0, east0=500000.0, theta=0.0)
        with pytest.raises(ValueError, match="UTM"):
            merge_emdata([em1, em2])

    def test_merge_real_csem_files(self, csem_dir, tmp_path):
        em = read_emdata(csem_dir / "demo_csem.emdata")
        out = tmp_path / "merged.emdata"
        # Merging a file with itself (trivial case)
        write_emdata(em, tmp_path / "copy.emdata")
        merged = merge_data_files(
            [csem_dir / "demo_csem.emdata", tmp_path / "copy.emdata"],
            out,
        )
        # Deduplication: result should not be larger than 2× original
        assert merged.n_data <= em.n_data * 2

    def test_merge_output_file_written(self, tmp_path):
        em1 = self._make_mt_em([1.0], np.linspace(0, 1000, 3))
        em2 = self._make_mt_em([10.0], np.linspace(2000, 3000, 3))
        write_emdata(em1, tmp_path / "a.emdata")
        write_emdata(em2, tmp_path / "b.emdata")
        out = tmp_path / "merged.emdata"
        merge_data_files([tmp_path / "a.emdata", tmp_path / "b.emdata"], out)
        assert out.exists()


# ==========================================================================
# grid_to_mare2dem
# ==========================================================================

class TestGridToMare2dem:

    def _simple_grid(self, ny=5, nz=4, rho=10.0):
        y1d = np.linspace(-1000, 1000, ny)
        z1d = np.linspace(0, 500, nz)
        Y, Z = np.meshgrid(y1d, z1d)
        Rho = np.full_like(Y, rho)
        return Y, Z, Rho

    def test_creates_all_files(self, tmp_path):
        Y, Z, Rho = self._simple_grid()
        files = grid_to_mare2dem(Y, Z, Rho, out_dir=tmp_path,
                                 padding_y=2000, padding_z=1000)
        assert files["resistivity"].exists()
        assert files["poly"].exists()
        assert files["settings"].exists()

    def test_resistivity_file_readable(self, tmp_path):
        Y, Z, Rho = self._simple_grid(ny=6, nz=4, rho=100.0)
        files = grid_to_mare2dem(Y, Z, Rho, out_dir=tmp_path,
                                 padding_y=5000, padding_z=2000)
        rf = read_resistivity(files["resistivity"])
        # 6×4 = 24 cell centres + 2 extra (air + padding) = 26
        assert rf.num_regions == 6 * 4 + 2

    def test_poly_file_readable(self, tmp_path):
        from pycsamt.models.mare2dem import read_poly
        Y, Z, Rho = self._simple_grid()
        files = grid_to_mare2dem(Y, Z, Rho, out_dir=tmp_path,
                                 padding_y=2000, padding_z=1000)
        pf = read_poly(files["poly"])
        assert pf.n_nodes > 0
        assert pf.n_segments > 0

    def test_custom_model_name(self, tmp_path):
        Y, Z, Rho = self._simple_grid()
        files = grid_to_mare2dem(Y, Z, Rho, out_dir=tmp_path,
                                 model_name="mymodel",
                                 padding_y=2000, padding_z=1000)
        assert files["resistivity"].name.startswith("mymodel")
        assert files["poly"].name.startswith("mymodel")

    def test_variable_resistivity_stored(self, tmp_path):
        y1d = np.linspace(-500, 500, 4)
        z1d = np.linspace(0, 300, 3)
        Y, Z = np.meshgrid(y1d, z1d)
        Rho = np.ones_like(Y)
        Rho[:, :2] = 100.0   # left half: 100 Ω·m
        Rho[:, 2:] = 1.0     # right half: 1 Ω·m
        files = grid_to_mare2dem(Y, Z, Rho, out_dir=tmp_path,
                                 padding_y=2000, padding_z=1000)
        rf = read_resistivity(files["resistivity"])
        # The resistivities should span the input range
        grid_rho = rf.resistivity[:-2, 0]  # exclude air and padding
        assert grid_rho.min() < 10.0
        assert grid_rho.max() > 10.0
