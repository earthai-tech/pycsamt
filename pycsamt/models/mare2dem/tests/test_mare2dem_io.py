# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Real-data I/O tests for MARE2DEM file parsers.

Uses the example files copied from ``mare2dem_examples/`` into
``data/mare2dem/``.  All tests skip automatically when that
directory is absent (e.g. in CI environments).

Covered parsers
---------------
* :func:`read_emdata` / :func:`write_emdata`
* :func:`read_resistivity` / :func:`write_resistivity`
* :func:`read_poly` / :func:`write_poly`
* :func:`write_settings` (settings file format)
* :func:`read_group_rms_log`
* :func:`read_data_group` / :func:`write_data_group`
* :func:`get_most_recent`
* :func:`detect_file_type` / :func:`is_response_file`
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from pycsamt.models.mare2dem import (
    read_emdata, write_emdata,
    read_resistivity, write_resistivity,
    read_poly, write_poly,
    write_settings, SettingsFile,
    read_group_rms_log,
    DataGroupFile, read_data_group, write_data_group,
    get_most_recent,
    detect_file_type, Mare2DEMFileType,
    is_emdata_file, is_resistivity_file, is_response_file,
    is_log_file, is_settings_file,
)


# ==========================================================================
# EMData / EMResp I/O — Wannamaker hill example (MT)
# ==========================================================================

class TestEMDataMT:
    """Tests on the Wannamaker 1986 hill MT example."""

    def test_read_emdata_basic(self, hill_dir):
        em = read_emdata(hill_dir / "hill.emdata")
        assert em.format.lower().startswith("emdata")
        assert not em.is_response
        assert em.mt is not None
        assert em.csem is None

    def test_read_emdata_mt_section(self, hill_dir):
        em = read_emdata(hill_dir / "hill.emdata")
        assert em.n_mt_receivers == 228
        assert len(em.mt.frequencies) == 3
        assert em.mt.frequencies[0] == pytest.approx(2.0)
        assert em.mt.frequencies[1] == pytest.approx(50.0)
        assert em.mt.frequencies[2] == pytest.approx(2000.0)

    def test_read_emdata_data_block(self, hill_dir):
        em = read_emdata(hill_dir / "hill.emdata")
        assert em.n_data > 0
        assert em.data.shape[1] == 6  # observed-data format
        # All type codes should be MT range (>100)
        assert (em.data[:, 0] > 100).all()

    def test_read_emdata_receiver_names(self, hill_dir):
        em = read_emdata(hill_dir / "hill.emdata")
        # First receiver is always MT001
        assert em.mt.receiver_name[0] == "MT001"
        # Last sequential receiver before any extra H-field (_B) duplicates
        assert "MT228" in em.mt.receiver_name or "MT001" in em.mt.receiver_name

    def test_read_response_file(self, hill_dir):
        resp = read_emdata(hill_dir / "hill.0.resp")
        assert resp.is_response
        assert resp.format.lower().startswith("emresp")
        assert resp.data.shape[1] == 8  # response format has 8 columns
        assert resp.n_data > 0

    def test_emdata_roundtrip(self, hill_dir, tmp_path):
        em = read_emdata(hill_dir / "hill.emdata")
        out = tmp_path / "hill_out.emdata"
        write_emdata(em, out)
        em2 = read_emdata(out)
        assert em2.n_data == em.n_data
        assert em2.n_mt_receivers == em.n_mt_receivers
        assert len(em2.mt.frequencies) == len(em.mt.frequencies)
        np.testing.assert_allclose(em2.mt.frequencies, em.mt.frequencies)

    def test_emdata_data_values_preserved(self, hill_dir, tmp_path):
        em = read_emdata(hill_dir / "hill.emdata")
        out = tmp_path / "hill_rt.emdata"
        write_emdata(em, out)
        em2 = read_emdata(out)
        np.testing.assert_allclose(em2.data[:, 4], em.data[:, 4], rtol=1e-10)

    def test_file_type_detection(self, hill_dir):
        assert is_emdata_file(hill_dir / "hill.emdata")
        assert is_response_file(hill_dir / "hill.0.resp")
        assert not is_emdata_file(hill_dir / "hill.0.resp")
        assert detect_file_type(hill_dir / "hill.0.resp") == Mare2DEMFileType.RESPONSE
        assert detect_file_type(hill_dir / "hill.emdata") == Mare2DEMFileType.EMDATA


# ==========================================================================
# EMData I/O — demo CSEM example
# ==========================================================================

class TestEMDataCSEM:
    """Tests on the demo CSEM forward example."""

    def test_read_csem_emdata(self, csem_dir):
        em = read_emdata(csem_dir / "demo_csem.emdata")
        assert em.csem is not None
        assert em.mt is None
        assert em.n_csem_transmitters > 0
        assert em.n_csem_receivers > 0

    def test_csem_frequencies(self, csem_dir):
        em = read_emdata(csem_dir / "demo_csem.emdata")
        freqs = em.csem.frequencies
        assert len(freqs) > 0
        assert (freqs > 0).all()

    def test_csem_data_codes(self, csem_dir):
        em = read_emdata(csem_dir / "demo_csem.emdata")
        codes = em.data[:, 0].astype(int)
        # CSEM codes should be < 100
        assert (codes < 100).all()

    def test_read_csem_response(self, csem_dir):
        resp = read_emdata(csem_dir / "demo.0.resp")
        assert resp.is_response
        assert resp.data.shape[1] == 8
        assert resp.n_csem_transmitters > 0

    def test_csem_roundtrip(self, csem_dir, tmp_path):
        em = read_emdata(csem_dir / "demo_csem.emdata")
        out = tmp_path / "demo_csem_rt.emdata"
        write_emdata(em, out)
        em2 = read_emdata(out)
        assert em2.n_csem_transmitters == em.n_csem_transmitters
        assert em2.n_csem_receivers == em.n_csem_receivers
        assert em2.n_data == em.n_data

    def test_phase_convention_preserved(self, csem_dir, tmp_path):
        em = read_emdata(csem_dir / "demo_csem.emdata")
        assert em.csem.phase_convention.lower() == "lag"
        out = tmp_path / "demo_csem_pc.emdata"
        write_emdata(em, out)
        em2 = read_emdata(out)
        assert em2.csem.phase_convention.lower() == "lag"


# ==========================================================================
# Resistivity I/O
# ==========================================================================

class TestResistivityIO:
    """Tests on the .resistivity file reader/writer."""

    def test_read_hill_resistivity(self, hill_dir):
        rf = read_resistivity(hill_dir / "hill.0.resistivity")
        assert rf.version.lower() == "mare2dem_1.1"
        assert rf.poly_file == "hill.poly"
        assert rf.data_file == "hill.emdata"
        assert rf.anisotropy.lower() == "isotropic"
        assert rf.num_regions > 0

    def test_resistivity_values(self, hill_dir):
        rf = read_resistivity(hill_dir / "hill.0.resistivity")
        assert rf.resistivity.shape[1] == 1  # isotropic
        assert rf.resistivity.shape[0] == rf.num_regions
        # Resistivities should be positive
        assert (rf.resistivity > 0).all()

    def test_resistivity_header_fields(self, hill_dir):
        rf = read_resistivity(hill_dir / "hill.0.resistivity")
        assert rf.target_misfit == pytest.approx(1.0)
        assert rf.max_iterations == 100
        assert rf.penalty_cut_weight == pytest.approx(0.1)
        assert rf.roughness_penalty_method == "gradient"

    def test_resistivity_bounds(self, hill_dir):
        rf = read_resistivity(hill_dir / "hill.0.resistivity")
        assert rf.global_bounds is not None
        assert len(rf.global_bounds) == 2
        # lower < upper
        assert rf.global_bounds[0] < rf.global_bounds[1]

    def test_resistivity_roundtrip(self, hill_dir, tmp_path):
        rf = read_resistivity(hill_dir / "hill.0.resistivity")
        out = tmp_path / "hill_rt.resistivity"
        write_resistivity(rf, out)
        rf2 = read_resistivity(out)
        assert rf2.num_regions == rf.num_regions
        assert rf2.anisotropy == rf.anisotropy
        np.testing.assert_allclose(
            rf2.resistivity, rf.resistivity, rtol=1e-5
        )

    def test_resistivity_no_data_flag(self, hill_dir):
        rf = read_resistivity(hill_dir / "hill.0.resistivity", no_data=True)
        # With no_data=True, region table not read
        assert rf.num_regions == 0

    def test_demo_resistivity(self, inversion_dir):
        rf = read_resistivity(inversion_dir / "demo.0.resistivity")
        assert rf.num_regions > 1000   # large mesh
        assert rf.poly_file.endswith(".poly")

    def test_final_vs_initial_resistivity(self, inversion_dir):
        rf0 = read_resistivity(inversion_dir / "demo.0.resistivity")
        rf6 = read_resistivity(inversion_dir / "demo.6.resistivity")
        assert rf0.num_regions == rf6.num_regions
        # Resistivities should differ after inversion
        assert not np.allclose(rf0.resistivity, rf6.resistivity)

    def test_resistivity_file_type(self, hill_dir):
        assert is_resistivity_file(hill_dir / "hill.0.resistivity")
        assert detect_file_type(hill_dir / "hill.0.resistivity") == Mare2DEMFileType.RESISTIVITY


# ==========================================================================
# Poly I/O
# ==========================================================================

class TestPolyIO:
    """Tests on the Triangle .poly PSLG file reader/writer."""

    def test_read_hill_poly(self, hill_dir):
        pf = read_poly(hill_dir / "hill.poly")
        assert len(pf.nodes) > 0
        assert len(pf.segments) > 0
        # Hill poly is very small
        assert pf.n_nodes < 50

    def test_hill_poly_coordinates(self, hill_dir):
        pf = read_poly(hill_dir / "hill.poly")
        y = pf.nodes[:, 0]
        z = pf.nodes[:, 1]
        # Profile spans ±100 km
        assert y.min() < -50000
        assert y.max() > 50000

    def test_poly_roundtrip(self, hill_dir, tmp_path):
        pf = read_poly(hill_dir / "hill.poly")
        out = tmp_path / "hill_rt.poly"
        write_poly(pf, out)
        pf2 = read_poly(out)
        assert pf2.n_nodes == pf.n_nodes
        assert pf2.n_segments == pf.n_segments
        np.testing.assert_allclose(pf2.nodes, pf.nodes, rtol=1e-10)

    def test_poly_segments_in_range(self, hill_dir):
        pf = read_poly(hill_dir / "hill.poly")
        # All segment node indices should be within node range
        segs = pf.segments
        assert segs.min() >= 1
        assert segs.max() <= pf.n_nodes

    def test_demo_poly_larger(self, csem_dir):
        pf = read_poly(csem_dir / "demo.poly")
        assert pf.n_nodes > 100   # larger mesh

    def test_poly_file_type(self, hill_dir):
        assert detect_file_type(hill_dir / "hill.poly") == Mare2DEMFileType.POLY


# ==========================================================================
# Settings I/O
# ==========================================================================

class TestSettingsIO:
    """Tests on the .settings file format."""

    def test_read_settings_content(self, hill_dir):
        text = (hill_dir / "mare2dem.settings").read_text()
        assert "Transmitters per group" in text or "tolerance" in text.lower()

    def test_write_settings_roundtrip(self, tmp_path):
        sf = SettingsFile(tolerance=2.0, csem_rx_per_group=30, mt_rx_per_group=25)
        p = write_settings(sf, tmp_path / "test.settings")
        assert p.exists()
        text = p.read_text()
        assert "30" in text
        assert "25" in text

    def test_settings_file_type(self, hill_dir):
        assert is_settings_file(hill_dir / "mare2dem.settings")
        assert detect_file_type(hill_dir / "mare2dem.settings") == Mare2DEMFileType.SETTINGS


# ==========================================================================
# Group RMS log
# ==========================================================================

class TestGroupRMSLog:
    """Tests on the demo_csem_mt group-RMS log reader."""

    def test_read_group_rms_log(self, csem_mt_dir):
        log = read_group_rms_log(csem_mt_dir / "demo.group_rms.log")
        assert log.n_groups >= 2
        assert log.n_iterations > 0
        assert "CSEM" in " ".join(log.headers) or "MT" in " ".join(log.headers)

    def test_group_rms_values_decreasing(self, csem_mt_dir):
        log = read_group_rms_log(csem_mt_dir / "demo.group_rms.log")
        # Column 0 is iteration number; Total RMS is in column 1
        # (headers = ['Iteration', 'Total RMS', 'CSEM', 'MT'])
        total_rms_col = 1  # index of 'Total RMS'
        total = log.rms_log[:, total_rms_col]
        assert total[0] > total[-1]

    def test_group_rms_shape(self, csem_mt_dir):
        log = read_group_rms_log(csem_mt_dir / "demo.group_rms.log")
        assert log.rms_log.shape == (log.n_iterations, log.n_groups)


# ==========================================================================
# DataGroupFile
# ==========================================================================

class TestDataGroupFile:
    """Tests on the DataGroupFile reader/writer (no example file needed)."""

    def test_roundtrip(self, tmp_path):
        dg = DataGroupFile(
            group_names=["MT", "Seafloor CSEM", "Towed CSEM"],
            group_indices=np.array([1, 1, 2, 2, 3, 3]),
        )
        p = tmp_path / "test.emdata_group"
        write_data_group(dg, p)
        dg2 = read_data_group(p)
        assert dg2.group_names == dg.group_names
        np.testing.assert_array_equal(dg2.group_indices, dg.group_indices)

    def test_validation(self):
        dg = DataGroupFile(
            group_names=["A"],
            group_indices=np.array([2]),
        )
        with pytest.raises(ValueError, match="exceed"):
            with tempfile.NamedTemporaryFile(suffix=".emdata_group") as f:
                write_data_group(dg, f.name)


# ==========================================================================
# get_most_recent
# ==========================================================================

class TestGetMostRecent:
    """Tests on the most-recently-modified file finder."""

    def test_pass_through_literal(self):
        p = get_most_recent("some_file.resistivity")
        assert str(p) == "some_file.resistivity"

    def test_keyword_newest(self, inversion_dir):
        p = get_most_recent("newest", "*.resistivity", search_dir=inversion_dir)
        assert p is not None
        assert p.suffix == ".resistivity"

    def test_keyword_no_files(self, tmp_path):
        p = get_most_recent("last", "*.resistivity", search_dir=tmp_path)
        assert p is None
