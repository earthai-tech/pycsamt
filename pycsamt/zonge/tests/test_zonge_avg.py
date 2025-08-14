# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the BaseAVG and AVG classes.
"""
import pytest
import pandas as pd
import numpy as np

# Assuming the package structure is set up correctly
from pycsamt.zonge.avg import BaseAVG, AVG
from pycsamt.zonge.info import DataInfo

# --- Test Data Fixtures -------------------------------------------

@pytest.fixture
def modern_avg_content():
    """
    Provides comprehensive content for a modern (kind-2) AVG file.
    """
    return (
        '\\ AMTAVG 7.76: "LCS01.fld", Dated 99-01-01\n'
        "\\ ASTATIC v3.60d updated data on 22/07/16\n"
        "$Survey.Type=CSAMT\n"
        "$Unit.E=nV/Am\n"
        "$Unit.B=pT/A\n"
        "$Unit.Phase=mrad\n"
        "\n"
        "$Rx.Stn=25\n"
        "$Rx.Cmp=ExHy\n"
        "Z.mwgt,Freq,Tx.Amp,ARes.mag,Z.phz,E.mag,B.mag,"
        "E.%err,B.%err,ARes.%err,E.perr,B.perr,Z.perr\n"
        "1,1024,12,1220.7,618.3,223.12,0.0892,3.3,3.7,1.5,"
        "22.4,22.8,16.1\n"
        "1,512,13,1454,651.1,260.67,0.1351,2.0,24.6,120.9,"
        "20.0,45.8,803\n"
        "\n"
        "$Rx.Stn=75\n"
        "$Rx.Cmp=ExHy\n"
        "Z.mwgt,Freq,Tx.Amp,ARes.mag,Z.phz,E.mag,B.mag,"
        "E.%err,B.%err,ARes.%err,E.perr,B.perr,Z.perr\n"
        "1,1024,12,876.62,618.3,150.0,0.0912,1.8,3.5,7.7,"
        "21.5,83.5,148.9\n"
        "1,512,13,990.45,365,177.88,0.1351,2.6,24.6,120.4,"
        "24.8,45.8,878\n"
    )

@pytest.fixture
def legacy_avg_content():
    """Provides content for a legacy (kind-1) AVG file."""
    return (
        "\\ AMTAVG 7.20: \"SAMCSAM.FLD\", Processed 29 Jul 93\n"
        "$ ASPACE= 183.0m\n"
        "$ XMTR  =    1.\n"
        " skp Station Freq  Comp Amps     Emag     Ephz      "
        "Hmag     Hphz  Resistivity   Phase\n"
        "\\-++------++----++---++----++---------++------++"
        "---------++------++---------++------+\n"
        "  2    100.0  8192 ExHy  4.5   1.1731e+3  1491.0  "
        "3.5150e-1   722.7  2.7195e+2   768.3\n"
        "  2    100.0  4096 ExHy  8.0   8.5835e+2  2087.0  "
        "3.9003e-1  1244.7  2.3648e+2   842.4\n"
        "  2    200.0  8192 ExHy  4.5   1.2513e+3 -2111.5  "
        "3.6216e-1 -3036.5  2.9142e+2   925.0\n"
    )

@pytest.fixture
def modern_avg_file(tmp_path, modern_avg_content):
    """Creates a temporary modern AVG file and returns its path."""
    f = tmp_path / "modern.avg"
    f.write_text(modern_avg_content)
    return f

@pytest.fixture
def legacy_avg_file(tmp_path, legacy_avg_content):
    """Creates a temporary legacy AVG file and returns its path."""
    f = tmp_path / "legacy.avg"
    f.write_text(legacy_avg_content)
    return f

# --- Tests for BaseAVG --------------------------------------------

class TestBaseAVG:
    def test_initialization(self):
        """Test that BaseAVG initializes correctly."""
        base = BaseAVG(verbose=True)
        assert isinstance(base.info, DataInfo)
        assert base.verbose is True
        assert base._kind is None
        assert base._source_path is None

    def test_read_from_modern_file(self, modern_avg_file):
        """Test reading a modern AVG file."""
        base = BaseAVG()
        base.read(modern_avg_file)
        assert base._kind == 2
        assert base.info.df is not None
        assert len(base.info.df) == 4
        assert "rho" in base.info.df.columns
        assert "station" in base.info.df.columns

    def test_read_from_legacy_file(self, legacy_avg_file):
        """Test reading and transforming a legacy AVG file."""
        base = BaseAVG(verbose=True)
        base.read(legacy_avg_file)
        assert base._kind == 1
        assert base.info.df is not None
        assert len(base.info.df) == 3
        # Check for a transformed column
        assert "rho" in base.info.df.columns
        # Check for an injected modern column
        assert "z.mwgt" in base.info.df.columns

    def test_read_from_dataframe(self):
        """Test reading directly from a pandas DataFrame."""
        df = pd.DataFrame({
            "station": [100, 100],
            "freq": [1024, 512],
            "rho": [50, 60],
            "phase": [45, 55],
            "comp": ["ExHy", "ExHy"]
        })
        meta = {"Survey.Type": "CSAMT"}
        base = BaseAVG()
        base.read(df, meta=meta)
        assert base._kind == 2  # Assumes modern
        pd.testing.assert_frame_equal(base.info.df, df)
        assert base.info.meta["Survey.Type"] == "CSAMT"

    def test_write_modern(self, modern_avg_file, tmp_path):
        """Test writing to a modern AVG file."""
        base = BaseAVG()
        base.read(modern_avg_file)
        out_path = tmp_path / "output_modern.avg"
        base.to_modern(out_path)
        assert out_path.exists()

        # Verify content by reading it back
        reloaded_base = BaseAVG()
        reloaded_base.read(out_path)
        assert len(reloaded_base.info.df) == 4
        assert reloaded_base.info.header.config.survey_type == "CSAMT"

    def test_write_legacy(self, modern_avg_file, tmp_path):
        """Test writing to a legacy AVG file."""
        base = BaseAVG()
        base.read(modern_avg_file)
        out_path = tmp_path / "output_legacy.avg"
        base.to_legacy(out_path)
        assert out_path.exists()

        # Verify content by reading it back
        reloaded_base = BaseAVG()
        reloaded_base.read(out_path)
        assert reloaded_base._kind == 1
        assert len(reloaded_base.info.df) == 4

    def test_write_dispatcher(self, modern_avg_file, tmp_path):
        """Test the main `write` method dispatch logic."""
        base = BaseAVG()
        base.read(modern_avg_file)
        
        # Test modern dispatch
        path_modern = tmp_path / "disp_modern.avg"
        base.write(path_modern, fmt="modern")
        assert path_modern.exists()

        # Test legacy dispatch
        path_legacy = tmp_path / "disp_legacy.avg"
        base.write(path_legacy, fmt="legacy")
        assert path_legacy.exists()

        # Test auto (default to modern)
        path_auto = tmp_path / "disp_auto.avg"
        base.write(path_auto)
        assert path_auto.exists()
        
        # Check auto content
        reloaded_base = BaseAVG()
        reloaded_base.read(path_auto)
        assert reloaded_base._kind == 2


class TestAVG:
    def test_inheritance(self):
        """Test that AVG inherits from BaseAVG."""
        assert issubclass(AVG, BaseAVG)

    def test_from_file_factory(self, modern_avg_file):
        """Test the primary `from_file` classmethod."""
        avg = AVG.from_file(modern_avg_file, verbose=True)
        assert isinstance(avg, AVG)
        assert avg.info.df is not None
        assert len(avg.info.df) == 4
        assert avg.header.config.survey_type == "CSAMT"
        assert avg.station.n_unique == 2

    def test_properties_are_populated(self, modern_avg_file):
        """Test that component properties are accessible."""
        avg = AVG.from_file(modern_avg_file)
        assert avg.df is not None
        assert avg.header is not None
        assert avg.station is not None
        assert avg.z is not None
        assert avg.resistivity is not None
        assert avg.phase is not None
        assert avg.frequency is not None

    def test_to_xarray(self, modern_avg_file):
        """Test exporting data to an xarray.Dataset."""
        avg = AVG.from_file(modern_avg_file)
        ds = avg.to_xarray(include_qc=False, include_z=False)
        
        assert "xarray" in str(type(ds))
        assert "rho" in ds.data_vars
        assert "phase" in ds.data_vars
        assert "station" in ds.coords
        assert "freq" in ds.coords
        assert ds.dims["station"] == 2
        assert ds.dims["freq"] == 2

    def test_to_tensor(self, modern_avg_file):
        """Test exporting a variable to a numpy tensor."""
        avg = AVG.from_file(modern_avg_file)
        # Test for resistivity
        T, freqs, stations = avg.to_tensor(var="rho")
        
        assert isinstance(T, np.ndarray)
        assert T.shape == (2, 2, 2, 2)  # (st, f, e, h)
        assert freqs.shape == (2,)
        assert stations.shape == (2,)

if __name__=='__main__': # pragma: no-cover 
   pytest.main( [__file__])