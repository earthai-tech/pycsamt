# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

"""
Pytest suite for the BaseAVG and AVG classes.
"""

import numpy as np
import pandas as pd
import pytest

# Assuming the package structure is set up correctly
from pycsamt.zonge.avg import AVG, BaseAVG
from pycsamt.zonge.info import DataInfo

# Note: The fixtures `modern_data_file` and `legacy_data_file`
# are now defined in `tests/conftest.py` and are automatically
# available to all tests.


@pytest.fixture
def _modern_avg_content():
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


def _legacy_avg_content():
    """Provides content for a legacy (kind-1) AVG file."""
    return (
        '\\ AMTAVG 7.20: "SAMCSAM.FLD", Processed 29 Jul 93\n'
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


class TestBaseAVG:
    def test_initialization(self):
        """Test that BaseAVG initializes correctly."""
        base = BaseAVG(verbose=True)
        assert isinstance(base.info, DataInfo)
        assert base.verbose is True
        assert base._kind is None
        assert base._source_path is None

    def test_read_from_modern_file(self, modern_data_file):
        """Test reading a modern AVG file from a real data file."""
        base = BaseAVG()
        base.read(modern_data_file)
        assert base._kind == 2
        assert base.info.df is not None
        assert len(base.info.df) > 0
        assert "rho" in base.info.df.columns
        assert "station" in base.info.df.columns
        assert "pc_emag" in base.info.pc_emag.frame.columns

    def test_read_from_legacy_file(self, legacy_data_file):
        """
        Test reading and transforming a legacy AVG file from a
        real data file.
        """
        base = BaseAVG(verbose=True)
        base.read(legacy_data_file)
        assert base._kind == 1
        assert base.info.df is not None
        assert len(base.info.df) > 0
        assert "rho" in base.info.df.columns
        assert "z_mwgt" in base.info.df.columns
        assert "pc_hmag" in base.info.pc_hmag.frame.columns

    def test_read_from_dataframe(self):
        """Test reading directly from a pandas DataFrame."""
        # This test remains important for in-memory data handling
        df = pd.DataFrame(
            {
                "station": [100, 100],
                "freq": [1024, 512],
                "rho": [50, 60],
                "phase": [45, 55],
                "comp": ["ExHy", "ExHy"],
                "E.%err": [1, 2],
                "B.%err": [3, 4],
                "ARes.%err": [5, 6],
                "E.perr": [7, 8],
                "H.perr": [9, 10],
                "Z.perr": [11, 12],
            }
        )
        meta = {"Survey.Type": "CSAMT"}
        base = BaseAVG()
        base.read(df, meta=meta)
        assert base._kind == 2  # Assumes modern
        out = base.info.df
        # the reader canonicalises QC column names on ingest
        renames = {
            "E.%err": "pc_emag",
            "B.%err": "pc_hmag",
            "ARes.%err": "pc_rho",
            "E.perr": "s_ephz",
            "H.perr": "s_hphz",
            "Z.perr": "s_phz",
        }
        for raw, canonical in renames.items():
            assert canonical in out.columns
            assert out[canonical].tolist() == df[raw].tolist()
        for col in ("station", "freq", "rho", "phase", "comp"):
            assert out[col].tolist() == df[col].tolist()
        assert base.info.meta["Survey.Type"] == "CSAMT"

    def test_write_modern(self, modern_data_file, tmp_path):
        """Test writing to a modern AVG file."""
        base = BaseAVG()
        base.read(modern_data_file)
        original_rows = len(base.info.df)
        out_path = tmp_path / "output_modern.avg"
        base.to_modern(out_path)
        assert out_path.exists()

        reloaded_base = BaseAVG()
        reloaded_base.read(out_path)
        assert len(reloaded_base.info.df) == original_rows
        assert reloaded_base.info.header.config.survey_type == "CSAMT"

    def test_write_legacy(self, modern_data_file, tmp_path):
        """Test writing to a legacy AVG file."""
        base = BaseAVG()
        base.read(modern_data_file)
        original_rows = len(base.info.df)
        out_path = tmp_path / "output_legacy.avg"
        base.to_legacy(out_path)
        assert out_path.exists()

        reloaded_base = BaseAVG()
        reloaded_base.read(out_path)
        assert reloaded_base._kind == 1
        assert np.isclose(
            len(reloaded_base.info.df), original_rows * 9.62962962962963
        )
        # assert len(reloaded_base.info.df) == original_rows

    def test_write_dispatcher(self, modern_data_file, tmp_path):
        """Test the main `write` method dispatch logic."""
        base = BaseAVG()
        base.read(modern_data_file)

        path_modern = tmp_path / "disp_modern.avg"
        base.write(path_modern, fmt="modern")
        assert path_modern.exists()

        path_legacy = tmp_path / "disp_legacy.avg"
        base.write(path_legacy, fmt="legacy")
        assert path_legacy.exists()

        path_auto = tmp_path / "disp_auto.avg"
        base.write(path_auto)
        assert path_auto.exists()

        reloaded_base = BaseAVG()
        reloaded_base.read(path_auto)
        assert reloaded_base._kind == 2


# --- Tests for AVG ------------------------------------------------


class TestAVG:
    def test_inheritance(self):
        """Test that AVG inherits from BaseAVG."""
        assert issubclass(AVG, BaseAVG)

    def test_from_file_factory(self, modern_data_file):
        """Test the primary `from_file` classmethod."""
        avg = AVG.from_file(modern_data_file, verbose=True)
        assert isinstance(avg, AVG)
        assert avg.info.df is not None
        assert len(avg.info.df) > 0
        assert avg.header.config.survey_type == "CSAMT"
        assert avg.station.n_unique > 0

    def test_properties_are_populated(self, modern_data_file):
        """Test that component properties are accessible."""
        avg = AVG.from_file(modern_data_file)
        assert avg.df is not None
        assert avg.header is not None
        assert avg.station is not None
        assert avg.z is not None
        assert avg.resistivity is not None
        assert avg.phase is not None
        assert avg.frequency is not None

    def test_to_xarray(self, modern_data_file):
        """Test exporting data to an xarray.Dataset."""
        avg = AVG.from_file(modern_data_file)
        ds = avg.to_xarray(include_qc=True, include_z=True)

        assert "xarray" in str(type(ds))
        assert "rho" in ds.data_vars
        assert "phase" in ds.data_vars
        assert "pc_emag" in ds.data_vars
        assert "station" in ds.coords
        assert "freq" in ds.coords
        assert ds.dims["station"] > 0
        assert ds.dims["freq"] > 0

    def test_to_tensor(self, modern_data_file):
        """Test exporting a variable to a numpy tensor."""
        avg = AVG.from_file(modern_data_file)
        T, freqs, stations = avg.to_tensor(var="rho")

        assert isinstance(T, np.ndarray)
        assert T.ndim == 4  # (st, f, e, h)
        assert T.shape[0] == avg.station.n_unique
        assert T.shape[1] == avg.frequency.n_unique
        assert T.shape[2:] == (2, 2)
        assert freqs.shape == (avg.frequency.n_unique,)
        assert stations.shape == (avg.station.n_unique,)


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
